from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
RESULTS_DIR = ROOT / "results"


@dataclass
class LinearFit:
    intercept: float
    slope: float
    intercept_q16: float
    intercept_q84: float
    slope_q16: float
    slope_q84: float


def ensure_dirs() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def percentile_rank(series: pd.Series, higher_better: bool) -> pd.Series:
    ranked = series.rank(pct=True, ascending=higher_better)
    return ranked.clip(1e-3, 1 - 1e-3)


def valley_radius_from_series(
    radius_series: pd.Series,
    *,
    min_points: int = 60,
    bandwidth: float = 0.18,
) -> dict[str, float] | None:
    radius = radius_series.dropna()
    radius = radius[(radius > 0.9) & (radius < 4.2)]
    if len(radius) < min_points:
        return None

    log_radius = np.log10(radius.to_numpy())
    grid = np.linspace(np.log10(1.0), np.log10(4.0), 512)
    density = gaussian_kde(log_radius, bw_method=bandwidth)(grid)
    peaks, _ = find_peaks(density)

    left_band = np.where(grid < np.log10(1.9))[0]
    right_band = np.where(grid > np.log10(1.9))[0]

    left_candidates = peaks[np.isin(peaks, left_band)]
    right_candidates = peaks[np.isin(peaks, right_band)]

    if len(left_candidates):
        left_peak = left_candidates[np.argmax(density[left_candidates])]
    else:
        left_peak = left_band[np.argmax(density[left_band])]

    if len(right_candidates):
        right_peak = right_candidates[np.argmax(density[right_candidates])]
    else:
        right_peak = right_band[np.argmax(density[right_band])]

    if left_peak > right_peak:
        left_peak, right_peak = right_peak, left_peak

    valley_index = left_peak + np.argmin(density[left_peak : right_peak + 1])
    return {
        "valley_radius": float(10 ** grid[valley_index]),
        "left_peak_radius": float(10 ** grid[left_peak]),
        "right_peak_radius": float(10 ** grid[right_peak]),
        "sample_size": int(len(radius)),
    }


def clean_confirmed(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df.copy()
    cleaned = cleaned[(cleaned["st_rad"] < 1.8) & cleaned["st_teff"].between(3200, 7200)]
    cleaned = cleaned[cleaned["pl_rade"].between(0.8, 4.5)]
    cleaned = cleaned[cleaned["pl_orbper"].between(0.3, 100)]
    cleaned = cleaned[cleaned["pl_insol"].between(0.5, 3000)]
    cleaned = cleaned.dropna(subset=["pl_name", "hostname"])
    cleaned["log_period"] = np.log10(cleaned["pl_orbper"])
    cleaned["log_radius"] = np.log10(cleaned["pl_rade"])
    cleaned["log_insol"] = np.log10(cleaned["pl_insol"])
    cleaned["log_teff"] = np.log10(cleaned["st_teff"])
    return cleaned.reset_index(drop=True)


def clean_toi(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df.copy()
    cleaned = cleaned[(cleaned["st_rad"] < 1.8) & cleaned["st_teff"].between(3200, 7200)]
    cleaned = cleaned[cleaned["pl_rade"].between(0.8, 4.5)]
    cleaned = cleaned[cleaned["pl_orbper"].between(0.3, 100)]
    cleaned = cleaned[cleaned["pl_insol"].between(0.5, 3000)]
    cleaned["log_period"] = np.log10(cleaned["pl_orbper"])
    cleaned["log_radius"] = np.log10(cleaned["pl_rade"])
    cleaned["log_insol"] = np.log10(cleaned["pl_insol"])
    cleaned["log_teff"] = np.log10(cleaned["st_teff"])
    return cleaned.reset_index(drop=True)


def estimate_binned_valleys(
    df: pd.DataFrame,
    *,
    feature_column: str,
    edges: np.ndarray,
    label_prefix: str,
    min_points: int = 60,
) -> pd.DataFrame:
    rows: list[dict[str, float | str | int]] = []
    for idx in range(len(edges) - 1):
        low = edges[idx]
        high = edges[idx + 1]
        if idx < len(edges) - 2:
            mask = df[feature_column].ge(low) & df[feature_column].lt(high)
        else:
            mask = df[feature_column].ge(low) & df[feature_column].le(high)
        subset = df.loc[mask]
        valley = valley_radius_from_series(subset["pl_rade"], min_points=min_points)
        row: dict[str, float | str | int] = {
            "bin_label": f"{label_prefix}_{idx + 1}",
            "bin_index": idx + 1,
            "n_planets": int(len(subset)),
            f"{feature_column}_low": float(low),
            f"{feature_column}_high": float(high),
            f"{feature_column}_median": float(subset[feature_column].median()),
        }
        if valley is None:
            row |= {
                "valley_radius": np.nan,
                "left_peak_radius": np.nan,
                "right_peak_radius": np.nan,
            }
        else:
            row |= valley
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap_flux_fit(
    df: pd.DataFrame,
    edges: np.ndarray,
    *,
    n_bootstrap: int,
    seed: int,
) -> LinearFit:
    base_table = estimate_binned_valleys(
        df,
        feature_column="log_insol",
        edges=edges,
        label_prefix="flux",
    ).dropna(subset=["valley_radius"])
    base_coeff = np.polyfit(base_table["log_insol_median"], np.log10(base_table["valley_radius"]), 1)

    rng = np.random.default_rng(seed)
    coeffs: list[np.ndarray] = []
    for _ in range(n_bootstrap):
        sample = df.iloc[rng.integers(0, len(df), len(df))]
        fit_table = estimate_binned_valleys(
            sample,
            feature_column="log_insol",
            edges=edges,
            label_prefix="flux",
        ).dropna(subset=["valley_radius"])
        if len(fit_table) < 4:
            continue
        coeffs.append(
            np.polyfit(fit_table["log_insol_median"], np.log10(fit_table["valley_radius"]), 1)
        )

    coeff_array = np.vstack(coeffs)
    slope_q16, slope_q84 = np.quantile(coeff_array[:, 0], [0.16, 0.84])
    intercept_q16, intercept_q84 = np.quantile(coeff_array[:, 1], [0.16, 0.84])
    return LinearFit(
        intercept=float(base_coeff[1]),
        slope=float(base_coeff[0]),
        intercept_q16=float(intercept_q16),
        intercept_q84=float(intercept_q84),
        slope_q16=float(slope_q16),
        slope_q84=float(slope_q84),
    )


def bootstrap_metallicity_delta(
    df: pd.DataFrame,
    edges: np.ndarray,
    *,
    n_bootstrap: int,
    seed: int,
) -> dict[str, float]:
    valid = df.dropna(subset=["st_met"]).copy()
    tertiles = estimate_binned_valleys(
        valid,
        feature_column="st_met",
        edges=edges,
        label_prefix="metallicity",
        min_points=120,
    ).dropna(subset=["valley_radius"])
    base_delta = float(tertiles.iloc[-1]["valley_radius"] - tertiles.iloc[0]["valley_radius"])

    rng = np.random.default_rng(seed)
    deltas: list[float] = []
    for _ in range(n_bootstrap):
        sample = valid.iloc[rng.integers(0, len(valid), len(valid))]
        boot = estimate_binned_valleys(
            sample,
            feature_column="st_met",
            edges=edges,
            label_prefix="metallicity",
            min_points=120,
        ).dropna(subset=["valley_radius"])
        if len(boot) < 3:
            continue
        deltas.append(float(boot.iloc[-1]["valley_radius"] - boot.iloc[0]["valley_radius"]))

    q16, median, q84 = np.quantile(deltas, [0.16, 0.5, 0.84])
    return {
        "base_delta": base_delta,
        "q16": float(q16),
        "median": float(median),
        "q84": float(q84),
    }


def science_case(row: pd.Series) -> str:
    if row["valley_distance_dex"] < 0.045:
        return "radius_valley_probe"
    if row["pl_rade"] < 1.8 and row["pl_orbper"] < 1.2:
        return "ultra_short_period_rocky"
    if 1.8 <= row["pl_rade"] <= 3.2 and 10 <= row["pl_insol"] <= 200:
        return "temperate_sub_neptune"
    return "general_demographic_gain"


def rank_toi_candidates(
    confirmed: pd.DataFrame,
    toi: pd.DataFrame,
    fit: LinearFit,
) -> pd.DataFrame:
    feature_columns = ["log_period", "log_radius", "log_insol", "log_teff"]
    feature_mean = confirmed[feature_columns].mean()
    feature_std = confirmed[feature_columns].std()

    confirmed_scaled = (confirmed[feature_columns] - feature_mean) / feature_std
    toi_scaled = (toi[feature_columns] - feature_mean) / feature_std

    tree = cKDTree(confirmed_scaled.to_numpy())
    neighbors = min(20, len(confirmed_scaled))
    distances, _ = tree.query(toi_scaled.to_numpy(), k=neighbors)
    if distances.ndim == 1:
        distances = distances[:, None]

    ranked = toi.copy()
    ranked["predicted_valley_radius"] = 10 ** (
        fit.intercept + fit.slope * ranked["log_insol"]
    )
    ranked["valley_distance_dex"] = np.abs(
        np.log10(ranked["pl_rade"]) - np.log10(ranked["predicted_valley_radius"])
    )
    ranked["novelty_raw"] = distances.mean(axis=1)
    ranked["novelty_pct"] = percentile_rank(ranked["novelty_raw"], higher_better=True)
    ranked["brightness_pct"] = percentile_rank(ranked["st_tmag"], higher_better=False)
    ranked["valley_pct"] = percentile_rank(ranked["valley_distance_dex"], higher_better=False)
    ranked["depth_pct"] = percentile_rank(ranked["pl_trandep"], higher_better=True)
    ranked["priority_score"] = 100 * (
        ranked["novelty_pct"] ** 0.35
        * ranked["brightness_pct"] ** 0.30
        * ranked["valley_pct"] ** 0.25
        * ranked["depth_pct"] ** 0.10
    )
    ranked["science_case"] = ranked.apply(science_case, axis=1)
    return ranked.sort_values("priority_score", ascending=False).reset_index(drop=True)


def save_candidate_views(ranked: pd.DataFrame) -> None:
    columns = [
        "toi",
        "toipfx",
        "st_tmag",
        "ra",
        "dec",
        "pl_rade",
        "pl_orbper",
        "pl_insol",
        "pl_eqt",
        "pl_trandep",
        "predicted_valley_radius",
        "valley_distance_dex",
        "novelty_raw",
        "priority_score",
        "science_case",
    ]
    ranked[columns].to_csv(RESULTS_DIR / "top_candidates_overall.csv", index=False)
    ranked.loc[ranked["science_case"] == "radius_valley_probe", columns].head(25).to_csv(
        RESULTS_DIR / "top_radius_valley_probes.csv",
        index=False,
    )
    ranked.loc[ranked["science_case"] == "temperate_sub_neptune", columns].head(25).to_csv(
        RESULTS_DIR / "top_temperate_sub_neptunes.csv",
        index=False,
    )


def plot_valley_overview(
    confirmed: pd.DataFrame,
    flux_table: pd.DataFrame,
    metallicity_table: pd.DataFrame,
    fit: LinearFit,
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax = axes[0]
    metal = confirmed["st_met"].fillna(confirmed["st_met"].median())
    scatter = ax.scatter(
        confirmed["pl_insol"],
        confirmed["pl_rade"],
        c=metal,
        cmap="coolwarm",
        alpha=0.22,
        s=18,
        edgecolors="none",
    )
    flux_valid = flux_table.dropna(subset=["valley_radius"]).copy()
    ax.plot(
        10 ** flux_valid["log_insol_median"],
        flux_valid["valley_radius"],
        "o",
        color="black",
        label="Flux-bin valley estimates",
    )
    log_grid = np.linspace(confirmed["log_insol"].min(), confirmed["log_insol"].max(), 200)
    ax.plot(
        10 ** log_grid,
        10 ** (fit.intercept + fit.slope * log_grid),
        color="#c43d2f",
        linewidth=2.2,
        label="Fitted valley track",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Insolation Flux [S_earth]")
    ax.set_ylabel("Planet Radius [R_earth]")
    ax.set_title("Radius Valley vs Insolation")
    ax.legend(loc="upper left")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Host [Fe/H]")

    ax = axes[1]
    metallicity_valid = metallicity_table.dropna(subset=["valley_radius"]).copy()
    ax.plot(
        metallicity_valid["st_met_median"],
        metallicity_valid["valley_radius"],
        "-o",
        color="#1f5f8b",
        linewidth=2.2,
    )
    ax.set_xlabel("Host [Fe/H]")
    ax.set_ylabel("Valley Radius [R_earth]")
    ax.set_title("Metallicity Tertiles")
    ax.axhline(1.8, color="#999999", linestyle="--", linewidth=1)

    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "radius_valley_overview.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_toi_priority_map(ranked: pd.DataFrame) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    scatter = ax.scatter(
        ranked["pl_orbper"],
        ranked["pl_rade"],
        c=ranked["priority_score"],
        cmap="viridis",
        s=30 + 0.02 * ranked["pl_trandep"].clip(upper=4000),
        alpha=0.75,
        edgecolors="none",
    )
    for _, row in ranked.head(12).iterrows():
        ax.annotate(
            row["toi"],
            (row["pl_orbper"], row["pl_rade"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xscale("log")
    ax.set_xlabel("Orbital Period [days]")
    ax.set_ylabel("Planet Radius [R_earth]")
    ax.set_title("TOI Follow-up Priority Map")
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Priority Score")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "toi_priority_map.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def to_native(value: object) -> object:
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    return value


def write_summary(
    confirmed: pd.DataFrame,
    toi_ranked: pd.DataFrame,
    flux_table: pd.DataFrame,
    metallicity_table: pd.DataFrame,
    fit: LinearFit,
    metallicity_delta: dict[str, float],
) -> None:
    timestamp = datetime.now().isoformat(timespec="seconds")
    summary = {
        "generated_at": timestamp,
        "confirmed_sample_size": int(len(confirmed)),
        "candidate_sample_size": int(len(toi_ranked)),
        "flux_fit": {
            "intercept": fit.intercept,
            "slope": fit.slope,
            "intercept_q16": fit.intercept_q16,
            "intercept_q84": fit.intercept_q84,
            "slope_q16": fit.slope_q16,
            "slope_q84": fit.slope_q84,
        },
        "flux_binned_valleys": [
            {key: to_native(value) for key, value in row.items()}
            for row in flux_table.to_dict(orient="records")
        ],
        "metallicity_binned_valleys": [
            {key: to_native(value) for key, value in row.items()}
            for row in metallicity_table.to_dict(orient="records")
        ],
        "metallicity_delta_high_minus_low": metallicity_delta,
        "science_case_counts": {
            key: int(value) for key, value in toi_ranked["science_case"].value_counts().to_dict().items()
        },
        "top_candidates": [
            {key: to_native(value) for key, value in row.items()}
            for row in toi_ranked[
                [
                    "toi",
                    "toipfx",
                    "st_tmag",
                    "pl_rade",
                    "pl_orbper",
                    "pl_insol",
                    "pl_eqt",
                    "pl_trandep",
                    "priority_score",
                    "science_case",
                ]
            ]
            .head(10)
            .to_dict(orient="records")
        ],
    }
    (RESULTS_DIR / "research_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze NASA Exoplanet Archive snapshots for radius-valley science and TOI prioritization."
    )
    parser.add_argument(
        "--ps-file",
        type=Path,
        default=RAW_DIR / "pscomppars_small_transit.csv",
        help="Path to the confirmed-planet snapshot CSV.",
    )
    parser.add_argument(
        "--toi-file",
        type=Path,
        default=RAW_DIR / "toi_pc_small.csv",
        help="Path to the TOI snapshot CSV.",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=120,
        help="Number of bootstrap resamples for uncertainty estimates.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for bootstrap and ranking reproducibility.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dirs()

    confirmed_raw = pd.read_csv(args.ps_file)
    toi_raw = pd.read_csv(args.toi_file)

    confirmed = clean_confirmed(confirmed_raw)
    toi = clean_toi(toi_raw)

    confirmed.to_csv(PROCESSED_DIR / "confirmed_planets_clean.csv", index=False)
    toi.to_csv(PROCESSED_DIR / "toi_candidates_clean.csv", index=False)

    flux_edges = np.quantile(confirmed["log_insol"], [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    flux_table = estimate_binned_valleys(
        confirmed,
        feature_column="log_insol",
        edges=flux_edges,
        label_prefix="flux",
    )
    flux_table.to_csv(RESULTS_DIR / "flux_binned_valleys.csv", index=False)

    metallicity_confirmed = confirmed.dropna(subset=["st_met"]).copy()
    metallicity_edges = np.quantile(metallicity_confirmed["st_met"], [0.0, 1 / 3, 2 / 3, 1.0])
    metallicity_table = estimate_binned_valleys(
        metallicity_confirmed,
        feature_column="st_met",
        edges=metallicity_edges,
        label_prefix="metallicity",
        min_points=120,
    )
    metallicity_table.to_csv(RESULTS_DIR / "metallicity_binned_valleys.csv", index=False)

    flux_fit = bootstrap_flux_fit(
        confirmed,
        flux_edges,
        n_bootstrap=args.bootstrap,
        seed=args.seed,
    )
    metallicity_delta = bootstrap_metallicity_delta(
        confirmed,
        metallicity_edges,
        n_bootstrap=args.bootstrap,
        seed=args.seed + 1,
    )

    toi_ranked = rank_toi_candidates(confirmed, toi, flux_fit)
    save_candidate_views(toi_ranked)
    plot_valley_overview(confirmed, flux_table, metallicity_table, flux_fit)
    plot_toi_priority_map(toi_ranked)
    write_summary(
        confirmed,
        toi_ranked,
        flux_table,
        metallicity_table,
        flux_fit,
        metallicity_delta,
    )

    print(f"Confirmed sample: {len(confirmed)} planets")
    print(f"TOI sample: {len(toi_ranked)} candidates")
    print(
        "Flux relation: log10(R_valley/R_earth) = "
        f"{flux_fit.intercept:.3f} + {flux_fit.slope:.3f} * log10(S/S_earth)"
    )
    print(
        "Metallicity delta (high - low tertile): "
        f"{metallicity_delta['base_delta']:.3f} R_earth"
    )
    print("Top 5 candidates:")
    for _, row in toi_ranked.head(5).iterrows():
        print(
            f"  {row['toi']}: score={row['priority_score']:.1f}, "
            f"R={row['pl_rade']:.2f} R_earth, P={row['pl_orbper']:.2f} d, "
            f"case={row['science_case']}"
        )


if __name__ == "__main__":
    main()
