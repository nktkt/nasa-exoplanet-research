from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CATALOG = ROOT / "data" / "raw" / "intbsc_black_hole_candidates.psv"


def parse_angle(value: str | float | int, *, is_ra: bool) -> float:
    text = str(value).strip()
    if not text:
        raise ValueError("Empty coordinate value")

    if re.fullmatch(r"[+-]?\d+(?:\.\d+)?", text):
        return float(text)

    parts = [part for part in re.split(r"[:\s]+", text) if part]
    if len(parts) != 3:
        raise ValueError(f"Unsupported coordinate format: {value}")

    sign = -1.0 if parts[0].startswith("-") else 1.0
    first = abs(float(parts[0]))
    second = float(parts[1])
    third = float(parts[2])
    total = first + second / 60.0 + third / 3600.0
    if is_ra:
        return total * 15.0
    return sign * total


def angular_separation_deg(
    ra1_deg: float,
    dec1_deg: float,
    ra2_deg: pd.Series,
    dec2_deg: pd.Series,
) -> pd.Series:
    ra1 = math.radians(ra1_deg)
    dec1 = math.radians(dec1_deg)
    ra2 = ra2_deg.map(math.radians)
    dec2 = dec2_deg.map(math.radians)
    cos_sep = (
        math.sin(dec1) * dec2.map(math.sin)
        + math.cos(dec1) * dec2.map(math.cos) * (ra2 - ra1).map(math.cos)
    )
    cos_sep = cos_sep.clip(-1.0, 1.0)
    return cos_sep.map(math.acos).map(math.degrees)


def load_catalog(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Catalog not found: {path}. Run scripts/fetch_nasa_black_hole_data.sh first."
        )

    df = pd.read_csv(path, sep="|")
    for column in df.columns:
        if pd.api.types.is_object_dtype(df[column]):
            df[column] = df[column].astype(str).str.strip()

    df["ra_deg"] = df["ra"].map(lambda value: parse_angle(value, is_ra=True))
    df["dec_deg"] = df["dec"].map(lambda value: parse_angle(value, is_ra=False))
    df["lii"] = pd.to_numeric(df["lii"], errors="coerce")
    df["bii"] = pd.to_numeric(df["bii"], errors="coerce")
    df["class"] = pd.to_numeric(df["class"], errors="coerce").astype("Int64")
    return df


def filter_by_name(df: pd.DataFrame, query: str) -> pd.DataFrame:
    needle = query.casefold()
    mask = (
        df["name"].str.casefold().str.contains(needle, na=False)
        | df["source_type"].str.casefold().str.contains(needle, na=False)
    )
    return df.loc[mask].copy()


def filter_by_source_type(df: pd.DataFrame, query: str) -> pd.DataFrame:
    needle = query.casefold()
    mask = df["source_type"].str.casefold().str.contains(needle, na=False)
    return df.loc[mask].copy()


def cone_search(
    df: pd.DataFrame,
    *,
    ra: str | float,
    dec: str | float,
    radius_arcmin: float,
) -> pd.DataFrame:
    ra_deg = parse_angle(ra, is_ra=True)
    dec_deg = parse_angle(dec, is_ra=False)
    separated = df.copy()
    separated["offset_deg"] = angular_separation_deg(
        ra_deg,
        dec_deg,
        separated["ra_deg"],
        separated["dec_deg"],
    )
    separated["offset_arcmin"] = separated["offset_deg"] * 60.0
    matched = separated.loc[separated["offset_arcmin"] <= radius_arcmin].copy()
    return matched.sort_values(["offset_arcmin", "name"]).reset_index(drop=True)


def run_search(
    df: pd.DataFrame,
    *,
    name: str | None,
    source_type: str | None,
    ra: str | float | None,
    dec: str | float | None,
    radius_arcmin: float,
    limit: int | None,
) -> pd.DataFrame:
    result = df.copy()
    if name:
        result = filter_by_name(result, name)
    if source_type:
        result = filter_by_source_type(result, source_type)
    if (ra is None) ^ (dec is None):
        raise ValueError("ra and dec must be provided together")
    if ra is not None and dec is not None:
        result = cone_search(result, ra=ra, dec=dec, radius_arcmin=radius_arcmin)
    else:
        result = result.sort_values("name").reset_index(drop=True)
    if limit is not None:
        result = result.head(limit).copy()
    return result


def format_table(df: pd.DataFrame) -> str:
    display_columns = [
        "name",
        "source_type",
        "ra",
        "dec",
        "lii",
        "bii",
    ]
    if "offset_arcmin" in df.columns:
        display_columns.insert(4, "offset_arcmin")
    present = [column for column in display_columns if column in df.columns]
    return df[present].to_string(index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search a NASA/HEASARC black-hole candidate snapshot derived from the "
            "INTEGRAL IBIS soft gamma-ray catalog."
        )
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=DEFAULT_CATALOG,
        help="Path to the local PSV catalog snapshot.",
    )
    parser.add_argument(
        "--name",
        help="Case-insensitive substring match against source name or source type.",
    )
    parser.add_argument(
        "--source-type",
        help="Case-insensitive source-type filter, e.g. BHC or Microquasar.",
    )
    parser.add_argument(
        "--ra",
        help="Right ascension in decimal degrees or sexagesimal hours.",
    )
    parser.add_argument(
        "--dec",
        help="Declination in decimal degrees or sexagesimal degrees.",
    )
    parser.add_argument(
        "--radius-arcmin",
        type=float,
        default=60.0,
        help="Cone-search radius in arcminutes when ra/dec are given.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of rows to print.",
    )
    parser.add_argument(
        "--format",
        choices=["table", "json", "csv"],
        default="table",
        help="Output format.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    catalog = load_catalog(args.catalog)
    result = run_search(
        catalog,
        name=args.name,
        source_type=args.source_type,
        ra=args.ra,
        dec=args.dec,
        radius_arcmin=args.radius_arcmin,
        limit=args.limit,
    )

    if args.format == "json":
        print(result.to_json(orient="records", indent=2))
        return
    if args.format == "csv":
        print(result.to_csv(index=False))
        return

    print(f"Catalog rows loaded: {len(catalog)}")
    print(f"Matches returned: {len(result)}")
    if result.empty:
        print("No black-hole candidates matched the search.")
        return
    print(format_table(result))


if __name__ == "__main__":
    main()
