from __future__ import annotations

import math

import numpy as np
import pandas as pd

from src.nasa_exoplanet_research import LinearFit
from src.nasa_exoplanet_research import clean_confirmed
from src.nasa_exoplanet_research import clean_toi
from src.nasa_exoplanet_research import percentile_rank
from src.nasa_exoplanet_research import rank_toi_candidates
from src.nasa_exoplanet_research import science_case
from src.nasa_exoplanet_research import valley_radius_from_series


def test_clean_confirmed_filters_out_of_range_rows() -> None:
    df = pd.DataFrame(
        [
            {
                "pl_name": "Planet A",
                "hostname": "Star A",
                "st_rad": 1.0,
                "st_teff": 5000,
                "pl_rade": 1.8,
                "pl_orbper": 10.0,
                "pl_insol": 100.0,
            },
            {
                "pl_name": "Planet B",
                "hostname": "Star B",
                "st_rad": 2.2,
                "st_teff": 5000,
                "pl_rade": 1.8,
                "pl_orbper": 10.0,
                "pl_insol": 100.0,
            },
            {
                "pl_name": "Planet C",
                "hostname": "Star C",
                "st_rad": 1.0,
                "st_teff": 2900,
                "pl_rade": 1.8,
                "pl_orbper": 10.0,
                "pl_insol": 100.0,
            },
        ]
    )

    cleaned = clean_confirmed(df)

    assert list(cleaned["pl_name"]) == ["Planet A"]
    assert {"log_period", "log_radius", "log_insol", "log_teff"}.issubset(cleaned.columns)


def test_clean_toi_filters_and_adds_log_columns() -> None:
    df = pd.DataFrame(
        [
            {
                "toi": 100.01,
                "st_rad": 1.1,
                "st_teff": 5100,
                "pl_rade": 2.0,
                "pl_orbper": 8.0,
                "pl_insol": 40.0,
            },
            {
                "toi": 101.01,
                "st_rad": 1.1,
                "st_teff": 5100,
                "pl_rade": 5.2,
                "pl_orbper": 8.0,
                "pl_insol": 40.0,
            },
        ]
    )

    cleaned = clean_toi(df)

    assert list(cleaned["toi"]) == [100.01]
    assert math.isclose(cleaned.iloc[0]["log_radius"], math.log10(2.0))


def test_percentile_rank_respects_direction() -> None:
    series = pd.Series([1.0, 2.0, 3.0])

    higher_is_better = percentile_rank(series, higher_better=True)
    lower_is_better = percentile_rank(series, higher_better=False)

    assert higher_is_better.iloc[-1] > higher_is_better.iloc[0]
    assert lower_is_better.iloc[0] > lower_is_better.iloc[-1]


def test_valley_radius_from_series_finds_midpoint_between_two_peaks() -> None:
    rng = np.random.default_rng(7)
    radius = pd.Series(
        np.concatenate(
            [
                rng.normal(1.35, 0.05, 200),
                rng.normal(2.55, 0.08, 200),
            ]
        )
    )

    valley = valley_radius_from_series(radius)

    assert valley is not None
    assert 1.6 < valley["valley_radius"] < 2.1
    assert valley["left_peak_radius"] < valley["valley_radius"] < valley["right_peak_radius"]


def test_science_case_labels_expected_regions() -> None:
    assert science_case(
        pd.Series({"valley_distance_dex": 0.03, "pl_rade": 2.0, "pl_orbper": 10.0, "pl_insol": 100.0})
    ) == "radius_valley_probe"
    assert science_case(
        pd.Series({"valley_distance_dex": 0.10, "pl_rade": 1.4, "pl_orbper": 0.8, "pl_insol": 500.0})
    ) == "ultra_short_period_rocky"
    assert science_case(
        pd.Series({"valley_distance_dex": 0.10, "pl_rade": 2.4, "pl_orbper": 12.0, "pl_insol": 80.0})
    ) == "temperate_sub_neptune"


def test_rank_toi_candidates_returns_sorted_scores() -> None:
    confirmed = pd.DataFrame(
        [
            {
                "log_period": 0.00,
                "log_radius": 0.18,
                "log_insol": 1.00,
                "log_teff": 3.68,
                "pl_rade": 1.50,
                "pl_insol": 10.0,
                "st_teff": 4800.0,
                "pl_orbper": 1.0,
                "st_tmag": 10.0,
                "pl_trandep": 300.0,
            },
            {
                "log_period": 0.40,
                "log_radius": 0.26,
                "log_insol": 1.40,
                "log_teff": 3.72,
                "pl_rade": 1.82,
                "pl_insol": 25.0,
                "st_teff": 5200.0,
                "pl_orbper": 2.5,
                "st_tmag": 10.5,
                "pl_trandep": 500.0,
            },
            {
                "log_period": 0.78,
                "log_radius": 0.40,
                "log_insol": 2.00,
                "log_teff": 3.76,
                "pl_rade": 2.50,
                "pl_insol": 100.0,
                "st_teff": 5750.0,
                "pl_orbper": 6.0,
                "st_tmag": 11.0,
                "pl_trandep": 700.0,
            },
        ]
    )
    toi = pd.DataFrame(
        [
            {
                "toi": 200.01,
                "toipfx": 200,
                "st_tmag": 8.8,
                "ra": 10.0,
                "dec": -20.0,
                "pl_rade": 2.00,
                "pl_orbper": 3.0,
                "pl_insol": 50.0,
                "pl_eqt": 900.0,
                "pl_trandep": 1200.0,
                "log_period": math.log10(3.0),
                "log_radius": math.log10(2.0),
                "log_insol": math.log10(50.0),
                "log_teff": math.log10(5100.0),
                "st_teff": 5100.0,
            },
            {
                "toi": 201.01,
                "toipfx": 201,
                "st_tmag": 11.2,
                "ra": 11.0,
                "dec": -21.0,
                "pl_rade": 3.20,
                "pl_orbper": 9.0,
                "pl_insol": 180.0,
                "pl_eqt": 1200.0,
                "pl_trandep": 450.0,
                "log_period": math.log10(9.0),
                "log_radius": math.log10(3.2),
                "log_insol": math.log10(180.0),
                "log_teff": math.log10(5600.0),
                "st_teff": 5600.0,
            },
        ]
    )

    fit = LinearFit(
        intercept=0.03,
        slope=0.10,
        intercept_q16=0.02,
        intercept_q84=0.04,
        slope_q16=0.08,
        slope_q84=0.12,
    )

    ranked = rank_toi_candidates(confirmed, toi, fit)

    assert list(ranked["priority_score"]) == sorted(ranked["priority_score"], reverse=True)
    assert ranked["predicted_valley_radius"].gt(0).all()
    assert ranked["science_case"].isin(
        {
            "radius_valley_probe",
            "ultra_short_period_rocky",
            "temperate_sub_neptune",
            "general_demographic_gain",
        }
    ).all()
