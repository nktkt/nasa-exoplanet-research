from __future__ import annotations

import math

import pandas as pd

from src.nasa_black_hole_search import cone_search
from src.nasa_black_hole_search import filter_by_name
from src.nasa_black_hole_search import parse_angle
from src.nasa_black_hole_search import run_search


def sample_catalog() -> pd.DataFrame:
    df = pd.DataFrame(
        [
            {
                "name": "IGR J17456-2901",
                "source_type": "Sgr A* ?",
                "ra": "17 45 43.0",
                "dec": "-28 59 28",
                "lii": 359.9641,
                "bii": -0.0467,
            },
            {
                "name": "SS 433",
                "source_type": "HMXB / Microquasar",
                "ra": "19 11 50.0",
                "dec": "04 58 58",
                "lii": 39.6950,
                "bii": -2.2462,
            },
            {
                "name": "SWIFT J1753.5-0127",
                "source_type": "LMXB, BHC",
                "ra": "17 53 28.0",
                "dec": "-01 27 09",
                "lii": 24.8964,
                "bii": 12.1867,
            },
        ]
    )
    df["ra_deg"] = df["ra"].map(lambda value: parse_angle(value, is_ra=True))
    df["dec_deg"] = df["dec"].map(lambda value: parse_angle(value, is_ra=False))
    return df


def test_parse_angle_supports_decimal_and_sexagesimal() -> None:
    assert math.isclose(parse_angle("17 45 43.0", is_ra=True), 266.4291666666667)
    assert math.isclose(parse_angle("-28 59 28", is_ra=False), -28.99111111111111)
    assert math.isclose(parse_angle("266.4", is_ra=True), 266.4)


def test_filter_by_name_matches_name_or_source_type() -> None:
    df = sample_catalog()

    sgr = filter_by_name(df, "sgr a")
    bhc = filter_by_name(df, "bhc")

    assert list(sgr["name"]) == ["IGR J17456-2901"]
    assert list(bhc["name"]) == ["SWIFT J1753.5-0127"]


def test_cone_search_returns_nearby_matches_with_offsets() -> None:
    df = sample_catalog()

    result = cone_search(
        df,
        ra="17 45 40",
        dec="-29 00 00",
        radius_arcmin=10.0,
    )

    assert list(result["name"]) == ["IGR J17456-2901"]
    assert result.iloc[0]["offset_arcmin"] < 10.0


def test_run_search_combines_filters_and_limit() -> None:
    df = sample_catalog()

    result = run_search(
        df,
        name=None,
        source_type="bhc",
        ra=None,
        dec=None,
        radius_arcmin=60.0,
        limit=1,
    )

    assert len(result) == 1
    assert result.iloc[0]["name"] == "SWIFT J1753.5-0127"
