from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
from typing import Any, Mapping

from .hpc_artifacts import (
    CoordinateCubicTerms,
    coordinate_decode_terms,
    coordinate_terms_to_json,
)
from .intervals import cube_root_four_interval


@dataclass(frozen=True)
class ExactSeriesVerification:
    order: int
    source_commit: str
    degree_term_counts: tuple[int, ...]
    cell_l1_upper: Fraction
    site_l1_upper: Fraction
    highest_degree: CoordinateCubicTerms


def _canonical_json_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def read_portable_canonical_gzip(path: str | Path) -> dict[str, Any]:
    """Decode canonical JSON while allowing platform-specific gzip streams.

    The proof artifact is separately bound by SHA-256.  Here canonicality means
    that the uncompressed payload is the unique sorted compact JSON encoding;
    this avoids treating zlib's platform/version-specific DEFLATE choices as
    mathematical content.
    """

    raw = Path(path).read_bytes()
    try:
        decoded = gzip.decompress(raw)
        payload: Any = json.loads(decoded)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("exact-series artifact is not gzip JSON") from error
    if not isinstance(payload, dict):
        raise ValueError("exact-series artifact root must be an object")
    if decoded != _canonical_json_bytes(payload):
        raise ValueError("exact-series JSON payload is not canonical")
    return payload


def _fraction_pair(value: object, *, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or isinstance(value[0], bool)
        or isinstance(value[1], bool)
        or not isinstance(value[0], int)
        or not isinstance(value[1], int)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a rational pair")
    return Fraction(value[0], value[1])


def verify_exact_monolithic_series(
    payload: Mapping[str, object],
    *,
    expected_order: int,
    expected_source_commit: str,
    coefficient_interval_decimal_digits: int = 30,
) -> ExactSeriesVerification:
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported exact-series schema")
    if payload.get("kind") != "issue128_exact_right_generator_monolithic":
        raise ValueError("unexpected exact-series artifact kind")
    if payload.get("formula_id") != "five_copy_suzuki_fourth_order_exact_cubic":
        raise ValueError("exact-series formula mismatch")
    if payload.get("stage_count") != 31:
        raise ValueError("exact-series stage count mismatch")
    if payload.get("order") != expected_order:
        raise ValueError("exact-series order mismatch")
    if payload.get("source_commit") != expected_source_commit:
        raise ValueError("exact-series source commit mismatch")
    raw_series = payload.get("series")
    if not isinstance(raw_series, list) or len(raw_series) != expected_order + 1:
        raise ValueError("exact-series degree coverage mismatch")

    degree_term_counts: list[int] = []
    highest_degree: CoordinateCubicTerms | None = None
    for degree, raw_degree in enumerate(raw_series):
        decoded = coordinate_decode_terms(raw_degree)
        if coordinate_terms_to_json(decoded) != raw_degree:
            raise ValueError("exact-series terms are not canonically sorted")
        degree_term_counts.append(len(decoded))
        if degree == expected_order:
            highest_degree = decoded

    if highest_degree is None:
        raise ValueError("exact-series highest degree is missing")

    root = cube_root_four_interval(coefficient_interval_decimal_digits)
    cell_l1 = sum(
        (coefficient.enclose(root).abs_upper() for coefficient in highest_degree.values()),
        Fraction(),
    )
    if cell_l1 != _fraction_pair(
        payload.get("cell_pauli_l1_upper"), field="cell Pauli-l1 upper"
    ):
        raise ValueError("exact-series cell Pauli-l1 mismatch")
    if cell_l1 / 4 != _fraction_pair(
        payload.get("site_pauli_l1_upper"), field="site Pauli-l1 upper"
    ):
        raise ValueError("exact-series site Pauli-l1 mismatch")
    return ExactSeriesVerification(
        order=expected_order,
        source_commit=expected_source_commit,
        degree_term_counts=tuple(degree_term_counts),
        cell_l1_upper=cell_l1,
        site_l1_upper=cell_l1 / 4,
        highest_degree=highest_degree,
    )
