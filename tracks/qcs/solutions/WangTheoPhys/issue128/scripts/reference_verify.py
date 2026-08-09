#!/usr/bin/env python3
"""Standalone downstream verifier for Issue 128 schema-v3 certificates.

This program intentionally imports only the Python standard library.  It
independently checks D4/D5 sidecar coverage and anticommutation, rational norm
enclosures, the complete finite-step ledger including the generator tail,
adjacent-step rejection, and resource arithmetic.  It does not regenerate the
large D4/D5 coefficient maps from the Suzuki formula; that remains the primary
deep verifier's separate algebraic obligation.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


EXPECTED_NORMALIZATION = "(XX+YY+ZZ)/4"
Interval = tuple[Fraction, Fraction]


def _integer(
    value: object,
    field: str,
    *,
    minimum: int | None = None,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{field} must be at least {minimum}")
    return value


def _pair(value: object, field: str) -> Fraction:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{field} must be a rational pair")
    numerator = _integer(value[0], f"{field} numerator")
    denominator = _integer(value[1], f"{field} denominator", minimum=1)
    return Fraction(numerator, denominator)


def _canonical_json(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def _sidecar_payload(
    certificate_path: Path,
    metadata: dict[str, object],
    *,
    compressed: bool,
    label: str,
) -> dict[str, object]:
    root = certificate_path.resolve().parent
    sidecar = (root / str(metadata["path"])).resolve()
    if sidecar.parent != root:
        raise ValueError(f"{label} sidecar path escapes certificate directory")
    raw = sidecar.read_bytes()
    if hashlib.sha256(raw).hexdigest() != str(metadata["sha256"]):
        raise ValueError(f"{label} sidecar digest mismatch")
    try:
        decoded = gzip.decompress(raw) if compressed else raw
        payload: Any = json.loads(decoded)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} sidecar is not canonical JSON") from error
    if not isinstance(payload, dict) or decoded != _canonical_json(payload):
        raise ValueError(f"{label} sidecar JSON is not canonical")
    return payload


def _anticommutes(
    left: tuple[int, int],
    right: tuple[int, int],
) -> bool:
    return (
        (left[0] & right[1]).bit_count()
        + (left[1] & right[0]).bit_count()
    ) & 1 == 1


def _verify_partition(
    payload: dict[str, object],
    *,
    require_same_support: bool,
    label: str,
) -> dict[str, object]:
    if payload.get("schema_version") != 1:
        raise ValueError(f"unsupported {label} sidecar schema")
    coefficient_denominator = _integer(
        payload.get("coefficient_denominator"),
        f"{label} coefficient denominator",
        minimum=1,
    )
    sqrt_denominator = _integer(
        payload.get("sqrt_denominator"),
        f"{label} square-root denominator",
        minimum=1,
    )
    terms = payload.get("terms")
    groups = payload.get("groups")
    if not isinstance(terms, list) or not isinstance(groups, list):
        raise ValueError(f"{label} terms and groups must be lists")

    paulis: list[tuple[int, int]] = []
    magnitudes: list[int] = []
    for position, term in enumerate(terms):
        if not isinstance(term, list) or len(term) != 4:
            raise ValueError(f"{label} term {position} is malformed")
        x_mask = _integer(term[0], f"{label} x mask", minimum=0)
        z_mask = _integer(term[1], f"{label} z mask", minimum=0)
        lower = _integer(term[2], f"{label} lower endpoint")
        upper = _integer(term[3], f"{label} upper endpoint")
        if lower > upper:
            raise ValueError(f"{label} coefficient interval is reversed")
        paulis.append((x_mask, z_mask))
        magnitudes.append(max(abs(lower), abs(upper)))
    if paulis != sorted(paulis) or len(paulis) != len(set(paulis)):
        raise ValueError(f"{label} Pauli terms are not unique and canonical")

    coverage = bytearray(len(paulis))
    total_bound_numerator = 0
    maximum_group_size = 0
    for group_position, raw_group in enumerate(groups):
        if not isinstance(raw_group, list) or len(raw_group) != 2:
            raise ValueError(f"{label} group {group_position} is malformed")
        raw_indices = raw_group[0]
        if not isinstance(raw_indices, list) or not raw_indices:
            raise ValueError(f"{label} group indices must be a nonempty list")
        indices = [
            _integer(index, f"{label} group index", minimum=0)
            for index in raw_indices
        ]
        if any(index >= len(paulis) for index in indices):
            raise ValueError(f"{label} group index is out of range")
        for index in indices:
            if coverage[index]:
                raise ValueError(f"{label} partition contains a duplicate term")
            coverage[index] = 1
        group_paulis = [paulis[index] for index in indices]
        if require_same_support and len(
            {x_mask | z_mask for x_mask, z_mask in group_paulis}
        ) != 1:
            raise ValueError(f"{label} group mixes distinct supports")
        for left_position, left in enumerate(group_paulis):
            for right in group_paulis[left_position + 1 :]:
                if not _anticommutes(left, right):
                    raise ValueError(
                        f"{label} group contains a commuting Pauli pair"
                    )

        bound_numerator = _integer(
            raw_group[1],
            f"{label} group bound numerator",
            minimum=0,
        )
        squared_weight_numerator = sum(
            magnitudes[index] ** 2 for index in indices
        )
        if (
            bound_numerator**2 * coefficient_denominator**2
            < squared_weight_numerator * sqrt_denominator**2
        ):
            raise ValueError(f"{label} group square-root bound is too small")
        total_bound_numerator += bound_numerator
        maximum_group_size = max(maximum_group_size, len(indices))

    if not all(coverage):
        raise ValueError(f"{label} partition omits coefficient terms")
    cell_bound = _pair(payload.get("cell_bound"), f"{label} cell bound")
    if Fraction(total_bound_numerator, sqrt_denominator) != cell_bound:
        raise ValueError(f"{label} cell bound does not equal the group sum")
    if "term_count" in payload and _integer(
        payload["term_count"],
        f"{label} term count",
    ) != len(terms):
        raise ValueError(f"{label} term count mismatch")
    if "group_count" in payload and _integer(
        payload["group_count"],
        f"{label} group count",
    ) != len(groups):
        raise ValueError(f"{label} group count mismatch")
    if "site_bound" in payload and _pair(
        payload["site_bound"],
        f"{label} site bound",
    ) != cell_bound / 4:
        raise ValueError(f"{label} site bound mismatch")
    return {
        "site_bound": cell_bound / 4,
        "term_count": len(terms),
        "group_count": len(groups),
        "maximum_group_size": maximum_group_size,
    }


def _verify_partition_metadata(
    metadata: dict[str, object],
    partition: dict[str, object],
    *,
    label: str,
) -> None:
    fields = (
        ("term_count", "term_count", "term count"),
        ("group_count", "group_count", "group count"),
        ("max_group_size", "maximum_group_size", "maximum group size"),
    )
    for metadata_key, partition_key, description in fields:
        declared = _integer(
            metadata.get(metadata_key),
            f"{label} {description}",
            minimum=0,
        )
        if declared != partition[partition_key]:
            raise ValueError(f"{label} {description} metadata mismatch")


def _interval_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def _interval_neg(value: Interval) -> Interval:
    return -value[1], -value[0]


def _interval_sub(left: Interval, right: Interval) -> Interval:
    return _interval_add(left, _interval_neg(right))


def _interval_mul(left: Interval, right: Interval) -> Interval:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return min(products), max(products)


def _interval_reciprocal(value: Interval) -> Interval:
    if value[0] <= 0 <= value[1]:
        raise ValueError("interval reciprocal contains zero")
    return 1 / value[1], 1 / value[0]


def _interval_scale(value: Interval, scalar: Fraction) -> Interval:
    return _interval_mul(value, (scalar, scalar))


def _outward_quantize(value: Interval, denominator: int) -> Interval:
    lower_scaled = value[0] * denominator
    upper_scaled = value[1] * denominator
    lower = lower_scaled.numerator // lower_scaled.denominator
    upper = -((-upper_scaled.numerator) // upper_scaled.denominator)
    return Fraction(lower, denominator), Fraction(upper, denominator)


def _cube_root_four_interval(decimal_digits: int) -> Interval:
    denominator = 10**decimal_digits
    target = 4 * denominator**3
    low = denominator
    high = 4 * denominator
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**3 <= target:
            low = middle
        else:
            high = middle
    return Fraction(low, denominator), Fraction(high, denominator)


def _suzuki_stages(decimal_digits: int) -> tuple[tuple[int, Interval], ...]:
    root = _cube_root_four_interval(decimal_digits)
    denominator = 10**decimal_digits
    divisor = _interval_sub(
        (Fraction(4), Fraction(4)),
        root,
    )
    u = _outward_quantize(
        _interval_reciprocal(divisor),
        denominator,
    )
    scales = (
        u,
        u,
        _interval_sub(
            (Fraction(1), Fraction(1)),
            _interval_scale(u, Fraction(4)),
        ),
        u,
        u,
    )
    raw: list[tuple[int, Interval]] = []
    for scale in scales:
        half = _interval_scale(scale, Fraction(1, 2))
        raw.extend((index, half) for index in range(3))
        raw.append((3, scale))
        raw.extend((index, half) for index in reversed(range(3)))
    merged: list[tuple[int, Interval]] = []
    for fragment, coefficient in raw:
        if merged and merged[-1][0] == fragment:
            previous_fragment, previous_coefficient = merged.pop()
            merged.append(
                (
                    previous_fragment,
                    _interval_add(previous_coefficient, coefficient),
                )
            )
        else:
            merged.append((fragment, coefficient))
    return tuple(merged)


def _tail_site_bound(
    stages: tuple[tuple[int, Interval], ...],
    steps: int,
) -> Fraction:
    prefix = Fraction()
    total = Fraction()
    for _, coefficient in stages:
        coefficient_upper = max(abs(coefficient[0]), abs(coefficient[1]))
        ratio = prefix / steps
        if ratio >= 1:
            raise ValueError("step count is outside the tail convergence region")
        total += coefficient_upper * ratio**8 / (1 - ratio)
        prefix += coefficient_upper
    return Fraction(3, 8) * total


def _ceil_nth_root(value: Fraction, degree: int) -> int:
    low, high = 0, 1
    while high**degree * value.denominator < value.numerator:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**degree * value.denominator >= value.numerator:
            high = middle
        else:
            low = middle
    return high


def verify(path: str | Path) -> dict[str, object]:
    certificate_path = Path(path)
    data = json.loads(certificate_path.read_text())
    if not isinstance(data, dict) or data.get("schema_version") != 3:
        raise ValueError("reference verifier supports schema v3 only")
    benchmark = data["benchmark"]
    if benchmark["normalization"] != EXPECTED_NORMALIZATION:
        raise ValueError("Hamiltonian normalization mismatch")
    length = _integer(benchmark["length"], "benchmark length", minimum=1)
    if length % 2:
        raise ValueError("benchmark length must be even")
    n_sites = length * length
    time = _pair(benchmark["time"], "benchmark time")
    tolerance = _pair(benchmark["tolerance"], "benchmark tolerance")
    if time != 1:
        raise ValueError("reference verifier pins T=1")

    published = data["published_baseline"]
    published_density = _pair(
        published["site_density_upper"],
        "published site density",
    )
    published_steps = _ceil_nth_root(
        published_density * n_sites / tolerance,
        4,
    )
    if published_steps != _integer(
        published["steps"],
        "published steps",
        minimum=1,
    ):
        raise ValueError("published step count mismatch")
    published_groups = 30 * published_steps + 1
    if published_groups != _integer(
        published["group_exponentials"],
        "published groups",
        minimum=1,
    ):
        raise ValueError("published group count mismatch")

    candidate = data["candidate"]
    d4_payload = _sidecar_payload(
        certificate_path,
        candidate["d4_certificate"],
        compressed=False,
        label="D4",
    )
    d4 = _verify_partition(
        d4_payload,
        require_same_support=False,
        label="D4",
    )
    d4_metadata = candidate["d4_certificate"]
    _verify_partition_metadata(d4_metadata, d4, label="D4")
    if _pair(
        d4_metadata["cell_norm_upper"],
        "D4 metadata cell bound",
    ) != 4 * d4["site_bound"]:
        raise ValueError("D4 metadata bound mismatch")

    d5: dict[str, object] | None = None
    if "d5_certificate" in candidate:
        d5_payload = _sidecar_payload(
            certificate_path,
            candidate["d5_certificate"],
            compressed=True,
            label="D5",
        )
        d5 = _verify_partition(
            d5_payload,
            require_same_support=True,
            label="D5",
        )
        _verify_partition_metadata(
            candidate["d5_certificate"],
            d5,
            label="D5",
        )
        if _pair(
            candidate["d5_certificate"]["site_norm_upper"],
            "D5 metadata site bound",
        ) != d5["site_bound"]:
            raise ValueError("D5 metadata bound mismatch")
    if "d6_certificate" in candidate:
        raise ValueError(
            "standalone reference verifier does not yet support D6 sidecars"
        )

    steps = _integer(candidate["steps"], "candidate steps", minimum=2)
    e5_site = _pair(candidate["e5_site_l1_upper"], "E5 site bound")
    e7_site = _pair(candidate["e7_site_majorant"], "E7 site majorant")
    h_e5 = 24 * e5_site
    hh_e5 = 28 * h_e5
    hhh_e5 = 32 * hh_e5
    h_e7 = 32 * e7_site
    d5_site = d5["site_bound"] if d5 is not None else 2 * h_e5
    d6_site = 7 * e7_site + Fraction(2, 3) * hh_e5
    d7_site = 3 * h_e7 + Fraction(1, 6) * hhh_e5
    stages = _suzuki_stages(
        _integer(
            candidate["coefficient_interval_decimal_digits"],
            "candidate coefficient precision",
            minimum=1,
        )
    )

    def ledger(at_steps: int) -> dict[str, Fraction]:
        contributions = {
            "degree4": Fraction(n_sites) * d4["site_bound"]
            / (5 * at_steps**4),
            "degree5": Fraction(n_sites) * d5_site / (6 * at_steps**5),
            "degree6": Fraction(n_sites) * d6_site / (7 * at_steps**6),
            "degree7": Fraction(n_sites) * d7_site / (8 * at_steps**7),
            "tail": Fraction(n_sites) * _tail_site_bound(stages, at_steps),
        }
        contributions["total"] = sum(contributions.values(), Fraction())
        return contributions

    accepted = ledger(steps)
    previous = ledger(steps - 1)
    submitted = candidate["contributions"]
    for name in ("degree4", "degree5", "degree6", "degree7", "tail"):
        if _pair(submitted[name], f"candidate {name}") != accepted[name]:
            raise ValueError(f"candidate {name} contribution mismatch")
    if _pair(
        candidate["global_error_upper"],
        "candidate global error",
    ) != accepted["total"]:
        raise ValueError("candidate global error mismatch")
    if _pair(
        candidate["previous_step_error_upper"],
        "candidate previous-step error",
    ) != previous["total"]:
        raise ValueError("candidate previous-step error mismatch")
    if accepted["total"] > tolerance or previous["total"] <= tolerance:
        raise ValueError("candidate adjacent-step tolerance boundary is invalid")

    candidate_groups = 30 * steps + 1
    baseline_bonds = published_groups * n_sites // 2
    candidate_bonds = candidate_groups * n_sites // 2
    expected_resources = {
        "published_steps": published_steps,
        "candidate_steps": steps,
        "published_group_exponentials": published_groups,
        "candidate_group_exponentials": candidate_groups,
        "published_bond_propagators": baseline_bonds,
        "candidate_bond_propagators": candidate_bonds,
        "published_cnot_upper": 3 * baseline_bonds,
        "candidate_cnot_upper": 3 * candidate_bonds,
    }
    if data["claimed_resources"] != expected_resources:
        raise ValueError("resource summary mismatch")
    ratio = Fraction(published_groups, candidate_groups)
    if _pair(
        data["claims"]["exact_improvement_ratio"],
        "exact improvement ratio",
    ) != ratio:
        raise ValueError("exact improvement ratio mismatch")

    return {
        "valid": True,
        "verification_scope": (
            "independent_sidecars_finite_step_ledger_and_resources"
        ),
        "coefficient_generation_replayed": False,
        "published_steps": published_steps,
        "candidate_steps": steps,
        "baseline_centers_scanned": 0,
        "candidate_error_upper": str(accepted["total"]),
        "previous_step_error_upper": str(previous["total"]),
        "d4_term_count": d4["term_count"],
        "d4_group_count": d4["group_count"],
        "d5_term_count": d5["term_count"] if d5 is not None else 0,
        "d5_group_count": d5["group_count"] if d5 is not None else 0,
        "exact_improvement_ratio": str(ratio),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate")
    arguments = parser.parse_args()
    print(json.dumps(verify(arguments.certificate), indent=2))


if __name__ == "__main__":
    main()
