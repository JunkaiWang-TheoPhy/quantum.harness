#!/usr/bin/env python3
"""Build and verify the Issue-128 finite-step signed-margin certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

from scripts.certify_dual_e7_pairing import verify_reduced_payload as verify_e7
from scripts.certify_dual_e9_pairing import verify_reduced_payload as verify_e9
from scripts.certify_dual_log_tail import verify_payload as verify_tail
from scripts.certify_extensive_commutant_witness import verify_extensive_payload
from trottercert.cubic_field import Cubic
from trottercert.dual_word_manifest import (
    load_manifest_index,
    load_word_manifest,
    verify_manifest_index,
)
from trottercert.finite_step_obstruction import (
    decide_affine_spectral_obstruction,
    decide_finite_step,
)
from trottercert.intervals import RationalInterval, cube_root_four_interval


ISSUE_ROOT = Path(__file__).resolve().parents[1]
PROCESSOR_ROOT = ISSUE_ROOT / "docs/experiments/processor-obstruction"
DEFAULT_E5 = PROCESSOR_ROOT / "extensive-commutant-witness.json"
DEFAULT_E7 = PROCESSOR_ROOT / "dual-e7-pairing.json"
DEFAULT_E9 = PROCESSOR_ROOT / "dual-e9-pairing.json"
DEFAULT_TAIL = PROCESSOR_ROOT / "dual-log-tail.json"
DEFAULT_OUTPUT = PROCESSOR_ROOT / "finite-step-obstruction.json"
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/finite_step_obstruction.py",
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/intervals.py",
    Path(__file__).resolve(),
)
ROOT_DECIMAL_DIGITS = 48


def _canonical_bytes(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _digest(payload: object) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path, field: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_bytes())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read {field}: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{field} must contain a JSON object")
    return payload


def _source_label(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ISSUE_ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def _rational_json(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _parse_rational(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(not isinstance(entry, int) or isinstance(entry, bool) for entry in value)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if _rational_json(result) != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [
        _rational_json(value.a0),
        _rational_json(value.a1),
        _rational_json(value.a2),
    ]


def _parse_cubic(value: object, field: str) -> Cubic:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must contain three cubic coordinates")
    return Cubic(*(_parse_rational(entry, field) for entry in value))


def _interval_json(value: RationalInterval) -> dict[str, list[int]]:
    return {
        "lower": _rational_json(value.lower),
        "upper": _rational_json(value.upper),
    }


def _verify_reduced_e9_digest(
    payload: Mapping[str, object],
    *,
    e5_sha256: str,
    e7_sha256: str,
    q5: Cubic,
    q7: Cubic,
) -> Cubic:
    if (
        payload.get("schema_version") != 1
        or payload.get("kind") != "issue128_dual_e9_pairing"
        or payload.get("degree") != 9
        or payload.get("length") != 12
        or payload.get("shard_count") != 64
    ):
        raise ValueError("reduced E9 configuration mismatch")
    claim = payload.get("claim")
    if not isinstance(claim, Mapping) or claim != {
        "full_e9_operator": "not_computed",
        "dual_e9_pairing": "exact",
        "finite_step_status": "inconclusive",
        "missing": "E11-and-higher dual tail",
    }:
        raise ValueError("reduced E9 claim mismatch")
    inputs = payload.get("inputs")
    if not isinstance(inputs, Mapping):
        raise ValueError("reduced E9 inputs are missing")
    if inputs.get("extensive_witness_sha256") != e5_sha256:
        raise ValueError("reduced E9 extensive-witness digest mismatch")
    if inputs.get("dual_e7_sha256") != e7_sha256:
        raise ValueError("reduced E9 E7 digest mismatch")
    if _parse_cubic(inputs.get("q5_per_cell"), "E9 q5 input") != q5:
        raise ValueError("reduced E9 q5 input mismatch")
    if _parse_cubic(inputs.get("q7_per_cell"), "E9 q7 input") != q7:
        raise ValueError("reduced E9 q7 input mismatch")
    pairings = payload.get("pairings")
    if not isinstance(pairings, Mapping):
        raise ValueError("reduced E9 pairings are missing")
    q9 = _parse_cubic(pairings.get("tau_w_per_cell"), "E9 q9")
    if _parse_cubic(pairings.get("q9_over_97_pow_8"), "scaled E9 q9") != q9 / 97**8:
        raise ValueError("reduced E9 scaled q9 mismatch")
    expected = q5 / 97**4 + q7 / 97**6 + q9 / 97**8
    if (
        _parse_cubic(
            pairings.get("exact_e5_e7_e9_per_cell_at_r97"),
            "combined E5/E7/E9 pairing",
        )
        != expected
    ):
        raise ValueError("reduced E9 combined pairing mismatch")
    mathematical_digest = payload.get("mathematical_payload_sha256")
    if not isinstance(mathematical_digest, str) or mathematical_digest != _digest(
        {key: value for key, value in payload.items() if key != "mathematical_payload_sha256"}
    ):
        raise ValueError("reduced E9 mathematical payload digest mismatch")
    return q9


def _full_verify_e9(
    payload: Mapping[str, object],
    index_path: Path,
    *,
    e5_path: Path,
    e7_path: Path,
) -> None:
    index = load_manifest_index(index_path)
    verify_manifest_index(index, index_path.parent)
    manifests = tuple(
        load_word_manifest(index_path.parent / record.path, index)
        for record in index.shards
    )
    verify_e9(
        payload,
        index=index,
        manifests=manifests,
        extensive_path=e5_path,
        dual_e7_path=e7_path,
    )


def build_payload(
    *,
    e5_path: Path = DEFAULT_E5,
    e7_path: Path = DEFAULT_E7,
    e9_path: Path = DEFAULT_E9,
    tail_path: Path = DEFAULT_TAIL,
    index_path: Path | None,
    require_full_e9: bool = True,
) -> dict[str, Any]:
    if require_full_e9 and index_path is None:
        raise ValueError("a complete E9 manifest index is required to build the certificate")
    e5 = _load(e5_path, "E5 artifact")
    e7 = _load(e7_path, "E7 artifact")
    e9 = _load(e9_path, "E9 artifact")
    tail = _load(tail_path, "dual-log tail artifact")
    verify_extensive_payload(e5)
    verify_e7(e7)
    verify_tail(tail)
    if index_path is not None:
        _full_verify_e9(
            e9,
            index_path,
            e5_path=e5_path,
            e7_path=e7_path,
        )

    stable = e5.get("stable_relations")
    e7_pairings = e7.get("pairings")
    if not isinstance(stable, Mapping) or not isinstance(e7_pairings, Mapping):
        raise ValueError("exact E5/E7 pairing records are missing")
    q5 = _parse_cubic(stable.get("tau_w_e5_per_cell"), "E5 q5")
    q7 = _parse_cubic(e7_pairings.get("tau_w_per_cell"), "E7 q7")
    e5_sha256 = _sha(e5_path)
    e7_sha256 = _sha(e7_path)
    q9 = _verify_reduced_e9_digest(
        e9,
        e5_sha256=e5_sha256,
        e7_sha256=e7_sha256,
        q5=q5,
        q7=q7,
    )
    instance = tail.get("instance")
    bounds = tail.get("bounds")
    moments = tail.get("moments")
    envelope = tail.get("envelope")
    if (
        not isinstance(instance, Mapping)
        or not isinstance(bounds, Mapping)
        or not isinstance(moments, Mapping)
        or not isinstance(envelope, Mapping)
    ):
        raise ValueError("dual-log tail instance, moments, envelope, or bounds are missing")
    steps = instance.get("steps")
    if not isinstance(steps, int) or isinstance(steps, bool) or steps != 97:
        raise ValueError("finite-step certificate requires exactly 97 steps")
    cells = instance.get("cells")
    if not isinstance(cells, int) or isinstance(cells, bool) or cells != 36:
        raise ValueError("finite-step certificate requires exactly 36 cells")
    h_hs_squared = _parse_rational(moments.get("h_hs_squared"), "H HS squared")
    h_hs_cap = _parse_rational(moments.get("h_hs_cap"), "H HS cap")
    if h_hs_squared != 54 or h_hs_cap != Fraction(15, 2):
        raise ValueError("finite-step spectral moment geometry mismatch")
    one_step_log_defect = _parse_rational(
        envelope.get("log_defect"), "one-step log defect"
    )
    tail_bound = _parse_rational(bounds.get("total_log_tail"), "total log tail")
    root = cube_root_four_interval(ROOT_DECIMAL_DIGITS)
    decision = decide_finite_step(
        q5,
        q7,
        q9,
        root,
        RationalInterval.point(tail_bound),
        Fraction(1, steps),
    )
    dual_margin_lower = (
        decision.signed_margin_interval.lower
        if decision.signed_margin_interval.lower > 0
        else Fraction()
    )
    effective_log_defect = steps * one_step_log_defect
    affine = decide_affine_spectral_obstruction(
        dual_margin_lower=dual_margin_lower,
        cells=cells,
        m2=h_hs_squared,
        m3=Fraction(-27),
        h_hs_cap=h_hs_cap,
        effective_log_defect=effective_log_defect,
    )
    payload: dict[str, Any] = {
        "schema_version": 1,
        "kind": "issue128_finite_step_obstruction",
        "formula": "five_copy_fourth_order_suzuki_four_matchings",
        "instance": dict(instance),
        "normalization": "dual local-log pairing per 2x2 cell",
        "series": "h^4 q5 + h^6 q7 + h^8 q9 + R_{>=11}^{dual}",
        "step_size": _rational_json(decision.h),
        "root_interval_decimal_digits": ROOT_DECIMAL_DIGITS,
        "root_interval": _interval_json(decision.root_interval),
        "pairings": {
            "q5_per_cell": _cubic_json(q5),
            "q7_per_cell": _cubic_json(q7),
            "q9_per_cell": _cubic_json(q9),
        },
        "intervals": {
            "leading": _interval_json(decision.leading_interval),
            "absolute_leading": _interval_json(decision.absolute_leading_interval),
            "tail_bound": _interval_json(decision.tail_interval),
            "signed_margin": _interval_json(decision.signed_margin_interval),
        },
        "affine_spectral_invariant": {
            "definition": "m2(H)^3*m3(A_centered)^2-m3(H)^2*m2(A_centered)^3",
            "orbit": "A=a*I+b*U*H*U^dagger",
            "m2_h": _rational_json(h_hs_squared),
            "m3_h": _rational_json(Fraction(-27)),
            "h_hs_cap": _rational_json(h_hs_cap),
            "one_step_log_defect_bound": _rational_json(one_step_log_defect),
            "effective_log_defect_bound": _rational_json(effective_log_defect),
            "centered_defect_cap": _rational_json(affine.centered_defect_cap),
            "linear_invariant_lower": _rational_json(
                affine.linear_invariant_lower
            ),
            "nonlinear_remainder_upper": _rational_json(
                affine.nonlinear_remainder_upper
            ),
            "invariant_margin": _rational_json(affine.invariant_margin),
            "status": affine.status,
        },
        "sources": {
            "e5": {"path": _source_label(e5_path), "sha256": e5_sha256},
            "e7": {"path": _source_label(e7_path), "sha256": e7_sha256},
            "e9": {
                "path": _source_label(e9_path),
                "sha256": _sha(e9_path),
                "manifest_index_sha256": e9.get("manifest_index_sha256"),
                "mathematical_payload_sha256": e9.get("mathematical_payload_sha256"),
            },
            "tail": {"path": _source_label(tail_path), "sha256": _sha(tail_path)},
        },
        "implementation_sources": {
            _source_label(path): _sha(path) for path in SOURCE_PATHS
        },
        "claim": {
            "local_log_status": (
                "certified_local_log_obstruction"
                if decision.status == "certified_obstruction"
                else decision.status
            ),
            "finite_step_effective_spectrum_status": affine.status,
            "promotion_rule": (
                "signed_margin.lower > 0 and affine_invariant_margin > 0"
            ),
            "promoted": (
                affine.status
                == "certified_affine_effective_spectral_obstruction"
            ),
        },
    }
    payload["mathematical_payload_sha256"] = _digest(
        {key: value for key, value in payload.items() if key != "mathematical_payload_sha256"}
    )
    return payload


def verify_payload(
    payload: Mapping[str, object],
    *,
    e5_path: Path = DEFAULT_E5,
    e7_path: Path = DEFAULT_E7,
    e9_path: Path = DEFAULT_E9,
    tail_path: Path = DEFAULT_TAIL,
    index_path: Path | None = None,
) -> None:
    expected = build_payload(
        e5_path=e5_path,
        e7_path=e7_path,
        e9_path=e9_path,
        tail_path=tail_path,
        index_path=index_path,
        require_full_e9=False,
    )
    if dict(payload) != expected:
        raise ValueError("finite-step obstruction payload regeneration mismatch")


def _write_payload(path: Path, payload: Mapping[str, object], *, replace: bool) -> None:
    if path.exists() and not replace:
        raise FileExistsError(f"refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.pending.{os.getpid()}")
    if temporary.exists():
        raise FileExistsError(f"temporary artifact already exists: {temporary}")
    try:
        temporary.write_bytes(_canonical_bytes(payload))
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--build", type=Path)
    action.add_argument("--verify", type=Path)
    parser.add_argument("--e5", type=Path, default=DEFAULT_E5)
    parser.add_argument("--e7", type=Path, default=DEFAULT_E7)
    parser.add_argument("--e9", type=Path, default=DEFAULT_E9)
    parser.add_argument("--tail", type=Path, default=DEFAULT_TAIL)
    parser.add_argument("--index", type=Path)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    if args.build is not None:
        payload = build_payload(
            e5_path=args.e5,
            e7_path=args.e7,
            e9_path=args.e9,
            tail_path=args.tail,
            index_path=args.index,
        )
        _write_payload(args.build, payload, replace=args.replace)
        print(f"wrote finite-step obstruction certificate: {args.build}", flush=True)
        return
    if args.replace:
        parser.error("--replace cannot be combined with --verify")
    payload = _load(args.verify, "finite-step obstruction certificate")
    verify_payload(
        payload,
        e5_path=args.e5,
        e7_path=args.e7,
        e9_path=args.e9,
        tail_path=args.tail,
        index_path=args.index,
    )
    print(f"finite-step obstruction certificate valid: {args.verify}", flush=True)


if __name__ == "__main__":
    main()
