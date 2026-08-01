"""Build or verify the exact scoped PF4--TFIM family certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

from trottercert.cubic_field import Cubic
from trottercert.pf4_bch_mapping import (
    pf4_suzuki_gamma,
    verify_pf4_suzuki_trace_identity,
)
from trottercert.tfim_obstruction import (
    CUBIC_RESULT_KEYS,
    RATIONAL_RESULT_KEYS,
    tfim_trace_moments,
)
from trottercert.trace_obstruction import (
    pf4_trace_identity_record,
    verify_identity_record,
)

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/tfim-family-obstruction.json"
)
CHECKED_LENGTHS = (4, 6, 8)
CHECK_H = Fraction(2)
CHECK_J = Fraction(3)
SUZUKI_U = Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/algebra.py",
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/cubic_local.py",
    ISSUE_ROOT / "src/trottercert/intervals.py",
    ISSUE_ROOT / "src/trottercert/pf4_bch_mapping.py",
    ISSUE_ROOT / "src/trottercert/tfim_obstruction.py",
    ISSUE_ROOT / "src/trottercert/trace_obstruction.py",
    Path(__file__).resolve(),
)
FAMILY_DOMAIN = (
    "even length >= 4; exact rational h,j; strict endpoint-conjugation "
    "obstruction iff h*j != 0; commuting exception otherwise"
)
OPERATOR_LOWER_BOUND_FORMULA = (
    "0 if h=j=0; otherwise gamma*abs(quadratic_core_over_d)/"
    "(length*(abs(h)+abs(j)))"
)
EXPECTED_CLAIM = {
    "trace_moment_status": "certified_exact_pauli_counting",
    "pf4_bch_mapping_status": "proved_exact_suzuki_pf4_free_trace_identity",
    "endpoint_conjugation_status": (
        "certified_leading_order_fixed_time_fixed_normalization"
    ),
    "operator_lower_bound_status": (
        "certified_exact_algebraic_leading_coefficient"
    ),
    "affine_gauge_status": "not_claimed",
    "finite_step_no_go": "not_claimed",
    "total_time_eigenphase_status": "not_claimed",
    "promotion_status": "promoted_scoped_family_theorem",
}
EXPECTED_FORMULAS = {
    "trace_c2_over_d": "128*length*h^4*j^2",
    "trace_d2_over_d": "128*length*h^2*j^4",
    "trace_cd_over_d": "0",
    "trace_h2_over_d": "length*(h^2+j^2)",
    "quadratic_core_over_d": (
        "trace_c2_over_d - 4*trace_cd_over_d + 8*trace_d2_over_d/3"
    ),
    "trace_h_e5_over_d": "gamma*quadratic_core_over_d",
    "operator_lower_bound": OPERATOR_LOWER_BOUND_FORMULA,
}
EXPECTED_QUADRATIC_COEFFICIENTS = {
    "trace_c2": [1, 1],
    "trace_cd": [-4, 1],
    "trace_d2": [8, 3],
}


def _pair(value: Fraction) -> list[int]:
    exact = Fraction(value)
    return [exact.numerator, exact.denominator]


def _cubic(value: Cubic) -> list[list[int]]:
    if not isinstance(value, Cubic):
        raise TypeError("cubic result must be an exact Cubic")
    return [_pair(value.a0), _pair(value.a1), _pair(value.a2)]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_label(path: Path) -> str:
    return path.resolve().relative_to(ISSUE_ROOT.resolve()).as_posix()


def canonical_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def _payload_digest(payload: Mapping[str, object]) -> str:
    unsigned = {
        key: value for key, value in payload.items() if key != "payload_sha256"
    }
    return hashlib.sha256(canonical_bytes(unsigned)).hexdigest()


def _result_json(
    result: Mapping[str, Fraction | Cubic],
) -> dict[str, list[int] | list[list[int]]]:
    expected_keys = RATIONAL_RESULT_KEYS + CUBIC_RESULT_KEYS
    if tuple(result) != expected_keys:
        raise ValueError("TFIM exact-result schema mismatch")
    encoded: dict[str, list[int] | list[list[int]]] = {}
    for key in RATIONAL_RESULT_KEYS:
        value = result[key]
        if not isinstance(value, Fraction):
            raise TypeError(f"{key} must be an exact Fraction")
        encoded[key] = _pair(value)
    for key in CUBIC_RESULT_KEYS:
        value = result[key]
        if not isinstance(value, Cubic):
            raise TypeError(f"{key} must be an exact Cubic")
        encoded[key] = _cubic(value)
    return encoded


def _suzuki_field_json() -> dict[str, object]:
    alpha = Cubic(0, 1, 0)
    if (Cubic(4, 0, 0) - alpha) * SUZUKI_U != Cubic.one():
        raise ArithmeticError("exact Suzuki coefficient identity failed")
    gamma = pf4_suzuki_gamma()
    if not all(coordinate > 0 for coordinate in (gamma.a0, gamma.a1, gamma.a2)):
        raise ArithmeticError("exact PF4 gamma positivity failed")
    return {
        "basis": ["1", "alpha", "alpha^2"],
        "minimal_polynomial": "alpha^3 - 4",
        "positive_real_embedding": "alpha > 0",
        "u_identity": "u = 1/(4-alpha)",
        "u_coordinates": _cubic(SUZUKI_U),
        "gamma_coordinates": _cubic(gamma),
        "gamma_status": "strictly_positive_by_positive_coordinates",
    }


def build_payload() -> dict[str, Any]:
    verify_pf4_suzuki_trace_identity()
    pf4_record = pf4_trace_identity_record()
    if not verify_identity_record(pf4_record, "pf4"):
        raise ValueError("canonical PF4 identity record is invalid")
    if pf4_record.identity_status != EXPECTED_CLAIM["pf4_bch_mapping_status"]:
        raise ValueError("PF4 identity gate changed unexpectedly")

    checked = {
        str(length): _result_json(
            tfim_trace_moments(length, CHECK_H, CHECK_J, periodic=True)
        )
        for length in CHECKED_LENGTHS
    }
    payload: dict[str, Any] = {
        "schema_version": 2,
        "kind": "tfim_exact_pf4_endpoint_conjugation_family",
        "family": {
            "hamiltonian": "H=A+B, A=h sum_i X_i, B=j sum_i Z_i Z_{i+1}",
            "boundary": "periodic",
            "domain": FAMILY_DOMAIN,
            "nested_commutators": {
                "C": "[A,[A,B]]",
                "D": "[B,[B,A]]",
            },
            "theorem_scope": (
                "leading order; pure endpoint conjugation; fixed time and "
                "fixed Hamiltonian normalization"
            ),
        },
        "checked_lengths": list(CHECKED_LENGTHS),
        "checked_parameters": {"h": _pair(CHECK_H), "j": _pair(CHECK_J)},
        "formulas": dict(EXPECTED_FORMULAS),
        "exact_coefficients": {
            "trace_c2_over_d": {
                "multiplicity": 128,
                "length_power": 1,
                "h_power": 4,
                "j_power": 2,
            },
            "trace_d2_over_d": {
                "multiplicity": 128,
                "length_power": 1,
                "h_power": 2,
                "j_power": 4,
            },
            "pf4_quadratic_core": dict(EXPECTED_QUADRATIC_COEFFICIENTS),
            "closed_family_pairing": (
                "128*length*gamma*h^2*j^2*(h^2 + 8*j^2/3)"
            ),
        },
        "suzuki_coefficient_field": _suzuki_field_json(),
        "checked_instances": checked,
        "assumptions": [
            "finite-dimensional periodic even TFIM chain",
            "normalized trace is evaluated by exact Pauli orthogonality",
            (
                "closed coefficients follow from local Pauli multiplicities, "
                "not finite-size fitting"
            ),
            "the exact PF4 mapping is verified as a cyclic free-word identity",
            (
                "the endpoint-conjugation lower bound uses the Pauli-l1 upper "
                "bound on the normalized trace norm of H"
            ),
            (
                "affine rescaling, finite-step TFIM, and total-time eigenphase "
                "claims are outside this certificate"
            ),
        ],
        "pf4_identity": {
            "kind": pf4_record.identity_kind,
            "status": pf4_record.identity_status,
            "formula": pf4_record.formula,
            "identity_digest": pf4_record.identity_digest,
        },
        "implementation_sources": {
            _source_label(path): _sha(path) for path in SOURCE_PATHS
        },
        "claim": dict(EXPECTED_CLAIM),
    }
    payload["payload_sha256"] = _payload_digest(payload)
    return payload


def _mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{field} must be a mapping")
    return value


def verify_payload(payload: Mapping[str, object]) -> None:
    if not isinstance(payload, Mapping):
        raise TypeError("TFIM family payload must be a mapping")
    if payload.get("schema_version") != 2:
        raise ValueError("schema version mismatch")
    if payload.get("kind") != "tfim_exact_pf4_endpoint_conjugation_family":
        raise ValueError("artifact kind mismatch")
    if payload.get("checked_lengths") != list(CHECKED_LENGTHS):
        raise ValueError("checked lengths mismatch")
    if payload.get("checked_parameters") != {
        "h": _pair(CHECK_H),
        "j": _pair(CHECK_J),
    }:
        raise ValueError("checked parameters mismatch")

    family = _mapping(payload.get("family"), "family")
    formulas = _mapping(payload.get("formulas"), "formulas")
    if (
        family.get("domain") != FAMILY_DOMAIN
        or formulas.get("operator_lower_bound") != OPERATOR_LOWER_BOUND_FORMULA
    ):
        raise ValueError("operator lower-bound piecewise rule mismatch")
    if formulas != EXPECTED_FORMULAS:
        raise ValueError("closed family formulas mismatch")

    exact_coefficients = _mapping(
        payload.get("exact_coefficients"),
        "exact coefficients",
    )
    quadratic = _mapping(
        exact_coefficients.get("pf4_quadratic_core"),
        "PF4 quadratic coefficients",
    )
    if quadratic != EXPECTED_QUADRATIC_COEFFICIENTS:
        raise ValueError("PF4 quadratic coefficients mismatch")

    suzuki = _mapping(payload.get("suzuki_coefficient_field"), "Suzuki field")
    coordinates = suzuki.get("u_coordinates")
    if not isinstance(coordinates, list) or len(coordinates) != 3:
        raise ValueError("Suzuki coordinate schema mismatch")
    if coordinates[2] != _pair(SUZUKI_U.a2):
        raise ValueError("Suzuki alpha^2 coefficient mismatch")
    if suzuki.get("gamma_coordinates") != _cubic(pf4_suzuki_gamma()):
        raise ValueError("gamma coordinate mismatch")

    verify_pf4_suzuki_trace_identity()
    record = pf4_trace_identity_record()
    if payload.get("pf4_identity") != {
        "kind": record.identity_kind,
        "status": record.identity_status,
        "formula": record.formula,
        "identity_digest": record.identity_digest,
    }:
        raise ValueError("PF4 identity record mismatch")

    instances = _mapping(payload.get("checked_instances"), "checked instances")
    for length in CHECKED_LENGTHS:
        instance = _mapping(instances.get(str(length)), f"length {length} instance")
        expected = _result_json(
            tfim_trace_moments(length, CHECK_H, CHECK_J, periodic=True)
        )
        if instance.get("trace_cd_over_d") != expected["trace_cd_over_d"]:
            raise ValueError(f"trace_cd_over_d mismatch at length {length}")
        if instance.get("trace_h_e5_over_d") != expected["trace_h_e5_over_d"]:
            raise ValueError(f"trace_h_e5_over_d mismatch at length {length}")
        if instance != expected:
            raise ValueError(f"exact TFIM result mismatch at length {length}")

    if payload.get("claim") != EXPECTED_CLAIM:
        raise ValueError("claim gate mismatch")
    expected_sources = {_source_label(path): _sha(path) for path in SOURCE_PATHS}
    if payload.get("implementation_sources") != expected_sources:
        raise ValueError("source hash mismatch")
    digest = payload.get("payload_sha256")
    if not isinstance(digest, str) or digest != _payload_digest(payload):
        raise ValueError("payload digest mismatch")
    if dict(payload) != build_payload():
        raise ValueError("payload regeneration mismatch")


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_payload(path: Path) -> dict[str, Any]:
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read TFIM family artifact: {path}") from error
    try:
        payload = json.loads(raw, object_pairs_hook=_unique_object)
    except json.JSONDecodeError as error:
        raise ValueError(f"cannot parse TFIM family artifact: {path}") from error
    if not isinstance(payload, dict):
        raise TypeError("TFIM family artifact must contain a JSON object")
    if raw != canonical_bytes(payload):
        raise ValueError("TFIM family artifact is not canonical JSON")
    return payload


def _write_payload(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.pending.{os.getpid()}")
    try:
        temporary.write_bytes(canonical_bytes(build_payload()))
        verify_payload(load_payload(temporary))
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--build", type=Path)
    action.add_argument("--verify", type=Path)
    arguments = parser.parse_args()
    try:
        if arguments.build is not None:
            _write_payload(arguments.build)
            print(f"built={arguments.build}")
        else:
            verify_payload(load_payload(arguments.verify))
            print(f"verified={arguments.verify}")
    except (ArithmeticError, OSError, TypeError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
