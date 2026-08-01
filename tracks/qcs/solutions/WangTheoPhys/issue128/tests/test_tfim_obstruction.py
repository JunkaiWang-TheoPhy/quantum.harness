from __future__ import annotations

import copy
import hashlib
import json
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from scripts.certify_tfim_obstruction import (
    DEFAULT_OUTPUT,
    EXPECTED_CLAIM,
    FAMILY_DOMAIN,
    OPERATOR_LOWER_BOUND_FORMULA,
    build_payload,
    canonical_bytes,
    load_payload,
    verify_payload,
)
from trottercert.cubic_field import Cubic
from trottercert.pf4_bch_mapping import pf4_suzuki_gamma
from trottercert.tfim_obstruction import (
    CUBIC_RESULT_KEYS,
    RATIONAL_RESULT_KEYS,
    tfim_trace_moments,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_tfim_obstruction.py"
EXPECTED_KEYS = set(RATIONAL_RESULT_KEYS) | set(CUBIC_RESULT_KEYS)


def _reseal(payload: dict[str, object]) -> None:
    unsigned = {
        key: value for key, value in payload.items() if key != "payload_sha256"
    }
    payload["payload_sha256"] = hashlib.sha256(canonical_bytes(unsigned)).hexdigest()


def _without_pythonpath() -> dict[str, str]:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment["PYTHONHASHSEED"] = "0"
    return environment


def _positive_cubic(value: Cubic) -> bool:
    return value.a0 > 0 and value.a1 > 0 and value.a2 > 0


@pytest.mark.parametrize("length", (4, 6, 8))
def test_tfim_even_periodic_trace_formulas(length: int) -> None:
    h = Fraction(2)
    j = Fraction(3)
    moments = tfim_trace_moments(length, h, j, periodic=True)

    assert set(moments) == EXPECTED_KEYS
    assert all(
        isinstance(moments[key], Fraction) for key in RATIONAL_RESULT_KEYS
    )
    assert all(isinstance(moments[key], Cubic) for key in CUBIC_RESULT_KEYS)
    assert moments["trace_c2_over_d"] == 128 * length * h**4 * j**2
    assert moments["trace_d2_over_d"] == 128 * length * h**2 * j**4
    assert moments["trace_cd_over_d"] == 0
    assert moments["trace_h2_over_d"] == length * (h**2 + j**2)


def test_quadratic_core_pairing_and_endpoint_bound_are_exact() -> None:
    length = 4
    h = Fraction(2)
    j = Fraction(3)
    moments = tfim_trace_moments(length, h, j, periodic=True)
    core = (
        moments["trace_c2_over_d"]
        - 4 * moments["trace_cd_over_d"]
        + Fraction(8, 3) * moments["trace_d2_over_d"]
    )
    norm_upper = length * (abs(h) + abs(j))

    assert moments["quadratic_core_over_d"] == core
    assert moments["trace_h_e5_over_d"] == pf4_suzuki_gamma() * core
    assert moments["hamiltonian_norm_upper"] == norm_upper
    assert moments["endpoint_conjugation_lower_bound"] == (
        pf4_suzuki_gamma() * core / norm_upper
    )


@pytest.mark.parametrize(
    ("h", "j", "strict"),
    (
        (Fraction(2), Fraction(3), True),
        (Fraction(-2), Fraction(3), True),
        (Fraction(2), Fraction(-3), True),
        (Fraction(0), Fraction(3), False),
        (Fraction(2), Fraction(0), False),
        (Fraction(0), Fraction(0), False),
        (Fraction(2, 3), Fraction(-5, 7), True),
    ),
)
def test_strict_obstruction_exactly_matches_noncommuting_domain(
    h: Fraction,
    j: Fraction,
    strict: bool,
) -> None:
    moments = tfim_trace_moments(4, h, j, periodic=True)
    pairing = moments["trace_h_e5_over_d"]
    lower_bound = moments["endpoint_conjugation_lower_bound"]

    assert isinstance(pairing, Cubic)
    assert isinstance(lower_bound, Cubic)
    assert _positive_cubic(pairing) is strict
    assert _positive_cubic(lower_bound) is strict
    if not strict:
        assert pairing == Cubic.zero()
        assert lower_bound == Cubic.zero()


def test_zero_hamiltonian_has_zero_exact_quantities() -> None:
    moments = tfim_trace_moments(4, Fraction(0), Fraction(0), periodic=True)

    assert all(moments[key] == 0 for key in RATIONAL_RESULT_KEYS)
    assert all(moments[key] == Cubic.zero() for key in CUBIC_RESULT_KEYS)


def test_family_domain_and_piecewise_bound_are_certificate_bound() -> None:
    payload = build_payload()

    assert payload["schema_version"] == 3
    assert payload["formulas"]["operator_lower_bound"] == (
        OPERATOR_LOWER_BOUND_FORMULA
    )
    assert payload["family"]["domain"] == FAMILY_DOMAIN


def test_implementation_source_map_is_exact_and_complete() -> None:
    sources = build_payload()["implementation_sources"]
    expected_paths = {
        "pyproject.toml",
        "requirements-reproducibility.txt",
        "scripts/__init__.py",
        "scripts/certify_tfim_obstruction.py",
        "scripts/reference_tfim_gauge_witness.py",
        "src/trottercert/__init__.py",
        "src/trottercert/algebra.py",
        "src/trottercert/cubic_field.py",
        "src/trottercert/intervals.py",
        "src/trottercert/pf4_bch_mapping.py",
        "src/trottercert/tfim_obstruction.py",
        "src/trottercert/trace_obstruction.py",
    }

    assert set(sources) == expected_paths
    assert sources == {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in expected_paths
    }


def test_gauge_witness_binding_is_exact_and_scope_limited() -> None:
    gauge_path = (
        ROOT / "docs/experiments/processor-obstruction/tfim-gauge-witness.json"
    )
    gauge_payload = json.loads(gauge_path.read_bytes())
    binding = build_payload()["gauge_witness_binding"]

    assert binding == {
        "path": "docs/experiments/processor-obstruction/tfim-gauge-witness.json",
        "kind": "tfim_pf4_centered_polynomial_gauge_witness",
        "schema_version": 1,
        "file_sha256": hashlib.sha256(gauge_path.read_bytes()).hexdigest(),
        "payload_sha256": gauge_payload["payload_sha256"],
        "qualitative_status": "certified_on_checked_instances",
        "scope": "periodic TFIM h=j=1 at L=4,6,8,10",
    }


@pytest.mark.parametrize(
    ("length", "h", "j", "message"),
    [
        (True, Fraction(1), Fraction(1), "length"),
        (3, Fraction(1), Fraction(1), "even"),
        (4, 1, Fraction(1), "Fraction"),
        (4, Fraction(1), 1, "Fraction"),
    ],
)
def test_family_domain_is_fail_closed(length, h, j, message: str) -> None:
    with pytest.raises((TypeError, ValueError), match=message):
        tfim_trace_moments(length, h, j, periodic=True)


def test_payload_promotes_only_the_scoped_family_theorem() -> None:
    payload = build_payload()
    claim = payload["claim"]

    assert claim == EXPECTED_CLAIM
    assert claim["pf4_bch_mapping_status"] == (
        "proved_exact_suzuki_pf4_free_trace_identity"
    )
    assert claim["endpoint_conjugation_status"] == (
        "certified_leading_order_fixed_time_fixed_normalization"
    )
    assert claim["operator_lower_bound_status"] == (
        "certified_exact_algebraic_leading_coefficient"
    )
    assert claim["promotion_status"] == "promoted_scoped_family_theorem"
    assert claim["affine_gauge_status"] == (
        "qualitative_certified_periodic_tfim_h_eq_j_eq_1_checked_"
        "L_4_6_8_10_only"
    )
    assert claim["affine_gauge_universal_induction_status"] == "not_claimed"
    assert claim["affine_gauge_quantitative_distance_status"] == "not_claimed"
    assert claim["finite_step_no_go"] == "not_claimed"
    assert claim["total_time_eigenphase_status"] == "not_claimed"
    verify_payload(payload)


def test_frozen_payload_is_canonical_and_verifies() -> None:
    raw = DEFAULT_OUTPUT.read_bytes()
    payload = load_payload(DEFAULT_OUTPUT)

    assert raw == canonical_bytes(payload)
    verify_payload(payload)


def test_cli_build_and_verify(tmp_path: Path) -> None:
    output = tmp_path / "tfim-family.json"
    clean_environment = _without_pythonpath()
    built = subprocess.run(
        [sys.executable, str(SCRIPT), "--build", str(output)],
        cwd=ROOT,
        env=clean_environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert built.returncode == 0, built.stderr
    verified = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify", str(output)],
        cwd=ROOT,
        env=clean_environment,
        check=False,
        capture_output=True,
        text=True,
    )

    assert verified.returncode == 0, verified.stderr
    assert output.read_bytes() == canonical_bytes(load_payload(output))
    assert output.read_bytes() == DEFAULT_OUTPUT.read_bytes()


@pytest.mark.parametrize("serialization", ("pretty", "duplicate"))
def test_cli_verify_rejects_noncanonical_json(
    tmp_path: Path,
    serialization: str,
) -> None:
    payload = build_payload()
    if serialization == "pretty":
        raw = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
        expected_error = "canonical"
    else:
        canonical = canonical_bytes(payload)
        raw = b'{"schema_version":2,' + canonical[1:]
        expected_error = "duplicate JSON key"
    path = tmp_path / f"{serialization}.json"
    path.write_bytes(raw)
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify", str(path)],
        cwd=ROOT,
        env=_without_pythonpath(),
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode != 0
    assert expected_error in completed.stderr


def test_mutated_schema_version_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["schema_version"] = 1
    with pytest.raises(ValueError, match="schema version mismatch"):
        verify_payload(forged)


def test_old_quadratic_coefficient_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["exact_coefficients"]["pf4_quadratic_core"]["trace_cd"] = [14, 3]
    with pytest.raises(ValueError, match="PF4 quadratic coefficients mismatch"):
        verify_payload(forged)


def test_mutated_gamma_coordinate_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["suzuki_coefficient_field"]["gamma_coordinates"][2] = [1, 30]
    with pytest.raises(ValueError, match="gamma coordinate mismatch"):
        verify_payload(forged)


def test_mutated_trace_pairing_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["checked_instances"]["4"]["trace_h_e5_over_d"][0][0] += 1
    with pytest.raises(ValueError, match="trace_h_e5_over_d mismatch"):
        verify_payload(forged)


def test_mutated_suzuki_alpha2_coefficient_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["suzuki_coefficient_field"]["u_coordinates"][2] = [1, 30]
    with pytest.raises(ValueError, match=r"Suzuki alpha\^2 coefficient mismatch"):
        verify_payload(forged)


def test_mutated_checked_length_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["checked_lengths"][0] = 5
    with pytest.raises(ValueError, match="checked lengths mismatch"):
        verify_payload(forged)


def test_mutated_source_hash_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    source = next(iter(forged["implementation_sources"]))
    forged["implementation_sources"][source] = "0" * 64
    with pytest.raises(ValueError, match="source hash mismatch"):
        verify_payload(forged)


@pytest.mark.parametrize(
    ("section", "field"),
    (("formulas", "operator_lower_bound"), ("family", "domain")),
)
def test_mutated_piecewise_rule_is_rejected(
    section: str,
    field: str,
) -> None:
    forged = copy.deepcopy(build_payload())
    forged[section][field] = "unconditional division"
    with pytest.raises(ValueError, match="operator lower-bound piecewise rule mismatch"):
        verify_payload(forged)


@pytest.mark.parametrize(
    "field",
    (
        "affine_gauge_universal_induction_status",
        "affine_gauge_quantitative_distance_status",
        "finite_step_no_go",
        "total_time_eigenphase_status",
    ),
)
def test_unproved_claim_cannot_be_promoted_by_rehashing(field: str) -> None:
    forged = copy.deepcopy(build_payload())
    forged["claim"][field] = "proved"
    forged["payload_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="claim gate mismatch"):
        verify_payload(forged)


@pytest.mark.parametrize("field", ("file_sha256", "payload_sha256", "scope"))
def test_resealed_gauge_binding_mutation_is_rejected(field: str) -> None:
    forged = copy.deepcopy(build_payload())
    forged["gauge_witness_binding"][field] = (
        "0" * 64 if field.endswith("sha256") else "all even lengths"
    )
    _reseal(forged)
    with pytest.raises(ValueError, match="gauge witness artifact binding mismatch"):
        verify_payload(forged)


def test_unknown_field_and_boolean_numeric_alias_are_rejected() -> None:
    unknown = copy.deepcopy(build_payload())
    unknown["unexpected"] = "resealed"
    _reseal(unknown)
    with pytest.raises(ValueError, match="field set mismatch"):
        verify_payload(unknown)

    boolean_alias = copy.deepcopy(build_payload())
    boolean_alias["checked_parameters"]["h"][0] = True
    _reseal(boolean_alias)
    with pytest.raises(TypeError, match="numeric alias"):
        verify_payload(boolean_alias)
