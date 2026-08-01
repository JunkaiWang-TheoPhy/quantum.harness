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
    build_payload,
    canonical_bytes,
    load_payload,
    verify_payload,
)
from trottercert.tfim_obstruction import tfim_trace_moments

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_tfim_obstruction.py"
EXPECTED_KEYS = {
    "trace_c2_over_d",
    "trace_d2_over_d",
    "trace_cd_over_d",
    "trace_h2_over_d",
    "obstruction_over_d",
    "operator_lower_bound",
}


@pytest.mark.parametrize("length", (4, 6, 8))
def test_tfim_even_periodic_trace_formulas(length: int) -> None:
    h = Fraction(2)
    j = Fraction(3)
    moments = tfim_trace_moments(length, h, j, periodic=True)

    assert set(moments) == EXPECTED_KEYS
    assert all(isinstance(value, Fraction) for value in moments.values())
    assert moments["trace_c2_over_d"] == 128 * length * h**4 * j**2
    assert moments["trace_d2_over_d"] == 128 * length * h**2 * j**4
    assert moments["trace_cd_over_d"] == 0
    assert moments["trace_h2_over_d"] == length * (h**2 + j**2)


def test_quadratic_core_and_conditional_bound_are_exact() -> None:
    moments = tfim_trace_moments(4, Fraction(2), Fraction(3), periodic=True)
    expected_core = (
        Fraction(1, 2) * moments["trace_c2_over_d"]
        + Fraction(14, 3) * moments["trace_cd_over_d"]
        + Fraction(4, 3) * moments["trace_d2_over_d"]
    )
    assert moments["obstruction_over_d"] == expected_core
    assert moments["operator_lower_bound"] == expected_core / (4 * (2 + 3))


def test_zero_hamiltonian_has_zero_conditional_quantities() -> None:
    moments = tfim_trace_moments(4, Fraction(0), Fraction(0), periodic=True)
    assert all(value == 0 for value in moments.values())


def test_zero_hamiltonian_piecewise_domain_is_certificate_bound() -> None:
    payload = build_payload()
    assert payload["formulas"]["operator_lower_bound"] == (
        "0 if h=j=0; otherwise "
        "abs(obstruction_over_d)/(length*(abs(h)+abs(j)))"
    )
    assert payload["family"]["domain"] == (
        "even length >= 4; exact rational h,j; "
        "operator_lower_bound=0 at h=j=0"
    )


def test_implementation_source_map_is_exact_and_complete() -> None:
    sources = build_payload()["implementation_sources"]
    expected_paths = {
        "scripts/certify_tfim_obstruction.py",
        "src/trottercert/algebra.py",
        "src/trottercert/cubic_field.py",
        "src/trottercert/intervals.py",
        "src/trottercert/tfim_obstruction.py",
        "src/trottercert/trace_obstruction.py",
    }
    assert set(sources) == expected_paths
    assert sources == {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in expected_paths
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


def test_payload_inherits_unverified_pf4_gate() -> None:
    payload = build_payload()
    claim = payload["claim"]

    assert claim["trace_moment_status"] == "certified_exact_pauli_counting"
    assert (
        claim["pf4_bch_mapping_status"]
        == "algebraic_form_only_unverified_bch_mapping"
    )
    assert claim["obstruction_status"] == "conditional_not_certified"
    assert claim["operator_lower_bound_status"] == "conditional_not_certified"
    assert claim["finite_step_no_go"] == "not_claimed"
    assert claim["promotion_status"] == "blocked_pending_pf4_bch_mapping"
    verify_payload(payload)


def test_frozen_payload_is_canonical_and_verifies() -> None:
    raw = DEFAULT_OUTPUT.read_bytes()
    payload = load_payload(DEFAULT_OUTPUT)
    assert raw == canonical_bytes(payload)
    verify_payload(payload)


def test_cli_build_and_verify(tmp_path: Path) -> None:
    output = tmp_path / "tfim-family.json"
    built = subprocess.run(
        [sys.executable, str(SCRIPT), "--build", str(output)],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        check=False,
        capture_output=True,
        text=True,
    )
    assert built.returncode == 0, built.stderr
    verified = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify", str(output)],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        check=False,
        capture_output=True,
        text=True,
    )
    assert verified.returncode == 0, verified.stderr
    assert output.read_bytes() == canonical_bytes(load_payload(output))


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
        raw = b'{"schema_version":1,' + canonical[1:]
        expected_error = "duplicate JSON key"
    path = tmp_path / f"{serialization}.json"
    path.write_bytes(raw)
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify", str(path)],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0
    assert expected_error in completed.stderr


def test_mutated_trace_cd_is_rejected_distinctly() -> None:
    forged = copy.deepcopy(build_payload())
    forged["checked_instances"]["4"]["trace_cd_over_d"] = [1, 1]
    with pytest.raises(ValueError, match="trace_cd_over_d mismatch"):
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
def test_mutated_zero_hamiltonian_piecewise_rule_is_rejected(
    section: str,
    field: str,
) -> None:
    forged = copy.deepcopy(build_payload())
    forged[section][field] = "unconditional division"
    with pytest.raises(ValueError, match="operator lower-bound piecewise rule mismatch"):
        verify_payload(forged)


def test_gate_cannot_be_promoted_by_rehashing() -> None:
    forged = copy.deepcopy(build_payload())
    forged["claim"]["pf4_bch_mapping_status"] = "proved"
    forged["claim"]["obstruction_status"] = "proved"
    forged["claim"]["operator_lower_bound_status"] = "proved"
    forged["claim"]["finite_step_no_go"] = "proved"
    forged["claim"]["promotion_status"] = "promoted"
    forged["payload_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="claim gate mismatch"):
        verify_payload(forged)
