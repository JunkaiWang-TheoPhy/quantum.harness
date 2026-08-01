from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import pytest

import scripts.run_heldout_model as heldout
from scripts.run_heldout_model import (
    FROZEN_LENGTHS,
    FROZEN_STEPS,
    SCIENTIFIC_DEPENDENCY_PATHS,
    run_heldout,
    verify_heldout_payload,
)
from trottercert.hamiltonian import tfim_terms

ROOT = Path(__file__).resolve().parents[1]
HEAD = heldout._git_output("rev-parse", "HEAD").decode().strip()


def test_tfim_terms_have_exact_signs_and_canonical_keys() -> None:
    terms = tfim_terms(4, Fraction(3), Fraction(2), False)
    assert list(terms) == sorted(terms)
    assert terms == {
        (0, 1): Fraction(-3),
        (0, 2): Fraction(-3),
        (0, 4): Fraction(-3),
        (0, 8): Fraction(-3),
        (3, 0): Fraction(-2),
        (6, 0): Fraction(-2),
        (12, 0): Fraction(-2),
    }


def test_tfim_periodic_closure_is_explicit() -> None:
    opened = tfim_terms(4, Fraction(1), Fraction(1), False)
    periodic = tfim_terms(4, Fraction(1), Fraction(1), True)
    assert (9, 0) not in opened
    assert periodic[(9, 0)] == -1
    assert len(periodic) == len(opened) + 1


@pytest.mark.parametrize(
    "arguments",
    [
        (4, 1, Fraction(1), False),
        (4, Fraction(1), 1, False),
        (4, Fraction(1), Fraction(1), 1),
    ],
)
def test_tfim_constructor_requires_strict_exact_types(arguments) -> None:
    with pytest.raises(TypeError):
        tfim_terms(*arguments)


def test_heldout_run_requires_freeze_commit() -> None:
    with pytest.raises(ValueError, match="rules_frozen_at_commit is required"):
        run_heldout(
            field=Fraction(1),
            coupling=Fraction(1),
            rules_frozen_at_commit="",
            lengths=(4,),
            steps=(8,),
        )


def test_freeze_fails_closed_when_rule_content_differs(monkeypatch) -> None:
    original = heldout._current_rule_bytes

    def mutated(relative: str) -> bytes:
        result = original(relative)
        if relative.endswith("paper-a-claim-matrix.json"):
            return result + b"\n"
        return result

    monkeypatch.setattr(heldout, "_current_rule_bytes", mutated)
    with pytest.raises(ValueError, match="frozen rule mismatch"):
        run_heldout(
            field=Fraction(1),
            coupling=Fraction(1),
            rules_frozen_at_commit=HEAD,
            lengths=(4,),
            steps=(8,),
        )


def test_later_rule_evolution_does_not_break_frozen_artifact_verify(
    monkeypatch,
) -> None:
    payload = run_heldout(
        field=Fraction(1),
        coupling=Fraction(1),
        rules_frozen_at_commit=HEAD,
        lengths=(4,),
        steps=(8,),
    )
    original = heldout._current_rule_bytes

    def evolved(relative: str) -> bytes:
        result = original(relative)
        if relative.endswith("paper-a-file-ownership.json"):
            return result + b"\n"
        return result

    monkeypatch.setattr(heldout, "_current_rule_bytes", evolved)
    assert verify_heldout_payload(payload)
    with pytest.raises(ValueError, match="frozen rule mismatch"):
        run_heldout(
            field=Fraction(1),
            coupling=Fraction(1),
            rules_frozen_at_commit=HEAD,
            lengths=(4,),
            steps=(8,),
        )


def test_small_heldout_run_is_byte_deterministic_and_honest() -> None:
    arguments = {
        "field": Fraction(1),
        "coupling": Fraction(1),
        "rules_frozen_at_commit": HEAD,
        "lengths": (4,),
        "steps": (8, 12),
    }
    left = run_heldout(**arguments)
    right = run_heldout(**arguments)
    assert heldout._canonical_json_bytes(left) == heldout._canonical_json_bytes(right)
    assert [row["status"] for row in left["rows"]] == [
        "unsupported",
        "unsupported",
    ]
    assert all(row["candidate_certificate_upper"] is None for row in left["rows"])
    assert all(
        row["status_reason"] == "candidate_compiler_model_family_unsupported"
        for row in left["rows"]
    )
    assert all(row["published_bound_status"] == "certified" for row in left["rows"])
    assert all(
        "diagnostic_not_a_certificate" in row["actual_operator_error_status"]
        for row in left["rows"]
    )


@pytest.mark.parametrize(
    "mutation",
    (
        "rules",
        "status",
        "scope",
        "grouped_transfer",
        "top_actual_status",
        "row_actual_status",
        "published_status",
        "actual_error",
        "row_step_float",
        "row_fraction_float",
        "bound",
        "certificate",
        "certificate_fraction_float",
        "provenance",
        "runtime_self_digest",
    ),
)
def test_heldout_verifier_rejects_mutations(mutation: str) -> None:
    payload = run_heldout(
        field=Fraction(1),
        coupling=Fraction(1),
        rules_frozen_at_commit=HEAD,
        lengths=(4,),
        steps=(8,),
    )
    changed = deepcopy(payload)
    if mutation == "rules":
        changed["frozen_rules"]["rules_digest"] = "0" * 64
    elif mutation == "status":
        changed["rows"][0]["status"] = "certified"
    elif mutation == "scope":
        changed["candidate_compiler_scope"] = "universal"
    elif mutation == "grouped_transfer":
        changed["candidate_grouped_transfer_claimed"] = True
    elif mutation == "top_actual_status":
        changed["actual_error_status"] = "certified"
    elif mutation == "row_actual_status":
        changed["rows"][0]["actual_operator_error_status"] = "certified"
    elif mutation == "published_status":
        changed["rows"][0]["published_bound_status"] = "diagnostic"
    elif mutation == "actual_error":
        changed["rows"][0]["actual_operator_error"] = "0.00000000e+00"
    elif mutation == "row_step_float":
        changed["rows"][0]["steps"] = 8.0
    elif mutation == "row_fraction_float":
        changed["rows"][0]["field"] = [1.0, 1]
    elif mutation == "bound":
        changed["rows"][0]["published_upper_exact"] = [0, 1]
    elif mutation == "certificate":
        changed["certificates"][0]["canonical_terms"][0][2] -= 1
    elif mutation == "certificate_fraction_float":
        changed["certificates"][0]["field"] = [1.0, 1]
    elif mutation == "provenance":
        path = next(iter(changed["provenance"]["scientific_dependency_sha256"]))
        changed["provenance"]["scientific_dependency_sha256"][path] = "0" * 64
    else:
        changed["provenance"]["runtime_environment"]["numpy"] = "attacker"
        changed["provenance_digest"] = heldout._canonical_digest(changed["provenance"])
    with pytest.raises((TypeError, ValueError)):
        verify_heldout_payload(changed)


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("schema_version", True),
        ("requested_lengths", [4.0]),
        ("requested_lengths", [True]),
        ("requested_steps", [8.0]),
        ("requested_steps", [True]),
        ("field", [True, 1]),
        ("field", [1.0, 1]),
        ("field", [2, 2]),
        ("field", [1, 0]),
    ),
)
def test_heldout_verifier_rejects_noncanonical_numeric_types(
    field: str,
    value: object,
) -> None:
    payload = run_heldout(
        field=Fraction(1),
        coupling=Fraction(1),
        rules_frozen_at_commit=HEAD,
        lengths=(4,),
        steps=(8,),
    )
    payload[field] = value
    with pytest.raises((TypeError, ValueError)):
        verify_heldout_payload(payload)


@pytest.mark.parametrize("variant", ("leading_space", "extra_digit"))
def test_actual_diagnostic_requires_canonical_scientific_string(
    variant: str,
) -> None:
    payload = run_heldout(
        field=Fraction(1),
        coupling=Fraction(1),
        rules_frozen_at_commit=HEAD,
        lengths=(4,),
        steps=(8,),
    )
    original = payload["rows"][0]["actual_operator_error"]
    payload["rows"][0]["actual_operator_error"] = (
        f" {original}" if variant == "leading_space" else original.replace("e", "0e")
    )
    with pytest.raises(ValueError, match="actual diagnostic mismatch"):
        verify_heldout_payload(payload)


def test_frozen_artifact_is_complete_canonical_and_fail_closed() -> None:
    path = ROOT / "benchmarks/paper-a/heldout-tfim.json"
    raw = path.read_bytes()
    payload = json.loads(raw)
    assert raw == heldout._canonical_json_bytes(payload)
    assert payload["requested_lengths"] == list(FROZEN_LENGTHS)
    assert payload["requested_steps"] == list(FROZEN_STEPS)
    assert len(payload["rows"]) == 18
    assert verify_heldout_payload(payload)
    assert all(
        row["status"] in {"unsupported", "inconclusive"} for row in payload["rows"]
    )
    assert all(row["status"] == "unsupported" for row in payload["rows"])
    for row in payload["rows"]:
        exact = Fraction(*row["published_upper_exact"])
        assert Fraction(Decimal(row["published_upper"])) >= exact
        assert Decimal(row["actual_operator_error"]) <= Decimal(row["published_upper"])
    assert payload["provenance"]["scientific_dependency_sha256"] == {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in SCIENTIFIC_DEPENDENCY_PATHS
    }
