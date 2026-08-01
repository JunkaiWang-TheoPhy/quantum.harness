from __future__ import annotations

import json
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import pairwise
from pathlib import Path

import pytest

from scripts.build_paper_a_ablation import (
    MONOTONIC_CHAIN,
    RECOMPUTE_FUNCTION_ALLOWLIST,
    SCIENTIFIC_SOURCE_NAMES,
    STAGE_NAMES,
    build_ablation,
    canonical_json,
    compute_payload_digest,
    validate_ablation,
)
from trottercert.resources import fourth_order_four_matching_resources

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"
ARTIFACT = ROOT / "benchmarks" / "paper-a" / "ablation.json"


def _artifact() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def _by_stage(payload: dict[str, object]) -> dict[str, dict[str, object]]:
    return {row["stage"]: row for row in payload["stages"]}


def _fraction(value: object) -> Fraction:
    assert isinstance(value, list) and len(value) == 2
    return Fraction(value[0], value[1])


def _reseal(payload: dict[str, object]) -> None:
    payload["payload_sha256"] = compute_payload_digest(payload)


def test_frozen_ablation_is_canonical_and_valid() -> None:
    raw = ARTIFACT.read_text()
    payload = json.loads(raw)
    assert raw == canonical_json(payload)
    assert payload["payload_sha256"] == compute_payload_digest(payload)
    assert validate_ablation(payload, ROOT, raw.encode()) == []


def test_stage_order_and_missing_evidence_are_explicit() -> None:
    payload = _artifact()
    assert [row["stage"] for row in payload["stages"]] == list(STAGE_NAMES)
    complete_formula = _by_stage(payload)["complete_formula"]
    assert complete_formula["status"] == "missing_evidence"
    assert "no frozen pre-Pauli" in complete_formula["missing_evidence"]
    for field in ("bound", "steps", "groups", "leading_d4_site_bound"):
        assert complete_formula[field] is None


def test_final_endpoint_binds_certificate_and_recomputed_resources() -> None:
    payload = _artifact()
    certificate = json.loads(CERTIFICATE.read_text())
    final = _by_stage(payload)["finite_ledger"]
    assert payload["reference_steps"] == final["steps"] == 95
    assert final["groups"] == 2_851
    assert final["bound"] == certificate["candidate"]["global_error_upper"]
    assert final["steps"] == certificate["candidate"]["steps"]
    assert (
        final["groups"]
        == certificate["claimed_resources"]["candidate_group_exponentials"]
    )

    sites = certificate["benchmark"]["length"] ** 2
    for row in payload["stages"]:
        if row["status"] == "computed":
            expected = fourth_order_four_matching_resources(
                sites,
                row["steps"],
                stage_count=31,
            ).group_exponentials
            assert row["groups"] == expected


def test_every_frozen_source_hash_is_live() -> None:
    sources = _artifact()["sources"]
    assert set(SCIENTIFIC_SOURCE_NAMES) < set(sources)
    for source in sources.values():
        path = ROOT / source["path"]
        assert source["sha256"] == sha256(path.read_bytes()).hexdigest()


def test_recompute_function_lists_equal_the_canonical_allowlist() -> None:
    for row in _artifact()["stages"]:
        assert row["recompute_functions"] == list(
            RECOMPUTE_FUNCTION_ALLOWLIST[row["stage"]]
        )


def test_only_evidence_backed_chain_is_monotone() -> None:
    payload = _artifact()
    assert payload["monotonic_chain"] == list(MONOTONIC_CHAIN)
    rows = [_by_stage(payload)[name] for name in MONOTONIC_CHAIN]
    bounds = [_fraction(row["bound"]) for row in rows]
    steps = [row["steps"] for row in rows]
    assert all(right <= left for left, right in pairwise(bounds))
    assert all(right <= left for left, right in pairwise(steps))


@pytest.mark.parametrize(
    "mutation", [lambda rows: rows.reverse(), lambda rows: rows.pop(1)]
)
def test_validator_rejects_stage_structure_mutations(mutation) -> None:
    payload = deepcopy(_artifact())
    mutation(payload["stages"])
    _reseal(payload)
    assert "ablation stage order mismatch" in validate_ablation(payload, ROOT)


def test_attack_fake_complete_formula_computation_is_rejected() -> None:
    payload = deepcopy(_artifact())
    complete = _by_stage(payload)["complete_formula"]
    complete.update(
        {
            "status": "computed",
            "bound": [1, 10**6],
            "steps": 95,
            "groups": 2_851,
            "leading_d4_site_bound": [1, 1],
            "missing_evidence": None,
        }
    )
    _reseal(payload)
    errors = validate_ablation(payload, ROOT)
    assert any(
        "complete_formula must remain missing_evidence" in error for error in errors
    )
    assert any("authoritative canonical rebuild" in error for error in errors)


@pytest.mark.parametrize(
    ("stage", "field"),
    [
        (stage, field)
        for stage in ("pauli_aggregation", "translation_aggregation", "finite_ledger")
        for field in ("bound", "steps", "groups")
    ],
)
def test_attack_hand_filled_scientific_numbers_are_rejected(
    stage: str,
    field: str,
) -> None:
    payload = deepcopy(_artifact())
    row = _by_stage(payload)[stage]
    if field == "bound":
        value = _fraction(row[field])
        row[field] = [value.numerator + value.denominator, value.denominator]
    else:
        row[field] += 1
    _reseal(payload)
    errors = validate_ablation(payload, ROOT)
    assert "payload fields differ from authoritative canonical rebuild" in errors
    assert "payload bytes differ from authoritative canonical rebuild" in errors


@pytest.mark.parametrize(
    "stage", ["pauli_aggregation", "translation_aggregation", "finite_ledger"]
)
def test_attack_recompute_function_substitution_is_rejected(stage: str) -> None:
    payload = deepcopy(_artifact())
    _by_stage(payload)[stage]["recompute_functions"][-1] = "builtins.sum"
    _reseal(payload)
    errors = validate_ablation(payload, ROOT)
    assert f"{stage}: recompute-function allowlist mismatch" in errors
    assert "payload fields differ from authoritative canonical rebuild" in errors


def test_attack_scientific_dependency_substitution_is_rejected() -> None:
    payload = deepcopy(_artifact())
    payload["sources"]["refined_error"]["sha256"] = "0" * 64
    _reseal(payload)
    errors = validate_ablation(payload, ROOT)
    assert "authoritative source binding mismatch for refined_error" in errors
    assert "payload fields differ from authoritative canonical rebuild" in errors


def test_attack_payload_change_without_resealing_is_rejected_immediately() -> None:
    payload = deepcopy(_artifact())
    _by_stage(payload)["pauli_aggregation"]["steps"] += 1
    assert "payload digest mismatch" in validate_ablation(payload, ROOT)


def test_attack_noncanonical_serialized_bytes_are_rejected() -> None:
    raw = ARTIFACT.read_text()
    noncanonical = raw.replace("{\n", "{  \n", 1).encode()
    errors = validate_ablation(json.loads(noncanonical), ROOT, noncanonical)
    assert "payload bytes differ from authoritative canonical rebuild" in errors


def test_validator_rejects_source_path_escape(tmp_path: Path) -> None:
    outside = tmp_path / "outside.json"
    outside.write_text("{}\n")
    payload = deepcopy(_artifact())
    payload["sources"]["certificate"] = {
        "path": str(outside),
        "sha256": sha256(outside.read_bytes()).hexdigest(),
    }
    _reseal(payload)
    errors = validate_ablation(payload, ROOT)
    assert "authoritative source binding mismatch for certificate" in errors
    assert "payload fields differ from authoritative canonical rebuild" in errors


@pytest.mark.slow
def test_frozen_ablation_rebuilds_to_identical_bytes() -> None:
    rebuilt = canonical_json(build_ablation(CERTIFICATE))
    assert rebuilt.encode() == ARTIFACT.read_bytes()
