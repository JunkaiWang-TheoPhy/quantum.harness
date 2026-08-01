from __future__ import annotations

import json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import pytest

from scripts.run_xxz_transfer import (
    CERTIFICATION_ONLY_FIELDS,
    PREREGISTERED_DELTAS,
    canonical_json,
    run_xxz_transfer,
)
from trottercert.hamiltonian import xxz_bond

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks" / "paper-a" / "xxz-transfer.json"


def _fixture_metrics(delta: Fraction) -> dict[str, object]:
    certified = delta == 1
    return {
        "status": "certified" if certified else "unsupported",
        "reason": "three-term exact fixture" if certified else "fixture unsupported",
        "d4_term_count": 3 if certified else None,
        "d4_group_count": 1 if certified else None,
        "d4_norm_bound": [1, 8] if certified else None,
        "d5_term_count": 3 if certified else None,
        "d5_group_count": 1 if certified else None,
        "d5_norm_bound": [1, 16] if certified else None,
        "accepted_steps": 7 if certified else None,
        "builder_seconds": None,
        "verifier_seconds": None,
        "source_certificate": "fixture.json" if certified else None,
    }


def test_xxz_bond_reduces_to_heisenberg() -> None:
    assert xxz_bond(Fraction(1)) == {
        "XX": Fraction(1, 4),
        "YY": Fraction(1, 4),
        "ZZ": Fraction(1, 4),
    }


def test_xxz_bond_uses_exact_anisotropy_normalization() -> None:
    assert xxz_bond(Fraction(3, 2)) == {
        "XX": Fraction(1, 4),
        "YY": Fraction(1, 4),
        "ZZ": Fraction(3, 8),
    }
    with pytest.raises(TypeError, match="delta must be a Fraction"):
        xxz_bond(1)  # type: ignore[arg-type]


def test_preregistered_deltas_are_exact() -> None:
    assert PREREGISTERED_DELTAS == (
        Fraction(0),
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(1),
        Fraction(3, 2),
        Fraction(2),
        Fraction(4),
    )


def test_every_preregistered_delta_has_an_explicit_status_row() -> None:
    payload = run_xxz_transfer(evaluator=_fixture_metrics)
    assert [row["delta"] for row in payload["rows"]] == [
        [delta.numerator, delta.denominator] for delta in PREREGISTERED_DELTAS
    ]
    assert len(payload["rows"]) == 7
    assert all(
        row["status"] in {"certified", "unsupported", "inconclusive"}
        for row in payload["rows"]
    )
    assert sum(row["status"] == "unsupported" for row in payload["rows"]) == 6


def test_small_fixture_serialization_is_byte_deterministic() -> None:
    first = canonical_json(run_xxz_transfer(evaluator=_fixture_metrics))
    second = canonical_json(run_xxz_transfer(evaluator=_fixture_metrics))
    assert first.encode() == second.encode()


@pytest.mark.parametrize("field", CERTIFICATION_ONLY_FIELDS)
def test_uncertified_fixture_nulls_every_certification_field(field: str) -> None:
    metrics = _fixture_metrics(Fraction(0))
    metrics[field] = 7
    with pytest.raises(ValueError, match=rf"must set {field} to null"):
        run_xxz_transfer((Fraction(0),), evaluator=lambda _delta: metrics)


def test_frozen_payload_rejects_execution_timing() -> None:
    metrics = _fixture_metrics(Fraction(1))
    metrics["verifier_seconds"] = "0.250000"
    with pytest.raises(ValueError, match="timing must not enter"):
        run_xxz_transfer((Fraction(1),), evaluator=lambda _delta: metrics)


def test_custom_evaluator_cannot_certify_nonisotropic_delta() -> None:
    forged = _fixture_metrics(Fraction(1))
    with pytest.raises(ValueError, match="can certify only delta=1"):
        run_xxz_transfer((Fraction(1, 2),), evaluator=lambda _delta: forged)


def test_frozen_sweep_is_canonical_complete_and_honest() -> None:
    raw = ARTIFACT.read_text()
    payload = json.loads(raw)
    assert raw == canonical_json(payload)
    assert payload["execution_timing_included"] is False
    for source in payload["sources"].values():
        path = ROOT / source["path"]
        assert source["sha256"] == sha256(path.read_bytes()).hexdigest()
    assert [row["delta"] for row in payload["rows"]] == [
        [delta.numerator, delta.denominator] for delta in PREREGISTERED_DELTAS
    ]
    assert len(payload["rows"]) == 7
    assert all(
        row["builder_seconds"] is None and row["verifier_seconds"] is None
        for row in payload["rows"]
    )

    isotropic = next(row for row in payload["rows"] if row["delta"] == [1, 1])
    assert isotropic["status"] == "certified"
    certificate = json.loads(
        (ROOT / payload["sources"]["certificate"]["path"]).read_text()
    )
    candidate = certificate["candidate"]
    d4 = candidate["d4_certificate"]
    d5 = candidate["d5_certificate"]
    assert isotropic["accepted_steps"] == candidate["steps"] == 95
    assert isotropic["d4_term_count"] == d4["term_count"] == 75_324
    assert isotropic["d4_group_count"] == d4["group_count"]
    assert isotropic["d4_norm_bound"] == d4["cell_norm_upper"]
    assert isotropic["d5_term_count"] == d5["term_count"] == 605_832
    assert isotropic["d5_group_count"] == d5["group_count"]
    assert isotropic["d5_norm_bound"] == d5["site_norm_upper"]
    assert isotropic["source_certificate"] == payload["sources"]["certificate"]["path"]
    assert isotropic["builder_seconds"] is None
    assert isotropic["verifier_seconds"] is None

    unsupported = [row for row in payload["rows"] if row["status"] == "unsupported"]
    assert len(unsupported) == 6
    assert all(
        all(row[field] is None for field in CERTIFICATION_ONLY_FIELDS)
        for row in unsupported
    )
