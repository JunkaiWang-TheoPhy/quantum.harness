from __future__ import annotations

import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import pytest

from scripts.certify_finite_step_obstruction import (
    DEFAULT_E5,
    DEFAULT_E7,
    DEFAULT_TAIL,
    _cubic_json,
    _digest,
    _load,
    _parse_cubic,
    build_payload,
    verify_payload,
)
from trottercert.cubic_field import Cubic


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")


def _fake_e9(path: Path, *, q9: Cubic = Cubic.zero()) -> None:
    e5 = _load(DEFAULT_E5, "E5 artifact")
    e7 = _load(DEFAULT_E7, "E7 artifact")
    q5 = _parse_cubic(e5["stable_relations"]["tau_w_e5_per_cell"], "q5")
    q7 = _parse_cubic(e7["pairings"]["tau_w_per_cell"], "q7")
    payload = {
        "schema_version": 1,
        "kind": "issue128_dual_e9_pairing",
        "degree": 9,
        "length": 12,
        "shard_count": 64,
        "manifest_index_sha256": "1" * 64,
        "inputs": {
            "extensive_witness_sha256": _sha(DEFAULT_E5),
            "dual_e7_sha256": _sha(DEFAULT_E7),
            "q5_per_cell": _cubic_json(q5),
            "q7_per_cell": _cubic_json(q7),
        },
        "pairings": {
            "tau_w_per_cell": _cubic_json(q9),
            "q9_over_97_pow_8": _cubic_json(q9 / 97**8),
            "exact_e5_e7_e9_per_cell_at_r97": _cubic_json(
                q5 / 97**4 + q7 / 97**6 + q9 / 97**8
            ),
        },
        "claim": {
            "full_e9_operator": "not_computed",
            "dual_e9_pairing": "exact",
            "finite_step_status": "inconclusive",
            "missing": "E11-and-higher dual tail",
        },
    }
    payload["mathematical_payload_sha256"] = _digest(payload)
    _write_json(path, payload)


def test_certificate_builds_from_digest_bound_inputs(tmp_path: Path) -> None:
    e9 = tmp_path / "dual-e9.json"
    _fake_e9(e9)

    payload = build_payload(e9_path=e9, index_path=None, require_full_e9=False)

    assert payload["series"] == "h^4 q5 + h^6 q7 + h^8 q9 + R_{>=11}^{dual}"
    assert payload["step_size"] == [1, 97]
    assert payload["sources"]["e9"]["sha256"] == _sha(e9)
    assert payload["claim"]["promotion_rule"] == "signed_margin.lower > 0"
    verify_payload(payload, e9_path=e9)


def test_build_requires_full_manifest_replay_by_default(tmp_path: Path) -> None:
    e9 = tmp_path / "dual-e9.json"
    _fake_e9(e9)

    with pytest.raises(ValueError, match="manifest index is required"):
        build_payload(e9_path=e9, index_path=None)


def test_mutated_certificate_status_is_rejected(tmp_path: Path) -> None:
    e9 = tmp_path / "dual-e9.json"
    _fake_e9(e9)
    payload = build_payload(e9_path=e9, index_path=None, require_full_e9=False)
    forged = copy.deepcopy(payload)
    forged["claim"]["finite_step_status"] = "inconclusive"
    forged["claim"]["promoted"] = False

    with pytest.raises(ValueError, match="regeneration mismatch"):
        verify_payload(forged, e9_path=e9)


def test_mutated_e9_pairing_is_rejected_before_promotion(tmp_path: Path) -> None:
    e9 = tmp_path / "dual-e9.json"
    _fake_e9(e9)
    payload = _load(e9, "E9 artifact")
    payload["pairings"]["tau_w_per_cell"][0][0] += 1
    payload["mathematical_payload_sha256"] = _digest(
        {key: value for key, value in payload.items() if key != "mathematical_payload_sha256"}
    )
    _write_json(e9, payload)

    with pytest.raises(ValueError, match="scaled q9 mismatch"):
        build_payload(e9_path=e9, index_path=None, require_full_e9=False)


def test_e9_digest_drift_invalidates_existing_certificate(tmp_path: Path) -> None:
    e9 = tmp_path / "dual-e9.json"
    _fake_e9(e9)
    certificate = build_payload(e9_path=e9, index_path=None, require_full_e9=False)
    _fake_e9(e9, q9=Cubic(Fraction(1, 10**6), 0, 0))

    with pytest.raises(ValueError, match="regeneration mismatch"):
        verify_payload(certificate, e9_path=e9)


def test_tail_shrink_attack_fails_regeneration(tmp_path: Path) -> None:
    e9 = tmp_path / "dual-e9.json"
    _fake_e9(e9)
    tail = tmp_path / "dual-log-tail.json"
    tail.write_bytes(DEFAULT_TAIL.read_bytes())
    payload = build_payload(
        e9_path=e9,
        tail_path=tail,
        index_path=None,
        require_full_e9=False,
    )
    forged_tail = _load(tail, "tail")
    forged_tail["bounds"]["total_log_tail"][0] //= 2
    _write_json(tail, forged_tail)

    with pytest.raises(ValueError, match="payload regeneration mismatch"):
        verify_payload(payload, e9_path=e9, tail_path=tail)
