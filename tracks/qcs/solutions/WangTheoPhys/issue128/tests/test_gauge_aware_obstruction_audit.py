import copy
import json
from pathlib import Path

import pytest

from scripts.audit_gauge_aware_obstruction import build_payload, verify_payload

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/experiments/processor-obstruction/exact-obstruction.json"
WITNESS = (
    ROOT
    / "docs/experiments/processor-obstruction/quadratic-commutant-witness.json"
)
EXTENSIVE = (
    ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
DUAL_E7 = (
    ROOT
    / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
)


def test_current_obstruction_has_calibrated_leading_order_witness() -> None:
    payload = build_payload(SOURCE)
    assert payload["fixed_time_endpoint_processor"]["status"] == "no_go"
    calibrated = payload["calibrated_spectral_obstruction"]
    assert calibrated["leading_order_status"] == "no_go"
    assert calibrated["finite_step_status"] == "inconclusive"
    assert calibrated["witness"] == "H^2 - 54 I + H/2"
    assert calibrated["witness_artifact"]["sha256"]
    assert calibrated["extensive_witness_artifact"]["sha256"]
    assert calibrated["verified_lengths"] == [6, 8, 10, 12]
    assert calibrated["rejected_alias_lengths"] == [4]
    assert calibrated["finite_step_gate"] == "dual_pairing_remainder"
    assert calibrated["squared_normalized_pairing_l12"]
    assert calibrated["dual_e7_artifact"]["sha256"]
    assert calibrated["dual_e7_pairing_status"] == "exact"
    assert calibrated["truncated_log_gate"]["95"]["sign"] == "negative"
    assert calibrated["truncated_log_gate"]["97"]["branch_path_lt_pi"] is True
    assert calibrated["truncated_log_gate"]["96"]["branch_path_lt_pi"] is False
    assert "E9-and-higher dual tail" in calibrated["finite_step_missing"]
    assert "running E7" not in payload["next_gate"]
    assert (
        payload["reference_examples"]["hamiltonian_parallel"][
            "calibrated_status"
        ]
        == "removable"
    )
    assert (
        payload["reference_examples"]["independent_commutant"]["status"]
        == "spectral_obstruction"
    )
    verify_payload(payload, SOURCE)


def test_audit_rejects_forged_status_and_digest() -> None:
    payload = build_payload(SOURCE)
    forged = copy.deepcopy(payload)
    forged["calibrated_spectral_obstruction"]["finite_step_status"] = "no_go"
    with pytest.raises(ValueError, match="calibrated"):
        verify_payload(forged, SOURCE)
    forged = copy.deepcopy(payload)
    forged["source"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="digest"):
        verify_payload(forged, SOURCE)


def test_audit_rejects_mutated_or_missing_witness(tmp_path: Path) -> None:
    witness_payload = json.loads(WITNESS.read_text())
    witness_payload["pairings"]["tau_w_e5_nonzero"] = False
    forged_witness = tmp_path / "forged-witness.json"
    forged_witness.write_text(json.dumps(witness_payload))
    with pytest.raises(ValueError, match="witness"):
        build_payload(SOURCE, forged_witness)

    with pytest.raises(ValueError, match="witness"):
        build_payload(SOURCE, tmp_path / "missing-witness.json")


def test_audit_rejects_mutated_or_missing_extensive_witness(tmp_path: Path) -> None:
    extensive_payload = json.loads(EXTENSIVE.read_text())
    extensive_payload["claim"]["finite_step_status"] = "proved"
    forged = tmp_path / "forged-extensive.json"
    forged.write_text(json.dumps(extensive_payload))
    with pytest.raises(ValueError, match="extensive"):
        build_payload(SOURCE, WITNESS, forged)

    with pytest.raises(ValueError, match="extensive"):
        build_payload(SOURCE, WITNESS, tmp_path / "missing-extensive.json")


def test_audit_rejects_mutated_or_missing_dual_e7(tmp_path: Path) -> None:
    dual_payload = json.loads(DUAL_E7.read_text())
    dual_payload["claim"]["dual_e7_pairing"] = "forged"
    forged = tmp_path / "forged-dual-e7.json"
    forged.write_text(json.dumps(dual_payload))
    with pytest.raises(ValueError, match="dual E7"):
        build_payload(SOURCE, WITNESS, EXTENSIVE, forged)

    with pytest.raises(ValueError, match="dual E7"):
        build_payload(
            SOURCE,
            WITNESS,
            EXTENSIVE,
            tmp_path / "missing-dual-e7.json",
        )


def test_payload_round_trips_through_canonical_json() -> None:
    payload = build_payload(SOURCE)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    decoded = json.loads(encoded)
    assert decoded == payload
    verify_payload(decoded, SOURCE)


def test_missing_or_false_source_witness_is_rejected(tmp_path: Path) -> None:
    missing = tmp_path / "missing.json"
    with pytest.raises(ValueError, match="read source"):
        build_payload(missing)

    source_payload = json.loads(SOURCE.read_text())
    source_payload["e5_hilbert_schmidt_overlap_with_h_per_cell"]["nonzero"] = False
    false_source = tmp_path / "false-source.json"
    false_source.write_text(json.dumps(source_payload))
    with pytest.raises(ValueError, match="nonzero"):
        build_payload(false_source)
