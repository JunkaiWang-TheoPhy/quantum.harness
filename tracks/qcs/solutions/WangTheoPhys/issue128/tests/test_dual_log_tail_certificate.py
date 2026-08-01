from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.certify_dual_log_tail import build_payload, verify_payload

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "docs/experiments/processor-obstruction/dual-log-tail.json"
)


def _integer_text_lengths(value: object) -> list[int]:
    if isinstance(value, bool):
        return []
    if isinstance(value, int):
        return [len(str(abs(value)))]
    if isinstance(value, dict):
        return [
            length
            for child in value.values()
            for length in _integer_text_lengths(child)
        ]
    if isinstance(value, list):
        return [
            length
            for child in value
            for length in _integer_text_lengths(child)
        ]
    return []


def test_payload_round_trip_is_source_regenerated() -> None:
    payload = build_payload()
    verify_payload(payload)
    assert payload["kind"] == "issue128_dual_log_tail"
    assert payload["claim"] == {
        "dual_e11_plus_tail": "certified",
        "finite_step_status": "inconclusive",
        "missing": "exact dual E9 pairing",
    }
    assert payload["bounds"]["total_log_tail"] == build_payload()[
        "bounds"
    ]["total_log_tail"]
    assert max(_integer_text_lengths(payload)) > 4300


def test_checked_in_artifact_is_canonical_and_current() -> None:
    encoded = ARTIFACT.read_bytes()
    payload = json.loads(encoded)
    verify_payload(payload)
    assert encoded == (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    assert payload == build_payload()


@pytest.mark.parametrize(
    ("path", "message"),
    [
        (("instance", "steps"), "payload regeneration mismatch"),
        (("moments", "w_hs_squared"), "payload regeneration mismatch"),
        (("analytic_caps", "dexp"), "payload regeneration mismatch"),
        (("envelope", "log_defect"), "payload regeneration mismatch"),
        (
            ("bounds", "direct_even_generator_tail"),
            "payload regeneration mismatch",
        ),
        (
            ("bounds", "dexp_correction_tail"),
            "payload regeneration mismatch",
        ),
        (
            ("claim", "finite_step_status"),
            "payload regeneration mismatch",
        ),
    ],
)
def test_payload_mutations_fail_closed(
    path: tuple[str, str],
    message: str,
) -> None:
    forged = copy.deepcopy(build_payload())
    parent, child = path
    value = forged[parent][child]
    forged[parent][child] = (
        "forged" if isinstance(value, str) else [0, 1]
    )
    with pytest.raises(ValueError, match=message):
        verify_payload(forged)
