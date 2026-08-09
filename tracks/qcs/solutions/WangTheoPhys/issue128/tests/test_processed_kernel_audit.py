from pathlib import Path

import pytest

from scripts.build_processed_kernel_audit import build_payload, write_payload
from scripts.verify_processed_kernel_audit import verify_payload


def test_processed_kernel_payload_round_trip_and_corruption(
    tmp_path: Path,
) -> None:
    path = tmp_path / "audit.json"
    payload = build_payload(("s10", "s11"))
    write_payload(path, payload)
    verified = verify_payload(path)
    assert verified == ("s10", "s11")

    payload["kernels"][0]["stage_count"] = 60
    write_payload(path, payload)
    with pytest.raises(ValueError, match="stage count"):
        verify_payload(path)


def test_processed_kernel_payload_rejects_bound_and_hash_corruption(
    tmp_path: Path,
) -> None:
    path = tmp_path / "audit.json"
    payload = build_payload(("s10",))
    payload["kernels"][0]["word_l1"]["processed_degree7"] = [0, 1]
    write_payload(path, payload)
    with pytest.raises(ValueError, match="word-l1"):
        verify_payload(path)

    payload = build_payload(("s10",))
    payload["implementation_sha256"]["effective_processor.py"] = "0" * 64
    write_payload(path, payload)
    with pytest.raises(ValueError, match="implementation hash"):
        verify_payload(path)


def test_processed_kernel_payload_requires_canonical_json(tmp_path: Path) -> None:
    path = tmp_path / "audit.json"
    payload = build_payload(("s10",))
    path.write_text(str(payload))
    with pytest.raises(ValueError, match="canonical JSON"):
        verify_payload(path)
