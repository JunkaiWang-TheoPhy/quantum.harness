from pathlib import Path

import pytest

from scripts.report_processed_local_budget import (
    build_report,
    verify_report,
    write_report,
)


def test_pre_e7_report_round_trip_and_missing_proof_status(
    tmp_path: Path,
) -> None:
    path = tmp_path / "pre-e7.json"
    payload = build_report()
    write_report(path, payload)
    assert verify_report(path) == (
        ("s10", 39),
        ("s10", 47),
        ("s11", 35),
        ("s11", 43),
    )
    for point in payload["points"]:
        assert point["status"] == "awaiting_degree7_and_local_log_remainder"
        assert point["degree7_site_l1"] is None
        assert point["higher_order_remainder"] is None
        assert point["local_log_theorem"] is None
        assert point["degree3_residual_total"] != [0, 1]
        assert point["degree5_residual_total"] != [0, 1]


def test_pre_e7_report_rejects_false_acceptance(tmp_path: Path) -> None:
    path = tmp_path / "pre-e7.json"
    payload = build_report()
    payload["points"][0]["status"] = "accepted"
    write_report(path, payload)
    with pytest.raises(ValueError, match="status"):
        verify_report(path)


def test_pre_e7_report_rejects_processor_bound_corruption(
    tmp_path: Path,
) -> None:
    path = tmp_path / "pre-e7.json"
    payload = build_report()
    payload["points"][0]["processor_distance_upper"] = [0, 1]
    write_report(path, payload)
    with pytest.raises(ValueError, match="processor distance"):
        verify_report(path)
