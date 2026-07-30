from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from trottercert.verify import _verify_d5_sidecar, verify_certificate


ROOT = Path(__file__).resolve().parents[1]


def test_d5_integrated_certificate_closes_at_95_steps() -> None:
    path = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"
    result = verify_certificate(path)
    assert result["valid"] is True
    assert result["candidate_steps"] == 95
    assert result["candidate_group_exponentials"] == 2_851
    assert result["d5_term_count"] == 605_832
    assert result["d5_group_count"] == 123_106
    assert result["exact_improvement_ratio"] == "11791/2851"
    data = json.loads(path.read_text())
    tolerance = Fraction(*data["benchmark"]["tolerance"])
    assert Fraction(*data["candidate"]["global_error_upper"]) <= tolerance
    assert Fraction(*data["candidate"]["previous_step_error_upper"]) > tolerance
    assert data["claimed_resources"]["candidate_bond_propagators"] == 205_272
    assert data["claimed_resources"]["candidate_cnot_upper"] == 615_816


def test_d5_integrated_certificate_rejects_detached_proof_method(
    tmp_path: Path,
) -> None:
    source = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"
    data = json.loads(source.read_text())
    data["candidate"]["proof_method"] = (
        "local_log_E5_grouped_D4_plus_E7_majorant_plus_exact_generator_tail"
    )
    path = tmp_path / "certificate.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="proof method and sidecar disagree"):
        verify_certificate(path)


def test_d5_integrated_certificate_rejects_digest_corruption() -> None:
    path = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"
    candidate = json.loads(path.read_text())["candidate"]
    candidate["d5_certificate"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="D5 sidecar digest mismatch"):
        _verify_d5_sidecar(path, candidate)
