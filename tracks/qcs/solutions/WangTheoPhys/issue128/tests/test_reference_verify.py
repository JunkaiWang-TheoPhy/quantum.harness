from __future__ import annotations

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil

import pytest

from scripts.reference_verify import verify


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"


def _bundle(tmp_path: Path, data: dict[str, object]) -> Path:
    root = tmp_path / "bundle"
    root.mkdir()
    for key in ("d4_certificate", "d5_certificate"):
        metadata = data["candidate"].get(key)
        if not isinstance(metadata, dict):
            continue
        source = CERTIFICATE.parent / str(metadata["path"])
        destination = root / str(metadata["path"])
        try:
            os.link(source, destination)
        except OSError:
            shutil.copyfile(source, destination)
    path = root / CERTIFICATE.name
    path.write_text(json.dumps(data))
    return path


def test_reference_verifier_accepts_d5_integrated_certificate() -> None:
    result = verify(CERTIFICATE)
    assert result["valid"] is True
    assert result["candidate_steps"] == 95
    assert result["d4_term_count"] == 75_324
    assert result["d5_term_count"] == 605_832
    assert result["exact_improvement_ratio"] == "11791/2851"
    assert result["coefficient_generation_replayed"] is False


def test_reference_verifier_rejects_retotalled_tail_attack(
    tmp_path: Path,
) -> None:
    data = json.loads(CERTIFICATE.read_text())
    mutated = deepcopy(data)
    contributions = mutated["candidate"]["contributions"]
    contributions["tail"] = [0, 1]
    from fractions import Fraction

    total = sum(
        (
            Fraction(*contributions[name])
            for name in ("degree4", "degree5", "degree6", "degree7", "tail")
        ),
        Fraction(),
    )
    mutated["candidate"]["global_error_upper"] = [
        total.numerator,
        total.denominator,
    ]

    with pytest.raises(ValueError, match="candidate tail contribution mismatch"):
        verify(_bundle(tmp_path, mutated))


@pytest.mark.parametrize(
    ("certificate_key", "metadata_key", "bad_value", "message"),
    (
        (
            "d4_certificate",
            "term_count",
            -1,
            "D4 term count must be at least 0",
        ),
        (
            "d4_certificate",
            "group_count",
            7_577,
            "D4 group count metadata mismatch",
        ),
        (
            "d4_certificate",
            "max_group_size",
            True,
            "D4 maximum group size must be an integer",
        ),
        (
            "d5_certificate",
            "term_count",
            "605832",
            "D5 term count must be an integer",
        ),
        (
            "d5_certificate",
            "group_count",
            123_107,
            "D5 group count metadata mismatch",
        ),
        (
            "d5_certificate",
            "max_group_size",
            True,
            "D5 maximum group size must be an integer",
        ),
    ),
)
def test_reference_verifier_rejects_partition_metadata_mutations(
    tmp_path: Path,
    certificate_key: str,
    metadata_key: str,
    bad_value: object,
    message: str,
) -> None:
    data = json.loads(CERTIFICATE.read_text())
    mutated = deepcopy(data)
    mutated["candidate"][certificate_key][metadata_key] = bad_value

    with pytest.raises(ValueError, match=message):
        verify(_bundle(tmp_path, mutated))
