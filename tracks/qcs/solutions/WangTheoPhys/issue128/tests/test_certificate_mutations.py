from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
import os
from pathlib import Path
import shutil

import pytest

from trottercert.verify import verify_certificate


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_CERTIFICATE = ROOT / "certificates" / "issue128-certificate.json"
D5_CERTIFICATE = (
    ROOT / "certificates" / "issue128-d5-integrated-certificate.json"
)


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _link_or_copy(source: Path, destination: Path) -> None:
    try:
        os.link(source, destination)
    except OSError:
        shutil.copyfile(source, destination)


def _write_bundle(
    root: Path,
    source_certificate: Path,
    data: dict[str, object],
) -> Path:
    root.mkdir()
    candidate = data["candidate"]
    for key in ("d4_certificate", "d5_certificate", "d6_certificate"):
        metadata = candidate.get(key)
        if not isinstance(metadata, dict):
            continue
        for path_key in (
            "path",
            "parent_path",
            "comparison_path",
            "groups_path",
        ):
            relative = metadata.get(path_key)
            if relative is None:
                continue
            source = source_certificate.parent / str(relative)
            destination = root / str(relative)
            destination.parent.mkdir(parents=True, exist_ok=True)
            _link_or_copy(source, destination)
    path = root / source_certificate.name
    path.write_text(json.dumps(data))
    return path


def _lower_contribution_and_retotal(
    data: dict[str, object],
    name: str,
) -> None:
    contributions = data["candidate"]["contributions"]
    contributions[name] = [0, 1]
    total = sum(
        (
            Fraction(*contributions[key])
            for key in ("degree4", "degree5", "degree6", "degree7", "tail")
        ),
        Fraction(),
    )
    data["candidate"]["global_error_upper"] = _pair(total)


@pytest.mark.parametrize("name", ["degree6", "degree7", "tail"])
def test_fast_verifier_rejects_retotalled_ledger_attack(
    tmp_path: Path,
    name: str,
) -> None:
    data = json.loads(PUBLIC_CERTIFICATE.read_text())
    _lower_contribution_and_retotal(data, name)
    path = _write_bundle(tmp_path / name, PUBLIC_CERTIFICATE, data)

    with pytest.raises(
        ValueError,
        match=rf"candidate {name} contribution regeneration mismatch",
    ):
        verify_certificate(path)


def test_fast_verifier_rejects_false_adjacent_step_value(
    tmp_path: Path,
) -> None:
    data = json.loads(PUBLIC_CERTIFICATE.read_text())
    data["candidate"]["previous_step_error_upper"] = [1, 100_000]
    path = _write_bundle(tmp_path / "adjacent", PUBLIC_CERTIFICATE, data)

    with pytest.raises(
        ValueError,
        match="candidate adjacent-step regeneration mismatch",
    ):
        verify_certificate(path)


def test_fast_verifier_rejects_formula_constant_attack(
    tmp_path: Path,
) -> None:
    data = json.loads(PUBLIC_CERTIFICATE.read_text())
    value = Fraction(*data["candidate"]["e7_site_majorant"])
    data["candidate"]["e7_site_majorant"] = _pair(value / 2)
    path = _write_bundle(tmp_path / "e7", PUBLIC_CERTIFICATE, data)

    with pytest.raises(ValueError, match="E7 site majorant regeneration mismatch"):
        verify_certificate(path)


def test_fast_verifier_rejects_d5_tail_attack(
    tmp_path: Path,
) -> None:
    data = json.loads(D5_CERTIFICATE.read_text())
    _lower_contribution_and_retotal(data, "tail")
    path = _write_bundle(tmp_path / "d5-tail", D5_CERTIFICATE, data)

    with pytest.raises(
        ValueError,
        match="candidate tail contribution regeneration mismatch",
    ):
        verify_certificate(path)


def test_verifier_rejects_boolean_fraction_components(
    tmp_path: Path,
) -> None:
    data = json.loads(PUBLIC_CERTIFICATE.read_text())
    data["benchmark"]["tolerance"] = [True, 1_000_000]
    path = _write_bundle(tmp_path / "boolean", PUBLIC_CERTIFICATE, data)

    with pytest.raises(
        ValueError,
        match="fraction numerator and denominator must be integers",
    ):
        verify_certificate(path)


def test_verifier_rejects_noninteger_schema_and_step_fields(
    tmp_path: Path,
) -> None:
    original = json.loads(PUBLIC_CERTIFICATE.read_text())

    schema_mutation = deepcopy(original)
    schema_mutation["schema_version"] = True
    schema_path = _write_bundle(
        tmp_path / "schema",
        PUBLIC_CERTIFICATE,
        schema_mutation,
    )
    with pytest.raises(
        ValueError,
        match="certificate schema version must be an integer",
    ):
        verify_certificate(schema_path)

    step_mutation = deepcopy(original)
    step_mutation["candidate"]["steps"] = "97"
    step_path = _write_bundle(
        tmp_path / "steps",
        PUBLIC_CERTIFICATE,
        step_mutation,
    )
    with pytest.raises(ValueError, match="candidate step count must be an integer"):
        verify_certificate(step_path)


def test_verifier_rejects_coercible_integer_metadata(
    tmp_path: Path,
) -> None:
    original = json.loads(PUBLIC_CERTIFICATE.read_text())

    metadata_mutation = deepcopy(original)
    metadata_mutation["candidate"]["d4_certificate"]["term_count"] = "75324"
    metadata_path = _write_bundle(
        tmp_path / "metadata",
        PUBLIC_CERTIFICATE,
        metadata_mutation,
    )
    with pytest.raises(ValueError, match="D4 term count must be an integer"):
        verify_certificate(metadata_path)

    precision_mutation = deepcopy(original)
    precision_mutation["candidate"]["coefficient_interval_decimal_digits"] = "12"
    precision_path = _write_bundle(
        tmp_path / "precision",
        PUBLIC_CERTIFICATE,
        precision_mutation,
    )
    with pytest.raises(
        ValueError,
        match="candidate coefficient precision must be an integer",
    ):
        verify_certificate(precision_path)


def test_verifier_rejects_resource_and_ratio_drift(
    tmp_path: Path,
) -> None:
    original = json.loads(PUBLIC_CERTIFICATE.read_text())

    resource_mutation = deepcopy(original)
    resource_mutation["claimed_resources"]["candidate_cnot_upper"] -= 1
    resource_path = _write_bundle(
        tmp_path / "resources",
        PUBLIC_CERTIFICATE,
        resource_mutation,
    )
    with pytest.raises(ValueError, match="resource summary mismatch"):
        verify_certificate(resource_path)

    ratio_mutation = deepcopy(original)
    ratio_mutation["claims"]["exact_improvement_ratio"] = [4, 1]
    ratio_path = _write_bundle(
        tmp_path / "ratio",
        PUBLIC_CERTIFICATE,
        ratio_mutation,
    )
    with pytest.raises(ValueError, match="claimed improvement ratio mismatch"):
        verify_certificate(ratio_path)
