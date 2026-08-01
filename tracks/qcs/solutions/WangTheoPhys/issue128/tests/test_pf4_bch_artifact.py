from __future__ import annotations

import copy
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

from scripts import derive_pf4_bch_mapping as builder
from trottercert.pf4_bch_mapping import payload_digest as mapping_digest

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/derive_pf4_bch_mapping.py"
REFERENCE = ROOT / "scripts/reference_pf4_bch_mapping.py"
ARTIFACT = ROOT / "docs/experiments/processor-obstruction/pf4-bch-mapping.json"


def _load_reference() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "reference_pf4_bch_mapping_for_artifact_tests",
        REFERENCE,
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def reference() -> ModuleType:
    return _load_reference()


def _run(script: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *arguments],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        check=False,
        capture_output=True,
        text=True,
    )


def _write(path: Path, payload: object) -> None:
    path.write_bytes(builder.canonical_bytes(payload))


def _reseal(payload: dict[str, object]) -> None:
    payload["payload_sha256"] = builder.payload_digest(payload)


def test_wrapper_schema_and_independent_rebuild_agree(reference: ModuleType) -> None:
    primary = builder.build_payload()
    independent = reference.build_reference_wrapper()

    assert tuple(primary) == (
        "schema_version",
        "kind",
        "mapping",
        "implementation_sources",
        "reference_algorithm",
        "payload_sha256",
    )
    assert primary["schema_version"] == 2
    assert primary["reference_algorithm"] == "recursive_bch_hall_lyndon_v1"
    assert primary == independent
    assert primary["payload_sha256"] == builder.payload_digest(primary)


def test_source_closure_is_complete_recursive_and_exact(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expected = {
        "src/trottercert/pf4_bch_mapping.py",
        "src/trottercert/cubic_field.py",
        "src/trottercert/intervals.py",
        "scripts/derive_pf4_bch_mapping.py",
        "scripts/reference_pf4_bch_mapping.py",
        "pyproject.toml",
        "requirements-reproducibility.txt",
    }
    assert set(builder.source_closure()) == expected
    assert builder.implementation_sources() == {
        relative: builder.hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in expected
    }

    monkeypatch.setattr(
        builder,
        "SOURCE_ROOTS",
        tuple(
            path
            for path in builder.SOURCE_ROOTS
            if path != "src/trottercert/intervals.py"
        ),
    )
    assert "src/trottercert/intervals.py" in builder.source_closure()


def test_local_import_resolution_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    original = builder._module_source

    def hide_interval(module: str) -> str | None:
        if module == "trottercert.intervals":
            return None
        return original(module)

    monkeypatch.setattr(builder, "_module_source", hide_interval)
    with pytest.raises(ValueError, match="unresolved local import"):
        builder._local_imports("src/trottercert/cubic_field.py")


def test_frozen_artifact_is_canonical_current_and_independently_valid(
    reference: ModuleType,
) -> None:
    raw = ARTIFACT.read_bytes()
    payload = builder.load_payload(ARTIFACT)

    assert raw == builder.canonical_bytes(payload)
    assert payload == builder.build_payload()
    builder.verify_payload(payload)
    assert reference.verify_artifact_payload(payload) == []


def test_cli_fresh_build_is_byte_identical_and_both_verifiers_accept(
    tmp_path: Path,
) -> None:
    rebuilt = tmp_path / "pf4-bch-mapping.json"
    built = _run(BUILDER, "--build", str(rebuilt))
    primary_verified = _run(BUILDER, "--verify", str(rebuilt))
    reference_verified = _run(REFERENCE, "--verify", str(rebuilt))

    assert built.returncode == 0, built.stderr
    assert primary_verified.returncode == 0, primary_verified.stderr
    assert reference_verified.returncode == 0, reference_verified.stdout
    assert rebuilt.read_bytes() == ARTIFACT.read_bytes()


@pytest.mark.parametrize("serialization", ("pretty", "duplicate"))
def test_both_clis_reject_noncanonical_or_duplicate_json(
    serialization: str,
    tmp_path: Path,
) -> None:
    payload = builder.build_payload()
    if serialization == "pretty":
        raw = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("ascii")
        expected = "canonical"
    else:
        raw = b'{"schema_version":2,' + builder.canonical_bytes(payload)[1:]
        expected = "duplicate JSON key"
    artifact = tmp_path / f"{serialization}.json"
    artifact.write_bytes(raw)

    for script in (BUILDER, REFERENCE):
        completed = _run(script, "--verify", str(artifact))
        assert completed.returncode != 0
        assert expected in completed.stdout + completed.stderr


def test_digest_resealed_mapping_mutation_is_rejected_by_both(
    tmp_path: Path,
) -> None:
    forged = copy.deepcopy(builder.build_payload())
    forged["mapping"]["solved_coefficients"]["g_cd"][0][0] += 1
    forged["mapping"]["payload_sha256"] = mapping_digest(forged["mapping"])
    _reseal(forged)
    artifact = tmp_path / "resealed-mapping.json"
    _write(artifact, forged)

    for script in (BUILDER, REFERENCE):
        completed = _run(script, "--verify", str(artifact))
        assert completed.returncode != 0
        assert "authoritative" in completed.stdout + completed.stderr or (
            "semantic mismatch" in completed.stdout + completed.stderr
        )


@pytest.mark.parametrize("mutation", ("hash", "missing_transitive"))
def test_digest_resealed_source_manifest_mutation_is_rejected_by_both(
    mutation: str,
    tmp_path: Path,
) -> None:
    forged = copy.deepcopy(builder.build_payload())
    sources = forged["implementation_sources"]
    if mutation == "hash":
        sources["scripts/reference_pf4_bch_mapping.py"] = "0" * 64
    else:
        del sources["src/trottercert/intervals.py"]
    _reseal(forged)
    artifact = tmp_path / f"resealed-source-{mutation}.json"
    _write(artifact, forged)

    for script in (BUILDER, REFERENCE):
        completed = _run(script, "--verify", str(artifact))
        assert completed.returncode != 0
        assert "source closure mismatch" in completed.stdout + completed.stderr


def test_live_source_drift_is_detected_even_when_artifact_is_unchanged(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    payload = builder.build_payload()
    original = builder._checked_path
    drifted = tmp_path / "reference_pf4_bch_mapping.py"
    drifted.write_bytes(
        original("scripts/reference_pf4_bch_mapping.py").read_bytes()
        + b"\n# simulated source drift\n"
    )

    def drifted_source(relative: str) -> Path:
        if relative == "scripts/reference_pf4_bch_mapping.py":
            return drifted
        return original(relative)

    monkeypatch.setattr(builder, "_checked_path", drifted_source)
    with pytest.raises(ValueError, match="source closure mismatch"):
        builder.verify_payload(payload)
