from __future__ import annotations

import ast
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import certify_tfim_gauge_witness as builder

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/certify_tfim_gauge_witness.py"
REFERENCE = ROOT / "scripts/reference_tfim_gauge_witness.py"
ARTIFACT = (
    ROOT
    / "docs/experiments/processor-obstruction/tfim-gauge-witness.json"
)


def _run(script: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    clean_environment = dict(os.environ)
    clean_environment.pop("PYTHONPATH", None)
    return subprocess.run(
        [sys.executable, str(script), *arguments],
        cwd=ROOT,
        env=clean_environment,
        check=False,
        capture_output=True,
        text=True,
        timeout=300,
    )


def _write(path: Path, payload: object) -> None:
    path.write_bytes(builder.canonical_bytes(payload))


def _reseal(payload: dict[str, object]) -> None:
    payload["payload_sha256"] = builder.payload_digest(payload)


def test_checked_artifact_is_canonical_source_closed_and_independently_valid() -> None:
    raw = ARTIFACT.read_bytes()
    payload = builder.load_payload(ARTIFACT)
    assert raw == builder.canonical_bytes(payload)
    builder.verify_payload(payload)

    primary = _run(BUILDER, "--verify", str(ARTIFACT))
    reference = _run(REFERENCE, "--verify", str(ARTIFACT))
    assert primary.returncode == 0, primary.stderr
    assert reference.returncode == 0, reference.stderr
    assert payload["checked_lengths"] == [4, 6, 8, 10]
    assert payload["search"]["power_order"] == [2, 3, 4]
    assert payload["claim"]["finite_step_status"] == "not_claimed"
    assert payload["claim"]["quantitative_distance_status"] == "not_claimed"
    assert payload["claim"]["universal_even_l_ge_6_status"] == (
        "analytic_formula_recorded_induction_not_machine_proved"
    )
    assert "src/trottercert/__init__.py" in builder.source_closure()
    assert "src/trottercert/__init__.py" in payload["implementation_sources"]


def test_reference_verifier_is_standard_library_only() -> None:
    tree = ast.parse(REFERENCE.read_text())
    imported_roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_roots.add((node.module or "").split(".", 1)[0])
    assert imported_roots.isdisjoint({"trottercert", "numpy", "sympy", "scripts"})


def test_fresh_rebuild_is_byte_identical(tmp_path: Path) -> None:
    rebuilt = tmp_path / "tfim-gauge-witness.json"
    completed = _run(BUILDER, "--build", "--output", str(rebuilt))
    assert completed.returncode == 0, completed.stderr
    assert rebuilt.read_bytes() == ARTIFACT.read_bytes()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("unknown_top", "field"),
        ("schema_version_bool", "numeric alias"),
        ("mapping_file_hash", "mapping"),
        ("mapping_payload_hash", "mapping"),
        ("source_hash", "source"),
        ("package_init_source_hash", "source"),
        ("search_order", "search"),
        ("duplicate_length", "length"),
        ("orbit_coefficient", "orbit"),
        ("witness_b", "witness"),
        ("witness_pairing", "pairing"),
        ("commutator_flag", "commut"),
        ("universal_overclaim", "universal"),
        ("finite_step_overclaim", "finite-step"),
        ("distance_overclaim", "distance"),
    ],
)
def test_resealed_semantic_mutations_are_rejected(
    tmp_path: Path, mutation: str, message: str
) -> None:
    payload = copy.deepcopy(builder.load_payload(ARTIFACT))
    if mutation == "unknown_top":
        payload["extra"] = True
    elif mutation == "schema_version_bool":
        payload["schema_version"] = True
    elif mutation == "mapping_file_hash":
        payload["pf4_mapping_binding"]["file_sha256"] = "0" * 64
    elif mutation == "mapping_payload_hash":
        payload["pf4_mapping_binding"]["mapping_payload_sha256"] = "0" * 64
    elif mutation == "source_hash":
        source = next(iter(payload["implementation_sources"]))
        payload["implementation_sources"][source] = "0" * 64
    elif mutation == "package_init_source_hash":
        payload["implementation_sources"][
            "src/trottercert/__init__.py"
        ] = "0" * 64
    elif mutation == "search_order":
        payload["search"]["power_order"] = [3, 2, 4]
    elif mutation == "duplicate_length":
        payload["checked_lengths"][-1] = 8
    elif mutation == "orbit_coefficient":
        payload["checked_instances"][1]["orbit_polynomials"]["tau_h4"][0][
            "coefficient"
        ][0] += 1
    elif mutation == "witness_b":
        payload["checked_instances"][1]["witnesses"][1][
            "hamiltonian_coefficient"
        ][0] += 1
    elif mutation == "witness_pairing":
        payload["checked_instances"][1]["witnesses"][1]["tau_w_e5"][0][0] += 1
    elif mutation == "commutator_flag":
        payload["checked_instances"][1]["witnesses"][1][
            "commutator_is_zero"
        ] = False
    elif mutation == "universal_overclaim":
        payload["claim"]["universal_even_l_ge_6_status"] = "machine_proved"
    elif mutation == "finite_step_overclaim":
        payload["claim"]["finite_step_status"] = "certified"
    elif mutation == "distance_overclaim":
        payload["claim"]["quantitative_distance_status"] = "certified"
    else:
        raise AssertionError(mutation)
    _reseal(payload)
    forged = tmp_path / f"{mutation}.json"
    _write(forged, payload)

    completed = _run(REFERENCE, "--verify", str(forged))
    assert completed.returncode != 0
    assert message in completed.stderr.lower()


@pytest.mark.parametrize("serialization", ["pretty", "duplicate"])
def test_both_verifiers_reject_noncanonical_or_duplicate_json(
    tmp_path: Path, serialization: str
) -> None:
    payload = builder.load_payload(ARTIFACT)
    path = tmp_path / f"{serialization}.json"
    if serialization == "pretty":
        path.write_text(json.dumps(payload, indent=2) + "\n")
        message = "canonical"
    else:
        canonical = builder.canonical_bytes(payload)
        path.write_bytes(b'{"schema_version":1,' + canonical[1:])
        message = "duplicate"
    for script, arguments in (
        (BUILDER, ("--verify", str(path))),
        (REFERENCE, ("--verify", str(path))),
    ):
        completed = _run(script, *arguments)
        assert completed.returncode != 0
        assert message in completed.stderr.lower()
