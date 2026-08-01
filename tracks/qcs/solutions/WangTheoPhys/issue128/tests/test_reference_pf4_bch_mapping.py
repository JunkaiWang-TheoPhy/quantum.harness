from __future__ import annotations

import ast
import copy
import importlib.util
import json
import subprocess
import sys
from collections.abc import Callable
from fractions import Fraction
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/reference_pf4_bch_mapping.py"


def _load_reference() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "reference_pf4_bch_mapping_for_tests",
        SCRIPT,
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def reference() -> ModuleType:
    return _load_reference()


@pytest.fixture(scope="module")
def canonical_payload(reference: ModuleType) -> dict[str, object]:
    return reference.build_reference_artifact()


def _write_payload(reference: ModuleType, path: Path, payload: object) -> None:
    path.write_bytes(reference._canonical_json_bytes(payload))


def _reseal(reference: ModuleType, payload: dict[str, object]) -> None:
    payload["payload_sha256"] = reference._payload_digest(payload)


def _run(path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--verify", str(path)],
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


def test_reference_source_is_standard_library_and_algorithmically_separate() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported_roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported_roots.add(node.module.split(".", 1)[0])

    assert not imported_roots.intersection(
        {"trottercert", "numpy", "sympy", "derive_pf4_bch_mapping"}
    )
    assert "cubic_formula_log_series" not in source
    assert "recursive_bch_hall_lyndon_v1" in source


def test_reference_uses_the_frozen_literal_eleven_stages(reference: ModuleType) -> None:
    u = reference.Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
    v = reference.ONE - 4 * u
    assert reference._literal_stages() == (
        (0, u / 2),
        (1, u),
        (0, u),
        (1, u),
        (0, (u + v) / 2),
        (1, v),
        (0, (u + v) / 2),
        (1, u),
        (0, u),
        (1, u),
        (0, u / 2),
    )


def test_recursive_bch_has_the_universal_low_degree_terms(reference: ModuleType) -> None:
    a = {(0,): reference.ONE}
    b = {(1,): reference.ONE}
    logarithm = reference._bch(a, b)
    words = reference._expand_lie(logarithm)

    assert {word: value for word, value in words.items() if len(word) == 1} == {
        (0,): reference.ONE,
        (1,): reference.ONE,
    }
    assert {word: value for word, value in words.items() if len(word) == 2} == {
        (0, 1): reference.Cubic(Fraction(1, 2)),
        (1, 0): reference.Cubic(Fraction(-1, 2)),
    }


def test_reference_derives_order_conditions_and_general_cyclic_solution(
    canonical_payload: dict[str, object],
) -> None:
    log_terms = canonical_payload["log_terms"]
    solved = canonical_payload["solved_coefficients"]
    projective = canonical_payload["projective_comparisons"]

    assert [entry["degree"] for entry in log_terms] == list(range(6))
    assert len(log_terms[1]["terms"]) == 2
    assert log_terms[2]["terms"] == []
    assert log_terms[3]["terms"] == []
    assert log_terms[4]["terms"] == []
    assert len(log_terms[5]["terms"]) == 30
    assert len(canonical_payload["cyclic_map"]["target_h_l5"]) == 10
    assert canonical_payload["cyclic_map"]["target_h_l5"] == (
        canonical_payload["cyclic_map"]["reconstructed"]
    )
    assert canonical_payload["residual"] == []
    assert solved == {
        "g_c2": [[37, 900000], [313, 14400000], [29, 1800000]],
        "g_cd": [[-37, 225000], [-313, 3600000], [-29, 450000]],
        "g_d2": [[37, 337500], [313, 5400000], [29, 675000]],
    }
    assert projective["historical_candidate"]["status"] == "proportional"
    assert projective["historical_candidate"]["residual"] == [
        [[0, 1], [0, 1], [0, 1]]
    ] * 3
    assert projective["legacy"]["status"] == "not_proportional"
    assert any(
        value != [[0, 1], [0, 1], [0, 1]]
        for value in projective["legacy"]["residual"]
    )
    assert canonical_payload["physical_bridge"]["relation"] == "E5=L5"
    assert canonical_payload["conventions"]["convention_digest"] == (
        "1ef04d6fcaf15a4350fbebd1528b3e84426a21bd946e5fe9ee3e990b470d95ac"
    )


def test_reference_cli_accepts_its_canonical_fixture(
    reference: ModuleType,
    canonical_payload: dict[str, object],
    tmp_path: Path,
) -> None:
    artifact = tmp_path / "mapping.json"
    _write_payload(reference, artifact, canonical_payload)
    completed = _run(artifact)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    output = json.loads(completed.stdout)
    assert output["valid"] is True
    assert output["algorithm"] == "recursive_bch_hall_lyndon_v1"


Mutation = Callable[[dict[str, object]], None]


def _mutate_stage(payload: dict[str, object]) -> None:
    payload["stages"][0]["coefficient"][0][0] += 1


def _mutate_log_coordinate(payload: dict[str, object]) -> None:
    payload["log_terms"][5]["terms"][0]["coefficient"][1][0] += 1


def _mutate_cyclic_word(payload: dict[str, object]) -> None:
    payload["cyclic_map"]["target_h_l5"][0]["word"] = [1] * 6


def _mutate_solved_coefficient(payload: dict[str, object]) -> None:
    payload["solved_coefficients"]["g_cd"][2][0] += 1


def _mutate_residual(payload: dict[str, object]) -> None:
    payload["residual"] = [
        {"word": [0, 0, 0, 0, 1, 1], "coefficient": [[1, 1], [0, 1], [0, 1]]}
    ]


def _mutate_root(payload: dict[str, object]) -> None:
    payload["sign_certificate"]["alpha_interval"][0][0] += 1


def _mutate_physical_bridge(payload: dict[str, object]) -> None:
    payload["physical_bridge"]["relation"] = "E5=-L5"


def _mutate_identity_status(payload: dict[str, object]) -> None:
    payload["identity_status"] = "proved_by_submitted_string"


@pytest.mark.parametrize(
    "mutation",
    (
        _mutate_stage,
        _mutate_log_coordinate,
        _mutate_cyclic_word,
        _mutate_solved_coefficient,
        _mutate_residual,
        _mutate_root,
        _mutate_physical_bridge,
        _mutate_identity_status,
    ),
)
def test_reference_cli_rejects_digest_resealed_semantic_mutations(
    mutation: Mutation,
    reference: ModuleType,
    canonical_payload: dict[str, object],
    tmp_path: Path,
) -> None:
    payload = copy.deepcopy(canonical_payload)
    mutation(payload)
    _reseal(reference, payload)
    artifact = tmp_path / f"{mutation.__name__}.json"
    _write_payload(reference, artifact, payload)

    completed = _run(artifact)

    assert completed.returncode != 0
    assert "semantic mismatch" in completed.stdout


def test_reference_cli_rejects_duplicate_keys(
    canonical_payload: dict[str, object], tmp_path: Path
) -> None:
    artifact = tmp_path / "duplicate.json"
    encoded = json.dumps(canonical_payload, sort_keys=True, separators=(",", ":"))
    artifact.write_text(
        encoded[:-1] + ',"schema_version":1}\n',
        encoding="ascii",
    )

    completed = _run(artifact)

    assert completed.returncode != 0
    assert "duplicate JSON key" in completed.stdout


def test_reference_cli_rejects_noncanonical_json(
    canonical_payload: dict[str, object], tmp_path: Path
) -> None:
    artifact = tmp_path / "pretty.json"
    artifact.write_text(json.dumps(canonical_payload, indent=2) + "\n", encoding="ascii")

    completed = _run(artifact)

    assert completed.returncode != 0
    assert "not canonical JSON" in completed.stdout
