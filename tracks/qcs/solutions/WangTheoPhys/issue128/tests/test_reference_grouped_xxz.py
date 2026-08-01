from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from itertools import islice
from pathlib import Path

import pytest

import trottercert.grouped_xxz_compressed as compressed_module
from scripts.reference_verify_grouped_xxz import (
    SOURCE_PATHS,
    _deterministic_gzip,
    _payload_digest,
    verify_artifact_bytes,
)
from trottercert.grouped_xxz import XXZCompileSpec, canonical_json_bytes
from trottercert.grouped_xxz_artifact import (
    build_per_delta_artifact,
    source_closure,
)
from trottercert.grouped_xxz_compressed import (
    build_compressed_xxz_ledger,
    compile_compressed_grouped_xxz,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/reference_verify_grouped_xxz.py"


@pytest.fixture(scope="module")
def two_record_artifact():
    monkeypatch = pytest.MonkeyPatch()
    records = tuple(islice(compressed_module._iter_raw_theorem_records(), 2))
    monkeypatch.setattr(
        compressed_module,
        "_projected_raw_record_count",
        lambda: 2,
    )
    monkeypatch.setattr(
        compressed_module,
        "_iter_raw_theorem_records",
        lambda: iter(records),
    )
    try:
        spec = XXZCompileSpec.pilot(Fraction(1, 2))
        ledger = build_compressed_xxz_ledger(spec)
        certificate = compile_compressed_grouped_xxz(spec, ledger=ledger)
        return build_per_delta_artifact(
            ledger,
            certificate,
            source_closure(ROOT),
        )
    finally:
        monkeypatch.undo()


def _seal(payload: dict[str, object]) -> dict[str, object]:
    payload["payload_sha256"] = ""
    payload["payload_sha256"] = _payload_digest(payload)
    return payload


def _rebind(
    summary: dict[str, object],
    witness: dict[str, object],
) -> tuple[bytes, bytes]:
    witness_payload = canonical_json_bytes(_seal(witness))
    witness_gzip = _deterministic_gzip(witness_payload)
    binding = summary["witness"]
    assert isinstance(binding, dict)
    binding["file_sha256"] = hashlib.sha256(witness_gzip).hexdigest()
    binding["payload_sha256"] = hashlib.sha256(witness_payload).hexdigest()
    return canonical_json_bytes(_seal(summary)), witness_gzip


def test_reference_verifier_is_standard_library_plus_independent_oracle() -> None:
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".", 1)[0])
    assert imports <= {
        "__future__",
        "argparse",
        "collections",
        "fractions",
        "gzip",
        "hashlib",
        "io",
        "json",
        "math",
        "pathlib",
        "reference_grouped_xxz_theorem",
        "sys",
    }
    assert "trottercert" not in imports


def test_two_record_builder_artifact_passes_independent_replay(
    two_record_artifact,
) -> None:
    summary, witness = verify_artifact_bytes(
        two_record_artifact.summary_bytes,
        two_record_artifact.witness_gzip_bytes,
        root=ROOT,
        expected_raw_records=2,
    )
    assert summary == two_record_artifact.summary
    assert witness == two_record_artifact.witness
    assert tuple(row["path"] for row in summary["source_closure"]["files"]) == (
        SOURCE_PATHS
    )


def test_clean_shell_cli_rejects_small_fixture_even_with_old_override_variable(
    two_record_artifact,
    tmp_path: Path,
) -> None:
    summary = tmp_path / "summary.json"
    witness = tmp_path / "witness.json.gz"
    summary.write_bytes(two_record_artifact.summary_bytes)
    witness.write_bytes(two_record_artifact.witness_gzip_bytes)
    environment = {
        "GROUPED_XXZ_REFERENCE_TEST_RECORDS": "2",
        "PATH": "",
        "PYTHONHASHSEED": "0",
    }
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), str(summary), str(witness)],
        cwd=tmp_path,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0
    assert "ledger theorem metadata or full coverage mismatch" in completed.stderr


@pytest.mark.parametrize("digest_field", ("file_sha256", "payload_sha256"))
def test_summary_witness_dual_hash_binding_is_enforced(
    two_record_artifact,
    digest_field: str,
) -> None:
    summary = copy.deepcopy(two_record_artifact.summary)
    summary["witness"][digest_field] = "0" * 64
    with pytest.raises(ValueError, match="witness .* hash mismatch"):
        verify_artifact_bytes(
            canonical_json_bytes(_seal(summary)),
            two_record_artifact.witness_gzip_bytes,
            root=ROOT,
            expected_raw_records=2,
        )


def test_parser_rejects_noncanonical_duplicate_bool_and_float(
    two_record_artifact,
) -> None:
    summary = copy.deepcopy(two_record_artifact.summary)
    summary["counts"]["raw_records"] = True
    with pytest.raises((TypeError, ValueError)):
        verify_artifact_bytes(
            canonical_json_bytes(_seal(summary)),
            two_record_artifact.witness_gzip_bytes,
            root=ROOT,
            expected_raw_records=2,
        )
    summary = copy.deepcopy(two_record_artifact.summary)
    summary["counts"]["raw_records"] = 2.0
    raw = (json.dumps(_seal(summary), sort_keys=True, separators=(",", ":")) + "\n").encode()
    with pytest.raises((TypeError, ValueError)):
        verify_artifact_bytes(
            raw,
            two_record_artifact.witness_gzip_bytes,
            root=ROOT,
            expected_raw_records=2,
        )
    duplicate = b'{"schema":"shadow",' + two_record_artifact.summary_bytes[1:]
    with pytest.raises(ValueError, match="duplicate"):
        verify_artifact_bytes(
            duplicate,
            two_record_artifact.witness_gzip_bytes,
            root=ROOT,
            expected_raw_records=2,
        )
    pretty = (json.dumps(two_record_artifact.summary, indent=2) + "\n").encode()
    with pytest.raises(ValueError, match="canonical"):
        verify_artifact_bytes(
            pretty,
            two_record_artifact.witness_gzip_bytes,
            root=ROOT,
            expected_raw_records=2,
        )


def test_gzip_must_be_one_deterministic_member(two_record_artifact) -> None:
    doubled = two_record_artifact.witness_gzip_bytes * 2
    summary = copy.deepcopy(two_record_artifact.summary)
    binding = summary["witness"]
    binding["file_sha256"] = hashlib.sha256(doubled).hexdigest()
    binding["payload_sha256"] = hashlib.sha256(
        gzip.decompress(doubled)
    ).hexdigest()
    with pytest.raises(ValueError, match="deterministic gzip"):
        verify_artifact_bytes(
            canonical_json_bytes(_seal(summary)),
            doubled,
            root=ROOT,
            expected_raw_records=2,
        )


def _closure_reseal(closure: dict[str, object]) -> None:
    unsigned = {
        "algorithm": closure["algorithm"],
        "files": closure["files"],
        "source_commit": closure["source_commit"],
    }
    closure["closure_sha256"] = hashlib.sha256(
        canonical_json_bytes(unsigned)
    ).hexdigest()


@pytest.mark.parametrize(
    "mutation",
    (
        "delta",
        "axis_weight",
        "length",
        "boundary",
        "periodic_bond",
        "stage_order",
        "theorem_word",
        "coefficient_interval",
        "pauli_mask",
        "group_coverage",
        "anticommutation",
        "candidate_step",
        "previous_step_bound",
        "baseline",
        "resource_count",
        "source_hash",
        "payload_digest",
        "unknown_field",
        "noncanonical_fraction",
    ),
)
def test_fully_resealed_semantic_mutation_matrix_is_rejected(
    two_record_artifact,
    mutation: str,
) -> None:
    summary = copy.deepcopy(two_record_artifact.summary)
    witness = copy.deepcopy(two_record_artifact.witness)
    ledger = witness["ledger"]
    certificate = witness["certificate"]
    blocks = witness["representative_blocks"]
    rows = witness["word_weights"]

    if mutation == "delta":
        summary["spec"]["delta"] = [2, 1]
        ledger["spec"]["delta"] = [2, 1]
        certificate["spec"]["delta"] = [2, 1]
    elif mutation == "axis_weight":
        blocks[0]["terms"][0][2][0] += 1
    elif mutation == "length":
        summary["spec"]["length"] = 6
        ledger["spec"]["length"] = 6
        certificate["spec"]["length"] = 6
    elif mutation == "boundary":
        summary["spec"]["boundary"] = "open"
        ledger["spec"]["boundary"] = "open"
        certificate["spec"]["boundary"] = "open"
    elif mutation == "periodic_bond":
        blocks[0]["terms"][0][1] ^= 1 << 15
    elif mutation == "stage_order":
        ledger["raw_stream_sha256"] = "0" * 64
        certificate["raw_stream_sha256"] = "0" * 64
    elif mutation == "theorem_word":
        rows[0]["actual_word"] = [3, 3, 3, 3, 3]
    elif mutation == "coefficient_interval":
        rows[0]["raw_weight_upper"][0] += 1
    elif mutation == "pauli_mask":
        blocks[0]["terms"][0][0] ^= 1
    elif mutation == "group_coverage":
        group = next(
            group
            for block in blocks
            for group in block["groups"]
            if len(group["term_indices"]) == 2
        )
        group["term_indices"].pop()
    elif mutation == "anticommutation":
        block = blocks[0]
        terms = block["terms"]
        commuting = next(
            (left, right)
            for left in range(len(terms))
            for right in range(left + 1, len(terms))
            if (
                (
                    (terms[left][0] & terms[right][1]).bit_count()
                    + (terms[left][1] & terms[right][0]).bit_count()
                )
                & 1
            )
            == 0
        )
        block["groups"][0]["term_indices"] = list(commuting)
    elif mutation == "candidate_step":
        certificate["candidate_steps"] += 1
    elif mutation == "previous_step_bound":
        certificate["candidate_previous_error"][0] += 1
    elif mutation == "baseline":
        certificate["baseline_steps"] += 1
    elif mutation == "resource_count":
        certificate["candidate_resources"] += 30
    elif mutation == "source_hash":
        for closure in (summary["source_closure"], witness["source_closure"]):
            closure["files"][0]["sha256"] = "0" * 64
            _closure_reseal(closure)
    elif mutation == "payload_digest":
        witness["payload_sha256"] = "0" * 64
        witness_payload = canonical_json_bytes(witness)
        witness_gzip = _deterministic_gzip(witness_payload)
        binding = summary["witness"]
        binding["file_sha256"] = hashlib.sha256(witness_gzip).hexdigest()
        binding["payload_sha256"] = hashlib.sha256(witness_payload).hexdigest()
        summary_bytes = canonical_json_bytes(_seal(summary))
        with pytest.raises((TypeError, ValueError)):
            verify_artifact_bytes(
                summary_bytes,
                witness_gzip,
                root=ROOT,
                expected_raw_records=2,
            )
        return
    elif mutation == "unknown_field":
        ledger["unknown"] = 1
    elif mutation == "noncanonical_fraction":
        ledger["grouped_constant"] = [2, 2]
    else:
        raise AssertionError(mutation)

    summary_bytes, witness_gzip = _rebind(summary, witness)
    with pytest.raises((TypeError, ValueError)):
        verify_artifact_bytes(
            summary_bytes,
            witness_gzip,
            root=ROOT,
            expected_raw_records=2,
        )
