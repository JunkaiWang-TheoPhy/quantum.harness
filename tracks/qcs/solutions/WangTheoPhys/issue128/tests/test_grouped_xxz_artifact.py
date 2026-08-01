from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import sys
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import islice
from pathlib import Path

import pytest

import scripts.compile_grouped_xxz as cli
import trottercert.grouped_xxz_artifact as artifact_module
import trottercert.grouped_xxz_compressed as compressed_module
from trottercert.grouped_xxz import XXZCompileSpec, canonical_json_bytes
from trottercert.grouped_xxz_artifact import (
    FROZEN_LOCAL_SOURCE_ALLOWLIST,
    SOURCE_CLOSURE_PATHS,
    build_merged_pilot_artifact,
    build_per_delta_artifact,
    derive_local_source_closure,
    deterministic_gzip_bytes,
    payload_digest,
    source_closure,
    verify_merged_pilot_artifact,
    verify_per_delta_artifact,
    write_merged_pilot_artifact,
    write_per_delta_artifact,
)
from trottercert.grouped_xxz_compressed import (
    build_compressed_xxz_ledger,
    compile_compressed_grouped_xxz,
)

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "compile_grouped_xxz.py"
REQUIRED_CLOSURE_ADDITIONS = (
    "src/trottercert/__init__.py",
    "src/trottercert/algebra.py",
    "src/trottercert/local_commutators.py",
    "pyproject.toml",
    "requirements-reproducibility.txt",
)


def _two_record_certificate(monkeypatch, delta: Fraction = Fraction(1, 2)):
    records = tuple(islice(compressed_module._iter_raw_theorem_records(), 2))
    monkeypatch.setattr(
        compressed_module,
        "_projected_raw_record_count",
        lambda: len(records),
    )
    monkeypatch.setattr(
        compressed_module,
        "_iter_raw_theorem_records",
        lambda: iter(records),
    )
    spec = XXZCompileSpec.pilot(delta)
    ledger = build_compressed_xxz_ledger(spec)
    certificate = compile_compressed_grouped_xxz(spec, ledger=ledger)
    return ledger, certificate


def _reseal(payload: dict[str, object]) -> dict[str, object]:
    payload["payload_sha256"] = ""
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def _rebind_summary(
    summary: dict[str, object],
    witness: dict[str, object],
) -> tuple[bytes, bytes]:
    witness = _reseal(witness)
    witness_payload = canonical_json_bytes(witness)
    witness_gzip = deterministic_gzip_bytes(witness_payload)
    binding = summary["witness"]
    assert isinstance(binding, dict)
    binding["file_sha256"] = sha256(witness_gzip).hexdigest()
    binding["payload_sha256"] = sha256(witness_payload).hexdigest()
    summary = _reseal(summary)
    return canonical_json_bytes(summary), witness_gzip


def _copy_source_closure(tmp_path: Path) -> None:
    for relative in SOURCE_CLOSURE_PATHS:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, target)


def _two_point_artifacts(monkeypatch, *, source_commit="fixture-commit"):
    closure = source_closure(ROOT, source_commit=source_commit)
    artifacts = []
    for delta in (Fraction(1, 2), Fraction(2)):
        ledger, certificate = _two_record_certificate(monkeypatch, delta)
        artifacts.append(build_per_delta_artifact(ledger, certificate, closure))
    return tuple(artifacts), closure


def test_recursive_ast_source_closure_matches_the_independent_allowlist(
    tmp_path: Path,
) -> None:
    _copy_source_closure(tmp_path)
    derived = derive_local_source_closure(tmp_path)
    assert derived == FROZEN_LOCAL_SOURCE_ALLOWLIST
    assert set(REQUIRED_CLOSURE_ADDITIONS) <= set(SOURCE_CLOSURE_PATHS)
    assert "src/trottercert/__init__.py" in derived

    # Exercise ``from package import submodule`` explicitly.  The child is
    # already in the true closure, so the derived/frozen sets stay identical.
    script = tmp_path / "scripts" / "compile_grouped_xxz.py"
    script.write_text(
        script.read_text() + "\nfrom trottercert import grouped_xxz_compressed\n"
    )
    assert derive_local_source_closure(tmp_path) == FROZEN_LOCAL_SOURCE_ALLOWLIST


@pytest.mark.parametrize("relative", REQUIRED_CLOSURE_ADDITIONS)
def test_source_closure_digest_detects_required_file_drift(
    tmp_path: Path,
    relative: str,
) -> None:
    _copy_source_closure(tmp_path)
    baseline = source_closure(tmp_path, source_commit="fixture-commit")
    target = tmp_path / relative
    target.write_bytes(target.read_bytes() + b"\n# injected source drift\n")
    observed = source_closure(tmp_path, source_commit="fixture-commit")
    assert observed["closure_sha256"] != baseline["closure_sha256"]
    baseline_files = {row["path"]: row["sha256"] for row in baseline["files"]}
    observed_files = {row["path"]: row["sha256"] for row in observed["files"]}
    assert observed_files[relative] != baseline_files[relative]


def test_two_record_artifact_is_canonical_deterministic_and_semantic(
    monkeypatch,
) -> None:
    ledger, certificate = _two_record_certificate(monkeypatch)
    closure = source_closure(ROOT, source_commit="fixture-commit")
    first = build_per_delta_artifact(ledger, certificate, closure)
    second = build_per_delta_artifact(ledger, certificate, closure)

    assert first == second
    assert first.summary_bytes == canonical_json_bytes(first.summary)
    assert first.witness_payload_bytes == canonical_json_bytes(first.witness)
    assert gzip.decompress(first.witness_gzip_bytes) == first.witness_payload_bytes
    assert first.witness_gzip_bytes[3] & 0x08 == 0
    assert first.witness_gzip_bytes[4:8] == b"\x00\x00\x00\x00"

    summary = first.summary
    binding = summary["witness"]
    assert binding["file_sha256"] == sha256(first.witness_gzip_bytes).hexdigest()
    assert binding["payload_sha256"] == sha256(first.witness_payload_bytes).hexdigest()
    assert summary["spec"]["delta"] == [1, 2]
    assert summary["method"] == certificate.method
    assert summary["status"] == "certified"
    assert summary["counts"] == {
        "actual_words": len(ledger.word_weights),
        "pair_groups": sum(
            len(group.terms) == 2
            for block in ledger.representative_blocks
            for group in block.groups
        ),
        "projected_raw_records": 2,
        "raw_records": 2,
        "representative_blocks": len(ledger.representative_blocks),
        "representative_terms": sum(
            len(block.terms) for block in ledger.representative_blocks
        ),
        "singleton_groups": sum(
            len(group.terms) == 1
            for block in ledger.representative_blocks
            for group in block.groups
        ),
    }
    assert "elapsed_seconds" not in first.summary_bytes.decode()
    assert "peak_rss" not in first.witness_payload_bytes.decode()

    blocks = first.witness["representative_blocks"]
    assert blocks
    for block in blocks:
        assert "terms" in block
        for group in block["groups"]:
            assert set(group) == {
                "norm_interval",
                "squared_norm",
                "term_indices",
            }
            assert all(isinstance(index, int) for index in group["term_indices"])
    observed_ledger, observed_certificate = verify_per_delta_artifact(
        first.summary_bytes,
        first.witness_gzip_bytes,
        expected_source_closure=closure,
    )
    assert observed_ledger == ledger
    assert observed_certificate == certificate
    with pytest.raises(TypeError, match="expected_source_closure"):
        verify_per_delta_artifact(  # type: ignore[call-arg]
            first.summary_bytes,
            first.witness_gzip_bytes,
        )


@pytest.mark.parametrize(
    "mutation",
    ["summary_bound", "witness_hash", "term", "group_index", "source_closure"],
)
def test_primary_artifact_verifier_rejects_fully_resealed_mutations(
    monkeypatch,
    mutation: str,
) -> None:
    ledger, certificate = _two_record_certificate(monkeypatch, Fraction(2))
    closure = source_closure(ROOT, source_commit="fixture-commit")
    bundle = build_per_delta_artifact(ledger, certificate, closure)
    summary = deepcopy(bundle.summary)
    witness = deepcopy(bundle.witness)

    if mutation == "summary_bound":
        summary["bounds"]["grouped_constant"][0] += 1
        summary_bytes = canonical_json_bytes(_reseal(summary))
        witness_gzip = bundle.witness_gzip_bytes
    elif mutation == "witness_hash":
        summary["witness"]["file_sha256"] = "0" * 64
        summary_bytes = canonical_json_bytes(_reseal(summary))
        witness_gzip = bundle.witness_gzip_bytes
    elif mutation == "term":
        witness["representative_blocks"][0]["terms"][0][2][0] += 1
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    elif mutation == "group_index":
        witness["representative_blocks"][0]["groups"][0]["term_indices"][0] = 10**9
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    elif mutation == "source_closure":
        witness["source_closure"]["source_commit"] = "forged"
        summary["source_closure"]["source_commit"] = "forged"
        for submitted in (witness["source_closure"], summary["source_closure"]):
            unsigned = {
                "algorithm": submitted["algorithm"],
                "files": submitted["files"],
                "source_commit": submitted["source_commit"],
            }
            submitted["closure_sha256"] = sha256(
                canonical_json_bytes(unsigned)
            ).hexdigest()
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    else:
        raise AssertionError(mutation)

    with pytest.raises(ValueError):
        verify_per_delta_artifact(
            summary_bytes,
            witness_gzip,
            expected_source_closure=closure,
        )


def test_atomic_writer_refuses_overwrite(monkeypatch, tmp_path: Path) -> None:
    ledger, certificate = _two_record_certificate(monkeypatch)
    bundle = build_per_delta_artifact(
        ledger,
        certificate,
        source_closure(ROOT, source_commit="fixture-commit"),
    )
    summary = tmp_path / "result.json"
    witness = tmp_path / "result-witness.json.gz"
    write_per_delta_artifact(bundle, summary, witness)
    assert summary.read_bytes() == bundle.summary_bytes
    assert witness.read_bytes() == bundle.witness_gzip_bytes
    with pytest.raises(FileExistsError):
        write_per_delta_artifact(bundle, summary, witness)


def test_atomic_writer_rolls_back_witness_when_summary_link_fails(
    monkeypatch,
    tmp_path: Path,
) -> None:
    ledger, certificate = _two_record_certificate(monkeypatch)
    bundle = build_per_delta_artifact(
        ledger,
        certificate,
        source_closure(ROOT, source_commit="fixture-commit"),
    )
    summary = tmp_path / "result.json"
    witness = tmp_path / "result-witness.json.gz"
    real_link = artifact_module.os.link
    calls = 0

    def fail_second_link(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected summary link failure")
        return real_link(source, target)

    monkeypatch.setattr(artifact_module.os, "link", fail_second_link)
    with pytest.raises(OSError, match="summary link failure"):
        write_per_delta_artifact(bundle, summary, witness)
    assert calls == 2
    assert not summary.exists()
    assert not witness.exists()
    assert not tuple(tmp_path.glob(".*.tmp"))


def test_cli_bootstraps_src_and_runs_a_real_compressed_profile(
    tmp_path: Path,
) -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--profile-only",
            "--delta",
            "1/2",
            "--pipeline",
            "direct-theorem",
            "--max-records",
            "2",
        ],
        cwd=tmp_path,
        env={"PATH": ""},
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["mode"] == "profile-only"
    assert payload["raw_records"] == 2
    assert payload["complete"] is False
    assert payload["active_words"] == 2
    assert payload["representative_blocks"] > 0


def test_cli_two_record_build_self_verifies_and_rejects_unsafe_modes(
    monkeypatch,
    tmp_path: Path,
) -> None:
    _two_record_certificate(monkeypatch)
    summary = tmp_path / "summary.json"
    witness = tmp_path / "witness.json.gz"
    arguments = [
        "--build",
        "--delta",
        "1/2",
        "--pipeline",
        "direct-theorem",
        "--summary",
        str(summary),
        "--witness",
        str(witness),
    ]
    assert cli.main(arguments) == 0
    verify_per_delta_artifact(
        summary.read_bytes(),
        witness.read_bytes(),
        expected_source_closure=source_closure(ROOT),
    )
    with pytest.raises(SystemExit):
        cli.main(arguments)
    with pytest.raises(SystemExit):
        cli.main(arguments + ["--max-records", "2"])


def test_merged_pilot_is_canonical_fresh_and_positive(monkeypatch) -> None:
    artifacts, closure = _two_point_artifacts(monkeypatch)
    left, right = artifacts
    merged = build_merged_pilot_artifact(
        left.summary_bytes,
        left.witness_gzip_bytes,
        right.summary_bytes,
        right.witness_gzip_bytes,
        closure,
    )
    reversed_inputs = build_merged_pilot_artifact(
        right.summary_bytes,
        right.witness_gzip_bytes,
        left.summary_bytes,
        left.witness_gzip_bytes,
        closure,
    )
    assert merged == reversed_inputs
    assert merged.summary_bytes == canonical_json_bytes(merged.summary)
    assert gzip.decompress(merged.witness_gzip_bytes) == merged.witness_payload_bytes
    assert [row["delta"] for row in merged.summary["rows"]] == [[1, 2], [2, 1]]
    assert len(merged.witness["per_delta_artifacts"]) == 2
    for row in merged.summary["rows"]:
        assert row["candidate"]["resources"] < row["baseline"]["resources"]
    verified = verify_merged_pilot_artifact(
        merged.summary_bytes,
        merged.witness_gzip_bytes,
        expected_source_closure=closure,
    )
    assert tuple(certificate.spec.delta for _, certificate in verified) == (
        Fraction(1, 2),
        Fraction(2),
    )
    with pytest.raises(ValueError, match="duplicate|Delta|two-point"):
        build_merged_pilot_artifact(
            left.summary_bytes,
            left.witness_gzip_bytes,
            left.summary_bytes,
            left.witness_gzip_bytes,
            closure,
        )


@pytest.mark.parametrize(
    "mutation",
    ["binding", "missing_row", "row_order", "embedded_summary", "source_closure"],
)
def test_merged_primary_verifier_rejects_resealed_mutations(
    monkeypatch,
    mutation: str,
) -> None:
    artifacts, closure = _two_point_artifacts(monkeypatch)
    merged = build_merged_pilot_artifact(
        artifacts[0].summary_bytes,
        artifacts[0].witness_gzip_bytes,
        artifacts[1].summary_bytes,
        artifacts[1].witness_gzip_bytes,
        closure,
    )
    summary = deepcopy(merged.summary)
    witness = deepcopy(merged.witness)
    if mutation == "binding":
        summary["witness"]["file_sha256"] = "0" * 64
        summary_bytes = canonical_json_bytes(_reseal(summary))
        witness_gzip = merged.witness_gzip_bytes
    elif mutation == "missing_row":
        witness["per_delta_artifacts"].pop()
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    elif mutation == "row_order":
        witness["per_delta_artifacts"].reverse()
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    elif mutation == "embedded_summary":
        embedded = witness["per_delta_artifacts"][0]["summary"]
        embedded["resources"]["candidate"] += 30
        _reseal(embedded)
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    elif mutation == "source_closure":
        witness["source_closure"]["source_commit"] = "forged"
        summary["source_closure"]["source_commit"] = "forged"
        for submitted in (witness["source_closure"], summary["source_closure"]):
            unsigned = {
                "algorithm": submitted["algorithm"],
                "files": submitted["files"],
                "source_commit": submitted["source_commit"],
            }
            submitted["closure_sha256"] = sha256(
                canonical_json_bytes(unsigned)
            ).hexdigest()
        summary_bytes, witness_gzip = _rebind_summary(summary, witness)
    else:
        raise AssertionError(mutation)
    with pytest.raises(ValueError):
        verify_merged_pilot_artifact(
            summary_bytes,
            witness_gzip,
            expected_source_closure=closure,
        )


def test_merge_cli_verifies_inputs_and_atomically_refuses_overwrite(
    monkeypatch,
    tmp_path: Path,
) -> None:
    artifacts, closure = _two_point_artifacts(monkeypatch, source_commit=None)
    inputs = []
    for index, artifact in enumerate(artifacts):
        summary = tmp_path / f"input-{index}.json"
        witness = tmp_path / f"input-{index}.json.gz"
        write_per_delta_artifact(artifact, summary, witness)
        inputs.extend((str(summary), str(witness)))
    output_summary = tmp_path / "pilot.json"
    output_witness = tmp_path / "pilot-witness.json.gz"
    arguments = [
        "--merge-pilot",
        *inputs,
        "--summary",
        str(output_summary),
        "--witness",
        str(output_witness),
    ]
    assert cli.main(arguments) == 0
    verify_merged_pilot_artifact(
        output_summary.read_bytes(),
        output_witness.read_bytes(),
        expected_source_closure=closure,
    )
    with pytest.raises(SystemExit):
        cli.main(arguments)
    with pytest.raises(SystemExit):
        cli.main(arguments + ["--delta", "1/2"])

    # The direct writer uses the same summary-last/no-overwrite primitive.
    merged = build_merged_pilot_artifact(
        artifacts[0].summary_bytes,
        artifacts[0].witness_gzip_bytes,
        artifacts[1].summary_bytes,
        artifacts[1].witness_gzip_bytes,
        closure,
    )
    with pytest.raises(FileExistsError):
        write_merged_pilot_artifact(merged, output_summary, output_witness)
