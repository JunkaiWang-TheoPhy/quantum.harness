from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_dual_e9_pairing import (
    E9_CONTRACTION_LENGTH,
    build_worker_payload,
    reduce_worker_payloads,
)
from scripts.reduce_dual_e9_run import audit_e9_run, recommended_rerun_shards
from trottercert.cubic_field import Cubic
from trottercert.dual_manifest_pairing import ManifestPairingPartial
from trottercert.dual_word_manifest import (
    FORMULA_NAME,
    ManifestIndex,
    ManifestShardRecord,
    WordGroup,
    WordManifest,
    WordRecord,
    _canonical_bytes,
    _index_payload,
    _manifest_payload,
    _sha256,
    _word_set_digest,
    manifest_shard_index,
)


ROOT = Path(__file__).resolve().parents[1]
E5 = ROOT / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
E7 = ROOT / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
SOURCES = {"source.py": "a" * 64}


def _encoded(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _runtime() -> dict[str, object]:
    return {"wall_seconds": "1.0", "peak_rss_bytes": 1024, "scheduler": None}


def _build_run(root: Path) -> tuple[Path, ...]:
    manifests_root = root / "manifests"
    shards_root = root / "shards"
    reruns_root = root / "independent-rerun"
    manifests_root.mkdir(parents=True)
    shards_root.mkdir()
    reruns_root.mkdir()
    groups = tuple(
        WordGroup(
            ordinal=ordinal,
            suffix=(0, 0, 0, 0, 0, 0, 0, ordinal),
            records=(
                WordRecord(
                    word=(0, 0, 0, 0, 0, 0, 0, 0, ordinal),
                    coefficient=Cubic(ordinal + 1, ordinal + 2, ordinal + 3),
                ),
            ),
        )
        for ordinal in range(4)
    )
    digest = _word_set_digest(groups)
    manifests = tuple(
        WordManifest(
            degree=9,
            shard_index=shard_index,
            shard_count=2,
            total_groups=4,
            total_words=4,
            word_set_digest=digest,
            formula=FORMULA_NAME,
            groups=tuple(
                group
                for group in groups
                if manifest_shard_index(group.ordinal, 2) == shard_index
            ),
            implementation_sources=tuple(SOURCES.items()),
        )
        for shard_index in range(2)
    )
    records: list[ManifestShardRecord] = []
    for manifest in manifests:
        encoded = _canonical_bytes(_manifest_payload(manifest))
        path = manifests_root / f"manifest-{manifest.shard_index:03d}.json"
        path.write_bytes(encoded)
        records.append(
            ManifestShardRecord(
                shard_index=manifest.shard_index,
                path=path.name,
                sha256=_sha256(encoded),
                group_count=len(manifest.groups),
                word_count=sum(len(group.records) for group in manifest.groups),
            )
        )
    index = ManifestIndex(
        degree=9,
        shard_count=2,
        total_groups=4,
        total_words=4,
        word_set_digest=digest,
        formula=FORMULA_NAME,
        implementation_sources=tuple(SOURCES.items()),
        shards=tuple(records),
    )
    (manifests_root / "index.json").write_bytes(_canonical_bytes(_index_payload(index)))

    workers: list[dict[str, object]] = []
    worker_paths: list[Path] = []
    for shard_index, manifest in enumerate(manifests):
        tau_h = Cubic(shard_index + 1, shard_index + 2, shard_index + 3)
        tau_h2 = Cubic(shard_index + 4, shard_index + 5, shard_index + 6)
        partial = ManifestPairingPartial(
            degree=9,
            length=E9_CONTRACTION_LENGTH,
            shard_index=shard_index,
            shard_count=2,
            total_groups=4,
            group_indices=tuple(group.ordinal for group in manifest.groups),
            word_count=len(manifest.groups),
            nonzero_word_count=len(manifest.groups),
            retained_term_count=4 * len(manifest.groups),
            tau_h=tau_h,
            tau_h2=tau_h2,
            tau_w=tau_h2 + tau_h / 2,
        )
        worker = build_worker_payload(
            partial,
            index,
            manifest,
            records[shard_index].sha256,
            implementation_sources=SOURCES,
            runtime=_runtime(),
        )
        path = shards_root / f"shard-{shard_index:03d}.json"
        path.write_bytes(_encoded(worker))
        workers.append(worker)
        worker_paths.append(path)
    reduced = reduce_worker_payloads(
        tuple(workers),
        index=index,
        manifests=manifests,
        parent_sha256=tuple(_sha256(path.read_bytes()) for path in worker_paths),
        extensive_path=E5,
        dual_e7_path=E7,
    )
    (root / "dual-e9-pairing.json").write_bytes(_encoded(reduced))

    rerun_paths: list[Path] = []
    for shard_index in recommended_rerun_shards(index):
        rerun = json.loads(_encoded(workers[shard_index]))
        rerun["runtime"] = {
            "wall_seconds": "2.0",
            "peak_rss_bytes": 2048,
            "scheduler": None,
        }
        path = reruns_root / f"shard-{shard_index:03d}.json"
        path.write_bytes(_encoded(rerun))
        rerun_paths.append(path)
    return tuple(rerun_paths)


def test_complete_run_and_two_reruns_are_audited(tmp_path: Path) -> None:
    run = tmp_path / "run"
    reruns = _build_run(run)

    audit = audit_e9_run(run, rerun_paths=reruns, expected_shards=2)

    assert audit["coverage"]["group_count"] == 4
    assert audit["shard_count"] == 2
    assert [record["shard_index"] for record in audit["independent_reruns"]] == [0, 1]
    assert audit["claim"]["finite_step_status"] == "inconclusive"


def test_missing_shard_is_rejected(tmp_path: Path) -> None:
    run = tmp_path / "run"
    reruns = _build_run(run)
    (run / "shards/shard-001.json").unlink()

    with pytest.raises(ValueError, match="missing shard 1"):
        audit_e9_run(run, rerun_paths=reruns, expected_shards=2)


def test_reducer_parent_hash_must_match_archived_worker(tmp_path: Path) -> None:
    run = tmp_path / "run"
    reruns = _build_run(run)
    path = run / "shards/shard-001.json"
    worker = json.loads(path.read_bytes())
    worker["runtime"]["wall_seconds"] = "3.0"
    path.write_bytes(_encoded(worker))

    with pytest.raises(ValueError, match="parent digest"):
        audit_e9_run(run, rerun_paths=reruns, expected_shards=2)


def test_both_deterministic_reruns_are_required(tmp_path: Path) -> None:
    run = tmp_path / "run"
    reruns = _build_run(run)

    with pytest.raises(ValueError, match="median and heavy"):
        audit_e9_run(run, rerun_paths=reruns[:1], expected_shards=2)


def test_rerun_mathematical_drift_is_rejected(tmp_path: Path) -> None:
    run = tmp_path / "run"
    reruns = list(_build_run(run))
    payload = json.loads(reruns[0].read_bytes())
    payload["mathematical_payload_sha256"] = "0" * 64
    reruns[0].write_bytes(_encoded(payload))

    with pytest.raises(ValueError, match="digest|mathematical"):
        audit_e9_run(run, rerun_paths=reruns, expected_shards=2)
