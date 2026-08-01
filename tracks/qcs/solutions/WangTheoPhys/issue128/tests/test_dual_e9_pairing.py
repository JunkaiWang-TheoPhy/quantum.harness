from __future__ import annotations

import copy
import hashlib
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import pytest

from scripts.certify_dual_e9_pairing import (
    E9_CONTRACTION_CELLS,
    E9_CONTRACTION_LENGTH,
    build_worker_payload,
    reduce_worker_payloads,
    verify_reduced_payload,
    verify_worker_payload,
)
from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.dual_manifest_pairing import (
    ManifestPairingPartial,
    contract_word_manifest,
)
from trottercert.dual_word_manifest import (
    FORMULA_NAME,
    ManifestIndex,
    ManifestShardRecord,
    WordGroup,
    WordManifest,
    WordRecord,
    build_word_manifests,
    manifest_shard_index,
)

ROOT = Path(__file__).resolve().parents[1]
EXTENSIVE = (
    ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
DUAL_E7 = ROOT / "docs/experiments/processor-obstruction/dual-e7-pairing.json"
SOURCES = {"source.py": "a" * 64}


def _manifest(index: int) -> WordManifest:
    ordinals = tuple(
        ordinal
        for ordinal in range(4)
        if manifest_shard_index(ordinal, 2) == index
    )
    groups = tuple(
        WordGroup(
            ordinal=ordinal,
            suffix=(0,) * 8,
            records=(
                WordRecord(
                    word=(ordinal % 4,) + (0,) * 8,
                    coefficient=Cubic(ordinal + 1, ordinal + 2, ordinal + 3),
                ),
            ),
        )
        for ordinal in ordinals
    )
    return WordManifest(
        degree=9,
        shard_index=index,
        shard_count=2,
        total_groups=4,
        total_words=4,
        word_set_digest="b" * 64,
        formula=FORMULA_NAME,
        groups=groups,
        implementation_sources=tuple(SOURCES.items()),
    )


def _index() -> ManifestIndex:
    return ManifestIndex(
        degree=9,
        shard_count=2,
        total_groups=4,
        total_words=4,
        word_set_digest="b" * 64,
        formula=FORMULA_NAME,
        implementation_sources=tuple(SOURCES.items()),
        shards=tuple(
            ManifestShardRecord(
                shard_index=index,
                path=f"manifest-{index:03d}.json",
                sha256=str(index + 1) * 64,
                group_count=len(_manifest(index).groups),
                word_count=len(_manifest(index).groups),
            )
            for index in range(2)
        ),
    )


def _partial(index: int) -> ManifestPairingPartial:
    tau_h = Cubic(index + 1, index + 2, index + 3)
    tau_h2 = Cubic(index + 4, index + 5, index + 6)
    return ManifestPairingPartial(
        degree=9,
        length=E9_CONTRACTION_LENGTH,
        shard_index=index,
        shard_count=2,
        total_groups=4,
        group_indices=tuple(group.ordinal for group in _manifest(index).groups),
        word_count=len(_manifest(index).groups),
        nonzero_word_count=len(_manifest(index).groups),
        retained_term_count=4 * len(_manifest(index).groups),
        tau_h=tau_h,
        tau_h2=tau_h2,
        tau_w=tau_h2 + tau_h / 2,
    )


def _runtime() -> dict[str, object]:
    return {
        "wall_seconds": "1.25",
        "peak_rss_bytes": 1024,
        "scheduler": None,
    }


def _worker(index: int) -> dict[str, object]:
    manifest = _manifest(index)
    record = _index().shards[index]
    return build_worker_payload(
        _partial(index),
        _index(),
        manifest,
        record.sha256,
        implementation_sources=SOURCES,
        runtime=_runtime(),
    )


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_e9_worker_payload_round_trips_and_rejects_attacks() -> None:
    index = _index()
    manifest = _manifest(0)
    payload = _worker(0)
    verify_worker_payload(payload, index, manifest)

    assert payload["degree"] == 9
    assert payload["claim"]["full_e9_operator"] == "not_computed"
    assert payload["claim"]["dual_e9_pairing"] == "exact"
    assert payload["claim"]["finite_step_status"] == "inconclusive"

    forged = copy.deepcopy(payload)
    forged["manifest_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="manifest"):
        verify_worker_payload(forged, index, manifest)

    forged = copy.deepcopy(payload)
    forged["runtime"]["wall_seconds"] = "2.0"
    verify_worker_payload(forged, index, manifest)

    forged = copy.deepcopy(payload)
    forged["pairings"]["tau_w"][0] = [0, 1]
    with pytest.raises(ValueError, match="pairing|digest"):
        verify_worker_payload(forged, index, manifest)


def test_e9_reducer_requires_complete_manifest_bound_coverage() -> None:
    index = _index()
    manifests = tuple(_manifest(shard) for shard in range(2))
    workers = tuple(_worker(shard) for shard in range(2))
    reduced = reduce_worker_payloads(
        workers,
        index=index,
        manifests=manifests,
        parent_sha256=("c" * 64, "d" * 64),
        extensive_path=EXTENSIVE,
        dual_e7_path=DUAL_E7,
    )
    verify_reduced_payload(reduced, index=index, manifests=manifests)

    assert reduced["coverage"]["group_count"] == 4
    assert reduced["pairings"]["tau_w_per_cell"]
    assert reduced["pairings"]["q9_over_97_pow_8"]
    assert reduced["claim"]["finite_step_status"] == "inconclusive"
    assert reduced["claim"]["missing"] == "E11-and-higher dual tail"
    assert reduced["inputs"]["extensive_witness_sha256"] == _sha(EXTENSIVE)
    assert reduced["inputs"]["dual_e7_sha256"] == _sha(DUAL_E7)

    reversed_reduction = reduce_worker_payloads(
        tuple(reversed(workers)),
        index=index,
        manifests=manifests,
        parent_sha256=("d" * 64, "c" * 64),
        extensive_path=EXTENSIVE,
        dual_e7_path=DUAL_E7,
    )
    assert [parent["sha256"] for parent in reversed_reduction["parents"]] == [
        "c" * 64,
        "d" * 64,
    ]
    assert reversed_reduction["pairings"] == reduced["pairings"]

    with pytest.raises(ValueError, match="complete|coverage|shard"):
        reduce_worker_payloads(
            workers[:-1],
            index=index,
            manifests=manifests,
            parent_sha256=("c" * 64,),
            extensive_path=EXTENSIVE,
            dual_e7_path=DUAL_E7,
        )

    forged = copy.deepcopy(reduced)
    forged["claim"]["finite_step_status"] = "proved"
    with pytest.raises(ValueError, match="finite-step|evidence|digest"):
        verify_reduced_payload(forged, index=index, manifests=manifests)


def test_e9_worker_rejects_aliased_six_by_six_configuration() -> None:
    with pytest.raises(ValueError, match="alias-free 12x12"):
        build_worker_payload(
            replace(_partial(0), length=6),
            _index(),
            _manifest(0),
            _index().shards[0].sha256,
            implementation_sources=SOURCES,
            runtime=_runtime(),
        )


def test_e9_per_cell_normalization_uses_target_torus_cells() -> None:
    assert E9_CONTRACTION_LENGTH == 12
    assert E9_CONTRACTION_CELLS == 36


def test_production_alias_counterexample_contracts_on_target_torus() -> None:
    suffix = (0, 0, 0, 0, 1, 1, 1, 0)
    first = Cubic(
        Fraction(97210657, 19440000000000),
        Fraction(15162241, 6480000000000),
        Fraction(52637263, 77760000000000),
    )
    repeated = Cubic(
        Fraction(9654311, 1944000000000),
        Fraction(18044341, 7776000000000),
        Fraction(1355431, 1944000000000),
    )
    group = WordGroup(
        ordinal=84,
        suffix=suffix,
        records=tuple(
            WordRecord((prefix,) + suffix, first if prefix == 0 else repeated)
            for prefix in range(4)
        ),
    )
    manifest = WordManifest(
        degree=9,
        shard_index=manifest_shard_index(group.ordinal, 64),
        shard_count=64,
        total_groups=65536,
        total_words=262140,
        word_set_digest="b" * 64,
        formula=FORMULA_NAME,
        groups=(group,),
        implementation_sources=tuple(SOURCES.items()),
    )

    with pytest.raises(ValueError, match="alias"):
        contract_word_manifest(manifest, length=6)
    result = contract_word_manifest(manifest, length=E9_CONTRACTION_LENGTH)

    assert result.retained_term_count == 3618
    assert result.tau_w == Cubic(
        Fraction(200931647, 103680000000000),
        Fraction(226274041, 207360000000000),
        Fraction(-52841453, 46080000000000),
    )


@pytest.mark.slow
def test_actual_degree_nine_one_group_worker_smoke() -> None:
    production = build_word_manifests(
        fourth_order_suzuki_cubic_stages(4),
        degree=9,
        shard_count=64,
        implementation_sources=SOURCES,
    )
    heavy_group = next(
        group for group in production[0].groups if group.ordinal == 4818
    )
    test_group = replace(heavy_group, ordinal=0)
    manifest = replace(
        production[0],
        shard_count=1,
        total_groups=1,
        total_words=len(test_group.records),
        groups=(test_group,),
    )
    index = ManifestIndex(
        degree=9,
        shard_count=1,
        total_groups=1,
        total_words=len(test_group.records),
        word_set_digest=manifest.word_set_digest,
        formula=FORMULA_NAME,
        implementation_sources=manifest.implementation_sources,
        shards=(
            ManifestShardRecord(
                shard_index=0,
                path="manifest-000.json",
                sha256="e" * 64,
                group_count=1,
                word_count=len(test_group.records),
            ),
        ),
    )
    partial = contract_word_manifest(manifest, length=E9_CONTRACTION_LENGTH)
    payload = build_worker_payload(
        partial,
        index,
        manifest,
        "e" * 64,
        implementation_sources=SOURCES,
        runtime=_runtime(),
    )

    verify_worker_payload(payload, index, manifest)
    assert payload["word_count"] == 4
    assert payload["group_indices"] == [0]
    assert partial.retained_term_count == 112299
    assert partial.tau_w == 4 * Cubic(
        Fraction(-8181997, 4608000000000),
        Fraction(-64514921, 55296000000000),
        Fraction(-39898969, 55296000000000),
    )
