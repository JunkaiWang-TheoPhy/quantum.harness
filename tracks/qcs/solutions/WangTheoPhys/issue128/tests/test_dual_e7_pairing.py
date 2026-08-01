import copy

import pytest

from scripts.certify_dual_e7_pairing import (
    build_shard_payload,
    reduce_shard_payloads,
    verify_reduced_payload,
    verify_shard_payload,
)
from trottercert.cubic_field import Cubic
from trottercert.dual_log_pairing import DualPairingPartial


def _partial(index: int) -> DualPairingPartial:
    tau_h = Cubic(index + 1, index + 2, index + 3)
    tau_h2 = Cubic(index + 4, index + 5, index + 6)
    return DualPairingPartial(
        degree=7,
        length=6,
        shard_index=index,
        shard_count=2,
        total_groups=4,
        group_indices=(index, index + 2),
        word_count=8,
        nonzero_word_count=5,
        retained_term_count=20,
        tau_h=tau_h,
        tau_h2=tau_h2,
        tau_w=tau_h2 + tau_h / 2,
    )


def _sources():
    return {"dual.py": "a" * 64}


def test_shard_payload_round_trip_and_mutation_rejection() -> None:
    payload = build_shard_payload(_partial(0), _sources())
    verify_shard_payload(payload)
    forged = copy.deepcopy(payload)
    forged["group_indices"] = [0, 1]
    with pytest.raises(ValueError, match="group"):
        verify_shard_payload(forged)
    forged = copy.deepcopy(payload)
    forged["pairings"]["tau_w"][0] = [0, 1]
    with pytest.raises(ValueError, match="pairing"):
        verify_shard_payload(forged)


def test_reducer_requires_complete_disjoint_shards() -> None:
    shards = [build_shard_payload(_partial(index), _sources()) for index in range(2)]
    reduced = reduce_shard_payloads(shards, parent_sha256=("1" * 64, "2" * 64))
    verify_reduced_payload(reduced)
    assert reduced["coverage"]["group_count"] == 4
    assert reduced["coverage"]["reduction_order_check"] == "forward_equals_reverse"
    assert reduced["claim"]["full_e7_operator"] == "not_computed"
    assert reduced["hpc_authorized"] is False

    with pytest.raises(ValueError, match="shard"):
        reduce_shard_payloads(shards[:1], parent_sha256=("1" * 64,))
    overlapping = copy.deepcopy(shards)
    overlapping[1]["shard_index"] = 0
    with pytest.raises(ValueError, match="shard"):
        reduce_shard_payloads(overlapping, parent_sha256=("1" * 64, "2" * 64))


def test_reducer_rejects_source_mismatch() -> None:
    first = build_shard_payload(_partial(0), _sources())
    second = build_shard_payload(_partial(1), {"dual.py": "b" * 64})
    with pytest.raises(ValueError, match="source"):
        reduce_shard_payloads(
            (first, second),
            parent_sha256=("1" * 64, "2" * 64),
        )
