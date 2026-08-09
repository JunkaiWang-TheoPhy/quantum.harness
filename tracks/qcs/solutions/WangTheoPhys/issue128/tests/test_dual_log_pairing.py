import json
from pathlib import Path

import pytest

from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.dual_log_pairing import (
    contract_log_degree_shard,
    suffix_group_ordinals,
)

ROOT = Path(__file__).resolve().parents[1]
EXTENSIVE = (
    ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)


def test_suffix_group_assignment_is_disjoint_and_complete() -> None:
    assignments = [suffix_group_ordinals(17, index, 4) for index in range(4)]
    assert set().union(*map(set, assignments)) == set(range(17))
    assert sum(map(len, assignments)) == 17
    assert all(set(left).isdisjoint(right) for left, right in zip(assignments, assignments[1:]))


def test_shard_arguments_fail_closed() -> None:
    stages = fourth_order_suzuki_cubic_stages(4)
    with pytest.raises(ValueError, match="shard"):
        contract_log_degree_shard(stages, 5, -1, 4)
    with pytest.raises(ValueError, match="odd"):
        contract_log_degree_shard(stages, 4, 0, 1)


@pytest.mark.slow
def test_degree_five_dual_contraction_matches_frozen_extensive_record() -> None:
    observed = contract_log_degree_shard(
        fourth_order_suzuki_cubic_stages(4),
        5,
        0,
        1,
    )
    payload = json.loads(EXTENSIVE.read_text())
    expected = next(record for record in payload["records"] if record["length"] == 6)

    def cubic(value):
        from fractions import Fraction

        from trottercert.cubic_field import Cubic

        return Cubic(*(Fraction(*coordinate) for coordinate in value))

    assert observed.tau_h == cubic(expected["moments"]["tau_h_e5"])
    assert observed.tau_h2 == cubic(expected["moments"]["tau_h2_e5"])
    assert observed.tau_w == cubic(expected["moments"]["tau_w_e5"])
    assert observed.group_indices == tuple(range(observed.total_groups))
