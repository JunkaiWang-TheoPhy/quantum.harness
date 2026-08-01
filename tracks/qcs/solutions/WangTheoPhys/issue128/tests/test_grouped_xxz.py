from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import pytest

from trottercert.algebra import PauliString
from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.grouped_xxz import (
    StageRecord,
    XXZCompileSpec,
    canonical_json_bytes,
    fraction_pair,
    replay_schedule_resources,
    strict_fraction_pair,
    verify_xxz_suzuki_schedule,
    xxz_suzuki_schedule,
)
from trottercert.hamiltonian import (
    finite_torus_xxz_hamiltonian,
    four_matching_fragments,
    four_matching_xxz_fragments,
    xxz_bond_operator,
)
from trottercert.lattice import SquareLattice


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "compile_grouped_xxz.py"


def test_pilot_spec_pins_the_physics() -> None:
    spec = XXZCompileSpec.pilot(Fraction(1, 2))
    assert spec.delta == Fraction(1, 2)
    assert spec.length == 4
    assert spec.boundary == "periodic"
    assert spec.time == Fraction(1)
    assert spec.tolerance == Fraction(1, 10**6)
    assert spec.normalization == "(XX+YY+delta*ZZ)/4"
    assert spec.formula_identifier == "five_copy_fourth_order_suzuki_four_matchings"
    assert spec.stage_count == 31
    assert spec.primary_metric == "merged_group_exponentials"


@pytest.mark.parametrize("delta", [True, 1, 0.5, "1/2"])
def test_spec_requires_exact_delta(delta: object) -> None:
    with pytest.raises(TypeError):
        XXZCompileSpec.pilot(delta)  # type: ignore[arg-type]


def test_compile_spec_is_immutable_and_rejects_manual_invalid_state() -> None:
    spec = XXZCompileSpec.pilot(Fraction(2))
    with pytest.raises(AttributeError):
        spec.length = 6  # type: ignore[misc]
    with pytest.raises(ValueError, match="pilot length"):
        replace(spec, length=6)
    with pytest.raises(ValueError, match="stage count"):
        replace(spec, stage_count=30)


def test_nonpilot_periodic_spec_fails_closed_until_scaling_verifier() -> None:
    with pytest.raises(NotImplementedError, match="Task 10"):
        XXZCompileSpec.periodic(Fraction(1, 2), 12, {})
    with pytest.raises(ValueError, match="even"):
        XXZCompileSpec.periodic(Fraction(1, 2), 5, {})


@pytest.mark.parametrize(
    ("raw", "expected"),
    [([0, 1], Fraction(0)), ([-3, 7], Fraction(-3, 7)), ([5, 1], Fraction(5))],
)
def test_fraction_pair_round_trip(raw: list[int], expected: Fraction) -> None:
    value = strict_fraction_pair(raw, field="delta")
    assert value == expected
    assert fraction_pair(value) == raw


@pytest.mark.parametrize(
    "raw",
    [
        (1, 2),
        [2, 4],
        [0, 2],
        [1, -2],
        [1, 0],
        [True, 2],
        [1.0, 2],
        [1, 2, 3],
    ],
)
def test_fraction_pair_rejects_noncanonical_values(raw: object) -> None:
    with pytest.raises((TypeError, ValueError), match="delta"):
        strict_fraction_pair(raw, field="delta")


def test_canonical_json_bytes_are_compact_sorted_and_newline_terminated() -> None:
    assert canonical_json_bytes({"z": [1, 2], "a": "x"}) == b'{"a":"x","z":[1,2]}\n'


def test_xxz_bond_operator_has_exact_axis_weights() -> None:
    bond = xxz_bond_operator(0, 5, Fraction(3, 2))
    assert bond.terms[PauliString({0: "X", 5: "X"})].real == Fraction(1, 4)
    assert bond.terms[PauliString({0: "Y", 5: "Y"})].real == Fraction(1, 4)
    assert bond.terms[PauliString({0: "Z", 5: "Z"})].real == Fraction(3, 8)
    assert all(coefficient.imag == 0 for coefficient in bond.terms.values())


@pytest.mark.parametrize("delta", [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2)])
def test_four_matching_xxz_fragments_cover_finite_torus(delta: Fraction) -> None:
    lattice = SquareLattice(4)
    matchings = lattice.four_matchings()
    assert len(matchings) == 4
    assert all(len(group) == 8 for group in matchings)
    assert sum(len(group) for group in matchings) == 32
    assert all(left != right for group in matchings for left, right in group)
    assert all(
        set(matchings[left]).isdisjoint(matchings[right])
        for left in range(4)
        for right in range(left + 1, 4)
    )

    fragments = four_matching_xxz_fragments(lattice, delta)
    assert len(fragments) == 4
    total = sum(fragments[1:], fragments[0])
    assert total == finite_torus_xxz_hamiltonian(lattice, delta)


def test_xxz_isotropic_fragments_equal_existing_fragments() -> None:
    lattice = SquareLattice(4)
    assert four_matching_xxz_fragments(lattice, Fraction(1)) == four_matching_fragments(
        lattice
    )


def test_xxz_suzuki_schedule_is_exact_and_replays_resources() -> None:
    schedule = xxz_suzuki_schedule()
    source = fourth_order_suzuki_cubic_stages(4)
    assert len(schedule) == 31
    assert all(0 <= stage.fragment_index <= 3 for stage in schedule)
    assert schedule == tuple(reversed(schedule))
    assert tuple(
        (stage.fragment_index, stage.coefficient_coordinates) for stage in schedule
    ) == tuple(
        (
            stage.fragment_index,
            (stage.coefficient.a0, stage.coefficient.a1, stage.coefficient.a2),
        )
        for stage in source
    )
    for fragment in range(4):
        total = sum(
            (
                Cubic(*stage.coefficient_coordinates)
                for stage in schedule
                if stage.fragment_index == fragment
            ),
            Cubic.zero(),
        )
        assert total == Cubic.one()
    for steps in (1, 2, 95):
        assert replay_schedule_resources(schedule, steps) == 30 * steps + 1


def test_schedule_verifier_rejects_order_coordinate_duplicate_and_fragment_mutations() -> None:
    schedule = xxz_suzuki_schedule()

    swapped = list(schedule)
    index = next(
        i
        for i in range(len(swapped) - 1)
        if swapped[i].fragment_index != swapped[i + 1].fragment_index
    )
    swapped[index], swapped[index + 1] = swapped[index + 1], swapped[index]

    coordinate = schedule[0].coefficient_coordinates
    flipped = list(schedule)
    flipped[0] = replace(
        flipped[0],
        coefficient_coordinates=(-coordinate[0], coordinate[1], coordinate[2]),
    )

    duplicated = list(schedule)
    duplicated.insert(1, duplicated[0])

    replaced_fragment = list(schedule)
    replaced_fragment[0] = replace(
        replaced_fragment[0], fragment_index=(replaced_fragment[0].fragment_index + 1) % 4
    )

    for mutated in (swapped, flipped, duplicated, replaced_fragment):
        with pytest.raises(ValueError, match="schedule"):
            verify_xxz_suzuki_schedule(tuple(mutated))


def test_stage_record_rejects_nonexact_or_out_of_range_fields() -> None:
    with pytest.raises(TypeError):
        StageRecord(0, (Fraction(1), Fraction(0), 0.0))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="fragment"):
        StageRecord(4, (Fraction(1), Fraction(0), Fraction(0)))


def test_profile_cli_validates_setup_without_building_artifacts() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(CLI),
            "--profile-only",
            "--delta",
            "1/2",
            "--pipeline",
            "direct-theorem",
            "--max-records",
            "200",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload == {
        "delta": [1, 2],
        "max_records": 200,
        "mode": "profile-only",
        "pipeline": "direct-theorem",
        "setup": "validated",
        "stage_count": 31,
    }


def test_build_cli_fails_closed_until_task_8() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(CLI),
            "--build",
            "--delta",
            "1/2",
            "--pipeline",
            "direct-theorem",
            "--summary",
            "unused.json",
            "--witness",
            "unused.json.gz",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert completed.returncode != 0
    assert "Task 8" in completed.stderr
