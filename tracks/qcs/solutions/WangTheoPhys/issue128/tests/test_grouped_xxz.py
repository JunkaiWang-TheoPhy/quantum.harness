from __future__ import annotations

import json
import os
import subprocess
import sys
from dataclasses import replace
from fractions import Fraction
from itertools import product
from pathlib import Path

import pytest

from trottercert.algebra import PauliString, PauliSum
from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.grouped_xxz import (
    FINITE_STEP_ERROR_FORMULA,
    THEOREM_CENTER,
    THEOREM_DUHAMEL_CONVENTION,
    THEOREM_FACTORIAL_DENOMINATOR,
    THEOREM_IDENTIFIER,
    THEOREM_ORDER,
    AnticommutingGroupRecord,
    StageRecord,
    SymplecticCoefficient,
    XXZCompileSpec,
    build_finite_xxz_ledger,
    canonical_json_bytes,
    discover_anticommuting_groups,
    discover_xxz_groups,
    fraction_pair,
    replay_schedule_resources,
    sqrt_fraction_interval,
    strict_fraction_pair,
    verify_anticommuting_groups,
    verify_finite_xxz_ledger,
    verify_xxz_groups,
    verify_xxz_suzuki_schedule,
    weighted_symplectic_nested_commutator,
    weighted_xxz_fragment,
    xxz_suzuki_schedule,
    xxz_triangle_baseline,
)
from trottercert.hamiltonian import (
    finite_torus_xxz_hamiltonian,
    four_matching_fragments,
    four_matching_xxz_fragments,
    xxz_bond_operator,
)
from trottercert.higher_order import nested_commutator
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


def _open_fragment_operators(
    fragment_bonds: tuple[tuple[tuple[int, int], ...], ...],
    delta: Fraction,
) -> tuple[PauliSum, ...]:
    fragments = []
    for bonds in fragment_bonds:
        operator = PauliSum.zero()
        for left, right in bonds:
            operator += xxz_bond_operator(left, right, delta)
        fragments.append(operator)
    return tuple(fragments)


def _pauli_masks(operator: PauliSum) -> dict[tuple[int, int], tuple[Fraction, Fraction]]:
    result: dict[tuple[int, int], tuple[Fraction, Fraction]] = {}
    for pauli, coefficient in operator.terms.items():
        x_mask = z_mask = 0
        for site, axis in pauli.ops:
            bit = 1 << site
            if axis in {"X", "Y"}:
                x_mask |= bit
            if axis in {"Y", "Z"}:
                z_mask |= bit
        result[(x_mask, z_mask)] = (coefficient.real, coefficient.imag)
    return result


@pytest.mark.parametrize("delta", [Fraction(1, 2), Fraction(2)])
def test_weighted_symplectic_backend_matches_every_reduced_word(
    delta: Fraction,
) -> None:
    fragment_bonds = (((0, 1),), ((1, 2),))
    weighted = weighted_xxz_fragment(fragment_bonds[0], delta)
    assert weighted.common_denominator == 4 * delta.denominator
    assert weighted.axis_numerators == (
        delta.denominator,
        delta.denominator,
        delta.numerator,
    )

    references = _open_fragment_operators(fragment_bonds, delta)
    cache: dict[tuple[int, ...], PauliSum] = {}
    for degree in (3, 5):
        for key in product(range(2), repeat=degree):
            observed = {
                (term.x_mask, term.z_mask): (term.real, term.imag)
                for term in weighted_symplectic_nested_commutator(
                    fragment_bonds,
                    delta,
                    key,
                )
            }
            assert observed == _pauli_masks(nested_commutator(references, key, cache))


def test_bounded_ledger_is_raw_prefix_with_exact_projection_and_blocks() -> None:
    progress: list[tuple[int, int]] = []
    ledger = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=40,
        progress=lambda completed, projected: progress.append((completed, projected)),
    )
    verify_finite_xxz_ledger(ledger)

    assert ledger.theorem_identifier == THEOREM_IDENTIFIER
    assert ledger.order == THEOREM_ORDER == 4
    assert ledger.center == THEOREM_CENTER == 20
    assert ledger.factorial_denominator == THEOREM_FACTORIAL_DENOMINATOR == 120
    assert ledger.duhamel_convention == THEOREM_DUHAMEL_CONVENTION
    assert ledger.finite_step_error_formula == FINITE_STEP_ERROR_FORMULA
    assert ledger.projected_record_count == 61_677
    assert len(ledger.raw_records) == 40
    assert not ledger.complete
    assert progress[-1] == (40, 61_677)
    assert len({record.record_id for record in ledger.raw_records}) == 40
    assert all(record.weight_interval.lower > 0 for record in ledger.raw_records)
    assert all(len(record.block_key) == 5 for record in ledger.raw_records)

    weights: dict[tuple[int, ...], Fraction] = {}
    for record in ledger.raw_records:
        weights[record.block_key] = (
            weights.get(record.block_key, Fraction()) + record.weight_interval.upper
        )
    assert {block.block_key: block.raw_weight_upper for block in ledger.blocks} == weights
    for block in ledger.blocks:
        assert block.terms == tuple(sorted(block.terms, key=lambda term: term.mask))
        assert all(not term.is_zero for term in block.terms)


def test_ledger_verifier_rejects_raw_weight_and_preweighted_block_mutations() -> None:
    ledger = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(2)),
        max_records=5,
    )
    first_record = ledger.raw_records[0]
    bad_record = replace(
        first_record,
        weight_interval=first_record.weight_interval * 2,
    )
    with pytest.raises(ValueError):
        verify_finite_xxz_ledger(
            replace(ledger, raw_records=(bad_record,) + ledger.raw_records[1:])
        )

    nonzero_index = next(index for index, block in enumerate(ledger.blocks) if block.terms)
    block = ledger.blocks[nonzero_index]
    preweighted = replace(
        block,
        terms=tuple(
            replace(
                term,
                real=term.real * block.raw_weight_upper,
                imag=term.imag * block.raw_weight_upper,
            )
            for term in block.terms
        ),
    )
    blocks = list(ledger.blocks)
    blocks[nonzero_index] = preweighted
    with pytest.raises(ValueError, match="unweighted"):
        verify_finite_xxz_ledger(replace(ledger, blocks=tuple(blocks)))


def _four_term_fixture() -> tuple[SymplecticCoefficient, ...]:
    return (
        SymplecticCoefficient(1, 0, Fraction(1), Fraction()),  # XI
        SymplecticCoefficient(0, 1, Fraction(1), Fraction()),  # ZI
        SymplecticCoefficient(2, 0, Fraction(1), Fraction()),  # IX
        SymplecticCoefficient(0, 2, Fraction(1), Fraction()),  # IZ
    )


def test_exact_group_verifier_checks_coverage_commutation_and_sqrt_interval() -> None:
    terms = _four_term_fixture()
    groups = discover_anticommuting_groups(terms)
    bound = verify_anticommuting_groups(terms, groups)
    assert len(groups) == 2
    assert bound == sum((group.norm_interval.upper for group in groups), Fraction())
    for group in groups:
        assert group.squared_norm == 2
        assert group.norm_interval.lower**2 <= 2 <= group.norm_interval.upper**2
        assert group.norm_interval.upper - group.norm_interval.lower <= Fraction(1, 10**30)

    duplicate = replace(groups[0], terms=groups[0].terms + (groups[0].terms[0],))
    with pytest.raises(ValueError, match="duplicate"):
        verify_anticommuting_groups(terms, (duplicate,) + groups[1:])
    with pytest.raises(ValueError, match="coverage"):
        verify_anticommuting_groups(terms, (groups[0],))

    commuting = AnticommutingGroupRecord(
        terms=(terms[0], terms[2]),
        squared_norm=Fraction(2),
        norm_interval=sqrt_fraction_interval(Fraction(2)),
    )
    remainder = AnticommutingGroupRecord(
        terms=(terms[1], terms[3]),
        squared_norm=Fraction(2),
        norm_interval=sqrt_fraction_interval(Fraction(2)),
    )
    with pytest.raises(ValueError, match="anticommute"):
        verify_anticommuting_groups(terms, (commuting, remainder))


def test_theorem_group_witness_is_block_local_delta_bound_and_unweighted() -> None:
    ledger = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=1,
    )
    groups = discover_xxz_groups(ledger)
    grouped = verify_xxz_groups(ledger, groups)
    baseline = xxz_triangle_baseline(ledger)
    assert grouped <= baseline
    assert all(record.delta == Fraction(1, 2) for record in groups)
    assert all(record.ledger_digest == ledger.ledger_digest for record in groups)
    assert {record.block_key for record in groups} == {
        block.block_key for block in ledger.blocks
    }

    other_delta = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(2)),
        max_records=1,
    )
    with pytest.raises(ValueError, match="Delta|digest"):
        verify_xxz_groups(other_delta, groups)

    nonempty_index = next(
        index
        for index, record in enumerate(groups)
        if record.groups and record.groups[0].terms
    )
    record = groups[nonempty_index]
    first_group = record.groups[0]
    doubled_norm = replace(
        first_group,
        norm_interval=first_group.norm_interval * 2,
    )
    mutated_groups = list(groups)
    mutated_groups[nonempty_index] = replace(
        record,
        groups=(doubled_norm,) + record.groups[1:],
    )
    with pytest.raises(ValueError, match="norm"):
        verify_xxz_groups(ledger, tuple(mutated_groups))
