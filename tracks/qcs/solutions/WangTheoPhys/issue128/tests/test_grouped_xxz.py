from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from dataclasses import replace
from fractions import Fraction
from itertools import product
from pathlib import Path
from time import monotonic

import pytest

from trottercert.algebra import PauliString, PauliSum
from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.grouped_xxz import (
    FINITE_STEP_ERROR_FORMULA,
    GROUPING_ALGORITHM_IDENTIFIER,
    THEOREM_CENTER,
    THEOREM_DUHAMEL_CONVENTION,
    THEOREM_FACTORIAL_DENOMINATOR,
    THEOREM_IDENTIFIER,
    THEOREM_ORDER,
    AnticommutingGroupRecord,
    StageRecord,
    SymplecticCoefficient,
    XXZCertificate,
    XXZCompileSpec,
    _close_xxz_certificate,
    build_finite_xxz_ledger,
    canonical_json_bytes,
    compile_grouped_xxz,
    discover_anticommuting_groups,
    discover_xxz_groups,
    fraction_pair,
    replay_schedule_resources,
    sqrt_fraction_interval,
    strict_fraction_pair,
    verify_anticommuting_groups,
    verify_finite_xxz_ledger,
    verify_xxz_finite_step_fields,
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


def test_profile_cli_runs_a_bounded_compressed_prefix_without_artifacts() -> None:
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
            "2",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["delta"] == [1, 2]
    assert payload["max_records"] == 2
    assert payload["mode"] == "profile-only"
    assert payload["pipeline"] == "direct-theorem"
    assert payload["raw_records"] == 2
    assert payload["projected_raw_records"] == 61_677
    assert payload["active_words"] == 2
    assert payload["representative_blocks"] > 0
    assert payload["representative_terms"] > 0
    assert payload["complete"] is False
    assert "summary" not in payload and "witness" not in payload


def test_build_cli_requires_outputs_and_rejects_max_records() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            str(CLI),
            "--build",
            "--delta",
            "1/2",
            "--pipeline",
            "direct-theorem",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert completed.returncode != 0
    assert "--summary and --witness" in completed.stderr

    bounded = subprocess.run(
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
            "--max-records",
            "2",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
    )
    assert bounded.returncode != 0
    assert "valid only with --profile-only" in bounded.stderr


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


def _test_anticommutes(
    left: tuple[int, int],
    right: tuple[int, int],
) -> bool:
    return bool(
        ((left[0] & right[1]).bit_count() + (left[1] & right[0]).bit_count())
        & 1
    )


def _naive_pair_only_groups(
    terms: tuple[SymplecticCoefficient, ...],
) -> tuple[tuple[SymplecticCoefficient, ...], ...]:
    ordered = sorted(terms, key=lambda term: (-abs(term.real), term.mask))
    unmatched = set(range(len(ordered)))
    groups: list[tuple[SymplecticCoefficient, ...]] = []
    for index, term in enumerate(ordered):
        if index not in unmatched:
            continue
        unmatched.remove(index)
        partner = next(
            (
                candidate
                for candidate in sorted(unmatched)
                if _test_anticommutes(term.mask, ordered[candidate].mask)
            ),
            None,
        )
        if partner is None:
            groups.append((term,))
        else:
            unmatched.remove(partner)
            groups.append((term, ordered[partner]))
    assert not unmatched
    return tuple(groups)


def _test_group_record(
    terms: tuple[SymplecticCoefficient, ...],
) -> AnticommutingGroupRecord:
    squared = sum((term.real * term.real for term in terms), Fraction())
    return AnticommutingGroupRecord(
        terms=terms,
        squared_norm=squared,
        norm_interval=sqrt_fraction_interval(squared),
    )


def test_pair_only_bitsets_choose_the_same_earliest_partners_as_naive_scan() -> None:
    fixtures = [_four_term_fixture()]
    generator = random.Random(20260801)
    for qubits, count in ((2, 9), (3, 20), (4, 40)):
        masks: set[tuple[int, int]] = set()
        while len(masks) < count:
            mask = (
                generator.randrange(1 << qubits),
                generator.randrange(1 << qubits),
            )
            if mask != (0, 0):
                masks.add(mask)
        fixtures.append(
            tuple(
                SymplecticCoefficient(
                    x_mask,
                    z_mask,
                    Fraction(generator.choice((-3, -2, -1, 1, 2, 3)), 7),
                    Fraction(),
                )
                for x_mask, z_mask in masks
            )
        )

    for terms in fixtures:
        observed = discover_anticommuting_groups(terms)
        expected = _naive_pair_only_groups(terms)
        assert tuple(group.terms for group in observed) == expected
        assert all(len(group.terms) in (1, 2) for group in observed)
        verify_anticommuting_groups(terms, observed)


def test_exact_group_verifier_checks_coverage_commutation_and_sqrt_interval() -> None:
    assert sqrt_fraction_interval.cache_info().maxsize == 8192
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


def test_pair_only_verifier_rejects_larger_group_even_if_pairwise_anticommuting() -> None:
    terms = (
        SymplecticCoefficient(1, 0, Fraction(1), Fraction()),
        SymplecticCoefficient(0, 1, Fraction(1), Fraction()),
        SymplecticCoefficient(1, 1, Fraction(1), Fraction()),
    )
    oversized = AnticommutingGroupRecord(
        terms=terms,
        squared_norm=Fraction(3),
        norm_interval=sqrt_fraction_interval(Fraction(3)),
    )
    with pytest.raises(ValueError, match="singleton or pair"):
        verify_anticommuting_groups(terms, (oversized,))


@pytest.mark.slow
def test_pair_only_200_record_profile_is_bounded_and_fully_verified() -> None:
    started = monotonic()
    ledger = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=200,
    )
    groups = discover_xxz_groups(ledger)
    grouped_bound = verify_xxz_groups(ledger, groups)
    elapsed = monotonic() - started
    group_sizes = [len(group.terms) for record in groups for group in record.groups]

    assert len(ledger.raw_records) == 200
    assert sum(len(block.terms) for block in ledger.blocks) == 402_720
    assert group_sizes.count(2) == 200_552
    assert group_sizes.count(1) == 1_616
    assert all(size in (1, 2) for size in group_sizes)
    assert grouped_bound <= xxz_triangle_baseline(ledger)
    assert elapsed < 180


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
    assert all(
        record.grouping_algorithm == GROUPING_ALGORITHM_IDENTIFIER
        for record in groups
    )
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

    wrong_algorithm = list(groups)
    wrong_algorithm[0] = replace(
        wrong_algorithm[0],
        grouping_algorithm="unbound_pairing_heuristic",
    )
    with pytest.raises(ValueError, match="algorithm identifier"):
        verify_xxz_groups(ledger, tuple(wrong_algorithm))


def test_group_verifier_rejects_alternative_valid_pairing() -> None:
    ledger = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=1,
    )
    records = list(discover_xxz_groups(ledger))
    record = records[0]
    alternative = list(record.groups)
    replacement: tuple[
        int,
        int,
        AnticommutingGroupRecord,
        AnticommutingGroupRecord,
    ] | None = None
    for left_index, left_group in enumerate(record.groups):
        if len(left_group.terms) != 2:
            continue
        for right_index in range(left_index + 1, len(record.groups)):
            right_group = record.groups[right_index]
            if len(right_group.terms) != 2:
                continue
            a, b = left_group.terms
            c, d = right_group.terms
            for first, second in (((a, c), (b, d)), ((a, d), (b, c))):
                if _test_anticommutes(first[0].mask, first[1].mask) and (
                    _test_anticommutes(second[0].mask, second[1].mask)
                ):
                    replacement = (
                        left_index,
                        right_index,
                        _test_group_record(first),
                        _test_group_record(second),
                    )
                    break
            if replacement is not None:
                break
        if replacement is not None:
            break
    assert replacement is not None
    left_index, right_index, first, second = replacement
    alternative[left_index] = first
    alternative[right_index] = second
    assert verify_anticommuting_groups(
        ledger.blocks[0].terms,
        tuple(alternative),
    ) > 0

    records[0] = replace(record, groups=tuple(alternative))
    with pytest.raises(ValueError, match="deterministic pair discovery"):
        verify_xxz_groups(ledger, tuple(records))


def test_fake_complete_ledger_closes_exact_adjacent_steps_and_resources() -> None:
    bounded = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=1,
    )
    fake_complete = replace(bounded, complete=True)
    certificate = _close_xxz_certificate(
        fake_complete,
        grouped_constant=Fraction(16, 10**6),
        triangle_constant=Fraction(81, 10**6),
    )
    assert isinstance(certificate, XXZCertificate)
    verify_xxz_finite_step_fields(certificate)

    assert certificate.method == "direct_finite_high_order_theorem_grouped_norm"
    assert certificate.grouping_algorithm == GROUPING_ALGORITHM_IDENTIFIER
    assert certificate.grouped_constant == Fraction(16, 10**6)
    assert certificate.triangle_constant == Fraction(81, 10**6)
    assert certificate.candidate_steps == 2
    assert certificate.candidate_error == Fraction(1, 10**6)
    assert certificate.candidate_previous_error == Fraction(16, 10**6)
    assert certificate.baseline_steps == 3
    assert certificate.baseline_error == Fraction(1, 10**6)
    assert certificate.baseline_previous_error == Fraction(81, 16 * 10**6)
    assert certificate.candidate_resources == 30 * 2 + 1
    assert certificate.baseline_resources == 30 * 3 + 1
    assert certificate.candidate_steps != 393
    assert certificate.bond_growth == 1
    assert certificate.cell_base == Fraction(5, 4)


def test_finite_step_verifier_rejects_nonminimal_step_and_resource_mutations() -> None:
    bounded = build_finite_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(2)),
        max_records=1,
    )
    certificate = _close_xxz_certificate(
        replace(bounded, complete=True),
        grouped_constant=Fraction(16, 10**6),
        triangle_constant=Fraction(81, 10**6),
    )
    attacks = (
        replace(certificate, candidate_steps=3),
        replace(certificate, candidate_error=Fraction()),
        replace(certificate, candidate_previous_error=Fraction()),
        replace(certificate, baseline_steps=4),
        replace(certificate, baseline_resources=certificate.baseline_resources + 30),
        replace(certificate, grouping_algorithm="unbound_pairing_heuristic"),
        replace(certificate, bond_growth=Fraction(1)),
        replace(certificate, cell_base=Fraction(1)),
    )
    for attacked in attacks:
        with pytest.raises(ValueError):
            verify_xxz_finite_step_fields(attacked)


def test_compile_grouped_xxz_rejects_real_bounded_ledger_before_grouping() -> None:
    spec = XXZCompileSpec.pilot(Fraction(1, 2))
    bounded = build_finite_xxz_ledger(spec, max_records=1)
    assert not bounded.complete
    with pytest.raises(ValueError, match="complete ledger"):
        compile_grouped_xxz(spec, ledger=bounded)
