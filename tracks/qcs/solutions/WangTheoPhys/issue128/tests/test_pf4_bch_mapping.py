from __future__ import annotations

import copy
import inspect
from collections.abc import Mapping
from dataclasses import replace
from fractions import Fraction

import numpy as np
import pytest

from trottercert.cubic_field import (
    Cubic,
    CubicStage,
    fourth_order_suzuki_cubic_stages,
)
from trottercert.pf4_bch_mapping import (
    CubicWordPolynomial,
    PF4ConventionRecord,
    Word,
    canonical_payload_bytes,
    compare_projectively,
    convention_record_digest,
    cyclic_trace_classes,
    derive_pf4_cyclic_mapping,
    parse_payload,
    payload_digest,
    pf4_convention_record,
    pf4_suzuki_gamma,
    pf4_suzuki_trace_polynomials,
    pf4_two_fragment_stages,
    physical_substitution_factor,
    to_payload,
    verify_pf4_convention_record,
    verify_pf4_suzuki_trace_identity,
)
from trottercert.trace_obstruction import (
    pf4_suzuki_trace_pairing,
    pf4_trace_identity_record,
    pf4_trace_quadratic_form,
    verify_identity_record,
)


def test_two_fragment_pf4_has_frozen_eleven_stage_table() -> None:
    u = Cubic(Fraction(4, 15), Fraction(1, 15), Fraction(1, 60))
    v = Cubic.one() - 4 * u
    expected = (
        CubicStage(0, u / 2),
        CubicStage(1, u),
        CubicStage(0, u),
        CubicStage(1, u),
        CubicStage(0, (u + v) / 2),
        CubicStage(1, v),
        CubicStage(0, (u + v) / 2),
        CubicStage(1, u),
        CubicStage(0, u),
        CubicStage(1, u),
        CubicStage(0, u / 2),
    )

    assert pf4_two_fragment_stages() == expected
    assert pf4_two_fragment_stages() == fourth_order_suzuki_cubic_stages(2)
    assert pf4_two_fragment_stages() == tuple(reversed(pf4_two_fragment_stages()))


def test_convention_record_binds_parity_and_physical_bridge() -> None:
    record = pf4_convention_record()

    assert isinstance(record, PF4ConventionRecord)
    assert record.stage_table == pf4_two_fragment_stages()
    assert record.parity_identity == "S4(-t)=S4(t)^(-1)"
    assert record.formal_log == "log S4(t)=t(A+B)+t^5*L5+O(t^7)"
    assert record.physical_log == "log U4(t)=-i*t*H-i*t^5*E5+O(t^7)"
    assert record.physical_bridge == "E5=L5"
    assert physical_substitution_factor(1) == (0, -1)
    assert physical_substitution_factor(5) == (0, -1)
    assert verify_pf4_convention_record(record)


def test_resealed_convention_mutations_are_rejected() -> None:
    canonical = pf4_convention_record()
    noncentral = list(canonical.stage_table)
    noncentral[1] = CubicStage(1, noncentral[1].coefficient + 1)
    central = list(canonical.stage_table)
    central[5] = CubicStage(1, central[5].coefficient + 1)
    adjacent_swap = list(canonical.stage_table)
    adjacent_swap[0], adjacent_swap[1] = adjacent_swap[1], adjacent_swap[0]
    global_swap = tuple(
        CubicStage(1 - stage.fragment_index, stage.coefficient)
        for stage in canonical.stage_table
    )
    mutations = (
        replace(canonical, stage_table=tuple(noncentral), convention_digest=""),
        replace(canonical, stage_table=tuple(central), convention_digest=""),
        replace(canonical, stage_table=tuple(adjacent_swap), convention_digest=""),
        replace(canonical, stage_table=global_swap, convention_digest=""),
        replace(canonical, multiplication_order="right_to_left", convention_digest=""),
        replace(canonical, commutator_d="[B,[A,B]]", convention_digest=""),
        replace(canonical, physical_bridge="E5=-L5", convention_digest=""),
    )

    for mutation in mutations:
        forged = replace(
            mutation,
            convention_digest=convention_record_digest(mutation),
        )
        assert not verify_pf4_convention_record(forged)


def test_mapping_is_solved_without_an_expected_vector() -> None:
    assert not inspect.signature(derive_pf4_cyclic_mapping).parameters
    mapping = derive_pf4_cyclic_mapping()

    assert mapping.log_terms[1] == {
        (0,): Cubic.one(),
        (1,): Cubic.one(),
    }
    assert all(mapping.log_terms[degree] == {} for degree in (2, 3, 4))
    assert mapping.log_terms[5]
    assert mapping.residual == {}
    assert mapping.target_cyclic == mapping.reconstructed_cyclic
    assert mapping.coefficients == (
        mapping.g_c2,
        mapping.g_cd,
        mapping.g_d2,
    )
    assert mapping.candidate_comparison == compare_projectively(
        mapping.coefficients,
        (Fraction(1), Fraction(-4), Fraction(8, 3)),
    )
    assert mapping.legacy_comparison == compare_projectively(
        mapping.coefficients,
        (Fraction(1, 2), Fraction(14, 3), Fraction(4, 3)),
    )
    assert len(mapping.payload_sha256) == 64
    assert derive_pf4_cyclic_mapping() is mapping


def test_projective_comparison_is_exact_and_post_solve() -> None:
    gamma = Cubic(Fraction(2, 7), Fraction(-3, 11), Fraction(5, 13))
    candidate = (Fraction(1), Fraction(-4), Fraction(8, 3))
    solved = tuple(gamma * value for value in candidate)
    proportional = compare_projectively(solved, candidate)
    nonproportional = compare_projectively(
        solved,
        (Fraction(1, 2), Fraction(14, 3), Fraction(4, 3)),
    )

    assert proportional.status == "proportional"
    assert proportional.scale == gamma
    assert proportional.residual == (Cubic.zero(),) * 3
    assert nonproportional.status == "not_proportional"
    assert any(value != Cubic.zero() for value in nonproportional.residual)


def test_general_mapping_payload_has_reference_compatible_schema() -> None:
    payload = to_payload()

    assert tuple(payload) == (
        "schema_version",
        "kind",
        "conventions",
        "stages",
        "log_terms",
        "cyclic_map",
        "solved_coefficients",
        "residual",
        "projective_comparisons",
        "optional_gamma",
        "sign_certificate",
        "physical_bridge",
        "identity_status",
        "payload_sha256",
    )
    assert set(payload["solved_coefficients"]) == {"g_c2", "g_cd", "g_d2"}
    assert payload["residual"] == []
    assert payload["optional_gamma"] is not None
    assert payload["payload_sha256"] == payload_digest(payload)
    assert canonical_payload_bytes(payload).endswith(b"\n")
    assert parse_payload(payload) is derive_pf4_cyclic_mapping()

    for value in payload["solved_coefficients"].values():
        assert len(value) == 3
        assert all(len(pair) == 2 and pair[1] > 0 for pair in value)
    for records in payload["cyclic_map"].values():
        assert records == sorted(records, key=lambda record: record["word"])


def test_resealed_mapping_payload_mutation_is_rejected() -> None:
    forged = copy.deepcopy(to_payload())
    forged["solved_coefficients"]["g_cd"][0][0] += 1
    forged["payload_sha256"] = payload_digest(forged)

    with pytest.raises(ValueError, match="authoritative"):
        parse_payload(forged)


def _add(
    left: Mapping[Word, Fraction], right: Mapping[Word, Fraction]
) -> dict[Word, Fraction]:
    result = dict(left)
    for word, coefficient in right.items():
        updated = result.get(word, Fraction()) + coefficient
        if updated:
            result[word] = updated
        else:
            result.pop(word, None)
    return result


def _scale(
    polynomial: Mapping[Word, Fraction], scalar: Fraction | int
) -> dict[Word, Fraction]:
    factor = Fraction(scalar)
    return {
        word: coefficient * factor
        for word, coefficient in polynomial.items()
        if coefficient * factor
    }


def _multiply(
    left: Mapping[Word, Fraction], right: Mapping[Word, Fraction]
) -> dict[Word, Fraction]:
    result: dict[Word, Fraction] = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            word = left_word + right_word
            result = _add(
                result,
                {word: left_coefficient * right_coefficient},
            )
    return result


def _commutator(
    left: Mapping[Word, Fraction], right: Mapping[Word, Fraction]
) -> dict[Word, Fraction]:
    return _add(_multiply(left, right), _scale(_multiply(right, left), -1))


def _mutated_quadratic(
    c2: Fraction, cd: Fraction, d2: Fraction
) -> CubicWordPolynomial:
    a = {(0,): Fraction(1)}
    b = {(1,): Fraction(1)}
    c = _commutator(a, _commutator(a, b))
    d = _commutator(b, _commutator(b, a))
    rational = _add(
        _scale(_multiply(c, c), c2),
        _add(
            _scale(_multiply(c, d), cd),
            _scale(_multiply(d, d), d2),
        ),
    )
    gamma = pf4_suzuki_gamma()
    return {word: gamma * coefficient for word, coefficient in rational.items()}


def _evaluate_trace(
    polynomial: Mapping[Word, Cubic],
    a: np.ndarray,
    b: np.ndarray,
) -> Cubic:
    dimension = a.shape[0]
    identity = np.eye(dimension, dtype=object)
    generators = (a, b)
    result = Cubic.zero()
    for word, coefficient in polynomial.items():
        term = identity
        for letter in word:
            term = term @ generators[letter]
        trace = sum((Fraction(term[i, i]) for i in range(dimension)), Fraction())
        result += coefficient * trace
    return result


def test_pf4_suzuki_trace_identity_matches_frozen_cyclic_classes() -> None:
    left, right = pf4_suzuki_trace_polynomials()
    left_classes = cyclic_trace_classes(left)
    right_classes = cyclic_trace_classes(right)
    gamma = pf4_suzuki_gamma()
    expected_multiples = {
        (0, 0, 0, 0, 1, 1): Fraction(2),
        (0, 0, 0, 1, 0, 1): Fraction(-8),
        (0, 0, 0, 1, 1, 1): Fraction(-8),
        (0, 0, 1, 0, 0, 1): Fraction(6),
        (0, 0, 1, 0, 1, 1): Fraction(12),
        (0, 0, 1, 1, 0, 1): Fraction(12),
        (0, 0, 1, 1, 1, 1): Fraction(16, 3),
        (0, 1, 0, 1, 0, 1): Fraction(-16),
        (0, 1, 0, 1, 1, 1): Fraction(-64, 3),
        (0, 1, 1, 0, 1, 1): Fraction(16),
    }

    assert pf4_suzuki_gamma() == Cubic(
        Fraction(37, 900000),
        Fraction(313, 14400000),
        Fraction(29, 1800000),
    )
    assert left_classes == {
        word: gamma * multiple for word, multiple in expected_multiples.items()
    }
    assert right_classes == left_classes
    verify_pf4_suzuki_trace_identity()


@pytest.mark.parametrize(
    ("c2", "cd", "d2"),
    (
        (Fraction(2), Fraction(-4), Fraction(8, 3)),
        (Fraction(1), Fraction(-3), Fraction(8, 3)),
        (Fraction(1), Fraction(-4), Fraction(7, 3)),
    ),
)
def test_pf4_trace_identity_rejects_mutated_quadratic_coefficients(
    c2: Fraction,
    cd: Fraction,
    d2: Fraction,
) -> None:
    left, _ = pf4_suzuki_trace_polynomials()

    assert cyclic_trace_classes(_mutated_quadratic(c2, cd, d2)) != (
        cyclic_trace_classes(left)
    )


@pytest.mark.parametrize("dimension", (2, 3))
def test_pf4_trace_identity_on_exact_symmetric_matrices(dimension: int) -> None:
    rng = np.random.default_rng(12840 + dimension)
    raw_a = rng.integers(-3, 4, size=(dimension, dimension))
    raw_b = rng.integers(-3, 4, size=(dimension, dimension))
    a = np.asarray(raw_a + raw_a.T, dtype=object)
    b = np.asarray(raw_b + raw_b.T, dtype=object)
    left, right = pf4_suzuki_trace_polynomials()

    assert _evaluate_trace(left, a, b) == _evaluate_trace(right, a, b)


def test_pf4_trace_record_binds_corrected_core_and_gamma() -> None:
    core = pf4_trace_quadratic_form(6, 5, 3, Fraction(7, 11))
    expected = Fraction(7, 11) * (
        Fraction(6) - 4 * 5 + Fraction(8, 3) * 3
    )
    record = pf4_trace_identity_record()

    assert core == expected
    assert pf4_suzuki_trace_pairing(6, 0, 3) == (
        pf4_suzuki_gamma()
        * pf4_trace_quadratic_form(6, 0, 3, Fraction(1))
    )
    assert record.identity_status == "proved_exact_suzuki_pf4_free_trace_identity"
    assert record.exact_coefficients == (
        ("trace_c2", Fraction(1)),
        ("trace_cd", Fraction(-4)),
        ("trace_d2", Fraction(8, 3)),
        ("gamma_a0", Fraction(37, 900000)),
        ("gamma_a1", Fraction(313, 14400000)),
        ("gamma_a2", Fraction(29, 1800000)),
    )
    assert verify_identity_record(record, "pf4")
