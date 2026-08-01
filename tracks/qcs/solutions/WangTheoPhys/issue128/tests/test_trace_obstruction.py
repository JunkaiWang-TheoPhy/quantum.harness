from __future__ import annotations

import json
from dataclasses import replace
from fractions import Fraction
from hashlib import sha256

import numpy as np
import pytest
from trottercert.processed_kernels import RationalStage
from trottercert.rational_words import (
    commutator_polynomial,
    multiply_polynomials,
    rational_formula_log_series,
    scale_polynomial,
)

from trottercert.pf4_bch_mapping import pf4_suzuki_gamma
from trottercert.trace_obstruction import (
    ObstructionException,
    classify_pf4_exception,
    pf2_trace_identity_record,
    pf2_trace_obstruction,
    pf4_suzuki_trace_pairing,
    pf4_trace_identity_record,
    pf4_trace_quadratic_form,
    verify_identity_digest,
    verify_identity_record,
)


def _strang_degree_three():
    stages = (
        RationalStage(0, Fraction(1, 2)),
        RationalStage(1, Fraction(1)),
        RationalStage(0, Fraction(1, 2)),
    )
    return rational_formula_log_series(stages, 3)[3]


def _cyclic_trace_classes(polynomial):
    classes = {}
    for word, coefficient in polynomial.items():
        representative = min(
            word[offset:] + word[:offset] for offset in range(len(word))
        )
        classes[representative] = (
            classes.get(representative, Fraction()) + coefficient
        )
    return {word: value for word, value in classes.items() if value}


def _evaluate_word_polynomial(polynomial, a, b):
    generators = (a, b)
    result = np.zeros_like(a)
    identity = np.eye(a.shape[0], dtype=a.dtype)
    for word, coefficient in polynomial.items():
        term = identity
        for letter in word:
            term = term @ generators[letter]
        result = result + coefficient * term
    return result


def test_pf2_free_trace_identity_is_exact() -> None:
    hamiltonian = {(0,): Fraction(1), (1,): Fraction(1)}
    degree_three = _strang_degree_three()
    bch_pairing = multiply_polynomials(hamiltonian, degree_three)

    a = {(0,): Fraction(1)}
    b = {(1,): Fraction(1)}
    commutator = commutator_polynomial(a, b)
    hilbert_schmidt_square = scale_polynomial(
        multiply_polynomials(commutator, commutator), -1
    )

    assert _cyclic_trace_classes(bch_pairing) == _cyclic_trace_classes(
        scale_polynomial(hilbert_schmidt_square, Fraction(1, 24))
    )


@pytest.mark.parametrize("dimension", [2, 3])
def test_pf2_seeded_hermitian_identity(dimension: int) -> None:
    rng = np.random.default_rng(1210 + dimension)
    raw_a = rng.integers(-3, 4, size=(dimension, dimension))
    raw_b = rng.integers(-3, 4, size=(dimension, dimension))
    a = np.asarray(raw_a + raw_a.T, dtype=object)
    b = np.asarray(raw_b + raw_b.T, dtype=object)

    degree_three = _strang_degree_three()
    defect = _evaluate_word_polynomial(degree_three, a, b)
    direct_pairing = sum(((a + b) @ defect)[index, index] for index in range(dimension))

    assert pf2_trace_obstruction(a, b) == direct_pairing
    assert isinstance(pf2_trace_obstruction(a, b), Fraction)


def test_pf2_numeric_complex_hermitian_identity() -> None:
    rng = np.random.default_rng(2402)
    matrices = []
    for _ in range(2):
        raw = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        matrices.append(raw + raw.conj().T)
    a, b = matrices

    defect = _evaluate_word_polynomial(_strang_degree_three(), a, b)
    direct_pairing = np.trace((a + b) @ defect).real

    assert pf2_trace_obstruction(a, b) == pytest.approx(
        direct_pairing, rel=2e-13, abs=2e-13
    )


def test_pf2_commuting_matrices_have_zero_obstruction() -> None:
    a = np.asarray([[Fraction(1), Fraction(0)], [Fraction(0), Fraction(2)]], dtype=object)
    b = np.asarray([[Fraction(3), Fraction(0)], [Fraction(0), Fraction(-1)]], dtype=object)

    assert pf2_trace_obstruction(a, b) == Fraction(0)


def test_pf2_rejects_nonhermitian_input() -> None:
    a = np.asarray([[0.0, 1.0], [0.0, 0.0]])
    b = np.eye(2)

    with pytest.raises(ValueError, match="Hermitian"):
        pf2_trace_obstruction(a, b)


@pytest.mark.parametrize(
    "a",
    [
        np.asarray([[True, False], [False, True]]),
        np.asarray([[True, Fraction(0)], [Fraction(0), Fraction(1)]], dtype=object),
        np.asarray([[np.bool_(False), 0.0], [0.0, 1.0]], dtype=object),
    ],
)
def test_pf2_rejects_any_boolean_matrix_entry(a: np.ndarray) -> None:
    with pytest.raises(TypeError, match="bool"):
        pf2_trace_obstruction(a, np.eye(2, dtype=int))


def test_pf4_quadratic_form_preserves_exact_derived_coefficients() -> None:
    value = pf4_trace_quadratic_form(
        trace_c2=Fraction(6),
        trace_cd=Fraction(5),
        trace_d2=Fraction(3),
        gamma=Fraction(7, 11),
    )
    expected = Fraction(7, 11) * (Fraction(6) - 4 * 5 + Fraction(8, 3) * 3)

    assert value == expected
    assert isinstance(value, Fraction)


def test_pf4_exact_suzuki_pairing_includes_cubic_gamma() -> None:
    core = pf4_trace_quadratic_form(6, 0, 3, 1)

    assert pf4_suzuki_trace_pairing(6, 0, 3) == pf4_suzuki_gamma() * core


def test_pf4_quadratic_form_refuses_inexact_inputs() -> None:
    with pytest.raises(TypeError, match="exact rational"):
        pf4_trace_quadratic_form(6.0, Fraction(5), Fraction(3), Fraction(7, 11))


def test_exception_categories_are_explicit_and_distinct() -> None:
    assert {member.value for member in ObstructionException} == {
        "commuting",
        "symmetry_protected_positive",
        "mixed_term_indefinite",
        "exactly_correctable",
    }
    assert classify_pf4_exception(0, 0, 0, 1, commuting=True) is ObstructionException.COMMUTING
    assert classify_pf4_exception(
        2, 0, 3, 1, symmetry_forces_mixed_zero=True
    ) is ObstructionException.SYMMETRY_PROTECTED_POSITIVE
    assert classify_pf4_exception(2, -5, 3, 1) is ObstructionException.MIXED_TERM_INDEFINITE
    assert classify_pf4_exception(
        0, 0, 0, 1, exactly_correctable=True
    ) is ObstructionException.EXACTLY_CORRECTABLE


def test_exact_correctability_is_not_inferred_from_a_zero_scalar() -> None:
    assert classify_pf4_exception(0, 0, 0, 1) is ObstructionException.MIXED_TERM_INDEFINITE
    with pytest.raises(ValueError, match="nonzero"):
        classify_pf4_exception(1, 0, 0, 1, exactly_correctable=True)


@pytest.mark.parametrize(
    ("commuting", "symmetry", "correctable"),
    [
        (True, True, False),
        (True, False, True),
        (False, True, True),
        (True, True, True),
    ],
)
def test_exception_premises_are_mutually_exclusive(
    commuting: bool,
    symmetry: bool,
    correctable: bool,
) -> None:
    with pytest.raises(ValueError, match="mutually exclusive"):
        classify_pf4_exception(
            0,
            0,
            0,
            1,
            commuting=commuting,
            symmetry_forces_mixed_zero=symmetry,
            exactly_correctable=correctable,
        )


def _attacker_digest(record) -> str:
    payload = {
        "identity_kind": record.identity_kind,
        "formula": record.formula,
        "assumptions": list(record.assumptions),
        "exact_coefficients": [
            [name, str(value)] for name, value in record.exact_coefficients
        ],
        "exception_conditions": [
            [name, condition] for name, condition in record.exception_conditions
        ],
        "identity_status": record.identity_status,
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return sha256(encoded).hexdigest()


def test_trace_identity_records_are_digest_bound_and_auditable() -> None:
    pf2 = pf2_trace_identity_record()
    pf4 = pf4_trace_identity_record()

    for record in (pf2, pf4):
        assert record.formula
        assert record.assumptions
        assert record.exact_coefficients
        assert record.exception_conditions
        assert len(record.identity_digest) == 64
        assert verify_identity_digest(record)

    assert verify_identity_record(pf2, "pf2")
    assert verify_identity_record(pf4, "pf4")
    assert not verify_identity_record(pf2, "pf4")
    assert pf2.identity_status == "proved_exact_free_trace_identity"
    assert pf4.identity_status == "proved_exact_suzuki_pf4_free_trace_identity"
    assert pf4.exact_coefficients == (
        ("trace_c2", Fraction(1)),
        ("trace_cd", Fraction(-4)),
        ("trace_d2", Fraction(8, 3)),
        ("gamma_a0", Fraction(37, 900000)),
        ("gamma_a1", Fraction(313, 14400000)),
        ("gamma_a2", Fraction(29, 1800000)),
    )
    assert not verify_identity_digest(replace(pf2, formula="tampered"))


def test_pf4_cannot_be_forged_as_proved_by_recomputing_digest() -> None:
    original = pf4_trace_identity_record()
    altered = replace(
        original,
        formula="gamma*((1/2)Tr(C^2) + (14/3)Tr(CD) + (4/3)Tr(D^2))",
        identity_status="algebraic_form_only_unverified_bch_mapping",
        identity_digest="",
    )
    forged = replace(altered, identity_digest=_attacker_digest(altered))

    assert not verify_identity_digest(forged)
    assert not verify_identity_record(forged, "pf4")
