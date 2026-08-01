"""Exact centered-polynomial PF4 gauge witnesses for the periodic TFIM.

The Hamiltonian convention is

``H = h * sum_i X_i + j * sum_i Z_i Z_(i+1)``.

For each frozen power ``p in (2, 3, 4)``, this module constructs

``W_p = H^p - a_p I - b_p H``

with exact normalized-trace constraints ``tau(W_p)=tau(W_p H)=0``.  Pairings
with the physical PF4 ``E5`` are evaluated by Pauli-coefficient overlap from
the authoritative word polynomial returned by
``derive_pf4_cyclic_mapping``.  No matrix materialization or norm claim is
made here.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product

from .algebra import PauliString, PauliSum, QComplex, commutator
from .cubic_field import Cubic
from .pf4_bch_mapping import derive_pf4_cyclic_mapping, pf4_suzuki_gamma


CENTERED_POLYNOMIAL_SEARCH_ORDER = (2, 3, 4)


@dataclass(frozen=True, slots=True)
class TFIMGaugeOrbitFormula:
    """Closed coefficients obtained from exact local-Pauli orbit counts."""

    length: int
    field: Fraction
    coupling: Fraction
    regime: str
    tau_h2: Fraction
    tau_h4: Fraction
    p3_hamiltonian_coefficient: Fraction
    tau_h_e5_over_gamma: Fraction
    tau_h3_e5_over_gamma: Fraction
    tau_w3_e5_over_gamma: Fraction
    nonzero_condition: str


@dataclass(frozen=True, slots=True)
class PauliOrbitPolynomialCertificate:
    """Parameter-free Pauli-word counts for one finite periodic ring.

    Monomial keys are ``(power_of_h, power_of_j)``.  Because the counter keeps
    these keys symbolic, equality with the closed expressions proves the
    identity for all exact ``h,j`` at this length; it is not a fit at selected
    coupling values.
    """

    length: int
    regime: str
    tau_h2_coefficients: tuple[tuple[tuple[int, int], Fraction], ...]
    tau_h4_coefficients: tuple[tuple[tuple[int, int], Fraction], ...]
    tau_h_e5_coefficients: tuple[tuple[tuple[int, int], Cubic], ...]
    tau_h3_e5_coefficients: tuple[tuple[tuple[int, int], Cubic], ...]
    matches_closed_formula: bool


@dataclass(frozen=True, slots=True)
class CenteredPolynomialWitness:
    """One exact ``H^p-aI-bH`` witness and its unnormalized dual pairing."""

    power: int
    identity_coefficient: Fraction
    hamiltonian_coefficient: Fraction
    operator: PauliSum
    tau_w: Fraction
    tau_wh: Fraction
    commutator_is_zero: bool
    tau_w_e5: Cubic


@dataclass(frozen=True, slots=True)
class TFIMGaugeWitnessResult:
    """Exact operator regression plus the separately derived orbit formula."""

    length: int
    field: Fraction
    coupling: Fraction
    tau_h2: Fraction
    tau_h4: Fraction
    tau_h_e5: Cubic
    tau_h3_e5: Cubic
    e5_term_count: int
    witnesses: tuple[
        CenteredPolynomialWitness,
        CenteredPolynomialWitness,
        CenteredPolynomialWitness,
    ]
    orbit_formula: TFIMGaugeOrbitFormula
    orbit_polynomial_certificate: PauliOrbitPolynomialCertificate

    @property
    def first_nonzero_power(self) -> int | None:
        return next(
            (
                witness.power
                for witness in self.witnesses
                if witness.tau_w_e5 != Cubic.zero()
            ),
            None,
        )


def _validate_inputs(length: int, h: Fraction, j: Fraction) -> None:
    if isinstance(length, bool) or not isinstance(length, int):
        raise TypeError("length must be an integer, not bool")
    if length < 4 or length % 2:
        raise ValueError("periodic TFIM length must be even and at least four")
    if not isinstance(h, Fraction) or not isinstance(j, Fraction):
        raise TypeError("h and j must be exact Fractions")
    if h == 0 and j == 0:
        raise ValueError("the zero Hamiltonian has no centered witness denominator")


def _tfim_split_operators(
    length: int,
    h: Fraction,
    j: Fraction,
) -> tuple[PauliSum, PauliSum]:
    field = PauliSum.zero()
    interaction = PauliSum.zero()
    for site in range(length):
        field += PauliSum.term(PauliString({site: "X"}), h)
        interaction += PauliSum.term(
            PauliString({site: "Z", (site + 1) % length: "Z"}), j
        )
    return field, interaction


def _normalized_trace(operator: PauliSum) -> Fraction:
    coefficient = operator.terms.get(PauliString(), QComplex())
    if coefficient.imag:
        raise ArithmeticError("normalized trace unexpectedly has an imaginary part")
    return coefficient.real


def _normalized_rational_overlap(left: PauliSum, right: PauliSum) -> Fraction:
    """Return ``tau(left^dagger right)`` by exact Pauli orthogonality."""

    if len(left.terms) <= len(right.terms):
        items = left.terms.items()
        lookup = right.terms
        conjugate_items = True
    else:
        items = right.terms.items()
        lookup = left.terms
        conjugate_items = False
    result = QComplex()
    for pauli, coefficient in items:
        if conjugate_items:
            result += coefficient.conjugate() * lookup.get(pauli, QComplex())
        else:
            result += lookup.get(pauli, QComplex()).conjugate() * coefficient
    if result.imag:
        raise ArithmeticError("rational Pauli overlap unexpectedly is not real")
    return result.real


def _normalized_product_trace(left: PauliSum, right: PauliSum) -> QComplex:
    """Return ``tau(left right)`` without materializing the product."""

    if len(left.terms) <= len(right.terms):
        items = left.terms.items()
        lookup = right.terms
    else:
        items = right.terms.items()
        lookup = left.terms
    return sum(
        (
            coefficient * lookup.get(pauli, QComplex())
            for pauli, coefficient in items
        ),
        QComplex(),
    )


def _word_operator(
    word: tuple[int, ...],
    field: PauliSum,
    interaction: PauliSum,
    cache: dict[tuple[int, ...], PauliSum],
) -> PauliSum:
    if word in cache:
        return cache[word]
    prefix = word[:-1]
    cache[word] = _word_operator(prefix, field, interaction, cache) * (
        field if word[-1] == 0 else interaction
    )
    return cache[word]


def _pf4_e5_pauli_coefficients(
    field: PauliSum,
    interaction: PauliSum,
) -> dict[PauliString, Cubic]:
    """Evaluate the exact PF4 degree-five word polynomial in the Pauli basis."""

    mapping = derive_pf4_cyclic_mapping()
    cache = {(): PauliSum.identity()}
    real_terms: dict[PauliString, Cubic] = {}
    imaginary_terms: dict[PauliString, Cubic] = {}
    for word, word_coefficient in mapping.log_terms[5].items():
        operator = _word_operator(word, field, interaction, cache)
        for pauli, coefficient in operator.terms.items():
            if coefficient.real:
                updated = real_terms.get(pauli, Cubic.zero()) + (
                    word_coefficient * coefficient.real
                )
                if updated == Cubic.zero():
                    real_terms.pop(pauli, None)
                else:
                    real_terms[pauli] = updated
            if coefficient.imag:
                updated = imaginary_terms.get(pauli, Cubic.zero()) + (
                    word_coefficient * coefficient.imag
                )
                if updated == Cubic.zero():
                    imaginary_terms.pop(pauli, None)
                else:
                    imaginary_terms[pauli] = updated
    if imaginary_terms:
        raise ArithmeticError("the physical PF4 E5 Pauli map is not Hermitian")
    return real_terms


def pauli_coefficient_overlap(
    rational_operator: PauliSum,
    cubic_operator: dict[PauliString, Cubic],
) -> Cubic:
    """Return normalized trace overlap in time linear in the smaller term map."""

    result = Cubic.zero()
    if len(rational_operator.terms) <= len(cubic_operator):
        for pauli, coefficient in rational_operator.terms.items():
            if coefficient.imag:
                raise ArithmeticError("witness Pauli coefficient is not real")
            result += cubic_operator.get(pauli, Cubic.zero()) * coefficient.real
    else:
        for pauli, coefficient in cubic_operator.items():
            rational = rational_operator.terms.get(pauli)
            if rational is None:
                continue
            if rational.imag:
                raise ArithmeticError("witness Pauli coefficient is not real")
            result += coefficient * rational.real
    return result


def _accumulate_rational_polynomial(
    polynomial: dict[tuple[int, int], Fraction],
    monomial: tuple[int, int],
    coefficient: Fraction,
) -> None:
    updated = polynomial.get(monomial, Fraction()) + coefficient
    if updated:
        polynomial[monomial] = updated
    else:
        polynomial.pop(monomial, None)


def _accumulate_cubic_polynomial(
    polynomial: dict[tuple[int, int], Cubic],
    monomial: tuple[int, int],
    coefficient: Cubic,
) -> None:
    updated = polynomial.get(monomial, Cubic.zero()) + coefficient
    if updated != Cubic.zero():
        polynomial[monomial] = updated
    else:
        polynomial.pop(monomial, None)


def _hamiltonian_trace_polynomial(
    degree: int,
    field: PauliSum,
    interaction: PauliSum,
    cache: dict[tuple[int, ...], PauliSum],
) -> dict[tuple[int, int], Fraction]:
    result: dict[tuple[int, int], Fraction] = {}
    imaginary: dict[tuple[int, int], Fraction] = {}
    for word in product((0, 1), repeat=degree):
        coefficient = _normalized_trace_qcomplex(
            _word_operator(word, field, interaction, cache)
        )
        monomial = (word.count(0), word.count(1))
        _accumulate_rational_polynomial(result, monomial, coefficient.real)
        _accumulate_rational_polynomial(imaginary, monomial, coefficient.imag)
    if imaginary:
        raise ArithmeticError("Hamiltonian trace polynomial is not real")
    return result


def _normalized_trace_qcomplex(operator: PauliSum) -> QComplex:
    return operator.terms.get(PauliString(), QComplex())


def _e5_trace_polynomial(
    prefix_degree: int,
    field: PauliSum,
    interaction: PauliSum,
    cache: dict[tuple[int, ...], PauliSum],
) -> dict[tuple[int, int], Cubic]:
    result: dict[tuple[int, int], Cubic] = {}
    imaginary: dict[tuple[int, int], Cubic] = {}
    log_five = derive_pf4_cyclic_mapping().log_terms[5]
    for prefix in product((0, 1), repeat=prefix_degree):
        prefix_operator = _word_operator(prefix, field, interaction, cache)
        for word, cubic_coefficient in log_five.items():
            trace = _normalized_product_trace(
                prefix_operator,
                _word_operator(word, field, interaction, cache),
            )
            monomial = (
                prefix.count(0) + word.count(0),
                prefix.count(1) + word.count(1),
            )
            _accumulate_cubic_polynomial(
                result, monomial, cubic_coefficient * trace.real
            )
            _accumulate_cubic_polynomial(
                imaginary, monomial, cubic_coefficient * trace.imag
            )
    if imaginary:
        raise ArithmeticError("PF4 trace polynomial is not real")
    return result


def _sorted_rational_polynomial(
    polynomial: dict[tuple[int, int], Fraction],
) -> tuple[tuple[tuple[int, int], Fraction], ...]:
    return tuple(sorted(polynomial.items()))


def _sorted_cubic_polynomial(
    polynomial: dict[tuple[int, int], Cubic],
) -> tuple[tuple[tuple[int, int], Cubic], ...]:
    return tuple(sorted(polynomial.items()))


@lru_cache(maxsize=None)
def count_tfim_pf4_pauli_orbits(length: int) -> PauliOrbitPolynomialCertificate:
    """Count the relevant normalized traces as exact bivariate polynomials.

    This enumerates Pauli coefficients on the requested finite ring while
    retaining the number of field and bond factors as symbolic monomial
    degrees.  Thus the result checks all ``h,j`` simultaneously at that L.
    """

    _validate_inputs(length, Fraction(1), Fraction(1))
    field, interaction = _tfim_split_operators(
        length, Fraction(1), Fraction(1)
    )
    cache: dict[tuple[int, ...], PauliSum] = {(): PauliSum.identity()}
    observed_h2 = _hamiltonian_trace_polynomial(2, field, interaction, cache)
    observed_h4 = _hamiltonian_trace_polynomial(4, field, interaction, cache)
    observed_h_e5 = _e5_trace_polynomial(1, field, interaction, cache)
    observed_h3_e5 = _e5_trace_polynomial(3, field, interaction, cache)
    gamma = pf4_suzuki_gamma()

    expected_h2 = {(2, 0): Fraction(length), (0, 2): Fraction(length)}
    expected_h_e5 = {
        (4, 2): gamma * (128 * length),
        (2, 4): gamma * Fraction(1024 * length, 3),
    }
    if length == 4:
        regime = "l4_wraparound"
        expected_h4 = {
            (4, 0): Fraction(40),
            (2, 2): Fraction(64),
            (0, 4): Fraction(64),
        }
        expected_h3_e5 = {
            (6, 2): gamma * 5120,
            (4, 4): gamma * Fraction(40960, 3),
            (2, 6): gamma * Fraction(65536, 3),
        }
    else:
        regime = "generic_even_l_ge_6"
        expected_h4 = {
            (4, 0): Fraction(length * (3 * length - 2)),
            (2, 2): Fraction(length * (6 * length - 8)),
            (0, 4): Fraction(length * (3 * length - 2)),
        }
        expected_h3_e5 = {
            (6, 2): gamma * (128 * length * (3 * length - 2)),
            (4, 4): gamma * Fraction(
                128 * length * (33 * length - 34), 3
            ),
            (2, 6): gamma
            * Fraction(1024 * length * (3 * length - 2), 3),
        }

    matches = (
        observed_h2 == expected_h2
        and observed_h4 == expected_h4
        and observed_h_e5 == expected_h_e5
        and observed_h3_e5 == expected_h3_e5
    )
    if not matches:
        raise ArithmeticError("finite-ring Pauli orbit count disagrees with formula")
    return PauliOrbitPolynomialCertificate(
        length=length,
        regime=regime,
        tau_h2_coefficients=_sorted_rational_polynomial(observed_h2),
        tau_h4_coefficients=_sorted_rational_polynomial(observed_h4),
        tau_h_e5_coefficients=_sorted_cubic_polynomial(observed_h_e5),
        tau_h3_e5_coefficients=_sorted_cubic_polynomial(observed_h3_e5),
        matches_closed_formula=True,
    )


def tfim_gauge_orbit_formula(
    length: int,
    h: Fraction,
    j: Fraction,
) -> TFIMGaugeOrbitFormula:
    """Return the separately derived finite-ring local-Pauli orbit identity.

    The coefficients count identity-producing ordered Pauli words after one
    local generator is fixed and its translations are restored.  The L=4
    ring is separate because opposite range-two supports wrap and merge.  For
    every even L>=6 the local orbit types are the generic ones.  These are
    algebraic counting identities, not an interpolation through L=6,8,10.
    ``verify_tfim_gauge_orbit_formula`` additionally compares them with a full
    exact finite-ring Pauli enumeration for each constructed instance.
    """

    _validate_inputs(length, h, j)
    x = h**2
    y = j**2
    tau_h2 = length * (x + y)
    tau_h_e5_over_gamma = 128 * length * x * y * (x + Fraction(8, 3) * y)
    if length == 4:
        tau_h4 = 40 * x**2 + 64 * x * y + 64 * y**2
        tau_h3_e5_over_gamma = (
            5120 * x**3 * y
            + Fraction(40960, 3) * x**2 * y**2
            + Fraction(65536, 3) * x * y**3
        )
        simplified = (
            Fraction(1024, 3)
            * h**4
            * j**4
            * (16 * j**2 - 9 * h**2)
            / (x + y)
        )
        regime = "l4_wraparound"
        nonzero_condition = "h*j!=0 and 9*h^2!=16*j^2"
    else:
        tau_h4 = length * (
            (3 * length - 2) * (x**2 + y**2)
            + (6 * length - 8) * x * y
        )
        tau_h3_e5_over_gamma = 128 * length * x * y * (
            (3 * length - 2) * x**2
            + Fraction(33 * length - 34, 3) * x * y
            + Fraction(8, 3) * (3 * length - 2) * y**2
        )
        simplified = (
            Fraction(2560, 3) * length * h**4 * j**6 / (x + y)
        )
        regime = "generic_even_l_ge_6"
        nonzero_condition = "h*j!=0"

    b3 = tau_h4 / tau_h2
    derived = tau_h3_e5_over_gamma - b3 * tau_h_e5_over_gamma
    if derived != simplified:
        raise ArithmeticError("centered orbit formula failed exact simplification")
    return TFIMGaugeOrbitFormula(
        length=length,
        field=h,
        coupling=j,
        regime=regime,
        tau_h2=tau_h2,
        tau_h4=tau_h4,
        p3_hamiltonian_coefficient=b3,
        tau_h_e5_over_gamma=tau_h_e5_over_gamma,
        tau_h3_e5_over_gamma=tau_h3_e5_over_gamma,
        tau_w3_e5_over_gamma=simplified,
        nonzero_condition=nonzero_condition,
    )


def _centered_witness(
    power: int,
    powers: dict[int, PauliSum],
    hamiltonian: PauliSum,
    e5_terms: dict[PauliString, Cubic],
    tau_h: Fraction,
    tau_h2: Fraction,
) -> CenteredPolynomialWitness:
    tau_hp = _normalized_trace(powers[power])
    tau_hp1 = (
        _normalized_trace(powers[power + 1])
        if power + 1 in powers
        else _normalized_rational_overlap(powers[2], powers[3])
    )
    denominator = tau_h2 - tau_h**2
    if denominator <= 0:
        raise ValueError("centered Hamiltonian variance must be positive")
    b = (tau_hp1 - tau_hp * tau_h) / denominator
    a = tau_hp - b * tau_h
    witness = (
        powers[power]
        - PauliSum.identity(a)
        - hamiltonian.scale(b)
    )
    tau_w = _normalized_trace(witness)
    tau_wh = _normalized_rational_overlap(witness, hamiltonian)
    commutator_is_zero = not commutator(witness, hamiltonian)
    if tau_w != 0 or tau_wh != 0 or not commutator_is_zero:
        raise ArithmeticError("centered polynomial witness identities failed")
    return CenteredPolynomialWitness(
        power=power,
        identity_coefficient=a,
        hamiltonian_coefficient=b,
        operator=witness,
        tau_w=tau_w,
        tau_wh=tau_wh,
        commutator_is_zero=commutator_is_zero,
        tau_w_e5=pauli_coefficient_overlap(witness, e5_terms),
    )


def build_tfim_centered_polynomial_witnesses(
    length: int,
    h: Fraction,
    j: Fraction,
) -> TFIMGaugeWitnessResult:
    """Build and operator-verify the frozen p=2,3,4 witness search."""

    _validate_inputs(length, h, j)
    field, interaction = _tfim_split_operators(length, h, j)
    hamiltonian = field + interaction
    powers = {0: PauliSum.identity(), 1: hamiltonian}
    for power in range(2, 5):
        powers[power] = powers[power - 1] * hamiltonian
    tau_h = _normalized_trace(hamiltonian)
    tau_h2 = _normalized_trace(powers[2])
    tau_h4 = _normalized_trace(powers[4])
    e5_terms = _pf4_e5_pauli_coefficients(field, interaction)
    witnesses = tuple(
        _centered_witness(
            power,
            powers,
            hamiltonian,
            e5_terms,
            tau_h,
            tau_h2,
        )
        for power in CENTERED_POLYNOMIAL_SEARCH_ORDER
    )
    result = TFIMGaugeWitnessResult(
        length=length,
        field=h,
        coupling=j,
        tau_h2=tau_h2,
        tau_h4=tau_h4,
        tau_h_e5=pauli_coefficient_overlap(hamiltonian, e5_terms),
        tau_h3_e5=pauli_coefficient_overlap(powers[3], e5_terms),
        e5_term_count=len(e5_terms),
        witnesses=(witnesses[0], witnesses[1], witnesses[2]),
        orbit_formula=tfim_gauge_orbit_formula(length, h, j),
        orbit_polynomial_certificate=count_tfim_pf4_pauli_orbits(length),
    )
    verify_tfim_gauge_orbit_formula(result)
    return result


def verify_tfim_gauge_orbit_formula(result: TFIMGaugeWitnessResult) -> None:
    """Cross-check derived orbit identities against the full exact operators."""

    if not isinstance(result, TFIMGaugeWitnessResult):
        raise TypeError("result must be a TFIMGaugeWitnessResult")
    formula = tfim_gauge_orbit_formula(result.length, result.field, result.coupling)
    if result.orbit_formula != formula:
        raise ValueError("stored TFIM orbit formula mismatch")
    orbit_certificate = count_tfim_pf4_pauli_orbits(result.length)
    if (
        result.orbit_polynomial_certificate != orbit_certificate
        or not orbit_certificate.matches_closed_formula
        or orbit_certificate.regime != formula.regime
    ):
        raise ValueError("Pauli orbit polynomial certificate mismatch")
    if tuple(witness.power for witness in result.witnesses) != (
        CENTERED_POLYNOMIAL_SEARCH_ORDER
    ):
        raise ValueError("centered-polynomial search order mismatch")
    if result.tau_h2 != formula.tau_h2 or result.tau_h4 != formula.tau_h4:
        raise ValueError("Hamiltonian moment orbit count mismatch")
    gamma = pf4_suzuki_gamma()
    if result.tau_h_e5 != gamma * formula.tau_h_e5_over_gamma:
        raise ValueError("tau(H E5) orbit count mismatch")
    if result.tau_h3_e5 != gamma * formula.tau_h3_e5_over_gamma:
        raise ValueError("tau(H^3 E5) orbit count mismatch")
    if any(
        witness.tau_w != 0
        or witness.tau_wh != 0
        or not witness.commutator_is_zero
        for witness in result.witnesses
    ):
        raise ValueError("centered witness operator identity mismatch")
    p2, p3, p4 = result.witnesses
    if p3.hamiltonian_coefficient != formula.p3_hamiltonian_coefficient:
        raise ValueError("cubic witness centering coefficient mismatch")
    if p2.tau_w_e5 != Cubic.zero() or p4.tau_w_e5 != Cubic.zero():
        raise ValueError("frozen p=2 or p=4 zero-pairing regression failed")
    if p3.tau_w_e5 != gamma * formula.tau_w3_e5_over_gamma:
        raise ValueError("cubic centered-witness pairing mismatch")
    expected_first = 3 if formula.tau_w3_e5_over_gamma else None
    if result.first_nonzero_power != expected_first:
        raise ValueError("first nonzero centered-polynomial witness mismatch")
