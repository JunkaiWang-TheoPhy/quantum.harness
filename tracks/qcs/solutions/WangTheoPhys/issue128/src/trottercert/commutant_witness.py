from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction

from .cubic_field import Cubic
from .lattice import SquareLattice
from .local_commutators import (
    CoordinateRegistry,
    SymplecticPauli,
    _iter_set_bits,
    _symplectic_anticommutes,
    _symplectic_product_phase,
)

RationalPauliTerms = dict[SymplecticPauli, Fraction]
CubicPauliTerms = dict[SymplecticPauli, Cubic]


def _density_coordinates(
    registry: CoordinateRegistry,
    pauli: SymplecticPauli,
) -> tuple[tuple[int, int, bool, bool], ...]:
    x_mask, z_mask = pauli
    result: list[tuple[int, int, bool, bool]] = []
    for site in _iter_set_bits(x_mask | z_mask):
        x, y = registry.coordinate(site)
        result.append((x, y, bool(x_mask & (1 << site)), bool(z_mask & (1 << site))))
    return tuple(result)


def _translate_density_pauli(
    coordinates: tuple[tuple[int, int, bool, bool], ...],
    lattice: SquareLattice,
    dx: int,
    dy: int,
) -> SymplecticPauli:
    x_mask = z_mask = 0
    occupied: set[int] = set()
    for x, y, has_x, has_z in coordinates:
        site = lattice.site(x + dx, y + dy)
        if site in occupied:
            raise ValueError("density coordinates alias on the periodic lattice")
        occupied.add(site)
        bit = 1 << site
        if has_x:
            x_mask |= bit
        if has_z:
            z_mask |= bit
    return x_mask, z_mask


def _translation_offsets(
    lattice: SquareLattice,
    unit_cell: tuple[int, int],
) -> tuple[tuple[int, int], ...]:
    step_x, step_y = unit_cell
    if step_x < 1 or step_y < 1:
        raise ValueError("unit-cell dimensions must be positive")
    if lattice.length % step_x or lattice.length % step_y:
        raise ValueError("unit cell must tile the periodic lattice")
    return tuple(
        (dx, dy)
        for dy in range(0, lattice.length, step_y)
        for dx in range(0, lattice.length, step_x)
    )


def lift_cubic_density(
    registry: CoordinateRegistry,
    density: Mapping[SymplecticPauli, Cubic],
    lattice: SquareLattice,
    *,
    unit_cell: tuple[int, int] = (2, 2),
) -> CubicPauliTerms:
    """Lift one canonical local density through all unit-cell translations."""

    result: CubicPauliTerms = {}
    offsets = _translation_offsets(lattice, unit_cell)
    for pauli, coefficient in density.items():
        if not isinstance(coefficient, Cubic):
            raise ValueError("density coefficients must be exact Cubic values")
        coordinates = _density_coordinates(registry, pauli)
        for dx, dy in offsets:
            translated = _translate_density_pauli(coordinates, lattice, dx, dy)
            updated = result.get(translated, Cubic.zero()) + coefficient
            if updated == Cubic.zero():
                result.pop(translated, None)
            else:
                result[translated] = updated
    return result


def pair_lifted_cubic_density(
    registry: CoordinateRegistry,
    density: Mapping[SymplecticPauli, Cubic],
    rational_operator: Mapping[SymplecticPauli, Fraction],
    lattice: SquareLattice,
    *,
    unit_cell: tuple[int, int] = (2, 2),
) -> Cubic:
    """Return normalized-trace pairing without materializing the lifted map."""

    result = Cubic.zero()
    offsets = _translation_offsets(lattice, unit_cell)
    for pauli, coefficient in density.items():
        if not isinstance(coefficient, Cubic):
            raise ValueError("density coefficients must be exact Cubic values")
        coordinates = _density_coordinates(registry, pauli)
        for dx, dy in offsets:
            translated = _translate_density_pauli(coordinates, lattice, dx, dy)
            rational = rational_operator.get(translated)
            if rational is not None:
                result += coefficient * rational
    return result


def heisenberg_symplectic_terms(lattice: SquareLattice) -> RationalPauliTerms:
    """Return H=sum_<uv>(XX+YY+ZZ)/4 in bit-packed Pauli coordinates."""

    result: RationalPauliTerms = {}
    for first, second in lattice.bonds():
        sites = (1 << first) | (1 << second)
        for pauli in ((sites, 0), (sites, sites), (0, sites)):
            result[pauli] = Fraction(1, 4)
    return result


def square_real_pauli_terms(
    operator: Mapping[SymplecticPauli, Fraction],
) -> RationalPauliTerms:
    """Square a real Hermitian Pauli map, cancelling anticommuting pairs."""

    items = list(operator.items())
    result: RationalPauliTerms = {}
    for index, (left, left_coefficient) in enumerate(items):
        if not isinstance(left_coefficient, Fraction):
            left_coefficient = Fraction(left_coefficient)
        result[(0, 0)] = result.get((0, 0), Fraction()) + left_coefficient**2
        for right, right_coefficient in items[index + 1 :]:
            if _symplectic_anticommutes(left, right):
                continue
            phase, product = _symplectic_product_phase(left, right)
            if phase not in (0, 2):
                raise ArithmeticError("commuting Hermitian Paulis had non-real phase")
            sign = 1 if phase == 0 else -1
            updated = result.get(product, Fraction()) + (
                2 * sign * left_coefficient * right_coefficient
            )
            if updated:
                result[product] = updated
            else:
                result.pop(product, None)
    return result


def rational_pauli_pairing(
    left: Mapping[SymplecticPauli, Fraction],
    right: Mapping[SymplecticPauli, Fraction],
) -> Fraction:
    """Return normalized Hilbert--Schmidt pairing by Pauli orthogonality."""

    if len(left) > len(right):
        left, right = right, left
    return sum(
        (coefficient * right.get(pauli, Fraction()) for pauli, coefficient in left.items()),
        Fraction(),
    )


def _rational_json(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [_rational_json(value.a0), _rational_json(value.a1), _rational_json(value.a2)]


def _parse_rational(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(not isinstance(entry, int) or isinstance(entry, bool) for entry in value)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if _rational_json(result) != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _parse_cubic(value: object, field: str) -> Cubic:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must contain three cubic coordinates")
    return Cubic(*(_parse_rational(entry, field) for entry in value))


def build_quadratic_witness_payload(
    *,
    length: int,
    tau_h2: Fraction,
    tau_h3: Fraction,
    tau_h_e5: Cubic,
    tau_h2_e5: Cubic,
    matching_reconstruction_terms: int,
) -> dict[str, object]:
    """Build the algebraic certificate from independently computed moments."""

    tau_h2 = Fraction(tau_h2)
    tau_h3 = Fraction(tau_h3)
    if tau_h2 == 0:
        raise ValueError("tau(H^2) must be nonzero")
    expected_terms = 6 * length * length
    if matching_reconstruction_terms != expected_terms:
        raise ValueError("matching densities did not reconstruct the Hamiltonian")
    h_coefficient = -tau_h3 / tau_h2
    identity_coefficient = -tau_h2
    tau_w_e5 = tau_h2_e5 + tau_h_e5 * h_coefficient
    if tau_w_e5 == Cubic.zero():
        raise ValueError("quadratic witness pairing must be nonzero")

    return {
        "schema_version": 1,
        "kind": "issue128_quadratic_commutant_witness",
        "model": {
            "hamiltonian": "sum_<uv>(XX+YY+ZZ)/4",
            "length": length,
            "boundary": "periodic",
            "normalized_trace": "tau(A)=Tr(A)/2^(L^2)",
        },
        "formula": {
            "id": "five_copy_fourth_order_suzuki",
            "fragment_order": [0, 1, 2, 3],
            "alpha_minimal_polynomial": "alpha^3-4",
        },
        "lifting": {
            "unit_cell": [2, 2],
            "translations": (length // 2) ** 2,
            "matching_reconstruction_terms": matching_reconstruction_terms,
            "matches_full_hamiltonian": True,
        },
        "moments": {
            "tau_h": [0, 1],
            "tau_h2": _rational_json(tau_h2),
            "tau_h3": _rational_json(tau_h3),
        },
        "witness": {
            "expression": "H^2 + identity_coefficient I + hamiltonian_coefficient H",
            "h2_coefficient": [1, 1],
            "identity_coefficient": _rational_json(identity_coefficient),
            "hamiltonian_coefficient": _rational_json(h_coefficient),
        },
        "proof": {
            "commutes_with_hamiltonian": True,
            "commutation_reason": "W is a polynomial in H",
            "orthogonal_to_identity": True,
            "orthogonal_to_hamiltonian": True,
        },
        "pairings": {
            "basis": ["1", "alpha", "alpha^2"],
            "tau_h_e5": _cubic_json(tau_h_e5),
            "tau_h2_e5": _cubic_json(tau_h2_e5),
            "tau_w_e5": _cubic_json(tau_w_e5),
            "tau_w_e5_nonzero": True,
        },
        "claim": {
            "quotient": "image(i ad_H) + span(I,H)",
            "leading_order_status": "no_go",
            "statement": "E5 is not in image(i ad_H) + span(I,H)",
            "finite_step_status": "inconclusive",
            "finite_step_missing": "certified local-log branch and all-order remainder",
        },
        "hpc_authorized": False,
    }


def verify_quadratic_witness_payload(payload: Mapping[str, object]) -> None:
    """Verify all lightweight algebraic obligations in a stored certificate."""

    try:
        model = payload["model"]
        lifting = payload["lifting"]
        moments = payload["moments"]
        pairings = payload["pairings"]
        claim = payload["claim"]
        witness = payload["witness"]
        proof = payload["proof"]
    except KeyError as exc:
        raise ValueError(f"certificate field is missing: {exc.args[0]}") from exc
    if not all(isinstance(section, Mapping) for section in (model, lifting, moments, pairings, claim, witness, proof)):
        raise ValueError("certificate sections must be mappings")
    length = model.get("length")
    if not isinstance(length, int) or isinstance(length, bool) or length < 4 or length % 2:
        raise ValueError("model length must be an even integer at least four")
    if lifting.get("matches_full_hamiltonian") is not True:
        raise ValueError("Hamiltonian reconstruction is not certified")
    if lifting.get("matching_reconstruction_terms") != 6 * length * length:
        raise ValueError("Hamiltonian reconstruction term count mismatch")
    if lifting.get("translations") != (length // 2) ** 2:
        raise ValueError("translation count mismatch")

    tau_h = _parse_rational(moments.get("tau_h"), "tau_h")
    tau_h2 = _parse_rational(moments.get("tau_h2"), "tau_h2")
    tau_h3 = _parse_rational(moments.get("tau_h3"), "tau_h3")
    identity_coefficient = _parse_rational(
        witness.get("identity_coefficient"), "identity coefficient"
    )
    h_coefficient = _parse_rational(
        witness.get("hamiltonian_coefficient"), "Hamiltonian coefficient"
    )
    if tau_h != 0 or identity_coefficient != -tau_h2:
        raise ValueError("witness is not orthogonal to identity")
    if tau_h3 + h_coefficient * tau_h2 != 0:
        raise ValueError("witness is not orthogonal to Hamiltonian")
    if proof.get("commutes_with_hamiltonian") is not True:
        raise ValueError("commutant proof is missing")
    if proof.get("orthogonal_to_identity") is not True:
        raise ValueError("identity orthogonality proof is missing")
    if proof.get("orthogonal_to_hamiltonian") is not True:
        raise ValueError("Hamiltonian orthogonality proof is missing")

    tau_h_e5 = _parse_cubic(pairings.get("tau_h_e5"), "tau_h_e5")
    tau_h2_e5 = _parse_cubic(pairings.get("tau_h2_e5"), "tau_h2_e5")
    submitted_pairing = _parse_cubic(pairings.get("tau_w_e5"), "tau_w_e5")
    expected_pairing = tau_h2_e5 + h_coefficient * tau_h_e5
    if submitted_pairing != expected_pairing:
        raise ValueError("quadratic witness pairing mismatch")
    if submitted_pairing == Cubic.zero() or pairings.get("tau_w_e5_nonzero") is not True:
        raise ValueError("quadratic witness pairing is not certified nonzero")
    if claim.get("leading_order_status") != "no_go":
        raise ValueError("leading-order claim mismatch")
    if claim.get("finite_step_status") != "inconclusive":
        raise ValueError("finite-step claim exceeds the certificate")
    if payload.get("hpc_authorized") is not False:
        raise ValueError("certificate does not authorize HPC")
