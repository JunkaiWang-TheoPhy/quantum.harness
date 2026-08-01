import copy
import json
from fractions import Fraction
from pathlib import Path

import pytest

from trottercert.commutant_witness import (
    build_quadratic_witness_payload,
    density_aliases_on_torus,
    heisenberg_symplectic_terms,
    lift_cubic_density,
    pair_lifted_cubic_density,
    quadratic_witness_moments,
    square_real_pauli_terms,
    verify_quadratic_witness_payload,
)
from trottercert.cubic_field import Cubic
from trottercert.cubic_local import exact_matching_density
from trottercert.lattice import SquareLattice
from trottercert.local_commutators import CoordinateRegistry

ROOT = Path(__file__).resolve().parents[1]
STORED_WITNESS = (
    ROOT
    / "docs/experiments/processor-obstruction/quadratic-commutant-witness.json"
)


def test_matching_densities_reconstruct_periodic_hamiltonian() -> None:
    lattice = SquareLattice(4)
    observed = {}
    for color in range(4):
        registry, density = exact_matching_density(color)
        lifted = lift_cubic_density(registry, density, lattice)
        for pauli, coefficient in lifted.items():
            assert coefficient.a1 == coefficient.a2 == 0
            observed[pauli] = observed.get(pauli, Fraction()) + coefficient.a0
    assert observed == heisenberg_symplectic_terms(lattice)


def test_square_and_trace_moments_are_exact_on_four_by_four() -> None:
    lattice = SquareLattice(4)
    hamiltonian = heisenberg_symplectic_terms(lattice)
    squared = square_real_pauli_terms(hamiltonian)
    assert squared[(0, 0)] == 6
    assert sum(
        coefficient * squared.get(pauli, Fraction())
        for pauli, coefficient in hamiltonian.items()
    ) == -3
    assert all(isinstance(value, Fraction) for value in squared.values())


def test_cubic_density_lift_and_sparse_pairing_agree() -> None:
    lattice = SquareLattice(4)
    registry = CoordinateRegistry()
    first = registry.site((0, 0))
    second = registry.site((1, 0))
    pauli = ((1 << first) | (1 << second), 0)
    density = {pauli: Cubic(1, 2, 3)}
    lifted = lift_cubic_density(registry, density, lattice)
    rational = {key: Fraction(index + 1, 7) for index, key in enumerate(lifted)}
    direct = Cubic.zero()
    for key, coefficient in lifted.items():
        direct += coefficient * rational[key]
    assert pair_lifted_cubic_density(
        registry,
        density,
        rational,
        lattice,
    ) == direct


def test_certificate_proves_calibrated_leading_obstruction() -> None:
    payload = build_quadratic_witness_payload(
        length=12,
        tau_h2=Fraction(54),
        tau_h3=Fraction(-27),
        tau_h_e5=Cubic(
            Fraction(10767, 200000),
            Fraction(91083, 3200000),
            Fraction(8439, 400000),
        ),
        tau_h2_e5=Cubic(
            Fraction(-8399, 80000),
            Fraction(-71051, 1280000),
            Fraction(-6583, 160000),
        ),
        matching_reconstruction_terms=864,
    )
    verify_quadratic_witness_payload(payload)
    assert payload["witness"]["identity_coefficient"] == [-54, 1]
    assert payload["witness"]["hamiltonian_coefficient"] == [1, 2]
    assert payload["proof"]["orthogonal_to_identity"] is True
    assert payload["proof"]["orthogonal_to_hamiltonian"] is True
    assert payload["pairings"]["tau_w_e5"][0] == [-7807, 100000]
    assert payload["claim"]["leading_order_status"] == "no_go"
    assert payload["claim"]["finite_step_status"] == "inconclusive"


def test_certificate_rejects_mutated_pairing_or_claim() -> None:
    payload = build_quadratic_witness_payload(
        length=12,
        tau_h2=Fraction(54),
        tau_h3=Fraction(-27),
        tau_h_e5=Cubic(1, 2, 3),
        tau_h2_e5=Cubic(4, 5, 6),
        matching_reconstruction_terms=864,
    )
    forged = copy.deepcopy(payload)
    forged["pairings"]["tau_w_e5"][0] = [0, 1]
    with pytest.raises(ValueError, match="pairing"):
        verify_quadratic_witness_payload(forged)
    forged = copy.deepcopy(payload)
    forged["claim"]["finite_step_status"] = "proved"
    with pytest.raises(ValueError, match="finite-step"):
        verify_quadratic_witness_payload(forged)


def test_stored_twelve_by_twelve_certificate_is_internally_valid() -> None:
    payload = json.loads(STORED_WITNESS.read_text())
    verify_quadratic_witness_payload(payload)
    assert payload["moments"]["tau_h2"] == [54, 1]
    assert payload["moments"]["tau_h3"] == [-27, 1]
    assert payload["pairings"]["tau_w_e5_nonzero"] is True


def test_lift_rejects_coordinate_aliasing() -> None:
    lattice = SquareLattice(4)
    registry = CoordinateRegistry()
    left = registry.site((0, 0))
    wrapped = registry.site((4, 0))
    pauli = ((1 << left) | (1 << wrapped), 0)
    with pytest.raises(ValueError, match="alias"):
        lift_cubic_density(registry, {pauli: Cubic.one()}, lattice)


def test_density_alias_detection_is_size_specific() -> None:
    registry = CoordinateRegistry()
    left = registry.site((0, 0))
    distant = registry.site((4, 0))
    pauli = ((1 << left) | (1 << distant), 0)
    density = {pauli: Cubic.one()}
    assert density_aliases_on_torus(registry, density, 4) is True
    assert density_aliases_on_torus(registry, density, 6) is False
    with pytest.raises(ValueError, match="even"):
        density_aliases_on_torus(registry, density, 5)


def test_quadratic_witness_moments_include_exact_normalization() -> None:
    lattice = SquareLattice(4)
    hamiltonian = heisenberg_symplectic_terms(lattice)
    squared = square_real_pauli_terms(hamiltonian)
    tau_h_e5 = Cubic(1, 2, 3)
    tau_h2_e5 = Cubic(4, 5, 6)
    moments = quadratic_witness_moments(
        hamiltonian,
        squared,
        tau_h_e5,
        tau_h2_e5,
    )
    assert moments.tau_h2 == Fraction(6)
    assert moments.tau_h3 == Fraction(-3)
    assert moments.h_coefficient == Fraction(1, 2)
    assert moments.tau_w2 > 0
    assert moments.tau_w_e5 == tau_h2_e5 + tau_h_e5 / 2
    assert (
        moments.squared_normalized_pairing
        == moments.tau_w_e5**2 / moments.tau_w2
    )
