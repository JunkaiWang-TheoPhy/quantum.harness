from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import groupby

from .commutant_witness import (
    heisenberg_symplectic_terms,
    square_real_pauli_terms,
)
from .cubic_field import Cubic, CubicStage
from .cubic_local import cubic_formula_log_series
from .lattice import SquareLattice
from .local_commutators import (
    SymplecticDyadicLocalDensityEvaluator,
    _iter_set_bits,
)


@dataclass(frozen=True, slots=True)
class DualPairingPartial:
    degree: int
    length: int
    shard_index: int
    shard_count: int
    total_groups: int
    group_indices: tuple[int, ...]
    word_count: int
    nonzero_word_count: int
    retained_term_count: int
    tau_h: Cubic
    tau_h2: Cubic
    tau_w: Cubic


def suffix_group_ordinals(
    total_groups: int,
    shard_index: int,
    shard_count: int,
) -> tuple[int, ...]:
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard index must lie in [0, shard_count)")
    if total_groups < 0:
        raise ValueError("total group count must be nonnegative")
    return tuple(range(shard_index, total_groups, shard_count))


def contract_log_degree_shard(
    stages: Sequence[CubicStage],
    degree: int,
    shard_index: int,
    shard_count: int,
    *,
    length: int = 6,
    progress: Callable[[int, int, DualPairingPartial], None] | None = None,
) -> DualPairingPartial:
    """Contract one exact odd logarithm degree with H, H^2, and W."""

    if degree < 3 or degree % 2 == 0:
        raise ValueError("logarithm degree must be an odd integer at least three")
    if length < 6 or length % 2:
        raise ValueError("contraction torus length must be even and at least six")
    suffix_group_ordinals(0, shard_index, shard_count)

    word_map = cubic_formula_log_series(stages, degree)[degree]
    ordered = sorted(
        word_map.items(),
        key=lambda item: (item[0][1:], item[0][0]),
    )
    groups = tuple(
        (suffix, tuple(items))
        for suffix, items in groupby(ordered, key=lambda item: item[0][1:])
    )
    selected = suffix_group_ordinals(len(groups), shard_index, shard_count)

    lattice = SquareLattice(length)
    cells = lattice.n_sites // 4
    hamiltonian = heisenberg_symplectic_terms(lattice)
    squared = square_real_pauli_terms(hamiltonian)
    evaluator = SymplecticDyadicLocalDensityEvaluator(shared_coordinates=True)
    registry = evaluator.registries[0]
    denominator = degree * (1 << evaluator.denominator_exponent((0,) * degree))

    tau_h = Cubic.zero()
    tau_h2 = Cubic.zero()
    word_count = nonzero_word_count = retained_term_count = 0

    for completed, group_index in enumerate(selected, start=1):
        _, items = groups[group_index]
        for word, word_coefficient in items:
            word_count += 1
            operator = evaluator.evaluate(word)
            if operator:
                nonzero_word_count += 1
            h_numerator = Fraction()
            h2_numerator = Fraction()
            for (x_mask, z_mask), numerator in operator.items():
                if (x_mask | z_mask).bit_count() > 4:
                    continue
                retained_term_count += 1
                global_x = global_z = 0
                occupied: set[int] = set()
                for site in _iter_set_bits(x_mask | z_mask):
                    x, y = registry.coordinate(site)
                    target_site = lattice.site(x, y)
                    if target_site in occupied:
                        raise ValueError(
                            "local Pauli coordinates alias on contraction torus"
                        )
                    occupied.add(target_site)
                    target_bit = 1 << target_site
                    source_bit = 1 << site
                    if x_mask & source_bit:
                        global_x |= target_bit
                    if z_mask & source_bit:
                        global_z |= target_bit
                pauli = (global_x, global_z)
                h_coefficient = hamiltonian.get(pauli)
                if h_coefficient is not None:
                    h_numerator += cells * numerator * h_coefficient
                h2_coefficient = squared.get(pauli)
                if h2_coefficient is not None:
                    h2_numerator += cells * numerator * h2_coefficient
            if h_numerator:
                tau_h += word_coefficient * (h_numerator / denominator)
            if h2_numerator:
                tau_h2 += word_coefficient * (h2_numerator / denominator)
            evaluator.cache.pop(word, None)
        evaluator.cache.clear()
        if progress is not None:
            progress(
                completed,
                len(selected),
                DualPairingPartial(
                    degree=degree,
                    length=length,
                    shard_index=shard_index,
                    shard_count=shard_count,
                    total_groups=len(groups),
                    group_indices=selected[:completed],
                    word_count=word_count,
                    nonzero_word_count=nonzero_word_count,
                    retained_term_count=retained_term_count,
                    tau_h=tau_h,
                    tau_h2=tau_h2,
                    tau_w=tau_h2 + tau_h / 2,
                ),
            )

    return DualPairingPartial(
        degree=degree,
        length=length,
        shard_index=shard_index,
        shard_count=shard_count,
        total_groups=len(groups),
        group_indices=selected,
        word_count=word_count,
        nonzero_word_count=nonzero_word_count,
        retained_term_count=retained_term_count,
        tau_h=tau_h,
        tau_h2=tau_h2,
        tau_w=tau_h2 + tau_h / 2,
    )
