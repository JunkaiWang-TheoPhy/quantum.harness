from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction

from .commutant_witness import (
    heisenberg_symplectic_terms,
    square_real_pauli_terms,
)
from .cubic_field import Cubic
from .dual_word_manifest import FORMULA_NAME, WordManifest
from .lattice import SquareLattice
from .local_commutators import (
    SymplecticDyadicLocalDensityEvaluator,
    _iter_set_bits,
)


@dataclass(frozen=True, slots=True)
class ManifestPairingPartial:
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


def _validate_manifest_for_contraction(manifest: WordManifest) -> None:
    if manifest.degree < 3 or manifest.degree % 2 == 0:
        raise ValueError("manifest degree must be odd and at least three")
    if manifest.formula != FORMULA_NAME:
        raise ValueError("manifest formula is unsupported")
    if manifest.shard_count < 1 or not 0 <= manifest.shard_index < manifest.shard_count:
        raise ValueError("manifest shard metadata is invalid")
    if manifest.total_groups < len(manifest.groups):
        raise ValueError("manifest group count exceeds global total")
    ordinals = tuple(group.ordinal for group in manifest.groups)
    if ordinals != tuple(sorted(set(ordinals))):
        raise ValueError("manifest groups are duplicated or unordered")
    for group in manifest.groups:
        if group.ordinal % manifest.shard_count != manifest.shard_index:
            raise ValueError("manifest group does not belong to its shard")
        if len(group.suffix) != manifest.degree - 1:
            raise ValueError("manifest suffix length does not match degree")
        if not 1 <= len(group.records) <= 4:
            raise ValueError("manifest group record count is invalid")
        words = tuple(record.word for record in group.records)
        if words != tuple(sorted(set(words))):
            raise ValueError("manifest words are duplicated or unordered")
        for record in group.records:
            if len(record.word) != manifest.degree:
                raise ValueError("manifest word length does not match degree")
            if record.word[1:] != group.suffix:
                raise ValueError("manifest word does not match suffix")
            if record.coefficient == Cubic.zero():
                raise ValueError("manifest word coefficient must be nonzero")


def contract_word_manifest(
    manifest: WordManifest,
    *,
    length: int = 6,
    progress: Callable[[int, int, ManifestPairingPartial], None] | None = None,
) -> ManifestPairingPartial:
    """Contract one verified exact word manifest with H, H^2, and W."""

    _validate_manifest_for_contraction(manifest)
    if length < 6 or length % 2:
        raise ValueError("contraction torus length must be even and at least six")

    lattice = SquareLattice(length)
    cells = lattice.n_sites // 4
    hamiltonian = heisenberg_symplectic_terms(lattice)
    squared = square_real_pauli_terms(hamiltonian)
    evaluator = SymplecticDyadicLocalDensityEvaluator(shared_coordinates=True)
    registry = evaluator.registries[0]
    denominator = manifest.degree * (
        1 << evaluator.denominator_exponent((0,) * manifest.degree)
    )

    tau_h = Cubic.zero()
    tau_h2 = Cubic.zero()
    word_count = nonzero_word_count = retained_term_count = 0
    completed_groups: list[int] = []

    for completed, group in enumerate(manifest.groups, start=1):
        for record in group.records:
            word_count += 1
            operator = evaluator.evaluate(record.word)
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
                tau_h += record.coefficient * (h_numerator / denominator)
            if h2_numerator:
                tau_h2 += record.coefficient * (h2_numerator / denominator)
            evaluator.cache.pop(record.word, None)
        evaluator.cache.clear()
        completed_groups.append(group.ordinal)
        if progress is not None:
            progress(
                completed,
                len(manifest.groups),
                ManifestPairingPartial(
                    degree=manifest.degree,
                    length=length,
                    shard_index=manifest.shard_index,
                    shard_count=manifest.shard_count,
                    total_groups=manifest.total_groups,
                    group_indices=tuple(completed_groups),
                    word_count=word_count,
                    nonzero_word_count=nonzero_word_count,
                    retained_term_count=retained_term_count,
                    tau_h=tau_h,
                    tau_h2=tau_h2,
                    tau_w=tau_h2 + tau_h / 2,
                ),
            )

    return ManifestPairingPartial(
        degree=manifest.degree,
        length=length,
        shard_index=manifest.shard_index,
        shard_count=manifest.shard_count,
        total_groups=manifest.total_groups,
        group_indices=tuple(completed_groups),
        word_count=word_count,
        nonzero_word_count=nonzero_word_count,
        retained_term_count=retained_term_count,
        tau_h=tau_h,
        tau_h2=tau_h2,
        tau_w=tau_h2 + tau_h / 2,
    )
