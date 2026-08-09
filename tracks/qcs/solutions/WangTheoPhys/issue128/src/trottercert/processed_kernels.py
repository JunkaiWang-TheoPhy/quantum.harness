from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


_SOURCE_DOI = "10.5281/zenodo.10814897"
_SOURCE_FILE = "metemprar.py"
_SOURCE_MD5 = "35a37af90096f371e7c82366cc7fad8b"


@dataclass(frozen=True)
class RationalStage:
    fragment_index: int
    coefficient: Fraction


@dataclass(frozen=True)
class ProcessedKernel:
    """Published effective-order kernel with exact decimal rationalization.

    ``source_strings`` contains the one-sided coefficients passed to the
    source archive's ``estendreX`` helper.  The missing center coefficient is
    reconstructed as ``1/2 - sum(source_strings)`` before reflection.
    """

    name: str
    source_strings: tuple[str, ...]
    seed_coefficients: tuple[Fraction, ...]
    source_doi: str = _SOURCE_DOI
    source_file: str = _SOURCE_FILE
    source_md5: str = _SOURCE_MD5

    def __post_init__(self) -> None:
        if not self.source_strings:
            raise ValueError("processed kernel needs source coefficients")
        if len(self.source_strings) != len(self.seed_coefficients):
            raise ValueError("source strings and exact coefficients disagree")
        if tuple(Fraction(value) for value in self.source_strings) != (
            self.seed_coefficients
        ):
            raise ValueError("kernel coefficients are not exact source decimals")

    @property
    def composition_coefficients(self) -> tuple[Fraction, ...]:
        center = Fraction(1, 2) - sum(self.seed_coefficients, Fraction())
        half = self.seed_coefficients + (center,)
        return half + tuple(reversed(half))

    @property
    def composition_stage_count(self) -> int:
        return len(self.composition_coefficients) // 2


def _kernel(name: str, source_strings: Sequence[str]) -> ProcessedKernel:
    frozen = tuple(source_strings)
    return ProcessedKernel(
        name=name,
        source_strings=frozen,
        seed_coefficients=tuple(Fraction(value) for value in frozen),
    )


_S8_STRINGS = (
    "0.1535",
    "0.146",
    "0.1535",
    "0.1564865138360775523331602",
    "0.1777546764340215024573463",
    "-0.3260392072026447259933467",
    "-0.3377852074639320941003920",
)
_S10_POSITIVE = "0.100838384835000970361478569216"
_S10_NEGATIVE = "-0.238737866770265639777320334936"
_S10_STRINGS = (_S10_POSITIVE,) * 7 + (_S10_NEGATIVE,) * 2
_S11_POSITIVE = "0.0852884432504611078507516523709"
_S11_NEGATIVE = "-0.211683070446329023994526003378"
_S11_STRINGS = (_S11_POSITIVE,) * 8 + (_S11_NEGATIVE,) * 2

_PUBLISHED_KERNELS = {
    kernel.name: kernel
    for kernel in (
        _kernel("s8", _S8_STRINGS),
        _kernel("s10", _S10_STRINGS),
        _kernel("s11", _S11_STRINGS),
    )
}


def published_effective_order_six_kernel(name: str) -> ProcessedKernel:
    try:
        return _PUBLISHED_KERNELS[name]
    except KeyError as error:
        raise ValueError(f"unknown processed kernel: {name}") from error


def alternating_kernel_stages(
    kernel: ProcessedKernel,
    n_fragments: int = 4,
) -> tuple[RationalStage, ...]:
    """Expand Eq. (1.8) in a fixed left-to-right product convention.

    The leftmost map is ``chi(alpha_2s)``.  With
    ``chi = phi_n ... phi_1``, its left-to-right exponential order is
    ``n-1,...,0``; the adjoint has the reverse order.  Equal neighboring
    fragment exponentials are merged exactly.
    """

    if n_fragments < 2:
        raise ValueError("alternating kernels require at least two fragments")
    merged: list[RationalStage] = []
    for position, coefficient in enumerate(
        reversed(kernel.composition_coefficients)
    ):
        fragment_order = (
            range(n_fragments - 1, -1, -1)
            if position % 2 == 0
            else range(n_fragments)
        )
        for fragment_index in fragment_order:
            stage = RationalStage(fragment_index, coefficient)
            if merged and merged[-1].fragment_index == fragment_index:
                previous = merged.pop()
                combined = previous.coefficient + coefficient
                if combined:
                    merged.append(RationalStage(fragment_index, combined))
            elif coefficient:
                merged.append(stage)
    return tuple(merged)


def merged_groups_per_step(
    kernel: ProcessedKernel,
    n_fragments: int = 4,
) -> int:
    return len(alternating_kernel_stages(kernel, n_fragments))


def groups_for_repetitions(
    kernel: ProcessedKernel,
    repetitions: int,
    n_fragments: int = 4,
) -> int:
    if (
        isinstance(repetitions, bool)
        or not isinstance(repetitions, int)
        or repetitions < 1
    ):
        raise ValueError("repetitions must be a positive integer")
    stages = alternating_kernel_stages(kernel, n_fragments)
    boundary_saving = int(stages[0].fragment_index == stages[-1].fragment_index)
    return repetitions * len(stages) - (repetitions - 1) * boundary_saving


def maximum_steps_for_group_budget(
    kernel: ProcessedKernel,
    group_budget: int,
    n_fragments: int = 4,
) -> int:
    if (
        isinstance(group_budget, bool)
        or not isinstance(group_budget, int)
        or group_budget < 1
    ):
        raise ValueError("group budget must be a positive integer")
    stages = alternating_kernel_stages(kernel, n_fragments)
    if stages[0].fragment_index != stages[-1].fragment_index:
        return group_budget // len(stages)
    return max(0, (group_budget - 1) // (len(stages) - 1))
