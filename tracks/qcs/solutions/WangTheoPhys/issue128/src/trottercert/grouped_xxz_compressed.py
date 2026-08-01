"""Exact D4-compressed production ledger for the grouped 4 by 4 XXZ pilot.

The raw theorem stream is authenticated record by record but is not retained.
Weights remain attached to actual five-letter words.  D4 representatives only
deduplicate the unweighted finite-torus Pauli maps and their deterministic
pair-only norm witnesses; theorem weights are never summed across an orbit.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from fractions import Fraction
from hashlib import sha256
from typing import Literal

from .grouped_xxz import (
    FINITE_STEP_ERROR_FORMULA,
    GROUPING_ALGORITHM_IDENTIFIER,
    THEOREM_CENTER,
    THEOREM_DUHAMEL_CONVENTION,
    THEOREM_FACTORIAL_DENOMINATOR,
    THEOREM_IDENTIFIER,
    THEOREM_ORDER,
    AnticommutingGroupRecord,
    RawTheoremRecord,
    SymplecticCoefficient,
    XXZCompileSpec,
    _finite_fragment_bonds,
    _FiniteXXZSymplecticEvaluator,
    _fourth_order_error,
    _iter_raw_theorem_records,
    _minimal_fourth_order_steps,
    _projected_raw_record_count,
    canonical_json_bytes,
    discover_anticommuting_groups,
    fraction_pair,
    replay_schedule_resources,
    verify_anticommuting_groups,
    xxz_suzuki_schedule,
)
from .grouped_xxz_orbits import (
    XXZ_4X4_SPATIAL_SYMMETRIES,
    FiveLetterWord,
    canonical_five_letter_word,
    transport_five_letter_word,
    verify_frozen_spatial_symmetries,
)

COMPRESSED_METHOD = "compressed_d4_direct_finite_high_order_theorem_grouped_norm"
CompressedStatus = Literal["certified"]


@dataclass(frozen=True, slots=True)
class ActualWordWeight:
    """One actual theorem word, its exact weight, and its D4 transport."""

    actual_word: FiveLetterWord
    representative: FiveLetterWord
    symmetry_name: str
    raw_weight_upper: Fraction


@dataclass(frozen=True, slots=True)
class RepresentativeNormBlock:
    """One unweighted canonical Pauli map and deterministic norm witness."""

    representative: FiveLetterWord
    terms: tuple[SymplecticCoefficient, ...]
    groups: tuple[AnticommutingGroupRecord, ...]
    grouped_norm: Fraction
    triangle_norm: Fraction


@dataclass(frozen=True, slots=True)
class CompressedXXZLedger:
    """Authenticated theorem stream with D4-deduplicated norm blocks."""

    spec: XXZCompileSpec
    theorem_identifier: str
    order: int
    center: int
    factorial_denominator: int
    duhamel_convention: str
    finite_step_error_formula: str
    grouping_algorithm: str
    max_degree: int
    projected_record_count: int
    record_limit: int | None
    complete: bool
    raw_record_count: int
    raw_stream_sha256: str
    word_weights: tuple[ActualWordWeight, ...]
    representative_blocks: tuple[RepresentativeNormBlock, ...]
    grouped_constant: Fraction
    triangle_constant: Fraction
    ledger_digest: str


@dataclass(frozen=True, slots=True)
class CompressedXXZCertificate:
    """Exact finite-step closure of one complete compressed ledger."""

    status: CompressedStatus
    method: str
    spec: XXZCompileSpec
    ledger_digest: str
    raw_stream_sha256: str
    raw_record_count: int
    theorem_identifier: str
    order: int
    center: int
    factorial_denominator: int
    duhamel_convention: str
    finite_step_error_formula: str
    grouping_algorithm: str
    grouped_constant: Fraction
    triangle_constant: Fraction
    candidate_steps: int
    candidate_error: Fraction
    candidate_previous_error: Fraction | None
    baseline_steps: int
    baseline_error: Fraction
    baseline_previous_error: Fraction | None
    candidate_resources: int
    baseline_resources: int
    bond_growth: Fraction
    cell_base: Fraction


def _interval_payload(record: RawTheoremRecord) -> list[list[int]]:
    return [
        fraction_pair(record.weight_interval.lower),
        fraction_pair(record.weight_interval.upper),
    ]


def canonical_raw_record_bytes(record: RawTheoremRecord) -> bytes:
    """Return the canonical, domain-labelled bytes for one raw stream item."""

    if not isinstance(record, RawTheoremRecord):
        raise TypeError("record must be a RawTheoremRecord")
    return canonical_json_bytes(
        {
            "adjoint_stage_indices": list(record.adjoint_stage_indices),
            "base_stage_index": record.base_stage_index,
            "block_key": list(record.block_key),
            "composition": list(record.composition),
            "fragment_word": list(record.fragment_word),
            "partial_sum_index": record.partial_sum_index,
            "record_id": record.record_id,
            "side": record.side,
            "stream_domain": "grouped_xxz_raw_theorem_record_v1",
            "weight_interval": _interval_payload(record),
        }
    )


def _spec_payload(spec: XXZCompileSpec) -> dict[str, object]:
    return {
        "boundary": spec.boundary,
        "delta": fraction_pair(spec.delta),
        "formula_identifier": spec.formula_identifier,
        "length": spec.length,
        "normalization": spec.normalization,
        "primary_metric": spec.primary_metric,
        "stage_count": spec.stage_count,
        "time": fraction_pair(spec.time),
        "tolerance": fraction_pair(spec.tolerance),
    }


def _term_payload(term: SymplecticCoefficient) -> list[object]:
    return [
        term.x_mask,
        term.z_mask,
        fraction_pair(term.real),
        fraction_pair(term.imag),
    ]


def compressed_xxz_ledger_digest(ledger: CompressedXXZLedger) -> str:
    """Hash all submitted compressed semantics without materializing one JSON blob."""

    if not isinstance(ledger, CompressedXXZLedger):
        raise TypeError("ledger must be a CompressedXXZLedger")
    digest = sha256()
    digest.update(
        canonical_json_bytes(
            {
                "complete": ledger.complete,
                "digest_domain": "grouped_xxz_compressed_ledger_v1",
                "duhamel_convention": ledger.duhamel_convention,
                "factorial_denominator": ledger.factorial_denominator,
                "finite_step_error_formula": ledger.finite_step_error_formula,
                "grouped_constant": fraction_pair(ledger.grouped_constant),
                "grouping_algorithm": ledger.grouping_algorithm,
                "max_degree": ledger.max_degree,
                "order": ledger.order,
                "center": ledger.center,
                "projected_record_count": ledger.projected_record_count,
                "raw_record_count": ledger.raw_record_count,
                "raw_stream_sha256": ledger.raw_stream_sha256,
                "record_limit": ledger.record_limit,
                "representative_block_count": len(ledger.representative_blocks),
                "spec": _spec_payload(ledger.spec),
                "theorem_identifier": ledger.theorem_identifier,
                "triangle_constant": fraction_pair(ledger.triangle_constant),
                "word_weight_count": len(ledger.word_weights),
            }
        )
    )
    for row in ledger.word_weights:
        digest.update(
            canonical_json_bytes(
                {
                    "actual_word": list(row.actual_word),
                    "raw_weight_upper": fraction_pair(row.raw_weight_upper),
                    "representative": list(row.representative),
                    "row_domain": "actual_word_weight_v1",
                    "symmetry_name": row.symmetry_name,
                }
            )
        )
    for block in ledger.representative_blocks:
        digest.update(
            canonical_json_bytes(
                {
                    "block_domain": "representative_norm_block_v1",
                    "group_count": len(block.groups),
                    "grouped_norm": fraction_pair(block.grouped_norm),
                    "representative": list(block.representative),
                    "term_count": len(block.terms),
                    "triangle_norm": fraction_pair(block.triangle_norm),
                }
            )
        )
        for term in block.terms:
            digest.update(canonical_json_bytes({"term": _term_payload(term)}))
        for group in block.groups:
            digest.update(
                canonical_json_bytes(
                    {
                        "group_masks": [list(term.mask) for term in group.terms],
                        "norm_interval": [
                            fraction_pair(group.norm_interval.lower),
                            fraction_pair(group.norm_interval.upper),
                        ],
                        "squared_norm": fraction_pair(group.squared_norm),
                    }
                )
            )
    return digest.hexdigest()


def _validate_build_inputs(
    spec: XXZCompileSpec,
    max_degree: int,
    max_records: int | None,
    progress: Callable[[int, int], None] | None,
) -> None:
    if not isinstance(spec, XXZCompileSpec):
        raise TypeError("spec must be an XXZCompileSpec")
    if isinstance(max_degree, bool) or not isinstance(max_degree, int):
        raise TypeError("max_degree must be an integer")
    if max_degree < THEOREM_ORDER + 1:
        raise ValueError("max_degree must permit the degree-five theorem words")
    if max_records is not None and (
        isinstance(max_records, bool)
        or not isinstance(max_records, int)
        or max_records < 1
    ):
        raise ValueError("max_records must be a positive integer")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable")


def _stream_prefix(
    limit: int,
    projected: int,
    progress: Callable[[int, int], None] | None,
    *,
    require_eof: bool,
) -> tuple[int, str, dict[FiveLetterWord, Fraction]]:
    digest = sha256()
    weights: dict[FiveLetterWord, Fraction] = {}
    count = 0
    next_progress = 1
    records = iter(_iter_raw_theorem_records())
    for _ in range(limit):
        try:
            record = next(records)
        except StopIteration as error:
            raise ArithmeticError(
                "raw theorem enumeration ended before its projection"
            ) from error
        if (
            record.fragment_word != record.block_key
            or len(record.fragment_word) != THEOREM_ORDER + 1
            or record.weight_interval.lower <= 0
        ):
            raise ArithmeticError("raw theorem stream contains an invalid record")
        digest.update(canonical_raw_record_bytes(record))
        word = record.fragment_word
        weights[word] = weights.get(word, Fraction()) + record.weight_interval.upper
        count += 1
        while next_progress <= 20 and count * 20 >= next_progress * limit:
            if progress is not None:
                progress(count, projected)
            next_progress += 1
    if require_eof:
        try:
            next(records)
        except StopIteration:
            pass
        else:
            raise ArithmeticError(
                "raw theorem enumeration continues beyond its projection"
            )
    return count, digest.hexdigest(), weights


def _binding_for_word(word: FiveLetterWord) -> tuple[FiveLetterWord, str]:
    representative = canonical_five_letter_word(word)
    symmetry = next(
        (
            symmetry
            for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
            if transport_five_letter_word(representative, symmetry) == word
        ),
        None,
    )
    if symmetry is None:
        raise ArithmeticError("actual word is absent from its frozen D4 orbit")
    return representative, symmetry.name


def _word_rows(
    weights: dict[FiveLetterWord, Fraction],
) -> tuple[ActualWordWeight, ...]:
    rows = []
    for word in sorted(weights):
        representative, symmetry_name = _binding_for_word(word)
        rows.append(
            ActualWordWeight(
                actual_word=word,
                representative=representative,
                symmetry_name=symmetry_name,
                raw_weight_upper=weights[word],
            )
        )
    return tuple(rows)


def _representative_block(
    representative: FiveLetterWord,
    evaluator: _FiniteXXZSymplecticEvaluator,
) -> RepresentativeNormBlock:
    terms = evaluator.evaluate(representative)
    groups = discover_anticommuting_groups(terms)
    grouped_norm = verify_anticommuting_groups(terms, groups)
    if any(term.imag for term in terms):
        raise ArithmeticError("degree-five representative is not Hermitian")
    triangle_norm = sum((abs(term.real) for term in terms), Fraction())
    return RepresentativeNormBlock(
        representative=representative,
        terms=terms,
        groups=groups,
        grouped_norm=grouped_norm,
        triangle_norm=triangle_norm,
    )


def _literal_constants(
    word_weights: tuple[ActualWordWeight, ...],
    representative_blocks: tuple[RepresentativeNormBlock, ...],
) -> tuple[Fraction, Fraction]:
    blocks = {block.representative: block for block in representative_blocks}
    if len(blocks) != len(representative_blocks):
        raise ValueError("compressed ledger contains duplicate representative blocks")
    grouped = Fraction()
    triangle = Fraction()
    # Deliberately sum one actual word at a time.  No representative-weight
    # subtotal exists in this implementation.
    for row in word_weights:
        try:
            block = blocks[row.representative]
        except KeyError as error:
            raise ValueError("actual word has no representative block") from error
        grouped += row.raw_weight_upper * block.grouped_norm
        triangle += row.raw_weight_upper * block.triangle_norm
    return (
        grouped / THEOREM_FACTORIAL_DENOMINATOR,
        triangle / THEOREM_FACTORIAL_DENOMINATOR,
    )


def build_compressed_xxz_ledger(
    spec: XXZCompileSpec,
    *,
    max_degree: int = 5,
    max_records: int | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> CompressedXXZLedger:
    """Build a full compressed ledger or a profiling-only bounded prefix."""

    _validate_build_inputs(spec, max_degree, max_records, progress)
    verify_frozen_spatial_symmetries()
    projected = _projected_raw_record_count()
    limit = projected if max_records is None else min(max_records, projected)
    raw_count, stream_digest, weights = _stream_prefix(
        limit,
        projected,
        progress,
        require_eof=max_records is None,
    )
    rows = _word_rows(weights)
    representatives = tuple(sorted({row.representative for row in rows}))
    evaluator = _FiniteXXZSymplecticEvaluator(
        _finite_fragment_bonds(spec),
        spec.delta,
    )
    blocks = tuple(
        _representative_block(representative, evaluator)
        for representative in representatives
    )
    grouped_constant, triangle_constant = _literal_constants(rows, blocks)
    unsealed = CompressedXXZLedger(
        spec=spec,
        theorem_identifier=THEOREM_IDENTIFIER,
        order=THEOREM_ORDER,
        center=THEOREM_CENTER,
        factorial_denominator=THEOREM_FACTORIAL_DENOMINATOR,
        duhamel_convention=THEOREM_DUHAMEL_CONVENTION,
        finite_step_error_formula=FINITE_STEP_ERROR_FORMULA,
        grouping_algorithm=GROUPING_ALGORITHM_IDENTIFIER,
        max_degree=max_degree,
        projected_record_count=projected,
        record_limit=max_records,
        complete=max_records is None and raw_count == projected,
        raw_record_count=raw_count,
        raw_stream_sha256=stream_digest,
        word_weights=rows,
        representative_blocks=blocks,
        grouped_constant=grouped_constant,
        triangle_constant=triangle_constant,
        ledger_digest="",
    )
    return replace(
        unsealed,
        ledger_digest=compressed_xxz_ledger_digest(unsealed),
    )


def verify_compressed_xxz_ledger(ledger: CompressedXXZLedger) -> None:
    """Primary semantic replay of the raw stream, D4 bindings, and rep blocks."""

    if not isinstance(ledger, CompressedXXZLedger):
        raise TypeError("ledger must be a CompressedXXZLedger")
    if (
        ledger.theorem_identifier != THEOREM_IDENTIFIER
        or ledger.order != THEOREM_ORDER
        or ledger.center != THEOREM_CENTER
        or ledger.factorial_denominator != THEOREM_FACTORIAL_DENOMINATOR
        or ledger.duhamel_convention != THEOREM_DUHAMEL_CONVENTION
        or ledger.finite_step_error_formula != FINITE_STEP_ERROR_FORMULA
        or ledger.grouping_algorithm != GROUPING_ALGORITHM_IDENTIFIER
    ):
        raise ValueError("compressed ledger metadata differs from the frozen theorem")
    if ledger.max_degree < THEOREM_ORDER + 1:
        raise ValueError("compressed ledger max_degree does not cover degree five")
    if ledger.record_limit is not None and (
        isinstance(ledger.record_limit, bool)
        or not isinstance(ledger.record_limit, int)
        or ledger.record_limit < 1
    ):
        raise ValueError("compressed ledger record limit is invalid")
    projected = _projected_raw_record_count()
    if ledger.projected_record_count != projected:
        raise ValueError("compressed ledger record projection is incorrect")
    expected_count = (
        projected
        if ledger.record_limit is None
        else min(ledger.record_limit, projected)
    )
    if ledger.raw_record_count != expected_count:
        raise ValueError("compressed ledger raw count is incorrect")
    expected_complete = ledger.record_limit is None and expected_count == projected
    if ledger.complete != expected_complete:
        raise ValueError("compressed ledger completeness flag is incorrect")

    verify_frozen_spatial_symmetries()
    count, stream_digest, weights = _stream_prefix(
        expected_count,
        projected,
        None,
        require_eof=ledger.record_limit is None,
    )
    if count != ledger.raw_record_count or stream_digest != ledger.raw_stream_sha256:
        raise ValueError("compressed ledger raw stream digest is incorrect")
    expected_rows = _word_rows(weights)
    if ledger.word_weights != expected_rows:
        raise ValueError(
            "compressed ledger actual-word weights or D4 bindings are incorrect"
        )

    expected_representatives = tuple(
        sorted({row.representative for row in expected_rows})
    )
    if tuple(block.representative for block in ledger.representative_blocks) != (
        expected_representatives
    ):
        raise ValueError("compressed ledger representative coverage is incorrect")
    evaluator = _FiniteXXZSymplecticEvaluator(
        _finite_fragment_bonds(ledger.spec),
        ledger.spec.delta,
    )
    for submitted in ledger.representative_blocks:
        expected = _representative_block(submitted.representative, evaluator)
        if submitted != expected:
            raise ValueError(
                "compressed representative map or deterministic groups are incorrect"
            )

    grouped_constant, triangle_constant = _literal_constants(
        ledger.word_weights,
        ledger.representative_blocks,
    )
    if (
        ledger.grouped_constant != grouped_constant
        or ledger.triangle_constant != triangle_constant
    ):
        raise ValueError("compressed theorem constants are not the literal word sum")
    if ledger.grouped_constant > ledger.triangle_constant:
        raise ValueError("grouped theorem constant exceeds its triangle baseline")
    if ledger.ledger_digest != compressed_xxz_ledger_digest(ledger):
        raise ValueError("compressed ledger outer digest is incorrect")


def _close_compressed_xxz_certificate(
    ledger: CompressedXXZLedger,
) -> CompressedXXZCertificate:
    if not isinstance(ledger, CompressedXXZLedger):
        raise TypeError("ledger must be a CompressedXXZLedger")
    if not ledger.complete:
        raise ValueError("compressed certificate closure requires a complete ledger")
    tolerance = ledger.spec.tolerance
    candidate_steps = _minimal_fourth_order_steps(
        ledger.grouped_constant,
        tolerance,
    )
    baseline_steps = _minimal_fourth_order_steps(
        ledger.triangle_constant,
        tolerance,
    )
    schedule = xxz_suzuki_schedule()
    delta_absolute = abs(ledger.spec.delta)
    return CompressedXXZCertificate(
        status="certified",
        method=COMPRESSED_METHOD,
        spec=ledger.spec,
        ledger_digest=ledger.ledger_digest,
        raw_stream_sha256=ledger.raw_stream_sha256,
        raw_record_count=ledger.raw_record_count,
        theorem_identifier=ledger.theorem_identifier,
        order=ledger.order,
        center=ledger.center,
        factorial_denominator=ledger.factorial_denominator,
        duhamel_convention=ledger.duhamel_convention,
        finite_step_error_formula=ledger.finite_step_error_formula,
        grouping_algorithm=ledger.grouping_algorithm,
        grouped_constant=ledger.grouped_constant,
        triangle_constant=ledger.triangle_constant,
        candidate_steps=candidate_steps,
        candidate_error=_fourth_order_error(
            ledger.grouped_constant,
            candidate_steps,
        ),
        candidate_previous_error=(
            None
            if candidate_steps == 1
            else _fourth_order_error(
                ledger.grouped_constant,
                candidate_steps - 1,
            )
        ),
        baseline_steps=baseline_steps,
        baseline_error=_fourth_order_error(
            ledger.triangle_constant,
            baseline_steps,
        ),
        baseline_previous_error=(
            None
            if baseline_steps == 1
            else _fourth_order_error(
                ledger.triangle_constant,
                baseline_steps - 1,
            )
        ),
        candidate_resources=replay_schedule_resources(schedule, candidate_steps),
        baseline_resources=replay_schedule_resources(schedule, baseline_steps),
        bond_growth=max(Fraction(1), (1 + delta_absolute) / 2),
        cell_base=(2 + delta_absolute) / 2,
    )


def verify_compressed_xxz_certificate(
    certificate: CompressedXXZCertificate,
    ledger: CompressedXXZLedger,
) -> None:
    """Replay a complete ledger and compare every finite-step field exactly."""

    if not isinstance(certificate, CompressedXXZCertificate):
        raise TypeError("certificate must be a CompressedXXZCertificate")
    if not isinstance(ledger, CompressedXXZLedger) or not ledger.complete:
        raise ValueError(
            "compressed certificate verification requires a complete ledger"
        )
    verify_compressed_xxz_ledger(ledger)
    _verify_compressed_xxz_certificate_fields(certificate, ledger)


def _verify_compressed_xxz_certificate_fields(
    certificate: CompressedXXZCertificate,
    ledger: CompressedXXZLedger,
) -> None:
    """Compare closure fields after the caller has authenticated the ledger."""

    expected = _close_compressed_xxz_certificate(ledger)
    if certificate != expected:
        raise ValueError("compressed certificate differs from semantic replay")


def compile_compressed_grouped_xxz(
    spec: XXZCompileSpec,
    *,
    ledger: CompressedXXZLedger | None = None,
) -> CompressedXXZCertificate:
    """Compile only a complete compressed stream; prefixes are profiles only."""

    if not isinstance(spec, XXZCompileSpec):
        raise TypeError("spec must be an XXZCompileSpec")
    selected = ledger if ledger is not None else build_compressed_xxz_ledger(spec)
    if selected.spec != spec:
        raise ValueError("supplied compressed ledger does not match the XXZ spec")
    if not selected.complete:
        raise ValueError("compile_compressed_grouped_xxz requires a complete ledger")
    verify_compressed_xxz_ledger(selected)
    certificate = _close_compressed_xxz_certificate(selected)
    _verify_compressed_xxz_certificate_fields(certificate, selected)
    return certificate
