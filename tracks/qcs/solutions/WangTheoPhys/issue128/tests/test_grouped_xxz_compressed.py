from __future__ import annotations

import json
import resource
from collections import defaultdict
from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
from itertools import islice
from time import monotonic

import pytest

import trottercert.grouped_xxz_compressed as compressed_module
from trottercert.grouped_xxz import (
    XXZCompileSpec,
    _finite_fragment_bonds,
    _FiniteXXZSymplecticEvaluator,
    _iter_raw_theorem_records,
)
from trottercert.grouped_xxz_compressed import (
    build_compressed_xxz_ledger,
    canonical_raw_record_bytes,
    compile_compressed_grouped_xxz,
    compressed_xxz_ledger_digest,
    verify_compressed_xxz_certificate,
    verify_compressed_xxz_ledger,
)
from trottercert.grouped_xxz_orbits import (
    XXZ_4X4_SPATIAL_SYMMETRIES,
    canonical_five_letter_word,
    transport_commutator_map,
    transport_five_letter_word,
)


def _reseal(ledger):
    unsealed = replace(ledger, ledger_digest="")
    return replace(unsealed, ledger_digest=compressed_xxz_ledger_digest(unsealed))


def _term_map(terms):
    return {term.mask: (term.real, term.imag) for term in terms}


def test_compressed_prefix_matches_the_canonical_raw_stream() -> None:
    spec = XXZCompileSpec.pilot(Fraction(1, 2))
    ledger = build_compressed_xxz_ledger(spec, max_records=8)
    expected_records = tuple(islice(_iter_raw_theorem_records(), 8))
    stream = sha256()
    expected_weights = defaultdict(Fraction)
    for record in expected_records:
        stream.update(canonical_raw_record_bytes(record))
        expected_weights[record.fragment_word] += record.weight_interval.upper

    assert ledger.raw_record_count == 8
    assert ledger.raw_stream_sha256 == stream.hexdigest()
    assert not ledger.complete
    assert not hasattr(ledger, "raw_records")
    assert {
        row.actual_word: row.raw_weight_upper for row in ledger.word_weights
    } == dict(expected_weights)
    assert tuple(row.actual_word for row in ledger.word_weights) == tuple(
        sorted(expected_weights)
    )
    verify_compressed_xxz_ledger(ledger)


def test_actual_words_bind_canonical_representatives_and_frozen_symmetries() -> None:
    ledger = build_compressed_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(2)),
        max_records=8,
    )
    symmetries = {symmetry.name: symmetry for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES}
    assert {row.representative for row in ledger.word_weights} == {
        block.representative for block in ledger.representative_blocks
    }
    for row in ledger.word_weights:
        assert row.representative == canonical_five_letter_word(row.actual_word)
        assert (
            transport_five_letter_word(
                row.representative,
                symmetries[row.symmetry_name],
            )
            == row.actual_word
        )


def test_exact_evaluator_maps_follow_all_frozen_transports() -> None:
    spec = XXZCompileSpec.pilot(Fraction(1, 2))
    evaluator = _FiniteXXZSymplecticEvaluator(
        _finite_fragment_bonds(spec),
        spec.delta,
    )
    representative = canonical_five_letter_word((0, 0, 0, 0, 1))
    source = _term_map(evaluator.evaluate(representative))
    assert source
    for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES:
        target_word = transport_five_letter_word(representative, symmetry)
        target = _term_map(evaluator.evaluate(target_word))
        assert target == transport_commutator_map(source, symmetry)


def test_primary_semantic_replay_rejects_resealed_mutations() -> None:
    ledger = build_compressed_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(1, 2)),
        max_records=4,
    )
    rows = list(ledger.word_weights)
    wrong = next(
        symmetry.name
        for symmetry in XXZ_4X4_SPATIAL_SYMMETRIES
        if symmetry.name != rows[0].symmetry_name
    )
    attacks = []
    changed = list(rows)
    changed[0] = replace(
        changed[0],
        raw_weight_upper=changed[0].raw_weight_upper + 1,
    )
    attacks.append(replace(ledger, word_weights=tuple(changed)))
    changed = list(rows)
    changed[0] = replace(changed[0], representative=(3, 3, 3, 3, 3))
    attacks.append(replace(ledger, word_weights=tuple(changed)))
    changed = list(rows)
    changed[0] = replace(changed[0], symmetry_name=wrong)
    attacks.append(replace(ledger, word_weights=tuple(changed)))
    attacks.append(replace(ledger, raw_stream_sha256="0" * 64))
    attacks.append(replace(ledger, word_weights=tuple(rows[1:])))
    for attacked in attacks:
        with pytest.raises(ValueError):
            verify_compressed_xxz_ledger(_reseal(attacked))


def test_bounded_profile_fails_closed_at_certificate_boundary() -> None:
    spec = XXZCompileSpec.pilot(Fraction(1, 2))
    ledger = build_compressed_xxz_ledger(spec, max_records=1)
    with pytest.raises(ValueError, match="complete"):
        compile_compressed_grouped_xxz(spec, ledger=ledger)


def test_full_stream_rejects_a_projection_that_is_one_too_small(monkeypatch) -> None:
    records = tuple(islice(_iter_raw_theorem_records(), 3))
    monkeypatch.setattr(compressed_module, "_projected_raw_record_count", lambda: 2)
    monkeypatch.setattr(
        compressed_module,
        "_iter_raw_theorem_records",
        lambda: iter(records),
    )
    with pytest.raises(ArithmeticError, match="projection|beyond"):
        build_compressed_xxz_ledger(XXZCompileSpec.pilot(Fraction(1, 2)))


def test_full_stream_rejects_an_extra_tail_but_bounded_prefix_allows_it(
    monkeypatch,
) -> None:
    expected = tuple(islice(_iter_raw_theorem_records(), 2))
    extra = replace(expected[-1], record_id="injected-extra-tail")
    stream = expected + (extra,)
    monkeypatch.setattr(compressed_module, "_projected_raw_record_count", lambda: 2)
    monkeypatch.setattr(
        compressed_module,
        "_iter_raw_theorem_records",
        lambda: iter(stream),
    )
    spec = XXZCompileSpec.pilot(Fraction(2))
    with pytest.raises(ArithmeticError, match="projection|beyond"):
        build_compressed_xxz_ledger(spec)
    bounded = build_compressed_xxz_ledger(spec, max_records=1)
    assert bounded.raw_record_count == 1
    assert not bounded.complete


def test_zero_representative_block_is_retained_and_required(monkeypatch) -> None:
    raw = next(_iter_raw_theorem_records())
    zero_word = (0, 0, 0, 0, 0)
    zero_record = replace(raw, fragment_word=zero_word, block_key=zero_word)
    monkeypatch.setattr(compressed_module, "_projected_raw_record_count", lambda: 1)
    monkeypatch.setattr(
        compressed_module,
        "_iter_raw_theorem_records",
        lambda: iter((zero_record,)),
    )
    ledger = build_compressed_xxz_ledger(XXZCompileSpec.pilot(Fraction(1, 2)))
    assert ledger.complete
    assert len(ledger.representative_blocks) == 1
    block = ledger.representative_blocks[0]
    assert block.representative == zero_word
    assert block.terms == ()
    assert block.groups == ()
    assert block.grouped_norm == block.triangle_norm == 0
    verify_compressed_xxz_ledger(ledger)

    deleted = _reseal(replace(ledger, representative_blocks=()))
    with pytest.raises(ValueError, match="representative coverage"):
        verify_compressed_xxz_ledger(deleted)


def test_primary_replay_rejects_resealed_rep_term_and_group_mutations() -> None:
    ledger = build_compressed_xxz_ledger(
        XXZCompileSpec.pilot(Fraction(2)),
        max_records=1,
    )
    block = ledger.representative_blocks[0]
    assert block.terms and block.groups

    terms = list(block.terms)
    terms[0] = replace(terms[0], real=terms[0].real + 1)
    term_attack = replace(block, terms=tuple(terms))
    attacked = _reseal(replace(ledger, representative_blocks=(term_attack,)))
    with pytest.raises(ValueError, match="representative map|groups"):
        verify_compressed_xxz_ledger(attacked)

    group_attack = replace(block, groups=block.groups[1:])
    attacked = _reseal(replace(ledger, representative_blocks=(group_attack,)))
    with pytest.raises(ValueError, match="representative map|groups"):
        verify_compressed_xxz_ledger(attacked)


def test_complete_semantic_replay_closes_exact_steps_and_resources(monkeypatch) -> None:
    records = tuple(islice(_iter_raw_theorem_records(), 2))
    monkeypatch.setattr(
        compressed_module,
        "_projected_raw_record_count",
        lambda: len(records),
    )
    monkeypatch.setattr(
        compressed_module,
        "_iter_raw_theorem_records",
        lambda: iter(records),
    )
    spec = XXZCompileSpec.pilot(Fraction(2))
    ledger = build_compressed_xxz_ledger(spec)
    assert ledger.complete
    assert ledger.record_limit is None
    limited_complete = _reseal(replace(ledger, record_limit=len(records)))
    with pytest.raises(ValueError, match="completeness|record limit"):
        verify_compressed_xxz_ledger(limited_complete)
    certificate = compile_compressed_grouped_xxz(spec, ledger=ledger)
    verify_compressed_xxz_certificate(certificate, ledger)
    assert certificate.candidate_error <= spec.tolerance
    assert certificate.baseline_error <= spec.tolerance
    if certificate.candidate_previous_error is not None:
        assert certificate.candidate_previous_error > spec.tolerance
    if certificate.baseline_previous_error is not None:
        assert certificate.baseline_previous_error > spec.tolerance
    assert certificate.candidate_resources == 30 * certificate.candidate_steps + 1
    assert certificate.baseline_resources == 30 * certificate.baseline_steps + 1

    attacks = (
        replace(certificate, candidate_steps=certificate.candidate_steps + 1),
        replace(
            certificate,
            candidate_resources=certificate.candidate_resources + 30,
        ),
        replace(certificate, grouped_constant=certificate.grouped_constant + 1),
    )
    monkeypatch.setattr(
        compressed_module,
        "verify_compressed_xxz_ledger",
        lambda submitted: None,
    )
    for attacked in attacks:
        with pytest.raises(ValueError):
            verify_compressed_xxz_certificate(attacked, ledger)


@pytest.mark.slow
@pytest.mark.parametrize("delta", [Fraction(1, 2), Fraction(2)])
def test_compressed_200_record_profile(delta: Fraction) -> None:
    started = monotonic()
    ledger = build_compressed_xxz_ledger(
        XXZCompileSpec.pilot(delta),
        max_records=200,
    )
    verify_compressed_xxz_ledger(ledger)
    elapsed = monotonic() - started
    pairs = sum(
        len(group.terms) == 2
        for block in ledger.representative_blocks
        for group in block.groups
    )
    singletons = sum(
        len(group.terms) == 1
        for block in ledger.representative_blocks
        for group in block.groups
    )
    representative_terms = sum(
        len(block.terms) for block in ledger.representative_blocks
    )
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    report = {
        "delta": str(delta),
        "raw_records": ledger.raw_record_count,
        "active_words": len(ledger.word_weights),
        "representatives": len(ledger.representative_blocks),
        "representative_terms": representative_terms,
        "pairs": pairs,
        "singletons": singletons,
        "grouped_triangle_ratio": float(
            ledger.grouped_constant / ledger.triangle_constant
        ),
        "elapsed_seconds": elapsed,
        "peak_rss_native": peak_rss,
    }
    print("COMPRESSED_PROFILE=" + json.dumps(report, sort_keys=True), flush=True)

    expected_records = tuple(islice(_iter_raw_theorem_records(), 200))
    stream = sha256()
    expected_weights = defaultdict(Fraction)
    for record in expected_records:
        stream.update(canonical_raw_record_bytes(record))
        expected_weights[record.fragment_word] += record.weight_interval.upper

    assert ledger.raw_record_count == 200
    assert ledger.raw_stream_sha256 == stream.hexdigest()
    assert {
        row.actual_word: row.raw_weight_upper for row in ledger.word_weights
    } == dict(expected_weights)
    assert len(ledger.word_weights) == 58
    assert len(ledger.representative_blocks) == 33
    assert representative_terms > 0
    assert pairs > 0
    assert singletons > 0
    assert ledger.grouped_constant <= ledger.triangle_constant
    assert elapsed < 180
