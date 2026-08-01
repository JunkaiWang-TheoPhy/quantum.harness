from __future__ import annotations

import ast
import copy
import hashlib
import json
import os
import subprocess
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import islice
from pathlib import Path

import pytest

from scripts.reference_grouped_xxz_theorem import (
    PROJECTED_RECORD_COUNT,
    Cubic,
    RawTheoremSummary,
    canonical_json_bytes,
    decode_record_payload,
    iter_raw_theorem_records,
    load_canonical_json,
    record_payload,
    strict_fraction_pair,
    summarize_raw_theorem,
    theorem_stages,
)
from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.grouped_xxz import (
    _iter_raw_theorem_records,
    _theorem_interval_stages,
)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "scripts/reference_grouped_xxz_theorem.py"
PREFIX_DIGEST = "e57f4e2146c192658d89ae295d0bad5963e30ddf51c69f392dcfb794a3df378c"
FULL_DIGEST = "c89a7708ced8cc2465a9f1dd61747488ff43850d3e7a07c7787c1e5d501a8e24"


def _primary_payload(record: object) -> dict[str, object]:
    weight = record.weight_interval.upper
    return {
        "record_id": record.record_id,
        "side": record.side,
        "partial_sum_index": record.partial_sum_index,
        "adjoint_stage_indices": list(record.adjoint_stage_indices),
        "composition": list(record.composition),
        "base_stage_index": record.base_stage_index,
        "fragment_word": list(record.fragment_word),
        "weight_upper": [weight.numerator, weight.denominator],
    }


@lru_cache(maxsize=2)
def _primary_summary(max_records: int | None) -> RawTheoremSummary:
    digest = hashlib.sha256()
    weights: dict[tuple[int, ...], Fraction] = {}
    count = 0
    records = _iter_raw_theorem_records()
    selected = records if max_records is None else islice(records, max_records)
    for record in selected:
        digest.update(canonical_json_bytes(_primary_payload(record)))
        weights[record.fragment_word] = (
            weights.get(record.fragment_word, Fraction())
            + record.weight_interval.upper
        )
        count += 1
    return RawTheoremSummary(
        record_count=count,
        projected_record_count=PROJECTED_RECORD_COUNT,
        complete=count == PROJECTED_RECORD_COUNT,
        stream_sha256=digest.hexdigest(),
        word_weights=tuple(sorted(weights.items())),
    )


def test_reference_oracle_imports_only_python_standard_library() -> None:
    tree = ast.parse(REFERENCE.read_text(encoding="utf-8"))
    imported_roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_roots.add(node.module.split(".", 1)[0])

    assert imported_roots <= {
        "__future__",
        "argparse",
        "collections",
        "dataclasses",
        "fractions",
        "hashlib",
        "json",
        "math",
        "pathlib",
    }


def test_q_alpha_and_31_stage_schedule_match_primary_exactly() -> None:
    alpha = Cubic(0, 1, 0)
    assert alpha * alpha * alpha == Cubic(4, 0, 0)

    reference = theorem_stages()
    primary_cubic = tuple(reversed(fourth_order_suzuki_cubic_stages()))
    primary_interval = _theorem_interval_stages()

    assert len(reference) == 31
    assert [stage.fragment_index for stage in reference] == [
        stage.fragment_index for stage in primary_cubic
    ]
    assert [
        (
            stage.algebraic_coefficient.a0,
            stage.algebraic_coefficient.a1,
            stage.algebraic_coefficient.a2,
        )
        for stage in reference
    ] == [
        (stage.coefficient.a0, stage.coefficient.a1, stage.coefficient.a2)
        for stage in primary_cubic
    ]
    assert [
        (stage.coefficient_interval.lower, stage.coefficient_interval.upper)
        for stage in reference
    ] == [
        (stage.coefficient.lower, stage.coefficient.upper)
        for stage in primary_interval
    ]


def test_200_record_prefix_matches_primary_record_by_record() -> None:
    reference = tuple(iter_raw_theorem_records(200))
    primary = tuple(islice(_iter_raw_theorem_records(), 200))

    assert len(reference) == len(primary) == 200
    assert [record_payload(record) for record in reference] == [
        _primary_payload(record) for record in primary
    ]
    assert summarize_raw_theorem(200) == _primary_summary(200)
    assert summarize_raw_theorem(200).stream_sha256 == PREFIX_DIGEST


@pytest.mark.slow
def test_full_raw_stream_count_word_weights_and_digest_match_primary() -> None:
    reference = summarize_raw_theorem()
    primary = _primary_summary(None)

    assert reference.record_count == primary.record_count == 61_677
    assert reference.projected_record_count == PROJECTED_RECORD_COUNT
    assert reference.complete and primary.complete
    assert len(reference.word_weights) == len(primary.word_weights) == 4**5
    assert reference.word_weights == primary.word_weights
    assert reference.stream_sha256 == primary.stream_sha256 == FULL_DIGEST


def test_record_decoder_replays_word_and_exact_positive_upper_weight() -> None:
    record = next(iter_raw_theorem_records())
    payload = record_payload(record)
    encoded = canonical_json_bytes(payload)

    assert load_canonical_json(encoded) == payload
    assert decode_record_payload(payload) == record
    assert record.weight_upper > 0


@pytest.mark.parametrize(
    "mutation",
    (
        lambda payload: payload.__setitem__("partial_sum_index", True),
        lambda payload: payload.__setitem__("partial_sum_index", 2.0),
        lambda payload: payload.__setitem__("unknown", 1),
        lambda payload: payload.__setitem__("weight_upper", [2, 2]),
    ),
)
def test_record_decoder_rejects_aliases_unknowns_and_noncanonical_fractions(
    mutation: object,
) -> None:
    payload = copy.deepcopy(record_payload(next(iter_raw_theorem_records())))
    mutation(payload)

    with pytest.raises((TypeError, ValueError)):
        decode_record_payload(payload)


def test_canonical_loader_rejects_bool_float_duplicate_and_pretty_json() -> None:
    payload = record_payload(next(iter_raw_theorem_records()))
    canonical = canonical_json_bytes(payload)
    duplicate = b'{"record_id":"shadow",' + canonical[1:]
    pretty = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()

    with pytest.raises(ValueError, match="duplicate JSON key"):
        load_canonical_json(duplicate)
    with pytest.raises(ValueError, match="canonical byte form"):
        load_canonical_json(pretty)
    for alias in (True, 2.0):
        forged = copy.deepcopy(payload)
        forged["partial_sum_index"] = alias
        with pytest.raises(TypeError, match="numeric aliases"):
            load_canonical_json(canonical_json_bytes(forged))


def test_strict_fraction_and_record_limit_reject_bool_aliases() -> None:
    with pytest.raises(TypeError):
        strict_fraction_pair([True, 1], "weight")
    with pytest.raises(ValueError):
        tuple(iter_raw_theorem_records(True))


def test_cli_runs_without_pythonpath_and_emits_canonical_prefix() -> None:
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, str(REFERENCE), "--max-records", "200"],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr.decode()
    payload = load_canonical_json(completed.stdout)
    assert payload["record_count"] == 200
    assert payload["coverage_status"] == "bounded_prefix"
    assert payload["stream_sha256"] == PREFIX_DIGEST
