"""Independent standard-library verifier for grouped-XXZ per-Delta artifacts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from fractions import Fraction
from io import BytesIO
from math import isqrt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from reference_grouped_xxz_theorem import (
    PROJECTED_RECORD_COUNT,
    canonical_json_bytes,
    fraction_pair,
    summarize_compressed_raw_theorem,
)

ROOT = SCRIPT_DIR.parent
SUMMARY_SCHEMA = "grouped_xxz_per_delta_summary_v1"
WITNESS_SCHEMA = "grouped_xxz_per_delta_witness_v1"
MERGED_SUMMARY_SCHEMA = "grouped_xxz_two_point_pilot_summary_v1"
MERGED_WITNESS_SCHEMA = "grouped_xxz_two_point_pilot_witness_v1"
MERGED_METHOD = "two_point_compressed_grouped_xxz_pilot"
THEOREM_IDENTIFIER = "published_high_order_triangle_v1"
DUHAMEL_CONVENTION = "published_triangle_duhamel_1_over_5_factorial"
ERROR_FORMULA = "E_r=K/r^4_for_T=1"
GROUPING_ALGORITHM = "deterministic_pair_only_bitset_v1"
METHOD = "compressed_d4_direct_finite_high_order_theorem_grouped_norm"
FORMULA = "five_copy_fourth_order_suzuki_four_matchings"
NORMALIZATION = "(XX+YY+delta*ZZ)/4"
PRIMARY_METRIC = "merged_group_exponentials"
SOURCE_PATHS = (
    "pyproject.toml",
    "requirements-reproducibility.txt",
    "scripts/compile_grouped_xxz.py",
    "src/trottercert/__init__.py",
    "src/trottercert/algebra.py",
    "src/trottercert/cubic_field.py",
    "src/trottercert/grouped_xxz.py",
    "src/trottercert/grouped_xxz_artifact.py",
    "src/trottercert/grouped_xxz_compressed.py",
    "src/trottercert/grouped_xxz_orbits.py",
    "src/trottercert/higher_order.py",
    "src/trottercert/intervals.py",
    "src/trottercert/lattice.py",
    "src/trottercert/local_commutators.py",
    "src/trottercert/rigorous_fourth.py",
)

Word = tuple[int, int, int, int, int]
Mask = tuple[int, int]
Term = tuple[int, int, Fraction, Fraction]


def _expect_keys(value: object, keys: set[str], label: str) -> dict[str, object]:
    if not isinstance(value, dict) or set(value) != keys:
        raise ValueError(f"{label} field set mismatch")
    return value


def _reject_floats(value: object, path: str = "payload") -> None:
    if isinstance(value, float):
        raise TypeError(f"{path}: floating-point values are forbidden")
    if isinstance(value, list):
        for index, item in enumerate(value):
            _reject_floats(item, f"{path}[{index}]")
    elif isinstance(value, Mapping):
        for key, item in value.items():
            _reject_floats(item, f"{path}.{key}")


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _load_canonical(raw: bytes, label: str) -> dict[str, object]:
    if not isinstance(raw, bytes):
        raise TypeError(f"{label} must be bytes")
    try:
        value = json.loads(
            raw,
            object_pairs_hook=_unique_object,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ValueError(f"invalid JSON constant: {token}")
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is invalid JSON") from error
    _reject_floats(value, label)
    if not isinstance(value, dict) or canonical_json_bytes(value) != raw:
        raise ValueError(f"{label} is not canonical JSON")
    return value


def _integer(value: object, label: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def _fraction(value: object, label: str) -> Fraction:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must be a fraction pair")
    numerator = _integer_signed(value[0], f"{label}.numerator")
    denominator = _integer(value[1], f"{label}.denominator", 1)
    result = Fraction(numerator, denominator)
    if value != fraction_pair(result):
        raise ValueError(f"{label} is not canonical and reduced")
    return result


def _integer_signed(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    return value


def _optional_fraction(value: object, label: str) -> Fraction | None:
    return None if value is None else _fraction(value, label)


def _word(value: object, label: str) -> Word:
    if not isinstance(value, list) or len(value) != 5:
        raise ValueError(f"{label} must be a five-letter word")
    decoded = tuple(_integer(item, f"{label} entry") for item in value)
    if any(item >= 4 for item in decoded):
        raise ValueError(f"{label} fragment index is outside 0..3")
    return decoded  # type: ignore[return-value]


def _sha(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{label} is not lowercase SHA-256")
    return value


def _payload_digest(payload: dict[str, object]) -> str:
    unsigned = dict(payload)
    unsigned["payload_sha256"] = ""
    return hashlib.sha256(canonical_json_bytes(unsigned)).hexdigest()


def _deterministic_gzip(payload: bytes) -> bytes:
    output = BytesIO()
    with gzip.GzipFile(
        filename="", mode="wb", compresslevel=9, fileobj=output, mtime=0
    ) as member:
        member.write(payload)
    return output.getvalue()


def _source_closure(value: object, root: Path) -> dict[str, object]:
    closure = _expect_keys(
        value,
        {"algorithm", "closure_sha256", "files", "source_commit"},
        "source closure",
    )
    if closure["algorithm"] != "sha256":
        raise ValueError("source closure algorithm mismatch")
    commit = closure["source_commit"]
    if not isinstance(commit, str) or not commit:
        raise ValueError("source closure commit is invalid")
    files = closure["files"]
    if not isinstance(files, list) or len(files) != len(SOURCE_PATHS):
        raise ValueError("source closure file coverage mismatch")
    observed_paths: list[str] = []
    for row_value in files:
        row = _expect_keys(row_value, {"path", "sha256"}, "source file")
        path = row["path"]
        if not isinstance(path, str):
            raise TypeError("source path is invalid")
        observed_paths.append(path)
        expected_hash = hashlib.sha256((root / path).read_bytes()).hexdigest()
        if _sha(row["sha256"], f"source hash {path}") != expected_hash:
            raise ValueError(f"source file hash differs from current bytes: {path}")
    if tuple(observed_paths) != SOURCE_PATHS:
        raise ValueError("source closure paths differ from frozen closure")
    unsigned = {"algorithm": "sha256", "files": files, "source_commit": commit}
    if closure["closure_sha256"] != hashlib.sha256(
        canonical_json_bytes(unsigned)
    ).hexdigest():
        raise ValueError("source closure digest mismatch")
    return closure


def _site(x: int, y: int) -> int:
    return (y % 4) * 4 + (x % 4)


def _bond(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _matchings() -> tuple[tuple[tuple[int, int], ...], ...]:
    groups: list[set[tuple[int, int]]] = [set(), set(), set(), set()]
    for y in range(4):
        for x in range(0, 4, 2):
            groups[0].add(_bond(_site(x, y), _site(x + 1, y)))
            groups[1].add(_bond(_site(x + 1, y), _site(x + 2, y)))
    for x in range(4):
        for y in range(0, 4, 2):
            groups[2].add(_bond(_site(x, y), _site(x, y + 1)))
            groups[3].add(_bond(_site(x, y + 1), _site(x, y + 2)))
    result = tuple(tuple(sorted(group)) for group in groups)
    all_bonds = {_bond(_site(x, y), _site(x + 1, y)) for y in range(4) for x in range(4)}
    all_bonds |= {_bond(_site(x, y), _site(x, y + 1)) for y in range(4) for x in range(4)}
    if any(len(group) != 8 for group in result):
        raise ArithmeticError("each matching must contain eight bonds")
    if any(len({site for bond in group for site in bond}) != 16 for group in result):
        raise ArithmeticError("matching bonds are not site-disjoint")
    if set().union(*(set(group) for group in result)) != all_bonds or len(all_bonds) != 32:
        raise ArithmeticError("four matchings do not cover the 32 periodic bonds")
    return result


def _symmetries() -> tuple[tuple[str, tuple[int, ...], tuple[int, ...]], ...]:
    actions = (
        ("identity", lambda x, y: (x, y)),
        ("rotate_90", lambda x, y: (-y, x)),
        ("rotate_180", lambda x, y: (-x, -y)),
        ("rotate_270", lambda x, y: (y, -x)),
        ("reflect_x_axis", lambda x, y: (x, -y)),
        ("reflect_y_axis", lambda x, y: (-x, y)),
        ("reflect_main_diagonal", lambda x, y: (y, x)),
        ("reflect_anti_diagonal", lambda x, y: (-y, -x)),
    )
    matchings = tuple(frozenset(group) for group in _matchings())
    result = []
    for name, action in actions:
        sites = tuple(_site(*action(site % 4, site // 4)) for site in range(16))
        fragments = []
        for group in matchings:
            transported = frozenset(_bond(sites[a], sites[b]) for a, b in group)
            if transported not in matchings:
                raise ArithmeticError("D4 action does not preserve matchings")
            fragments.append(matchings.index(transported))
        result.append((name, sites, tuple(fragments)))
    if len({sites for _, sites, _ in result}) != 8:
        raise ArithmeticError("D4 site actions are not distinct")
    return tuple(result)


def _transport_word(word: Word, fragments: tuple[int, ...]) -> Word:
    return tuple(fragments[index] for index in word)  # type: ignore[return-value]


def _word_binding(word: Word) -> tuple[Word, str]:
    symmetries = _symmetries()
    representative = min(_transport_word(word, fragments) for _, _, fragments in symmetries)
    for name, _, fragments in symmetries:
        if _transport_word(representative, fragments) == word:
            return representative, name
    raise ArithmeticError("word is missing from its D4 orbit")


def _anticommutes(left: Mask, right: Mask) -> bool:
    return bool(
        ((left[0] & right[1]).bit_count() + (left[1] & right[0]).bit_count()) & 1
    )


def _product_phase(left: Mask, right: Mask) -> tuple[int, Mask]:
    result = (left[0] ^ right[0], left[1] ^ right[1])
    exponent = (
        (left[0] & left[1]).bit_count()
        + (right[0] & right[1]).bit_count()
        + 2 * (left[1] & right[0]).bit_count()
        - (result[0] & result[1]).bit_count()
    ) % 4
    return exponent, result


def _times_i(value: tuple[int, int], exponent: int) -> tuple[int, int]:
    real, imag = value
    return ((real, imag), (-imag, real), (-real, -imag), (imag, -real))[exponent]


def _fragment_terms(
    bonds: Sequence[tuple[int, int]], delta: Fraction
) -> tuple[int, dict[Mask, int]]:
    denominator = 4 * delta.denominator
    weights = (delta.denominator, delta.denominator, delta.numerator)
    terms: dict[Mask, int] = {}
    for left, right in bonds:
        sites = (1 << left) | (1 << right)
        for mask, numerator in zip(
            ((sites, 0), (sites, sites), (0, sites)), weights, strict=True
        ):
            if numerator:
                terms[mask] = terms.get(mask, 0) + numerator
    return denominator, terms


def _commutator_terms(word: Word, delta: Fraction) -> tuple[Term, ...]:
    fragments = tuple(_fragment_terms(group, delta) for group in _matchings())
    denominators = {denominator for denominator, _ in fragments}
    if len(denominators) != 1:
        raise ArithmeticError("fragment denominators differ")
    denominator = denominators.pop()
    cache: dict[tuple[int, ...], dict[Mask, tuple[int, int]]] = {}

    def evaluate(key: tuple[int, ...]) -> dict[Mask, tuple[int, int]]:
        if key in cache:
            return cache[key]
        if len(key) == 1:
            result = {mask: (value, 0) for mask, value in fragments[key[0]][1].items()}
        else:
            result: dict[Mask, tuple[int, int]] = {}
            inner = evaluate(key[1:])
            for left, left_value in fragments[key[0]][1].items():
                for right, right_value in inner.items():
                    if not _anticommutes(left, right):
                        continue
                    exponent, product = _product_phase(left, right)
                    phased = _times_i(right_value, exponent)
                    contribution = (2 * left_value * phased[0], 2 * left_value * phased[1])
                    previous = result.get(product, (0, 0))
                    updated = (previous[0] + contribution[0], previous[1] + contribution[1])
                    if updated == (0, 0):
                        result.pop(product, None)
                    else:
                        result[product] = updated
        cache[key] = result
        return result

    scale = denominator ** len(word)
    return tuple(
        (mask[0], mask[1], Fraction(value[0], scale), Fraction(value[1], scale))
        for mask, value in sorted(evaluate(word).items())
        if value != (0, 0)
    )


def _sqrt_interval(value: Fraction) -> tuple[Fraction, Fraction]:
    if value < 0:
        raise ValueError("negative square-root radicand")
    if value == 0:
        return Fraction(), Fraction()
    scale = 10**30
    quotient = value.numerator * scale * scale // value.denominator
    lower_integer = isqrt(quotient)
    lower = Fraction(lower_integer, scale)
    upper_integer = lower_integer + (lower * lower < value)
    return lower, Fraction(upper_integer, scale)


def _groups(terms: tuple[Term, ...]) -> tuple[tuple[int, ...], ...]:
    ordered = sorted(range(len(terms)), key=lambda index: (-abs(terms[index][2]), terms[index][:2]))
    unmatched = set(ordered)
    groups = []
    for index in ordered:
        if index not in unmatched:
            continue
        unmatched.remove(index)
        partner = next(
            (
                candidate
                for candidate in ordered
                if candidate in unmatched
                and _anticommutes(terms[index][:2], terms[candidate][:2])
            ),
            None,
        )
        if partner is None:
            groups.append((index,))
        else:
            unmatched.remove(partner)
            groups.append((index, partner))
    return tuple(groups)


def _decode_spec(value: object) -> tuple[dict[str, object], Fraction]:
    spec = _expect_keys(
        value,
        {
            "boundary",
            "delta",
            "formula_identifier",
            "length",
            "normalization",
            "primary_metric",
            "stage_count",
            "time",
            "tolerance",
        },
        "spec",
    )
    delta = _fraction(spec["delta"], "spec.delta")
    expected = {
        "boundary": "periodic",
        "delta": fraction_pair(delta),
        "formula_identifier": FORMULA,
        "length": 4,
        "normalization": NORMALIZATION,
        "primary_metric": PRIMARY_METRIC,
        "stage_count": 31,
        "time": [1, 1],
        "tolerance": [1, 1_000_000],
    }
    if spec != expected:
        raise ValueError("spec differs from frozen 4x4 grouped XXZ setup")
    return expected, delta


def _decode_terms(value: object) -> tuple[Term, ...]:
    if not isinstance(value, list):
        raise TypeError("terms must be a list")
    terms = []
    for row in value:
        if not isinstance(row, list) or len(row) != 4:
            raise ValueError("term schema mismatch")
        x_mask = _integer(row[0], "term.x_mask")
        z_mask = _integer(row[1], "term.z_mask")
        if x_mask >= 1 << 16 or z_mask >= 1 << 16:
            raise ValueError("term mask escapes 16-site universe")
        real = _fraction(row[2], "term.real")
        imag = _fraction(row[3], "term.imag")
        if not real and not imag:
            raise ValueError("zero Pauli term is forbidden")
        terms.append((x_mask, z_mask, real, imag))
    result = tuple(terms)
    if result != tuple(sorted(result, key=lambda term: term[:2])):
        raise ValueError("terms are not canonically mask-sorted")
    if len({term[:2] for term in result}) != len(result):
        raise ValueError("duplicate Pauli mask")
    return result


def _verify_blocks(
    value: object, delta: Fraction, representatives: tuple[Word, ...]
) -> tuple[tuple[dict[str, object], tuple[Term, ...], Fraction, Fraction], ...]:
    if not isinstance(value, list) or len(value) != len(representatives):
        raise ValueError("representative block coverage mismatch")
    decoded = []
    for raw, representative in zip(value, representatives, strict=True):
        block = _expect_keys(
            raw,
            {"grouped_norm", "groups", "representative", "terms", "triangle_norm"},
            "representative block",
        )
        if _word(block["representative"], "representative") != representative:
            raise ValueError("representative ordering or identity mismatch")
        terms = _decode_terms(block["terms"])
        if terms != _commutator_terms(representative, delta):
            raise ValueError("representative Pauli map mismatch")
        expected_groups = _groups(terms)
        raw_groups = block["groups"]
        if not isinstance(raw_groups, list) or len(raw_groups) != len(expected_groups):
            raise ValueError("deterministic group coverage mismatch")
        grouped_norm = Fraction()
        for raw_group, indices in zip(raw_groups, expected_groups, strict=True):
            group = _expect_keys(
                raw_group,
                {"norm_interval", "squared_norm", "term_indices"},
                "group",
            )
            if group["term_indices"] != list(indices):
                raise ValueError("pair-only earliest matching differs")
            if len(indices) == 2 and not _anticommutes(terms[indices[0]][:2], terms[indices[1]][:2]):
                raise ValueError("submitted pair does not anticommute")
            squared = sum((terms[index][2] ** 2 for index in indices), Fraction())
            interval = _sqrt_interval(squared)
            if _fraction(group["squared_norm"], "group.squared_norm") != squared:
                raise ValueError("group squared norm mismatch")
            raw_interval = group["norm_interval"]
            if not isinstance(raw_interval, list) or len(raw_interval) != 2:
                raise ValueError("group norm interval schema mismatch")
            decoded_interval = (
                _fraction(raw_interval[0], "group.norm_interval.lower"),
                _fraction(raw_interval[1], "group.norm_interval.upper"),
            )
            if decoded_interval != interval:
                raise ValueError("group outward sqrt interval mismatch")
            grouped_norm += interval[1]
        triangle_norm = sum((abs(term[2]) for term in terms), Fraction())
        if _fraction(block["grouped_norm"], "block.grouped_norm") != grouped_norm:
            raise ValueError("block grouped norm mismatch")
        if _fraction(block["triangle_norm"], "block.triangle_norm") != triangle_norm:
            raise ValueError("block triangle norm mismatch")
        decoded.append((block, terms, grouped_norm, triangle_norm))
    return tuple(decoded)


def _ledger_digest(
    ledger: dict[str, object],
    spec: dict[str, object],
    rows: list[dict[str, object]],
    blocks: tuple[tuple[dict[str, object], tuple[Term, ...], Fraction, Fraction], ...],
) -> str:
    digest = hashlib.sha256()
    digest.update(
        canonical_json_bytes(
            {
                "complete": ledger["complete"],
                "digest_domain": "grouped_xxz_compressed_ledger_v1",
                "duhamel_convention": ledger["duhamel_convention"],
                "factorial_denominator": ledger["factorial_denominator"],
                "finite_step_error_formula": ledger["finite_step_error_formula"],
                "grouped_constant": ledger["grouped_constant"],
                "grouping_algorithm": ledger["grouping_algorithm"],
                "max_degree": ledger["max_degree"],
                "order": ledger["order"],
                "center": ledger["center"],
                "projected_record_count": ledger["projected_record_count"],
                "raw_record_count": ledger["raw_record_count"],
                "raw_stream_sha256": ledger["raw_stream_sha256"],
                "record_limit": ledger["record_limit"],
                "representative_block_count": len(blocks),
                "spec": spec,
                "theorem_identifier": ledger["theorem_identifier"],
                "triangle_constant": ledger["triangle_constant"],
                "word_weight_count": len(rows),
            }
        )
    )
    for row in rows:
        digest.update(
            canonical_json_bytes(
                {
                    "actual_word": row["actual_word"],
                    "raw_weight_upper": row["raw_weight_upper"],
                    "representative": row["representative"],
                    "row_domain": "actual_word_weight_v1",
                    "symmetry_name": row["symmetry_name"],
                }
            )
        )
    for block, terms, grouped_norm, triangle_norm in blocks:
        groups = block["groups"]
        assert isinstance(groups, list)
        digest.update(
            canonical_json_bytes(
                {
                    "block_domain": "representative_norm_block_v1",
                    "group_count": len(groups),
                    "grouped_norm": fraction_pair(grouped_norm),
                    "representative": block["representative"],
                    "term_count": len(terms),
                    "triangle_norm": fraction_pair(triangle_norm),
                }
            )
        )
        for term in terms:
            digest.update(
                canonical_json_bytes(
                    {
                        "term": [
                            term[0],
                            term[1],
                            fraction_pair(term[2]),
                            fraction_pair(term[3]),
                        ]
                    }
                )
            )
        for group in groups:
            assert isinstance(group, dict)
            indices = group["term_indices"]
            assert isinstance(indices, list)
            digest.update(
                canonical_json_bytes(
                    {
                        "group_masks": [[terms[index][0], terms[index][1]] for index in indices],
                        "norm_interval": group["norm_interval"],
                        "squared_norm": group["squared_norm"],
                    }
                )
            )
    return digest.hexdigest()


def _minimal_steps(constant: Fraction, tolerance: Fraction) -> int:
    if constant <= tolerance:
        return 1
    lower, upper = 1, 2
    while constant > tolerance * upper**4:
        lower, upper = upper, 2 * upper
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if constant <= tolerance * middle**4:
            upper = middle
        else:
            lower = middle
    return upper


def _verify_per_delta_artifact_bytes(
    summary_bytes: bytes,
    witness_gzip_bytes: bytes,
    *,
    root: Path = ROOT,
    expected_raw_records: int = PROJECTED_RECORD_COUNT,
) -> tuple[dict[str, object], dict[str, object]]:
    """Verify both artifact members and independently replay every claim."""

    expected_raw_records = _integer(expected_raw_records, "expected_raw_records", 1)
    if expected_raw_records > PROJECTED_RECORD_COUNT:
        raise ValueError("expected raw record count exceeds theorem projection")
    summary = _load_canonical(summary_bytes, "summary")
    _expect_keys(
        summary,
        {
            "adjacent_steps",
            "bounds",
            "counts",
            "method",
            "payload_sha256",
            "resources",
            "schema",
            "source_closure",
            "spec",
            "status",
            "theorem",
            "witness",
        },
        "summary",
    )
    if summary["schema"] != SUMMARY_SCHEMA or summary["payload_sha256"] != _payload_digest(summary):
        raise ValueError("summary schema or payload digest mismatch")
    binding = _expect_keys(
        summary["witness"],
        {"compression", "file_sha256", "payload_sha256"},
        "witness binding",
    )
    if binding["compression"] != "gzip-9-mtime0-empty-filename":
        raise ValueError("witness compression domain mismatch")
    if _sha(binding["file_sha256"], "witness file hash") != hashlib.sha256(witness_gzip_bytes).hexdigest():
        raise ValueError("witness file hash mismatch")
    try:
        witness_payload_bytes = gzip.decompress(witness_gzip_bytes)
    except (gzip.BadGzipFile, EOFError, OSError) as error:
        raise ValueError("invalid gzip witness") from error
    if _deterministic_gzip(witness_payload_bytes) != witness_gzip_bytes:
        raise ValueError("witness is not one deterministic gzip member")
    if _sha(binding["payload_sha256"], "witness payload hash") != hashlib.sha256(witness_payload_bytes).hexdigest():
        raise ValueError("witness payload hash mismatch")
    witness = _load_canonical(witness_payload_bytes, "witness")
    _expect_keys(
        witness,
        {
            "certificate",
            "ledger",
            "payload_sha256",
            "representative_blocks",
            "schema",
            "source_closure",
            "word_weights",
        },
        "witness",
    )
    if witness["schema"] != WITNESS_SCHEMA or witness["payload_sha256"] != _payload_digest(witness):
        raise ValueError("witness schema or payload digest mismatch")
    summary_closure = _source_closure(summary["source_closure"], root)
    witness_closure = _source_closure(witness["source_closure"], root)
    if summary_closure != witness_closure:
        raise ValueError("summary and witness source closures differ")

    spec, delta = _decode_spec(summary["spec"])
    ledger = _expect_keys(
        witness["ledger"],
        {
            "center",
            "complete",
            "duhamel_convention",
            "factorial_denominator",
            "finite_step_error_formula",
            "grouped_constant",
            "grouping_algorithm",
            "ledger_digest",
            "max_degree",
            "order",
            "projected_record_count",
            "raw_record_count",
            "raw_stream_sha256",
            "record_limit",
            "spec",
            "theorem_identifier",
            "triangle_constant",
        },
        "ledger",
    )
    ledger_spec, ledger_delta = _decode_spec(ledger["spec"])
    if ledger_spec != spec or ledger_delta != delta:
        raise ValueError("ledger and summary specs differ")
    expected_metadata = {
        "center": 20,
        "complete": True,
        "duhamel_convention": DUHAMEL_CONVENTION,
        "factorial_denominator": 120,
        "finite_step_error_formula": ERROR_FORMULA,
        "grouping_algorithm": GROUPING_ALGORITHM,
        "max_degree": 5,
        "order": 4,
        "projected_record_count": expected_raw_records,
        "raw_record_count": expected_raw_records,
        "record_limit": None,
        "theorem_identifier": THEOREM_IDENTIFIER,
    }
    if any(ledger[key] != value for key, value in expected_metadata.items()):
        raise ValueError("ledger theorem metadata or full coverage mismatch")
    oracle = summarize_compressed_raw_theorem(
        None if expected_raw_records == PROJECTED_RECORD_COUNT else expected_raw_records
    )
    if ledger["raw_stream_sha256"] != oracle.stream_sha256:
        raise ValueError("ledger raw theorem stream digest mismatch")
    rows_value = witness["word_weights"]
    if not isinstance(rows_value, list):
        raise TypeError("word weights must be a list")
    expected_weights = dict(oracle.word_weights)
    decoded_rows: list[dict[str, object]] = []
    observed_weights: dict[Word, Fraction] = {}
    representatives = []
    previous_word: Word | None = None
    for value in rows_value:
        row = _expect_keys(
            value,
            {"actual_word", "raw_weight_upper", "representative", "symmetry_name"},
            "word-weight row",
        )
        actual = _word(row["actual_word"], "actual_word")
        if previous_word is not None and actual <= previous_word:
            raise ValueError("actual words are not strictly sorted")
        previous_word = actual
        weight = _fraction(row["raw_weight_upper"], "raw_weight_upper")
        representative, symmetry_name = _word_binding(actual)
        if _word(row["representative"], "representative") != representative or row["symmetry_name"] != symmetry_name:
            raise ValueError("D4 representative or transport binding mismatch")
        observed_weights[actual] = weight
        representatives.append(representative)
        decoded_rows.append(row)
    if observed_weights != expected_weights:
        raise ValueError("actual-word theorem weights mismatch")
    expected_representatives = tuple(sorted(set(representatives)))
    blocks = _verify_blocks(witness["representative_blocks"], delta, expected_representatives)
    norm_by_rep = {
        representative: (block[2], block[3])
        for representative, block in zip(expected_representatives, blocks, strict=True)
    }
    grouped = sum(
        (weight * norm_by_rep[_word_binding(word)[0]][0] for word, weight in observed_weights.items()),
        Fraction(),
    ) / 120
    triangle = sum(
        (weight * norm_by_rep[_word_binding(word)[0]][1] for word, weight in observed_weights.items()),
        Fraction(),
    ) / 120
    if _fraction(ledger["grouped_constant"], "ledger.grouped_constant") != grouped:
        raise ValueError("ledger grouped constant is not literal actual-word sum")
    if _fraction(ledger["triangle_constant"], "ledger.triangle_constant") != triangle:
        raise ValueError("ledger triangle constant is not literal actual-word sum")
    expected_ledger_digest = _ledger_digest(ledger, spec, decoded_rows, blocks)
    if _sha(ledger["ledger_digest"], "ledger digest") != expected_ledger_digest:
        raise ValueError("compressed ledger digest mismatch")

    candidate_steps = _minimal_steps(grouped, Fraction(1, 1_000_000))
    baseline_steps = _minimal_steps(triangle, Fraction(1, 1_000_000))
    candidate_error = grouped / candidate_steps**4
    baseline_error = triangle / baseline_steps**4
    candidate_previous = None if candidate_steps == 1 else grouped / (candidate_steps - 1) ** 4
    baseline_previous = None if baseline_steps == 1 else triangle / (baseline_steps - 1) ** 4
    certificate = _expect_keys(
        witness["certificate"],
        {
            "baseline_error",
            "baseline_previous_error",
            "baseline_resources",
            "baseline_steps",
            "bond_growth",
            "candidate_error",
            "candidate_previous_error",
            "candidate_resources",
            "candidate_steps",
            "cell_base",
            "center",
            "duhamel_convention",
            "factorial_denominator",
            "finite_step_error_formula",
            "grouped_constant",
            "grouping_algorithm",
            "ledger_digest",
            "method",
            "order",
            "raw_record_count",
            "raw_stream_sha256",
            "spec",
            "status",
            "theorem_identifier",
            "triangle_constant",
        },
        "certificate",
    )
    delta_abs = abs(delta)
    expected_certificate = {
        "baseline_error": fraction_pair(baseline_error),
        "baseline_previous_error": None if baseline_previous is None else fraction_pair(baseline_previous),
        "baseline_resources": 30 * baseline_steps + 1,
        "baseline_steps": baseline_steps,
        "bond_growth": fraction_pair(max(Fraction(1), (1 + delta_abs) / 2)),
        "candidate_error": fraction_pair(candidate_error),
        "candidate_previous_error": None if candidate_previous is None else fraction_pair(candidate_previous),
        "candidate_resources": 30 * candidate_steps + 1,
        "candidate_steps": candidate_steps,
        "cell_base": fraction_pair((2 + delta_abs) / 2),
        "center": 20,
        "duhamel_convention": DUHAMEL_CONVENTION,
        "factorial_denominator": 120,
        "finite_step_error_formula": ERROR_FORMULA,
        "grouped_constant": fraction_pair(grouped),
        "grouping_algorithm": GROUPING_ALGORITHM,
        "ledger_digest": expected_ledger_digest,
        "method": METHOD,
        "order": 4,
        "raw_record_count": expected_raw_records,
        "raw_stream_sha256": oracle.stream_sha256,
        "spec": spec,
        "status": "certified",
        "theorem_identifier": THEOREM_IDENTIFIER,
        "triangle_constant": fraction_pair(triangle),
    }
    if certificate != expected_certificate:
        raise ValueError("certificate differs from independent finite-step replay")
    if candidate_error > Fraction(1, 1_000_000) or baseline_error > Fraction(1, 1_000_000):
        raise ValueError("accepted step does not meet tolerance")
    if candidate_previous is not None and candidate_previous <= Fraction(1, 1_000_000):
        raise ValueError("candidate previous-step minimality fails")
    if baseline_previous is not None and baseline_previous <= Fraction(1, 1_000_000):
        raise ValueError("baseline previous-step minimality fails")

    expected_counts = {
        "actual_words": len(decoded_rows),
        "pair_groups": sum(len(group["term_indices"]) == 2 for block, *_ in blocks for group in block["groups"]),
        "projected_raw_records": expected_raw_records,
        "raw_records": expected_raw_records,
        "representative_blocks": len(blocks),
        "representative_terms": sum(len(terms) for _, terms, _, _ in blocks),
        "singleton_groups": sum(len(group["term_indices"]) == 1 for block, *_ in blocks for group in block["groups"]),
    }
    expected_summary = {
        "adjacent_steps": {
            "baseline": {
                "accepted_error": fraction_pair(baseline_error),
                "previous_error": None if baseline_previous is None else fraction_pair(baseline_previous),
                "steps": baseline_steps,
            },
            "candidate": {
                "accepted_error": fraction_pair(candidate_error),
                "previous_error": None if candidate_previous is None else fraction_pair(candidate_previous),
                "steps": candidate_steps,
            },
        },
        "bounds": {"grouped_constant": fraction_pair(grouped), "triangle_constant": fraction_pair(triangle)},
        "counts": expected_counts,
        "method": METHOD,
        "payload_sha256": summary["payload_sha256"],
        "resources": {"baseline": 30 * baseline_steps + 1, "candidate": 30 * candidate_steps + 1, "metric": PRIMARY_METRIC},
        "schema": SUMMARY_SCHEMA,
        "source_closure": summary_closure,
        "spec": spec,
        "status": "certified",
        "theorem": {
            "center": 20,
            "duhamel_convention": DUHAMEL_CONVENTION,
            "factorial_denominator": 120,
            "finite_step_error_formula": ERROR_FORMULA,
            "grouping_algorithm": GROUPING_ALGORITHM,
            "identifier": THEOREM_IDENTIFIER,
            "order": 4,
        },
        "witness": binding,
    }
    if summary != expected_summary:
        raise ValueError("summary differs from independently verified witness")
    return summary, witness


def verify_merged_artifact_bytes(
    summary_bytes: bytes,
    witness_gzip_bytes: bytes,
    *,
    root: Path = ROOT,
    expected_raw_records: int = PROJECTED_RECORD_COUNT,
) -> tuple[tuple[dict[str, object], dict[str, object]], ...]:
    """Verify the canonical two-point pilot and both embedded artifacts."""

    expected_raw_records = _integer(expected_raw_records, "expected_raw_records", 1)
    if expected_raw_records > PROJECTED_RECORD_COUNT:
        raise ValueError("expected raw record count exceeds theorem projection")
    summary = _load_canonical(summary_bytes, "merged summary")
    _expect_keys(
        summary,
        {
            "method",
            "payload_sha256",
            "rows",
            "schema",
            "source_closure",
            "status",
            "witness",
        },
        "merged summary",
    )
    if (
        summary["schema"] != MERGED_SUMMARY_SCHEMA
        or summary["method"] != MERGED_METHOD
        or summary["status"] != "certified"
        or summary["payload_sha256"] != _payload_digest(summary)
    ):
        raise ValueError("merged summary identity or payload digest mismatch")
    closure = _source_closure(summary["source_closure"], root)
    binding = _expect_keys(
        summary["witness"],
        {"compression", "file_sha256", "payload_sha256"},
        "merged witness binding",
    )
    if binding["compression"] != "gzip-9-mtime0-empty-filename":
        raise ValueError("merged witness compression domain mismatch")
    if _sha(binding["file_sha256"], "merged witness file hash") != hashlib.sha256(
        witness_gzip_bytes
    ).hexdigest():
        raise ValueError("merged witness file hash mismatch")
    try:
        witness_payload_bytes = gzip.decompress(witness_gzip_bytes)
    except (gzip.BadGzipFile, EOFError, OSError) as error:
        raise ValueError("invalid merged gzip witness") from error
    if _deterministic_gzip(witness_payload_bytes) != witness_gzip_bytes:
        raise ValueError("merged witness is not one deterministic gzip member")
    if _sha(
        binding["payload_sha256"], "merged witness payload hash"
    ) != hashlib.sha256(witness_payload_bytes).hexdigest():
        raise ValueError("merged witness payload hash mismatch")
    witness = _load_canonical(witness_payload_bytes, "merged witness")
    _expect_keys(
        witness,
        {"payload_sha256", "per_delta_artifacts", "schema", "source_closure"},
        "merged witness",
    )
    if (
        witness["schema"] != MERGED_WITNESS_SCHEMA
        or witness["payload_sha256"] != _payload_digest(witness)
    ):
        raise ValueError("merged witness identity or payload digest mismatch")
    if _source_closure(witness["source_closure"], root) != closure:
        raise ValueError("merged source closures differ")
    embedded = witness["per_delta_artifacts"]
    if not isinstance(embedded, list) or len(embedded) != 2:
        raise ValueError("merged witness must contain exactly two artifacts")

    verified: list[tuple[dict[str, object], dict[str, object]]] = []
    expected_rows = []
    deltas = []
    expected_embedded = []
    for item_value in embedded:
        item = _expect_keys(
            item_value,
            {"summary", "witness_payload"},
            "embedded per-Delta artifact",
        )
        per_summary = item["summary"]
        per_witness = item["witness_payload"]
        if not isinstance(per_summary, dict) or not isinstance(per_witness, dict):
            raise TypeError("embedded per-Delta payloads must be mappings")
        per_summary_bytes = canonical_json_bytes(per_summary)
        per_witness_payload_bytes = canonical_json_bytes(per_witness)
        per_witness_gzip = _deterministic_gzip(per_witness_payload_bytes)
        verified_summary, verified_witness = _verify_per_delta_artifact_bytes(
            per_summary_bytes,
            per_witness_gzip,
            root=root,
            expected_raw_records=expected_raw_records,
        )
        if verified_summary != per_summary or verified_witness != per_witness:
            raise ValueError("embedded artifact differs from independent replay")
        per_closure = _source_closure(per_summary["source_closure"], root)
        if per_closure != closure:
            raise ValueError("embedded and merged source closures differ")
        per_binding = _expect_keys(
            per_summary["witness"],
            {"compression", "file_sha256", "payload_sha256"},
            "embedded witness binding",
        )
        expected_file_hash = hashlib.sha256(per_witness_gzip).hexdigest()
        expected_payload_hash = hashlib.sha256(per_witness_payload_bytes).hexdigest()
        if (
            per_binding["file_sha256"] != expected_file_hash
            or per_binding["payload_sha256"] != expected_payload_hash
        ):
            raise ValueError("embedded original witness hash semantics mismatch")
        spec = _expect_keys(
            per_summary["spec"],
            {
                "boundary",
                "delta",
                "formula_identifier",
                "length",
                "normalization",
                "primary_metric",
                "stage_count",
                "time",
                "tolerance",
            },
            "embedded spec",
        )
        delta = _fraction(spec["delta"], "embedded delta")
        deltas.append(delta)
        adjacent = _expect_keys(
            per_summary["adjacent_steps"], {"baseline", "candidate"}, "adjacent steps"
        )
        resources = _expect_keys(
            per_summary["resources"], {"baseline", "candidate", "metric"}, "resources"
        )
        baseline = _expect_keys(
            adjacent["baseline"], {"accepted_error", "previous_error", "steps"}, "baseline"
        )
        candidate = _expect_keys(
            adjacent["candidate"], {"accepted_error", "previous_error", "steps"}, "candidate"
        )
        candidate_resources = _integer(
            resources["candidate"], "candidate resources", 1
        )
        baseline_resources = _integer(resources["baseline"], "baseline resources", 1)
        if candidate_resources >= baseline_resources:
            raise ValueError("two-point row lacks positive resource transfer")
        expected_rows.append(
            {
                "baseline": {
                    "resources": baseline_resources,
                    "steps": _integer(baseline["steps"], "baseline steps", 1),
                },
                "candidate": {
                    "resources": candidate_resources,
                    "steps": _integer(candidate["steps"], "candidate steps", 1),
                },
                "delta": fraction_pair(delta),
                "method": per_summary["method"],
                "per_delta_bindings": {
                    "summary_sha256": hashlib.sha256(per_summary_bytes).hexdigest(),
                    "witness_file_sha256": expected_file_hash,
                    "witness_payload_sha256": expected_payload_hash,
                },
                "status": per_summary["status"],
            }
        )
        expected_embedded.append(
            {"summary": verified_summary, "witness_payload": verified_witness}
        )
        verified.append((verified_summary, verified_witness))
    if tuple(deltas) != (Fraction(1, 2), Fraction(2)):
        raise ValueError("merged rows must be canonical Delta order {1/2,2}")
    if any(
        row["method"] != METHOD or row["status"] != "certified"
        for row in expected_rows
    ):
        raise ValueError("merged row method or status mismatch")

    expected_witness = {
        "payload_sha256": witness["payload_sha256"],
        "per_delta_artifacts": expected_embedded,
        "schema": MERGED_WITNESS_SCHEMA,
        "source_closure": closure,
    }
    if witness != expected_witness:
        raise ValueError("merged witness differs from exact reconstruction")
    expected_summary = {
        "method": MERGED_METHOD,
        "payload_sha256": summary["payload_sha256"],
        "rows": expected_rows,
        "schema": MERGED_SUMMARY_SCHEMA,
        "source_closure": closure,
        "status": "certified",
        "witness": binding,
    }
    if summary != expected_summary:
        raise ValueError("merged summary differs from exact reconstruction")
    return tuple(verified)


def verify_artifact_bytes(
    summary_bytes: bytes,
    witness_gzip_bytes: bytes,
    *,
    root: Path = ROOT,
    expected_raw_records: int = PROJECTED_RECORD_COUNT,
) -> object:
    """Auto-dispatch one per-Delta or merged two-point artifact."""

    summary = _load_canonical(summary_bytes, "summary dispatch")
    schema = summary.get("schema")
    if schema == SUMMARY_SCHEMA:
        return _verify_per_delta_artifact_bytes(
            summary_bytes,
            witness_gzip_bytes,
            root=root,
            expected_raw_records=expected_raw_records,
        )
    if schema == MERGED_SUMMARY_SCHEMA:
        return verify_merged_artifact_bytes(
            summary_bytes,
            witness_gzip_bytes,
            root=root,
            expected_raw_records=expected_raw_records,
        )
    raise ValueError("unrecognized grouped-XXZ artifact summary schema")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("witness", type=Path)
    arguments = parser.parse_args(argv)
    verify_artifact_bytes(
        arguments.summary.read_bytes(),
        arguments.witness.read_bytes(),
        root=ROOT,
        expected_raw_records=PROJECTED_RECORD_COUNT,
    )
    print(f"verified={arguments.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
