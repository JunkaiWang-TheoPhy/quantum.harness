from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import groupby
from pathlib import Path
from typing import Any

from .cubic_field import Cubic, CubicStage
from .cubic_local import cubic_formula_log_series

FORMULA_NAME = "five_copy_fourth_order_suzuki_four_matchings"


@dataclass(frozen=True, slots=True)
class WordRecord:
    word: tuple[int, ...]
    coefficient: Cubic


@dataclass(frozen=True, slots=True)
class WordGroup:
    ordinal: int
    suffix: tuple[int, ...]
    records: tuple[WordRecord, ...]


@dataclass(frozen=True, slots=True)
class WordManifest:
    degree: int
    shard_index: int
    shard_count: int
    total_groups: int
    total_words: int
    word_set_digest: str
    formula: str
    groups: tuple[WordGroup, ...]
    implementation_sources: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class ManifestShardRecord:
    shard_index: int
    path: str
    sha256: str
    group_count: int
    word_count: int


@dataclass(frozen=True, slots=True)
class ManifestIndex:
    degree: int
    shard_count: int
    total_groups: int
    total_words: int
    word_set_digest: str
    formula: str
    implementation_sources: tuple[tuple[str, str], ...]
    shards: tuple[ManifestShardRecord, ...]


def _canonical_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()


def _sha256(encoded: bytes) -> str:
    return hashlib.sha256(encoded).hexdigest()


def _rational_json(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _parse_rational(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(
            not isinstance(entry, int) or isinstance(entry, bool)
            for entry in value
        )
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if _rational_json(result) != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [
        _rational_json(value.a0),
        _rational_json(value.a1),
        _rational_json(value.a2),
    ]


def _parse_cubic(value: object, field: str) -> Cubic:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must contain three cubic coordinates")
    return Cubic(*(_parse_rational(entry, field) for entry in value))


def _word_text(word: tuple[int, ...]) -> str:
    if any(letter not in range(4) for letter in word):
        raise ValueError("word letters must lie in 0..3")
    return "".join(str(letter) for letter in word)


def _parse_word(value: object, degree: int, field: str) -> tuple[int, ...]:
    if (
        not isinstance(value, str)
        or len(value) != degree
        or any(letter not in "0123" for letter in value)
    ):
        raise ValueError(f"{field} must be a length-{degree} base-four word")
    return tuple(int(letter) for letter in value)


def _validate_digest(value: object, field: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{field} must be a lowercase SHA-256 digest")
    return value


def _sources_tuple(
    implementation_sources: Mapping[str, str],
) -> tuple[tuple[str, str], ...]:
    if not implementation_sources:
        raise ValueError("implementation source manifest is empty")
    sources = tuple(sorted(implementation_sources.items()))
    for path, digest in sources:
        if not isinstance(path, str) or not path:
            raise ValueError("implementation source path is invalid")
        _validate_digest(digest, f"implementation source {path}")
    return sources


def _source_json(sources: tuple[tuple[str, str], ...]) -> dict[str, str]:
    return dict(sources)


def _group_json(group: WordGroup) -> dict[str, object]:
    return {
        "ordinal": group.ordinal,
        "suffix": _word_text(group.suffix),
        "records": [
            {
                "word": _word_text(record.word),
                "coefficient": _cubic_json(record.coefficient),
            }
            for record in group.records
        ],
    }


def _word_set_digest(groups: Sequence[WordGroup]) -> str:
    return _sha256(_canonical_bytes([_group_json(group) for group in groups]))


def _manifest_payload(manifest: WordManifest) -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "issue128_dual_word_manifest",
        "degree": manifest.degree,
        "shard_index": manifest.shard_index,
        "shard_count": manifest.shard_count,
        "total_groups": manifest.total_groups,
        "total_words": manifest.total_words,
        "word_set_digest": manifest.word_set_digest,
        "formula": manifest.formula,
        "groups": [_group_json(group) for group in manifest.groups],
        "implementation_sources": _source_json(
            manifest.implementation_sources
        ),
    }


def _parse_manifest_payload(payload: object) -> WordManifest:
    if not isinstance(payload, dict):
        raise ValueError("word manifest must be a JSON object")
    if payload.get("schema_version") != 1:
        raise ValueError("word manifest schema version mismatch")
    if payload.get("kind") != "issue128_dual_word_manifest":
        raise ValueError("word manifest kind mismatch")
    degree = payload.get("degree")
    shard_index = payload.get("shard_index")
    shard_count = payload.get("shard_count")
    total_groups = payload.get("total_groups")
    total_words = payload.get("total_words")
    if not isinstance(degree, int) or degree < 3 or degree % 2 == 0:
        raise ValueError("word manifest degree must be odd and at least three")
    if (
        not isinstance(shard_count, int)
        or isinstance(shard_count, bool)
        or shard_count < 1
        or not isinstance(shard_index, int)
        or isinstance(shard_index, bool)
        or not 0 <= shard_index < shard_count
    ):
        raise ValueError("word manifest shard metadata is invalid")
    if (
        not isinstance(total_groups, int)
        or isinstance(total_groups, bool)
        or total_groups < 1
        or not isinstance(total_words, int)
        or isinstance(total_words, bool)
        or total_words < 1
    ):
        raise ValueError("word manifest total counts are invalid")
    if payload.get("formula") != FORMULA_NAME:
        raise ValueError("word manifest formula mismatch")
    word_set_digest = _validate_digest(
        payload.get("word_set_digest"), "word-set digest"
    )
    raw_sources = payload.get("implementation_sources")
    if not isinstance(raw_sources, dict):
        raise ValueError("implementation source manifest is missing")
    sources = _sources_tuple(raw_sources)
    raw_groups = payload.get("groups")
    if not isinstance(raw_groups, list):
        raise ValueError("word manifest groups are missing")

    groups: list[WordGroup] = []
    for raw_group in raw_groups:
        if not isinstance(raw_group, dict):
            raise ValueError("word group must be a JSON object")
        ordinal = raw_group.get("ordinal")
        if (
            not isinstance(ordinal, int)
            or isinstance(ordinal, bool)
            or ordinal < 0
            or ordinal >= total_groups
            or ordinal % shard_count != shard_index
        ):
            raise ValueError("word group ordinal is inconsistent with shard")
        suffix = _parse_word(raw_group.get("suffix"), degree - 1, "suffix")
        raw_records = raw_group.get("records")
        if not isinstance(raw_records, list) or not 1 <= len(raw_records) <= 4:
            raise ValueError("word group must contain one to four records")
        records: list[WordRecord] = []
        for raw_record in raw_records:
            if not isinstance(raw_record, dict):
                raise ValueError("word record must be a JSON object")
            word = _parse_word(raw_record.get("word"), degree, "word")
            if word[1:] != suffix:
                raise ValueError("word does not match its suffix group")
            coefficient = _parse_cubic(
                raw_record.get("coefficient"), "word coefficient"
            )
            if coefficient == Cubic.zero():
                raise ValueError("word coefficient must be nonzero")
            records.append(WordRecord(word, coefficient))
        if records != sorted(records, key=lambda record: record.word):
            raise ValueError("word records are not canonically ordered")
        groups.append(WordGroup(ordinal, suffix, tuple(records)))
    if groups != sorted(groups, key=lambda group: group.ordinal):
        raise ValueError("word groups are not canonically ordered")
    if len({group.ordinal for group in groups}) != len(groups):
        raise ValueError("word group ordinals are duplicated")

    return WordManifest(
        degree=degree,
        shard_index=shard_index,
        shard_count=shard_count,
        total_groups=total_groups,
        total_words=total_words,
        word_set_digest=word_set_digest,
        formula=FORMULA_NAME,
        groups=tuple(groups),
        implementation_sources=sources,
    )


def build_word_manifests(
    stages: Sequence[CubicStage],
    degree: int,
    shard_count: int,
    *,
    implementation_sources: Mapping[str, str],
) -> tuple[WordManifest, ...]:
    if degree < 3 or degree % 2 == 0:
        raise ValueError("logarithm degree must be odd and at least three")
    if shard_count < 1:
        raise ValueError("shard count must be positive")
    sources = _sources_tuple(implementation_sources)
    word_map = cubic_formula_log_series(stages, degree)[degree]
    ordered = sorted(
        (
            (word, coefficient)
            for word, coefficient in word_map.items()
            if coefficient != Cubic.zero()
        ),
        key=lambda item: (item[0][1:], item[0][0]),
    )
    groups = tuple(
        WordGroup(
            ordinal=ordinal,
            suffix=suffix,
            records=tuple(
                sorted(
                    (
                        WordRecord(word, coefficient)
                        for word, coefficient in items
                    ),
                    key=lambda record: record.word,
                )
            ),
        )
        for ordinal, (suffix, items) in enumerate(
            (
                (suffix, tuple(items))
                for suffix, items in groupby(
                    ordered, key=lambda item: item[0][1:]
                )
            )
        )
    )
    if not groups:
        raise ValueError("logarithm degree has no nonzero words")
    if any(not 1 <= len(group.records) <= 4 for group in groups):
        raise ValueError("suffix group size exceeds the four-letter alphabet")
    words = [record.word for group in groups for record in group.records]
    if len(words) != len(set(words)):
        raise ValueError("logarithm word map contains duplicate words")
    digest = _word_set_digest(groups)
    total_words = len(words)
    return tuple(
        WordManifest(
            degree=degree,
            shard_index=shard_index,
            shard_count=shard_count,
            total_groups=len(groups),
            total_words=total_words,
            word_set_digest=digest,
            formula=FORMULA_NAME,
            groups=tuple(
                group
                for group in groups
                if group.ordinal % shard_count == shard_index
            ),
            implementation_sources=sources,
        )
        for shard_index in range(shard_count)
    )


def _index_payload(index: ManifestIndex) -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "issue128_dual_word_manifest_index",
        "degree": index.degree,
        "shard_count": index.shard_count,
        "total_groups": index.total_groups,
        "total_words": index.total_words,
        "word_set_digest": index.word_set_digest,
        "formula": index.formula,
        "implementation_sources": _source_json(
            index.implementation_sources
        ),
        "shards": [
            {
                "shard_index": record.shard_index,
                "path": record.path,
                "sha256": record.sha256,
                "group_count": record.group_count,
                "word_count": record.word_count,
            }
            for record in index.shards
        ],
    }


def _parse_index_payload(payload: object) -> ManifestIndex:
    if not isinstance(payload, dict):
        raise ValueError("manifest index must be a JSON object")
    if payload.get("schema_version") != 1:
        raise ValueError("manifest index schema version mismatch")
    if payload.get("kind") != "issue128_dual_word_manifest_index":
        raise ValueError("manifest index kind mismatch")
    degree = payload.get("degree")
    shard_count = payload.get("shard_count")
    total_groups = payload.get("total_groups")
    total_words = payload.get("total_words")
    if not isinstance(degree, int) or degree < 3 or degree % 2 == 0:
        raise ValueError("manifest index degree must be odd and at least three")
    if not isinstance(shard_count, int) or shard_count < 1:
        raise ValueError("manifest index shard count is invalid")
    if not isinstance(total_groups, int) or total_groups < 1:
        raise ValueError("manifest index group count is invalid")
    if not isinstance(total_words, int) or total_words < 1:
        raise ValueError("manifest index word count is invalid")
    if payload.get("formula") != FORMULA_NAME:
        raise ValueError("manifest index formula mismatch")
    word_set_digest = _validate_digest(
        payload.get("word_set_digest"), "word-set digest"
    )
    raw_sources = payload.get("implementation_sources")
    if not isinstance(raw_sources, dict):
        raise ValueError("manifest index sources are missing")
    sources = _sources_tuple(raw_sources)
    raw_shards = payload.get("shards")
    if not isinstance(raw_shards, list) or len(raw_shards) != shard_count:
        raise ValueError("manifest index shard list is incomplete")
    shards: list[ManifestShardRecord] = []
    for raw_shard in raw_shards:
        if not isinstance(raw_shard, dict):
            raise ValueError("manifest shard record must be a JSON object")
        shard_index = raw_shard.get("shard_index")
        path = raw_shard.get("path")
        group_count = raw_shard.get("group_count")
        word_count = raw_shard.get("word_count")
        if not isinstance(shard_index, int) or not 0 <= shard_index < shard_count:
            raise ValueError("manifest shard index is invalid")
        if path != f"manifest-{shard_index:03d}.json":
            raise ValueError("manifest shard path is not canonical")
        if not isinstance(group_count, int) or group_count < 0:
            raise ValueError("manifest shard group count is invalid")
        if not isinstance(word_count, int) or word_count < 0:
            raise ValueError("manifest shard word count is invalid")
        shards.append(
            ManifestShardRecord(
                shard_index=shard_index,
                path=path,
                sha256=_validate_digest(
                    raw_shard.get("sha256"), "manifest shard digest"
                ),
                group_count=group_count,
                word_count=word_count,
            )
        )
    if [record.shard_index for record in shards] != list(range(shard_count)):
        raise ValueError("manifest shard records are not complete and ordered")
    if sum(record.group_count for record in shards) != total_groups:
        raise ValueError("manifest index group total mismatch")
    if sum(record.word_count for record in shards) != total_words:
        raise ValueError("manifest index word total mismatch")
    return ManifestIndex(
        degree=degree,
        shard_count=shard_count,
        total_groups=total_groups,
        total_words=total_words,
        word_set_digest=word_set_digest,
        formula=FORMULA_NAME,
        implementation_sources=sources,
        shards=tuple(shards),
    )


def load_manifest_index(path: Path) -> ManifestIndex:
    try:
        payload = json.loads(path.read_bytes())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read manifest index: {path}") from exc
    return _parse_index_payload(payload)


def load_word_manifest(path: Path, index: ManifestIndex) -> WordManifest:
    try:
        encoded = path.read_bytes()
        payload = json.loads(encoded)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read word manifest: {path}") from exc
    matching = tuple(record for record in index.shards if record.path == path.name)
    if len(matching) != 1:
        raise ValueError("word manifest path is not present exactly once in index")
    if _sha256(encoded) != matching[0].sha256:
        raise ValueError("word manifest digest does not match index")
    manifest = _parse_manifest_payload(payload)
    if (
        manifest.degree != index.degree
        or manifest.shard_count != index.shard_count
        or manifest.total_groups != index.total_groups
        or manifest.total_words != index.total_words
        or manifest.word_set_digest != index.word_set_digest
        or manifest.formula != index.formula
        or manifest.implementation_sources != index.implementation_sources
    ):
        raise ValueError("word manifest configuration does not match index")
    return manifest


def verify_manifest_index(index: ManifestIndex, directory: Path) -> None:
    groups: list[WordGroup] = []
    for record in index.shards:
        path = directory / record.path
        try:
            encoded = path.read_bytes()
        except OSError as exc:
            raise ValueError(f"cannot read child manifest: {path}") from exc
        if _sha256(encoded) != record.sha256:
            raise ValueError(f"child manifest digest mismatch: {record.path}")
        manifest = load_word_manifest(path, index)
        if manifest.shard_index != record.shard_index:
            raise ValueError("child manifest shard index mismatch")
        if len(manifest.groups) != record.group_count:
            raise ValueError("child manifest group count mismatch")
        if sum(len(group.records) for group in manifest.groups) != record.word_count:
            raise ValueError("child manifest word count mismatch")
        groups.extend(manifest.groups)
    ordered = sorted(groups, key=lambda group: group.ordinal)
    if [group.ordinal for group in ordered] != list(range(index.total_groups)):
        raise ValueError("manifest group coverage is incomplete or overlapping")
    words = [record.word for group in ordered for record in group.records]
    if len(words) != index.total_words or len(set(words)) != index.total_words:
        raise ValueError("manifest word coverage is incomplete or duplicated")
    if _word_set_digest(ordered) != index.word_set_digest:
        raise ValueError("manifest word-set digest mismatch")


def write_manifest_set(
    output_directory: Path,
    stages: Sequence[CubicStage],
    degree: int,
    shard_count: int,
    *,
    implementation_sources: Mapping[str, str],
) -> Path:
    if output_directory.exists():
        raise FileExistsError(f"manifest output already exists: {output_directory}")
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}-",
            dir=output_directory.parent,
        )
    )
    try:
        manifests = build_word_manifests(
            stages,
            degree,
            shard_count,
            implementation_sources=implementation_sources,
        )
        records: list[ManifestShardRecord] = []
        for manifest in manifests:
            path = temporary / f"manifest-{manifest.shard_index:03d}.json"
            encoded = _canonical_bytes(_manifest_payload(manifest))
            path.write_bytes(encoded)
            parsed = _parse_manifest_payload(json.loads(encoded))
            if parsed != manifest:
                raise ArithmeticError("word manifest codec round trip failed")
            records.append(
                ManifestShardRecord(
                    shard_index=manifest.shard_index,
                    path=path.name,
                    sha256=_sha256(encoded),
                    group_count=len(manifest.groups),
                    word_count=sum(
                        len(group.records) for group in manifest.groups
                    ),
                )
            )
        first = manifests[0]
        index = ManifestIndex(
            degree=first.degree,
            shard_count=first.shard_count,
            total_groups=first.total_groups,
            total_words=first.total_words,
            word_set_digest=first.word_set_digest,
            formula=first.formula,
            implementation_sources=first.implementation_sources,
            shards=tuple(records),
        )
        index_path = temporary / "index.json"
        index_path.write_bytes(_canonical_bytes(_index_payload(index)))
        if load_manifest_index(index_path) != index:
            raise ArithmeticError("manifest index codec round trip failed")
        verify_manifest_index(index, temporary)
        temporary.replace(output_directory)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output_directory / "index.json"
