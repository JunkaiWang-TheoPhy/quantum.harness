"""Canonical per-Delta artifacts for the compressed grouped-XXZ compiler."""

# Malformed decoded JSON is a semantic artifact error even when the malformed
# field has the wrong Python type.
# ruff: noqa: TRY004

from __future__ import annotations

import ast
import gzip
import json
import os
import subprocess
import tempfile
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from pathlib import Path

from .grouped_xxz import (
    AnticommutingGroupRecord,
    SymplecticCoefficient,
    XXZCompileSpec,
    canonical_json_bytes,
    fraction_pair,
    strict_fraction_pair,
)
from .grouped_xxz_compressed import (
    ActualWordWeight,
    CompressedXXZCertificate,
    CompressedXXZLedger,
    RepresentativeNormBlock,
    verify_compressed_xxz_certificate,
)
from .intervals import RationalInterval

SUMMARY_SCHEMA = "grouped_xxz_per_delta_summary_v1"
WITNESS_SCHEMA = "grouped_xxz_per_delta_witness_v1"
SOURCE_ENTRYPOINTS = (
    "scripts/compile_grouped_xxz.py",
    "src/trottercert/grouped_xxz_artifact.py",
)
FROZEN_LOCAL_SOURCE_ALLOWLIST = (
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
ENVIRONMENT_CLOSURE_PATHS = (
    "pyproject.toml",
    "requirements-reproducibility.txt",
)
SOURCE_CLOSURE_PATHS = tuple(
    sorted((*FROZEN_LOCAL_SOURCE_ALLOWLIST, *ENVIRONMENT_CLOSURE_PATHS))
)


@dataclass(frozen=True, slots=True)
class PerDeltaArtifact:
    summary: dict[str, object]
    witness: dict[str, object]
    summary_bytes: bytes
    witness_payload_bytes: bytes
    witness_gzip_bytes: bytes


def payload_digest(payload: dict[str, object]) -> str:
    """Return the canonical digest with the payload digest field blanked."""

    if not isinstance(payload, dict):
        raise TypeError("payload must be a dictionary")
    unsigned = dict(payload)
    unsigned["payload_sha256"] = ""
    return sha256(canonical_json_bytes(unsigned)).hexdigest()


def _seal(payload: dict[str, object]) -> dict[str, object]:
    sealed = dict(payload)
    sealed["payload_sha256"] = ""
    sealed["payload_sha256"] = payload_digest(sealed)
    return sealed


def deterministic_gzip_bytes(payload: bytes) -> bytes:
    """Compress one payload with an empty filename and zero modification time."""

    if not isinstance(payload, bytes):
        raise TypeError("gzip payload must be bytes")
    output = BytesIO()
    with gzip.GzipFile(
        filename="",
        mode="wb",
        compresslevel=9,
        fileobj=output,
        mtime=0,
    ) as member:
        member.write(payload)
    return output.getvalue()


def _module_paths(root: Path, module: str) -> tuple[str, ...]:
    """Resolve one repo-local module and all executed parent package inits."""

    if module != "trottercert" and not module.startswith("trottercert."):
        return ()
    parts = module.split(".")
    resolved: list[str] = []
    for length in range(1, len(parts)):
        package = Path("src", *parts[:length], "__init__.py").as_posix()
        if (root / package).is_file():
            resolved.append(package)
    module_file = Path("src", *parts).with_suffix(".py").as_posix()
    package_init = Path("src", *parts, "__init__.py").as_posix()
    if (root / module_file).is_file():
        resolved.append(module_file)
    elif (root / package_init).is_file():
        resolved.append(package_init)
    return tuple(resolved)


def _module_context(relative: str) -> tuple[str | None, str | None]:
    path = Path(relative)
    if not relative.startswith("src/trottercert/"):
        return None, None
    parts = list(path.with_suffix("").parts[1:])
    if parts[-1] == "__init__":
        parts.pop()
        module = ".".join(parts)
        return module, module
    module = ".".join(parts)
    return module, ".".join(parts[:-1])


def _absolute_import_module(
    node: ast.ImportFrom,
    package: str | None,
) -> str | None:
    if node.level == 0:
        return node.module
    if package is None:
        return None
    package_parts = package.split(".")
    remove = node.level - 1
    if remove > len(package_parts):
        return None
    prefix = package_parts[: len(package_parts) - remove]
    if node.module:
        prefix.extend(node.module.split("."))
    return ".".join(prefix)


def derive_local_source_closure(root: str | Path) -> tuple[str, ...]:
    """Recursively derive every executed repo-local Python source by AST imports."""

    root_path = Path(root).resolve()
    pending = list(SOURCE_ENTRYPOINTS)
    visited: set[str] = set()
    while pending:
        relative = pending.pop()
        if relative in visited:
            continue
        path = root_path / relative
        if not path.is_file():
            raise FileNotFoundError(f"source closure file is missing: {relative}")
        visited.add(relative)
        try:
            tree = ast.parse(path.read_bytes(), filename=relative)
        except SyntaxError as error:
            raise ValueError(
                f"source closure file is not valid Python: {relative}"
            ) from error
        _, package = _module_context(relative)
        discovered: set[str] = set()
        for node in ast.walk(tree):
            modules: list[str] = []
            if isinstance(node, ast.Import):
                modules.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                base = _absolute_import_module(node, package)
                if base:
                    modules.append(base)
                    # ``from package import submodule`` executes both the
                    # parent package and the imported child when it exists.
                    modules.extend(
                        f"{base}.{alias.name}"
                        for alias in node.names
                        if alias.name != "*"
                    )
            for module in modules:
                discovered.update(_module_paths(root_path, module))
        pending.extend(sorted(discovered - visited, reverse=True))
    return tuple(sorted(visited))


def source_closure(
    root: str | Path,
    *,
    source_commit: str | None = None,
) -> dict[str, object]:
    """Hash the exact source files needed to generate and verify this artifact."""

    root_path = Path(root).resolve()
    if source_commit is None:
        completed = subprocess.run(
            ["git", "-C", str(root_path), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        source_commit = completed.stdout.strip()
    if not isinstance(source_commit, str) or not source_commit:
        raise ValueError("source commit must be a nonempty string")
    derived = derive_local_source_closure(root_path)
    if derived != FROZEN_LOCAL_SOURCE_ALLOWLIST:
        missing = sorted(set(FROZEN_LOCAL_SOURCE_ALLOWLIST) - set(derived))
        unexpected = sorted(set(derived) - set(FROZEN_LOCAL_SOURCE_ALLOWLIST))
        raise ValueError(
            "recursive source closure differs from frozen allowlist: "
            f"missing={missing}, unexpected={unexpected}"
        )
    files = []
    for relative in SOURCE_CLOSURE_PATHS:
        path = root_path / relative
        if not path.is_file():
            raise FileNotFoundError(f"source closure file is missing: {relative}")
        files.append(
            {
                "path": relative,
                "sha256": sha256(path.read_bytes()).hexdigest(),
            }
        )
    unsigned = {
        "algorithm": "sha256",
        "files": files,
        "source_commit": source_commit,
    }
    return {
        **unsigned,
        "closure_sha256": sha256(canonical_json_bytes(unsigned)).hexdigest(),
    }


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


def _optional_fraction(value: Fraction | None) -> list[int] | None:
    return None if value is None else fraction_pair(value)


def _word_weight_payload(row: ActualWordWeight) -> dict[str, object]:
    return {
        "actual_word": list(row.actual_word),
        "raw_weight_upper": fraction_pair(row.raw_weight_upper),
        "representative": list(row.representative),
        "symmetry_name": row.symmetry_name,
    }


def _representative_payload(block: RepresentativeNormBlock) -> dict[str, object]:
    term_indices = {term.mask: index for index, term in enumerate(block.terms)}
    if len(term_indices) != len(block.terms):
        raise ValueError("representative terms contain duplicate masks")
    return {
        "grouped_norm": fraction_pair(block.grouped_norm),
        "groups": [
            {
                "norm_interval": [
                    fraction_pair(group.norm_interval.lower),
                    fraction_pair(group.norm_interval.upper),
                ],
                "squared_norm": fraction_pair(group.squared_norm),
                "term_indices": [term_indices[term.mask] for term in group.terms],
            }
            for group in block.groups
        ],
        "representative": list(block.representative),
        "terms": [
            [
                term.x_mask,
                term.z_mask,
                fraction_pair(term.real),
                fraction_pair(term.imag),
            ]
            for term in block.terms
        ],
        "triangle_norm": fraction_pair(block.triangle_norm),
    }


def _ledger_payload(ledger: CompressedXXZLedger) -> dict[str, object]:
    return {
        "center": ledger.center,
        "complete": ledger.complete,
        "duhamel_convention": ledger.duhamel_convention,
        "factorial_denominator": ledger.factorial_denominator,
        "finite_step_error_formula": ledger.finite_step_error_formula,
        "grouped_constant": fraction_pair(ledger.grouped_constant),
        "grouping_algorithm": ledger.grouping_algorithm,
        "ledger_digest": ledger.ledger_digest,
        "max_degree": ledger.max_degree,
        "order": ledger.order,
        "projected_record_count": ledger.projected_record_count,
        "raw_record_count": ledger.raw_record_count,
        "raw_stream_sha256": ledger.raw_stream_sha256,
        "record_limit": ledger.record_limit,
        "spec": _spec_payload(ledger.spec),
        "theorem_identifier": ledger.theorem_identifier,
        "triangle_constant": fraction_pair(ledger.triangle_constant),
    }


def _certificate_payload(
    certificate: CompressedXXZCertificate,
) -> dict[str, object]:
    return {
        "baseline_error": fraction_pair(certificate.baseline_error),
        "baseline_previous_error": _optional_fraction(
            certificate.baseline_previous_error
        ),
        "baseline_resources": certificate.baseline_resources,
        "baseline_steps": certificate.baseline_steps,
        "bond_growth": fraction_pair(certificate.bond_growth),
        "candidate_error": fraction_pair(certificate.candidate_error),
        "candidate_previous_error": _optional_fraction(
            certificate.candidate_previous_error
        ),
        "candidate_resources": certificate.candidate_resources,
        "candidate_steps": certificate.candidate_steps,
        "cell_base": fraction_pair(certificate.cell_base),
        "center": certificate.center,
        "duhamel_convention": certificate.duhamel_convention,
        "factorial_denominator": certificate.factorial_denominator,
        "finite_step_error_formula": certificate.finite_step_error_formula,
        "grouped_constant": fraction_pair(certificate.grouped_constant),
        "grouping_algorithm": certificate.grouping_algorithm,
        "ledger_digest": certificate.ledger_digest,
        "method": certificate.method,
        "order": certificate.order,
        "raw_record_count": certificate.raw_record_count,
        "raw_stream_sha256": certificate.raw_stream_sha256,
        "spec": _spec_payload(certificate.spec),
        "status": certificate.status,
        "theorem_identifier": certificate.theorem_identifier,
        "triangle_constant": fraction_pair(certificate.triangle_constant),
    }


def _witness_payload(
    ledger: CompressedXXZLedger,
    certificate: CompressedXXZCertificate,
    closure: dict[str, object],
) -> dict[str, object]:
    return _seal(
        {
            "certificate": _certificate_payload(certificate),
            "ledger": _ledger_payload(ledger),
            "payload_sha256": "",
            "representative_blocks": [
                _representative_payload(block) for block in ledger.representative_blocks
            ],
            "schema": WITNESS_SCHEMA,
            "source_closure": closure,
            "word_weights": [_word_weight_payload(row) for row in ledger.word_weights],
        }
    )


def _summary_payload(
    ledger: CompressedXXZLedger,
    certificate: CompressedXXZCertificate,
    closure: dict[str, object],
    witness_payload_bytes: bytes,
    witness_gzip_bytes: bytes,
) -> dict[str, object]:
    return _seal(
        {
            "adjacent_steps": {
                "baseline": {
                    "accepted_error": fraction_pair(certificate.baseline_error),
                    "previous_error": _optional_fraction(
                        certificate.baseline_previous_error
                    ),
                    "steps": certificate.baseline_steps,
                },
                "candidate": {
                    "accepted_error": fraction_pair(certificate.candidate_error),
                    "previous_error": _optional_fraction(
                        certificate.candidate_previous_error
                    ),
                    "steps": certificate.candidate_steps,
                },
            },
            "bounds": {
                "grouped_constant": fraction_pair(certificate.grouped_constant),
                "triangle_constant": fraction_pair(certificate.triangle_constant),
            },
            "counts": {
                "actual_words": len(ledger.word_weights),
                "pair_groups": sum(
                    len(group.terms) == 2
                    for block in ledger.representative_blocks
                    for group in block.groups
                ),
                "projected_raw_records": ledger.projected_record_count,
                "raw_records": ledger.raw_record_count,
                "representative_blocks": len(ledger.representative_blocks),
                "representative_terms": sum(
                    len(block.terms) for block in ledger.representative_blocks
                ),
                "singleton_groups": sum(
                    len(group.terms) == 1
                    for block in ledger.representative_blocks
                    for group in block.groups
                ),
            },
            "method": certificate.method,
            "payload_sha256": "",
            "resources": {
                "baseline": certificate.baseline_resources,
                "candidate": certificate.candidate_resources,
                "metric": certificate.spec.primary_metric,
            },
            "schema": SUMMARY_SCHEMA,
            "source_closure": closure,
            "spec": _spec_payload(certificate.spec),
            "status": certificate.status,
            "theorem": {
                "center": certificate.center,
                "duhamel_convention": certificate.duhamel_convention,
                "factorial_denominator": certificate.factorial_denominator,
                "finite_step_error_formula": certificate.finite_step_error_formula,
                "grouping_algorithm": certificate.grouping_algorithm,
                "identifier": certificate.theorem_identifier,
                "order": certificate.order,
            },
            "witness": {
                "compression": "gzip-9-mtime0-empty-filename",
                "file_sha256": sha256(witness_gzip_bytes).hexdigest(),
                "payload_sha256": sha256(witness_payload_bytes).hexdigest(),
            },
        }
    )


def build_per_delta_artifact(
    ledger: CompressedXXZLedger,
    certificate: CompressedXXZCertificate,
    closure: dict[str, object],
) -> PerDeltaArtifact:
    """Encode one already-certified per-Delta ledger without timing metadata."""

    if not isinstance(ledger, CompressedXXZLedger):
        raise TypeError("ledger must be a CompressedXXZLedger")
    if not isinstance(certificate, CompressedXXZCertificate):
        raise TypeError("certificate must be a CompressedXXZCertificate")
    if (
        certificate.spec != ledger.spec
        or certificate.ledger_digest != ledger.ledger_digest
    ):
        raise ValueError("certificate does not bind the submitted compressed ledger")
    _validate_source_closure(closure)
    witness = _witness_payload(ledger, certificate, closure)
    witness_payload_bytes = canonical_json_bytes(witness)
    witness_gzip_bytes = deterministic_gzip_bytes(witness_payload_bytes)
    summary = _summary_payload(
        ledger,
        certificate,
        closure,
        witness_payload_bytes,
        witness_gzip_bytes,
    )
    return PerDeltaArtifact(
        summary=summary,
        witness=witness,
        summary_bytes=canonical_json_bytes(summary),
        witness_payload_bytes=witness_payload_bytes,
        witness_gzip_bytes=witness_gzip_bytes,
    )


def _validate_source_closure(closure: object) -> dict[str, object]:
    _expect_mapping_keys(
        closure,
        {"algorithm", "closure_sha256", "files", "source_commit"},
        "source_closure",
    )
    assert isinstance(closure, dict)
    if closure["algorithm"] != "sha256":
        raise ValueError("source closure algorithm is incorrect")
    if not isinstance(closure["source_commit"], str) or not closure["source_commit"]:
        raise ValueError("source closure commit is invalid")
    files = closure["files"]
    if not isinstance(files, list) or len(files) != len(SOURCE_CLOSURE_PATHS):
        raise ValueError("source closure file coverage is incorrect")
    expected_paths = []
    for entry in files:
        _expect_mapping_keys(entry, {"path", "sha256"}, "source closure file")
        assert isinstance(entry, dict)
        if not isinstance(entry["path"], str):
            raise ValueError("source closure path is invalid")
        _validate_sha(entry["sha256"], "source closure file digest")
        expected_paths.append(entry["path"])
    if tuple(expected_paths) != SOURCE_CLOSURE_PATHS:
        raise ValueError("source closure paths are not frozen")
    unsigned = {
        "algorithm": closure["algorithm"],
        "files": files,
        "source_commit": closure["source_commit"],
    }
    expected = sha256(canonical_json_bytes(unsigned)).hexdigest()
    if closure["closure_sha256"] != expected:
        raise ValueError("source closure digest is incorrect")
    return closure


def _expect_mapping_keys(value: object, keys: set[str], label: str) -> None:
    if not isinstance(value, dict) or set(value) != keys:
        raise ValueError(f"{label} schema keys are incorrect")


def _validate_sha(value: object, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{label} must be lowercase SHA-256 hex")
    return value


def _canonical_json_load(raw: bytes, label: str) -> dict[str, object]:
    if not isinstance(raw, bytes):
        raise TypeError(f"{label} must be bytes")

    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{label} contains a duplicate key")
            result[key] = value
        return result

    try:
        payload = json.loads(raw, object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid UTF-8 JSON") from error
    if not isinstance(payload, dict) or canonical_json_bytes(payload) != raw:
        raise ValueError(f"{label} is not canonical JSON")
    return payload


def _integer(value: object, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} is not a valid integer")
    return value


def _word(value: object, label: str) -> tuple[int, int, int, int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 5
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
        or any(not 0 <= item < 4 for item in value)
    ):
        raise ValueError(f"{label} is not a five-letter fragment word")
    return tuple(value)  # type: ignore[return-value]


def _fraction(value: object, label: str) -> Fraction:
    return strict_fraction_pair(value, field=label)


def _optional_fraction_decode(value: object, label: str) -> Fraction | None:
    return None if value is None else _fraction(value, label)


def _decode_spec(payload: object) -> XXZCompileSpec:
    keys = {
        "boundary",
        "delta",
        "formula_identifier",
        "length",
        "normalization",
        "primary_metric",
        "stage_count",
        "time",
        "tolerance",
    }
    _expect_mapping_keys(payload, keys, "XXZ spec")
    assert isinstance(payload, dict)
    return XXZCompileSpec(
        delta=_fraction(payload["delta"], "spec.delta"),
        length=_integer(payload["length"], "spec.length", minimum=1),
        boundary=payload["boundary"],  # type: ignore[arg-type]
        time=_fraction(payload["time"], "spec.time"),
        tolerance=_fraction(payload["tolerance"], "spec.tolerance"),
        normalization=payload["normalization"],  # type: ignore[arg-type]
        formula_identifier=payload["formula_identifier"],  # type: ignore[arg-type]
        stage_count=_integer(payload["stage_count"], "spec.stage_count", minimum=1),
        primary_metric=payload["primary_metric"],  # type: ignore[arg-type]
    )


def _decode_word_weights(payload: object) -> tuple[ActualWordWeight, ...]:
    if not isinstance(payload, list):
        raise ValueError("word_weights must be a list")
    rows = []
    for value in payload:
        _expect_mapping_keys(
            value,
            {"actual_word", "raw_weight_upper", "representative", "symmetry_name"},
            "word weight",
        )
        assert isinstance(value, dict)
        if not isinstance(value["symmetry_name"], str):
            raise ValueError("word-weight symmetry name is invalid")
        rows.append(
            ActualWordWeight(
                actual_word=_word(value["actual_word"], "actual_word"),
                representative=_word(value["representative"], "representative"),
                symmetry_name=value["symmetry_name"],
                raw_weight_upper=_fraction(
                    value["raw_weight_upper"],
                    "raw_weight_upper",
                ),
            )
        )
    return tuple(rows)


def _decode_representatives(payload: object) -> tuple[RepresentativeNormBlock, ...]:
    if not isinstance(payload, list):
        raise ValueError("representative_blocks must be a list")
    blocks = []
    for value in payload:
        _expect_mapping_keys(
            value,
            {"grouped_norm", "groups", "representative", "terms", "triangle_norm"},
            "representative block",
        )
        assert isinstance(value, dict)
        raw_terms = value["terms"]
        if not isinstance(raw_terms, list):
            raise ValueError("representative terms must be a list")
        terms = []
        for raw_term in raw_terms:
            if not isinstance(raw_term, list) or len(raw_term) != 4:
                raise ValueError("representative term schema is incorrect")
            terms.append(
                SymplecticCoefficient(
                    _integer(raw_term[0], "term.x_mask"),
                    _integer(raw_term[1], "term.z_mask"),
                    _fraction(raw_term[2], "term.real"),
                    _fraction(raw_term[3], "term.imag"),
                )
            )
        raw_groups = value["groups"]
        if not isinstance(raw_groups, list):
            raise ValueError("representative groups must be a list")
        groups = []
        for raw_group in raw_groups:
            _expect_mapping_keys(
                raw_group,
                {"norm_interval", "squared_norm", "term_indices"},
                "representative group",
            )
            assert isinstance(raw_group, dict)
            indices = raw_group["term_indices"]
            if not isinstance(indices, list) or not indices:
                raise ValueError("group term indices must be a nonempty list")
            decoded_indices = tuple(
                _integer(index, "group term index") for index in indices
            )
            if len(set(decoded_indices)) != len(decoded_indices) or any(
                index >= len(terms) for index in decoded_indices
            ):
                raise ValueError("group term index is invalid or duplicated")
            interval = raw_group["norm_interval"]
            if not isinstance(interval, list) or len(interval) != 2:
                raise ValueError("group norm interval schema is incorrect")
            groups.append(
                AnticommutingGroupRecord(
                    terms=tuple(terms[index] for index in decoded_indices),
                    squared_norm=_fraction(
                        raw_group["squared_norm"],
                        "group.squared_norm",
                    ),
                    norm_interval=RationalInterval(
                        _fraction(interval[0], "group.norm_interval.lower"),
                        _fraction(interval[1], "group.norm_interval.upper"),
                    ),
                )
            )
        blocks.append(
            RepresentativeNormBlock(
                representative=_word(value["representative"], "representative"),
                terms=tuple(terms),
                groups=tuple(groups),
                grouped_norm=_fraction(value["grouped_norm"], "grouped_norm"),
                triangle_norm=_fraction(value["triangle_norm"], "triangle_norm"),
            )
        )
    return tuple(blocks)


def _decode_ledger(
    payload: object,
    word_weights: tuple[ActualWordWeight, ...],
    blocks: tuple[RepresentativeNormBlock, ...],
) -> CompressedXXZLedger:
    keys = {
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
    }
    _expect_mapping_keys(payload, keys, "compressed ledger")
    assert isinstance(payload, dict)
    if not isinstance(payload["complete"], bool):
        raise ValueError("ledger complete flag is invalid")
    record_limit = payload["record_limit"]
    if record_limit is not None:
        record_limit = _integer(record_limit, "ledger.record_limit", minimum=1)
    for field in (
        "duhamel_convention",
        "finite_step_error_formula",
        "grouping_algorithm",
        "theorem_identifier",
    ):
        if not isinstance(payload[field], str):
            raise ValueError(f"ledger {field} is invalid")
    return CompressedXXZLedger(
        spec=_decode_spec(payload["spec"]),
        theorem_identifier=payload["theorem_identifier"],
        order=_integer(payload["order"], "ledger.order", minimum=1),
        center=_integer(payload["center"], "ledger.center", minimum=1),
        factorial_denominator=_integer(
            payload["factorial_denominator"],
            "ledger.factorial_denominator",
            minimum=1,
        ),
        duhamel_convention=payload["duhamel_convention"],
        finite_step_error_formula=payload["finite_step_error_formula"],
        grouping_algorithm=payload["grouping_algorithm"],
        max_degree=_integer(payload["max_degree"], "ledger.max_degree", minimum=1),
        projected_record_count=_integer(
            payload["projected_record_count"],
            "ledger.projected_record_count",
            minimum=1,
        ),
        record_limit=record_limit,
        complete=payload["complete"],
        raw_record_count=_integer(
            payload["raw_record_count"],
            "ledger.raw_record_count",
            minimum=1,
        ),
        raw_stream_sha256=_validate_sha(
            payload["raw_stream_sha256"],
            "ledger raw stream digest",
        ),
        word_weights=word_weights,
        representative_blocks=blocks,
        grouped_constant=_fraction(
            payload["grouped_constant"],
            "ledger.grouped_constant",
        ),
        triangle_constant=_fraction(
            payload["triangle_constant"],
            "ledger.triangle_constant",
        ),
        ledger_digest=_validate_sha(payload["ledger_digest"], "ledger digest"),
    )


def _decode_certificate(payload: object) -> CompressedXXZCertificate:
    keys = {
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
    }
    _expect_mapping_keys(payload, keys, "compressed certificate")
    assert isinstance(payload, dict)
    for field in (
        "duhamel_convention",
        "finite_step_error_formula",
        "grouping_algorithm",
        "method",
        "status",
        "theorem_identifier",
    ):
        if not isinstance(payload[field], str):
            raise ValueError(f"certificate {field} is invalid")
    return CompressedXXZCertificate(
        status=payload["status"],  # type: ignore[arg-type]
        method=payload["method"],
        spec=_decode_spec(payload["spec"]),
        ledger_digest=_validate_sha(
            payload["ledger_digest"], "certificate ledger digest"
        ),
        raw_stream_sha256=_validate_sha(
            payload["raw_stream_sha256"],
            "certificate raw stream digest",
        ),
        raw_record_count=_integer(
            payload["raw_record_count"],
            "certificate.raw_record_count",
            minimum=1,
        ),
        theorem_identifier=payload["theorem_identifier"],
        order=_integer(payload["order"], "certificate.order", minimum=1),
        center=_integer(payload["center"], "certificate.center", minimum=1),
        factorial_denominator=_integer(
            payload["factorial_denominator"],
            "certificate.factorial_denominator",
            minimum=1,
        ),
        duhamel_convention=payload["duhamel_convention"],
        finite_step_error_formula=payload["finite_step_error_formula"],
        grouping_algorithm=payload["grouping_algorithm"],
        grouped_constant=_fraction(
            payload["grouped_constant"],
            "certificate.grouped_constant",
        ),
        triangle_constant=_fraction(
            payload["triangle_constant"],
            "certificate.triangle_constant",
        ),
        candidate_steps=_integer(
            payload["candidate_steps"],
            "certificate.candidate_steps",
            minimum=1,
        ),
        candidate_error=_fraction(
            payload["candidate_error"],
            "certificate.candidate_error",
        ),
        candidate_previous_error=_optional_fraction_decode(
            payload["candidate_previous_error"],
            "certificate.candidate_previous_error",
        ),
        baseline_steps=_integer(
            payload["baseline_steps"],
            "certificate.baseline_steps",
            minimum=1,
        ),
        baseline_error=_fraction(
            payload["baseline_error"],
            "certificate.baseline_error",
        ),
        baseline_previous_error=_optional_fraction_decode(
            payload["baseline_previous_error"],
            "certificate.baseline_previous_error",
        ),
        candidate_resources=_integer(
            payload["candidate_resources"],
            "certificate.candidate_resources",
            minimum=1,
        ),
        baseline_resources=_integer(
            payload["baseline_resources"],
            "certificate.baseline_resources",
            minimum=1,
        ),
        bond_growth=_fraction(payload["bond_growth"], "certificate.bond_growth"),
        cell_base=_fraction(payload["cell_base"], "certificate.cell_base"),
    )


def verify_per_delta_artifact(
    summary_bytes: bytes,
    witness_gzip_bytes: bytes,
    *,
    expected_source_closure: dict[str, object],
) -> tuple[CompressedXXZLedger, CompressedXXZCertificate]:
    """Verify canonical encoding, both witness hashes, and all scientific fields."""

    summary = _canonical_json_load(summary_bytes, "summary")
    _expect_mapping_keys(
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
    if summary["schema"] != SUMMARY_SCHEMA:
        raise ValueError("summary schema identifier is incorrect")
    if summary["payload_sha256"] != payload_digest(summary):
        raise ValueError("summary payload digest is incorrect")
    _expect_mapping_keys(
        summary["witness"],
        {"compression", "file_sha256", "payload_sha256"},
        "summary witness binding",
    )
    binding = summary["witness"]
    assert isinstance(binding, dict)
    if binding["compression"] != "gzip-9-mtime0-empty-filename":
        raise ValueError("witness compression identifier is incorrect")
    if (
        _validate_sha(binding["file_sha256"], "witness file digest")
        != sha256(witness_gzip_bytes).hexdigest()
    ):
        raise ValueError("witness gzip file digest is incorrect")
    try:
        witness_payload_bytes = gzip.decompress(witness_gzip_bytes)
    except (gzip.BadGzipFile, EOFError, OSError) as error:
        raise ValueError("witness is not a valid gzip member") from error
    if deterministic_gzip_bytes(witness_payload_bytes) != witness_gzip_bytes:
        raise ValueError("witness gzip encoding is not deterministic")
    if (
        _validate_sha(binding["payload_sha256"], "witness payload digest")
        != sha256(witness_payload_bytes).hexdigest()
    ):
        raise ValueError("witness uncompressed payload digest is incorrect")

    witness = _canonical_json_load(witness_payload_bytes, "witness payload")
    _expect_mapping_keys(
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
    if witness["schema"] != WITNESS_SCHEMA:
        raise ValueError("witness schema identifier is incorrect")
    if witness["payload_sha256"] != payload_digest(witness):
        raise ValueError("witness payload digest is incorrect")
    summary_closure = _validate_source_closure(summary["source_closure"])
    witness_closure = _validate_source_closure(witness["source_closure"])
    if summary_closure != witness_closure:
        raise ValueError("summary and witness source closures differ")
    _validate_source_closure(expected_source_closure)
    if summary_closure != expected_source_closure:
        raise ValueError("artifact source closure differs from expected source")

    rows = _decode_word_weights(witness["word_weights"])
    blocks = _decode_representatives(witness["representative_blocks"])
    ledger = _decode_ledger(witness["ledger"], rows, blocks)
    certificate = _decode_certificate(witness["certificate"])
    verify_compressed_xxz_certificate(certificate, ledger)

    expected_witness = _witness_payload(ledger, certificate, summary_closure)
    if witness != expected_witness:
        raise ValueError("witness differs from canonical semantic reconstruction")
    expected_summary = _summary_payload(
        ledger,
        certificate,
        summary_closure,
        witness_payload_bytes,
        witness_gzip_bytes,
    )
    if summary != expected_summary:
        raise ValueError("summary differs from the verified witness")
    return ledger, certificate


def _write_temporary(path: Path, raw: bytes) -> Path:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return temporary


def write_per_delta_artifact(
    artifact: PerDeltaArtifact,
    summary_path: str | Path,
    witness_path: str | Path,
) -> None:
    """Publish witness first and the no-overwrite summary commit marker last.

    A process crash between the two atomic links may leave a witness-only
    orphan, but can never leave an apparently committed summary without its
    witness.  An ordinary second-link exception rolls the witness back.
    """

    if not isinstance(artifact, PerDeltaArtifact):
        raise TypeError("artifact must be a PerDeltaArtifact")
    summary = Path(summary_path)
    witness = Path(witness_path)
    if summary.resolve() == witness.resolve():
        raise ValueError("summary and witness paths must be distinct")
    if summary.exists() or witness.exists():
        raise FileExistsError("artifact output already exists")
    summary.parent.mkdir(parents=True, exist_ok=True)
    witness.parent.mkdir(parents=True, exist_ok=True)
    summary_temp = _write_temporary(summary, artifact.summary_bytes)
    witness_temp = _write_temporary(witness, artifact.witness_gzip_bytes)
    witness_published = False
    try:
        os.link(witness_temp, witness)
        witness_published = True
        os.link(summary_temp, summary)
    except BaseException:
        if witness_published:
            witness.unlink(missing_ok=True)
        raise
    finally:
        summary_temp.unlink(missing_ok=True)
        witness_temp.unlink(missing_ok=True)
