#!/usr/bin/env python3
"""Fail-closed audit of Paper A claims and their frozen evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path, PurePosixPath
from typing import Any

CERTIFICATE = "certificates/issue128-d5-integrated-certificate.json"
REFERENCE_CHECKER = "scripts/reference_verify.py"
TRANSCRIPT = "artifacts/d5-integrated/verification-transcript.txt"
MUTATION_CHECKER = "tests/test_certificate_mutations.py"
VERIFICATION_SECTION = "docs/manuscript/sections/verification.tex"
LIMITATIONS_SECTION = "docs/manuscript/sections/limitations.tex"
CLAIM_MANIFEST = "artifacts/publication/paper-a-claim-SHA256SUMS"

CERTIFICATE_CLAIMS = {
    "benchmark_model",
    "benchmark_normalization",
    "lattice_size",
    "simulation_time",
    "tolerance",
    "accepted_steps",
    "adjacent_rejected_steps",
    "published_steps",
    "published_groups",
    "candidate_groups",
    "published_ratio",
    "d4_term_count",
    "d4_group_count",
    "d5_term_count",
    "d5_group_count",
}
REQUIRED_CLAIM_NAMES = CERTIFICATE_CLAIMS | {
    "test_evidence",
    "trusted_computing_base",
}
REQUIRED_RECORD_FIELDS = {
    "value",
    "source",
    "checker",
    "manuscript_locations",
}
TOP_LEVEL_FIELDS = {
    "schema_version",
    "sha256_manifest",
    "claims",
    "excluded_pending_claims",
}
FORBIDDEN_PHRASES = (
    "physically minimal",
    "universal 4.1357",
    "certified fivefold",
)
FIVEFOLD_OPEN_MARKER = re.compile(
    r"\\subsection\{The fivefold boundary is open\}.*?"
    r"This is arithmetic, not a 78-step error certificate\..*?"
    r"\$r=78\$ is therefore not certified\.",
    re.DOTALL,
)
POSITIVE_FIVEFOLD_PATTERNS = (
    re.compile(
        r"\bfivefold\b[^.\n]{0,100}\b(?:is\s+|has\s+been\s+)?"
        r"(?<!not\s)(?:achieved|certified|met|attained|closed)\b",
    ),
    re.compile(
        r"\b(?<!not\s)(?:achieved|certified|met|attained)\b[^.\n]{0,100}"
        r"\bfivefold\b",
    ),
    re.compile(
        r"\br\s*=\s*\$?78\$?[^.\n]{0,40}\b(?:is|has\s+been)\s+"
        r"(?<!not\s)(?:certified|accepted)\b",
    ),
)

CLAIM_BINDINGS: dict[str, tuple[str, re.Pattern[str]]] = {
    "benchmark_model": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(r"periodic \$12\\times12\$ spin-\$1/2\$ isotropic Heisenberg model"),
    ),
    "benchmark_normalization": (
        "docs/manuscript/sections/problem.tex",
        re.compile(r"h_\{uv\}=\\frac\{X_uX_v\+Y_uY_v\+Z_uZ_v\}\{4\}"),
    ),
    "lattice_size": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(r"periodic \$12\\times12\$"),
    ),
    "simulation_time": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(r"\$T=1\$"),
    ),
    "tolerance": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(r"tolerance \$10\^\{-6\}\$"),
    ),
    "accepted_steps": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(r"accepts (?<!\d)95(?!\d) steps"),
    ),
    "adjacent_rejected_steps": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(r"rejects (?<!\d)94(?!\d)(?: steps)?"),
    ),
    "published_steps": (
        "docs/manuscript/sections/introduction.tex",
        re.compile(r"pinned (?<!\d)393(?!\d)-step instantiation"),
    ),
    "published_groups": (
        "docs/manuscript/sections/results.tex",
        re.compile(
            r"Pinned published-theorem control\s*&\s*(?<!\d)393(?!\d)\s*&\s*"
            r"(?<!\d)11,791(?!\d)"
        ),
    ),
    "candidate_groups": (
        "docs/manuscript/sections/results.tex",
        re.compile(
            r"Proof-carrying model-aware bound\s*&\s*(?<!\d)95(?!\d)\s*&\s*"
            r"(?<!\d)2,851(?!\d)"
        ),
    ),
    "published_ratio": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(
            r"(?<!\d)11791/2851=4\.135741844966678(?!\d)"
        ),
    ),
    "d4_term_count": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(
            r"D4 and D5\s+sidecars contain respectively "
            r"(?<!\d)75,324(?!\d) and (?<!\d)605,832(?!\d)",
            re.DOTALL,
        ),
    ),
    "d4_group_count": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(
            r"partitioned into (?<!\d)7,576(?!\d) and (?<!\d)123,106(?!\d)"
        ),
    ),
    "d5_term_count": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(
            r"D4 and D5\s+sidecars contain respectively "
            r"(?<!\d)75,324(?!\d) and (?<!\d)605,832(?!\d)",
            re.DOTALL,
        ),
    ),
    "d5_group_count": (
        "docs/manuscript/sections/abstract.tex",
        re.compile(
            r"partitioned into (?<!\d)7,576(?!\d) and (?<!\d)123,106(?!\d)"
        ),
    ),
    "test_evidence": (
        "docs/manuscript/sections/reproducibility.tex",
        re.compile(r"Verifier falsifiability.*?field-mutation tests.*?must fail", re.DOTALL),
    ),
    "trusted_computing_base": (
        VERIFICATION_SECTION,
        re.compile(r"untrusted producers.*?Python standard library", re.DOTALL),
    ),
}

SOURCE_CHECKER_ALLOWLIST = {
    **{
        name: (CERTIFICATE, REFERENCE_CHECKER)
        for name in CERTIFICATE_CLAIMS
    },
    "test_evidence": (TRANSCRIPT, MUTATION_CHECKER),
    "trusted_computing_base": (VERIFICATION_SECTION, REFERENCE_CHECKER),
}
CLAIM_MANIFEST_PATHS = {
    CERTIFICATE,
    REFERENCE_CHECKER,
    TRANSCRIPT,
    MUTATION_CHECKER,
    LIMITATIONS_SECTION,
    *(location for location, _marker in CLAIM_BINDINGS.values()),
}


def _load_json(path: Path) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key {key!r}")
            result[key] = value
        return result

    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=reject_duplicates)


def _canonical_relative(raw_path: object, field: str) -> tuple[str | None, str | None]:
    if not isinstance(raw_path, str) or not raw_path or "\\" in raw_path:
        return None, f"{field} must be a canonical repository-relative path"
    pure = PurePosixPath(raw_path)
    if pure.is_absolute() or any(part in {"", ".", ".."} for part in pure.parts):
        return None, f"{field} must be a canonical repository-relative path"
    if pure.as_posix() != raw_path:
        return None, f"{field} must be a canonical repository-relative path"
    return raw_path, None


def _path_in_root(root: Path, raw_path: object, field: str) -> tuple[Path | None, str | None]:
    relative, error = _canonical_relative(raw_path, field)
    if error:
        return None, error
    assert relative is not None
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return None, f"{field} escapes the Paper A root"
    if path != root / relative:
        return None, f"{field} resolves through a noncanonical path"
    return path, None


def _canonical_matrix_path(matrix_path: Path, root: Path) -> tuple[Path | None, str | None]:
    candidate = matrix_path if matrix_path.is_absolute() else root / matrix_path
    if any(part in {".", ".."} for part in candidate.parts):
        return None, "matrix_path must be canonical"
    resolved = candidate.resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        return None, "matrix_path must be inside the Paper A root"
    if resolved != candidate:
        return None, "matrix_path must be canonical"
    return resolved, None


def _load_manifest(path: Path, root: Path) -> tuple[dict[str, str], list[str]]:
    errors: list[str] = []
    entries: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        return {}, [f"cannot read SHA-256 manifest {path}: {error}"]
    for line_number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-fA-F]{64}", parts[0]):
            errors.append(f"malformed SHA-256 manifest line {line_number}")
            continue
        relative, path_error = _canonical_relative(parts[1].strip(), "manifest path")
        if path_error:
            errors.append(f"line {line_number}: {path_error}")
            continue
        assert relative is not None
        if relative in entries:
            errors.append(f"duplicate SHA-256 manifest path {relative}")
            continue
        entries[relative] = parts[0].lower()
    return entries, errors


def _certificate_expectations(certificate: dict[str, Any]) -> dict[str, object]:
    benchmark = certificate.get("benchmark", {})
    baseline = certificate.get("published_baseline", {})
    candidate = certificate.get("candidate", {})
    resources = certificate.get("claimed_resources", {})
    d4 = candidate.get("d4_certificate", {})
    d5 = candidate.get("d5_certificate", {})
    length = benchmark.get("length")
    steps = candidate.get("steps")
    return {
        "benchmark_model": benchmark.get("model"),
        "benchmark_normalization": benchmark.get("normalization"),
        "lattice_size": [length, length],
        "simulation_time": benchmark.get("time"),
        "tolerance": benchmark.get("tolerance"),
        "accepted_steps": steps,
        "adjacent_rejected_steps": steps - 1 if isinstance(steps, int) else None,
        "published_steps": baseline.get("steps"),
        "published_groups": resources.get("published_group_exponentials"),
        "candidate_groups": resources.get("candidate_group_exponentials"),
        "published_ratio": certificate.get("claims", {}).get("exact_improvement_ratio"),
        "d4_term_count": d4.get("term_count"),
        "d4_group_count": d4.get("group_count"),
        "d5_term_count": d5.get("term_count"),
        "d5_group_count": d5.get("group_count"),
    }


def _parse_test_evidence(text: str) -> dict[str, int] | None:
    normal = re.search(
        r"NORMAL TEST SUITE.*?Observed:\s*(\d+) passed, (\d+) deselected",
        text,
        re.DOTALL,
    )
    mutation = re.search(
        r"FOCUSED ADVERSARIAL SUITES.*?Observed:\s*(\d+) passed",
        text,
        re.DOTALL,
    )
    if normal is None or mutation is None:
        return None
    return {
        "normal_passed": int(normal.group(1)),
        "normal_deselected": int(normal.group(2)),
        "mutation_passed": int(mutation.group(1)),
    }


def _validate_excluded_pending(value: object) -> list[str]:
    if not isinstance(value, list):
        return ["excluded_pending_claims must be a list"]
    errors: list[str] = []
    names: set[str] = set()
    for position, item in enumerate(value):
        if not isinstance(item, dict) or set(item) != {"name", "reason"}:
            errors.append(
                f"excluded_pending_claims[{position}] must contain exactly name and reason"
            )
            continue
        name, reason = item["name"], item["reason"]
        if not isinstance(name, str) or not name or name in names:
            errors.append(f"excluded_pending_claims[{position}] has an invalid name")
        else:
            names.add(name)
        if not isinstance(reason, str) or not reason.strip():
            errors.append(f"excluded_pending_claims[{position}] has an invalid reason")
    return errors


def audit_claims(matrix_path: Path, root: Path) -> list[str]:
    """Return every publication-audit failure; return ``[]`` only on success."""

    root = root.resolve()
    canonical_matrix, matrix_error = _canonical_matrix_path(matrix_path, root)
    if matrix_error:
        return [matrix_error]
    assert canonical_matrix is not None
    try:
        matrix = _load_json(canonical_matrix)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        return [f"cannot load claim matrix: {error}"]
    if not isinstance(matrix, dict):
        return ["claim matrix must be a JSON object"]

    errors: list[str] = []
    missing_top = TOP_LEVEL_FIELDS - set(matrix)
    unknown_top = set(matrix) - TOP_LEVEL_FIELDS
    if missing_top:
        errors.append(f"claim matrix missing top-level fields {', '.join(sorted(missing_top))}")
    if unknown_top:
        errors.append(f"claim matrix has unknown top-level fields {', '.join(sorted(unknown_top))}")
    if matrix.get("schema_version") != 1:
        errors.append("claim matrix schema_version must equal 1")
    errors.extend(_validate_excluded_pending(matrix.get("excluded_pending_claims")))
    excluded_names = {
        item["name"]
        for item in matrix.get("excluded_pending_claims", [])
        if isinstance(item, dict) and isinstance(item.get("name"), str)
    }
    if "fivefold_global_certificate" not in excluded_names:
        errors.append("fivefold_global_certificate must remain an excluded pending claim")

    claims = matrix.get("claims")
    if not isinstance(claims, dict):
        return errors + ["claim matrix claims must be an object"]
    missing_claims = REQUIRED_CLAIM_NAMES - set(claims)
    unknown_claims = set(claims) - REQUIRED_CLAIM_NAMES
    if missing_claims:
        errors.append(f"claim matrix missing required claims {', '.join(sorted(missing_claims))}")
    if unknown_claims:
        errors.append(f"claim matrix has unknown claims {', '.join(sorted(unknown_claims))}")

    if matrix.get("sha256_manifest") != CLAIM_MANIFEST:
        errors.append(f"sha256_manifest must be {CLAIM_MANIFEST}")
    manifest_path, path_error = _path_in_root(
        root, matrix.get("sha256_manifest"), "sha256_manifest"
    )
    if path_error:
        return errors + [path_error]
    assert manifest_path is not None
    manifest, manifest_errors = _load_manifest(manifest_path, root)
    errors.extend(manifest_errors)

    referenced_paths: set[str] = set()
    manuscript_paths: set[Path] = set()
    for name, record in claims.items():
        if name not in REQUIRED_CLAIM_NAMES:
            continue
        if not isinstance(record, dict):
            errors.append(f"claim {name!r} must be an object")
            continue
        if set(record) != REQUIRED_RECORD_FIELDS:
            missing = REQUIRED_RECORD_FIELDS - set(record)
            unknown = set(record) - REQUIRED_RECORD_FIELDS
            if missing:
                errors.append(f"{name}: missing required fields {', '.join(sorted(missing))}")
            if unknown:
                errors.append(f"{name}: unknown fields {', '.join(sorted(unknown))}")

        expected_source, expected_checker = SOURCE_CHECKER_ALLOWLIST[name]
        if record.get("source") != expected_source:
            errors.append(f"{name}: source must be {expected_source}")
        if record.get("checker") != expected_checker:
            errors.append(f"{name}: checker must be {expected_checker}")
        expected_location, marker = CLAIM_BINDINGS[name]
        if record.get("manuscript_locations") != [expected_location]:
            errors.append(f"{name}: manuscript location must be {expected_location}")

        for field, relative in (
            ("source", record.get("source")),
            ("checker", record.get("checker")),
            ("manuscript location", expected_location),
        ):
            path, current_error = _path_in_root(root, relative, f"{name} {field}")
            if current_error:
                errors.append(current_error)
                continue
            assert path is not None and isinstance(relative, str)
            referenced_paths.add(relative)
            if not path.is_file():
                errors.append(f"{name}: missing {field} file {relative}")
                continue
            if field == "manuscript location":
                manuscript_paths.add(path)
                try:
                    text = path.read_text(encoding="utf-8")
                except (OSError, UnicodeDecodeError) as error:
                    errors.append(f"{name}: cannot read manuscript location: {error}")
                else:
                    if marker.search(text) is None:
                        errors.append(f"{name}: stable manuscript marker is absent")

    for relative in sorted(referenced_paths):
        expected = manifest.get(relative)
        if expected is None:
            errors.append(f"referenced file is not hash-bound: {relative}")
            continue
        path = root / relative
        if path.is_file():
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != expected:
                errors.append(f"SHA-256 mismatch for {relative}")
    if set(manifest) != CLAIM_MANIFEST_PATHS:
        missing_entries = CLAIM_MANIFEST_PATHS - set(manifest)
        extra_entries = set(manifest) - CLAIM_MANIFEST_PATHS
        if missing_entries:
            errors.append(
                "claim manifest missing exact entries "
                + ", ".join(sorted(missing_entries))
            )
        if extra_entries:
            errors.append(
                "claim manifest has non-claim entries "
                + ", ".join(sorted(extra_entries))
            )

    try:
        certificate = _load_json(root / CERTIFICATE)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        errors.append(f"cannot load frozen certificate: {error}")
    else:
        if not isinstance(certificate, dict):
            errors.append("frozen certificate must be an object")
        else:
            for name, expected in _certificate_expectations(certificate).items():
                record = claims.get(name)
                if isinstance(record, dict) and record.get("value") != expected:
                    errors.append(f"{name} does not match the frozen certificate")

    if isinstance(claims.get("published_ratio"), dict) and claims["published_ratio"].get(
        "value"
    ) != [11_791, 2_851]:
        errors.append("published ratio must equal 11791/2851")

    transcript_path = root / TRANSCRIPT
    try:
        transcript_evidence = _parse_test_evidence(
            transcript_path.read_text(encoding="utf-8")
        )
    except (OSError, UnicodeDecodeError) as error:
        errors.append(f"cannot read test transcript: {error}")
    else:
        if transcript_evidence is None:
            errors.append("test transcript lacks exact normal and mutation counts")
        else:
            test_record = claims.get("test_evidence")
            if not isinstance(test_record, dict) or test_record.get("value") != transcript_evidence:
                errors.append("test_evidence does not match the verification transcript")

    tcb_record = claims.get("trusted_computing_base")
    expected_tcb = {
        "discovery_trusted": False,
        "reference_verifier_imports_trottercert": False,
    }
    if not isinstance(tcb_record, dict) or tcb_record.get("value") != expected_tcb:
        errors.append("trusted_computing_base value is invalid")
    try:
        reference_source = (root / REFERENCE_CHECKER).read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        errors.append(f"cannot inspect reference verifier imports: {error}")
    else:
        if re.search(r"^\s*(?:from|import)\s+trottercert\b", reference_source, re.MULTILINE):
            errors.append("reference verifier imports trottercert")

    manuscript_tree = root / "docs/manuscript"
    if manuscript_tree.is_dir():
        manuscript_paths.update(manuscript_tree.rglob("*.tex"))
    limitations_path = root / LIMITATIONS_SECTION
    if not limitations_path.is_file():
        errors.append(f"missing fivefold limitations file {LIMITATIONS_SECTION}")
    else:
        try:
            limitations_text = limitations_path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError) as error:
            errors.append(f"cannot read fivefold limitations file: {error}")
        else:
            if FIVEFOLD_OPEN_MARKER.search(limitations_text) is None:
                errors.append("fivefold open/not-certified manuscript marker is absent")
    for path in manuscript_paths:
        try:
            lowered = path.read_text(encoding="utf-8").lower()
        except (OSError, UnicodeDecodeError):
            continue
        for phrase in FORBIDDEN_PHRASES:
            if phrase in lowered:
                errors.append(f"{path.relative_to(root)}: forbidden phrase {phrase!r}")
        if (
            "strengthened_control_ratio" in excluded_names
            and re.search(r"(?<!\d)10,?591(?!\d)", lowered)
        ):
            errors.append(
                f"{path.relative_to(root)}: excluded strengthened-control numeric token"
            )
        for pattern in POSITIVE_FIVEFOLD_PATTERNS:
            if pattern.search(lowered):
                errors.append(
                    f"{path.relative_to(root)}: forbidden positive fivefold claim"
                )
                break
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    arguments = parser.parse_args()
    errors = audit_claims(arguments.matrix, arguments.root)
    if errors:
        print("\n".join(errors))
        return 1
    print("paper_a_claims=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
