"""Audit Issue 128 ownership and the curated Paper A release boundary."""

from __future__ import annotations

import argparse
import ast
import json
import os
import re
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

ALLOWED_LABELS = {
    "paper-a",
    "paper-b",
    "shared",
    "exploratory",
    "obsolete",
}

_ENTRY_FIELDS = {"ownership", "release_included"}
_IGNORED_FALLBACK_PARTS = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
}
_PAPER_A_RELEASE_PATHS = {
    "CITATION.cff",
    "README.md",
    "requirements-reproducibility.txt",
    "artifacts/d5-integrated/verification-transcript.txt",
    "artifacts/publication/paper-a-claim-SHA256SUMS",
    "artifacts/publication/paper-a-claim-matrix.json",
    "artifacts/publication/paper-a-file-ownership.json",
    "certificates/issue128-certificate.json",
    "certificates/issue128-d4-groups.json",
    "certificates/issue128-d5-groups.json.gz",
    "certificates/issue128-d5-integrated-certificate.json",
    "scripts/audit_paper_a_claims.py",
    "scripts/audit_publication_scope.py",
    "scripts/reference_verify.py",
    "scripts/verify.py",
    "tests/test_certificate_mutations.py",
    "tests/test_paper_a_claims.py",
    "tests/test_publication_scope.py",
    "tests/test_reference_verify.py",
}
_PAPER_A_RELEASE_PREFIX = "docs/manuscript/"
_SELECTED_SHARED_RELEASE_PATHS = {
    "pyproject.toml",
    "src/trottercert/__init__.py",
    "src/trottercert/algebra.py",
    "src/trottercert/anticommuting.py",
    "src/trottercert/baseline.py",
    "src/trottercert/clusters.py",
    "src/trottercert/crosscheck.py",
    "src/trottercert/cubic_field.py",
    "src/trottercert/cubic_local.py",
    "src/trottercert/defect_series.py",
    "src/trottercert/d6_physical_channels.py",
    "src/trottercert/exact_series_certificate.py",
    "src/trottercert/formulas.py",
    "src/trottercert/hamiltonian.py",
    "src/trottercert/higher_order.py",
    "src/trottercert/hpc_artifacts.py",
    "src/trottercert/intervals.py",
    "src/trottercert/lattice.py",
    "src/trottercert/lie_series.py",
    "src/trottercert/local_commutators.py",
    "src/trottercert/orbits.py",
    "src/trottercert/pareto.py",
    "src/trottercert/refined_error.py",
    "src/trottercert/representation.py",
    "src/trottercert/resources.py",
    "src/trottercert/rigorous_fourth.py",
    "src/trottercert/series.py",
    "src/trottercert/su2clusters.py",
    "src/trottercert/support_groups.py",
    "src/trottercert/verify.py",
}
_D6_RUNTIME_EXCEPTION = "src/trottercert/d6_physical_channels.py"
_FIVEFOLD_LIMITATIONS_EXCEPTION = "docs/manuscript/figures/fivefold_gap.pdf"
_FORBIDDEN_RELEASE_TOKENS = (
    "d8",
    "e9",
    "processed",
    "processor-obstruction",
    "processor_obstruction",
)
_REQUIRED_OWNERSHIP = {
    "artifacts/d5-integrated/SHA256SUMS": "paper-a",
    "artifacts/publication/paper-a-claim-SHA256SUMS": "paper-a",
    "certificates/issue128-certificate.json": "paper-a",
    _D6_RUNTIME_EXCEPTION: "shared",
    _FIVEFOLD_LIMITATIONS_EXCEPTION: "paper-a",
    "scripts/package_delivery.py": "obsolete",
    "tests/test_delivery_package.py": "obsolete",
}
_CLAIM_SHA_MANIFEST = "artifacts/publication/paper-a-claim-SHA256SUMS"


@dataclass(frozen=True)
class PublicationEntry:
    ownership: str
    release_included: bool


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"{key}: duplicate JSON key")
        result[key] = value
    return result


def _validate_relative_path(name: object) -> str:
    if not isinstance(name, str) or not name:
        raise ValueError("ownership path must be a non-empty string")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or name != path.as_posix():
        raise ValueError(f"{name}: ownership path must be canonical and relative")
    return name


def load_publication_manifest(path: Path) -> dict[str, PublicationEntry]:
    """Load the strict ownership and release-inclusion manifest."""

    raw = json.loads(path.read_text(), object_pairs_hook=_unique_object)
    if not isinstance(raw, dict):
        raise ValueError("ownership manifest must be a JSON object")  # noqa: TRY004

    result: dict[str, PublicationEntry] = {}
    for raw_name, raw_entry in raw.items():
        name = _validate_relative_path(raw_name)
        if not isinstance(raw_entry, dict) or set(raw_entry) != _ENTRY_FIELDS:
            raise ValueError(
                f"{name}: entry must contain exactly ownership and release_included"
            )
        ownership = raw_entry["ownership"]
        release_included = raw_entry["release_included"]
        if not isinstance(ownership, str):
            raise ValueError(  # noqa: TRY004
                f"{name}: exactly one ownership label is required"
            )
        if ownership not in ALLOWED_LABELS:
            raise ValueError(f"{name}: unknown ownership label {ownership!r}")
        if type(release_included) is not bool:
            raise ValueError(
                f"{name}: release_included must be a boolean"
            )
        result[name] = PublicationEntry(ownership, release_included)
    return result


def load_ownership(path: Path) -> dict[str, str]:
    """Return the path-to-label view retained for Task 1 callers."""

    return {
        name: entry.ownership
        for name, entry in load_publication_manifest(path).items()
    }


def _git_repository_root(root: Path) -> Path | None:
    completed = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
        check=False,
        capture_output=True,
    )
    if completed.returncode:
        return None
    return Path(os.fsdecode(completed.stdout.rstrip(b"\0\r\n"))).resolve()


def _git_paths(root: Path, repository: Path) -> set[str]:
    relative_root = root.resolve().relative_to(repository).as_posix()
    pathspec = "." if relative_root == "." else relative_root
    completed = subprocess.run(
        [
            "git",
            "-C",
            str(repository),
            "ls-files",
            "-z",
            "--cached",
            "--others",
            "--exclude-standard",
            "--deleted",
            "--",
            pathspec,
        ],
        check=True,
        capture_output=True,
    )
    names = [os.fsdecode(item) for item in completed.stdout.split(b"\0") if item]
    if relative_root == ".":
        return set(names)
    prefix = relative_root + "/"
    return {name[len(prefix) :] for name in names if name.startswith(prefix)}


def _fallback_paths(root: Path) -> set[str]:
    return {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if (path.is_file() or path.is_symlink())
        and not any(part in _IGNORED_FALLBACK_PARTS for part in path.parts)
        and path.suffix != ".pyc"
    }


def repository_paths(root: Path) -> set[str]:
    """Return tracked/deleted and non-ignored untracked files below ``root``."""

    resolved = root.resolve()
    repository = _git_repository_root(resolved)
    if repository is not None and resolved.is_relative_to(repository):
        return _git_paths(resolved, repository)
    return _fallback_paths(resolved)


def expected_release_inclusion(name: str, ownership: str) -> bool:
    """Return the complete, path-based Paper A release policy decision."""

    lowered = name.lower()
    if "d6" in lowered and name != _D6_RUNTIME_EXCEPTION:
        return False
    if "fivefold" in lowered and name != _FIVEFOLD_LIMITATIONS_EXCEPTION:
        return False
    if any(token in lowered for token in _FORBIDDEN_RELEASE_TOKENS):
        return False
    if name.endswith("/.gitignore"):
        return False
    if ownership == "paper-a":
        return name in _PAPER_A_RELEASE_PATHS or name.startswith(
            _PAPER_A_RELEASE_PREFIX
        )
    if ownership == "shared":
        return name in _SELECTED_SHARED_RELEASE_PATHS
    return False


def audit_release_policy(entries: dict[str, PublicationEntry]) -> list[str]:
    """Validate every release bit against the curated Paper A policy."""

    errors: list[str] = []
    for name, expected_ownership in sorted(_REQUIRED_OWNERSHIP.items()):
        entry = entries.get(name)
        if entry is not None and entry.ownership != expected_ownership:
            errors.append(
                f"ownership policy mismatch for {name}: expected {expected_ownership}"
            )
    for name, entry in sorted(entries.items()):
        expected = expected_release_inclusion(name, entry.ownership)
        if entry.release_included != expected:
            errors.append(
                f"release_included policy mismatch for {name}: "
                f"expected {str(expected).lower()}"
            )
    return errors


def _python_dependencies(name: str, root: Path) -> set[str]:
    path = root / name
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, SyntaxError):
        return set()
    dependencies: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if node.level == 1 and name.startswith("src/trottercert/"):
                candidate = f"src/trottercert/{module.split('.')[0]}.py"
                if module and (root / candidate).is_file():
                    dependencies.add(candidate)
                continue
            if module == "trottercert":
                dependencies.add("src/trottercert/__init__.py")
                for alias in node.names:
                    candidate = f"src/trottercert/{alias.name}.py"
                    if (root / candidate).is_file():
                        dependencies.add(candidate)
            elif module.startswith("trottercert."):
                leaf = module.split(".", 1)[1].replace(".", "/")
                dependencies.add(f"src/trottercert/{leaf}.py")
            elif module.startswith("scripts."):
                leaf = module.split(".", 1)[1].replace(".", "/")
                dependencies.add(f"scripts/{leaf}.py")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "trottercert":
                    dependencies.add("src/trottercert/__init__.py")
                elif alias.name.startswith("trottercert."):
                    leaf = alias.name.split(".", 1)[1].replace(".", "/")
                    dependencies.add(f"src/trottercert/{leaf}.py")
                elif alias.name.startswith("scripts."):
                    leaf = alias.name.split(".", 1)[1].replace(".", "/")
                    dependencies.add(f"scripts/{leaf}.py")
    return dependencies


def _tex_dependencies(name: str, root: Path) -> set[str]:
    text = (root / name).read_text(encoding="utf-8")
    manuscript = PurePosixPath("docs/manuscript")
    dependencies: set[str] = set()
    for raw in re.findall(r"\\(?:input|include)\{([^}]+)\}", text):
        target = manuscript / raw
        if not target.suffix:
            target = target.with_suffix(".tex")
        dependencies.add(target.as_posix())
    for block in re.findall(r"\\bibliography\{([^}]+)\}", text):
        for raw in block.split(","):
            dependencies.add((manuscript / raw.strip()).with_suffix(".bib").as_posix())
    for raw in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", text):
        target = manuscript / raw
        if target.suffix:
            dependencies.add(target.as_posix())
            continue
        candidates = [target.with_suffix(suffix) for suffix in (".pdf", ".png", ".jpg")]
        existing = next((item for item in candidates if (root / item).is_file()), candidates[0])
        dependencies.add(existing.as_posix())
    return dependencies


def _claim_sha_paths(root: Path) -> tuple[set[str], list[str]]:
    paths: set[str] = set()
    errors: list[str] = []
    try:
        lines = (root / _CLAIM_SHA_MANIFEST).read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        return set(), [f"cannot read claim SHA manifest: {exc}"]
    for line_number, line in enumerate(lines, 1):
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-f]{64}", parts[0]):
            errors.append(f"malformed claim SHA entry on line {line_number}")
            continue
        try:
            paths.add(_validate_relative_path(parts[1].strip()))
        except ValueError as exc:
            errors.append(f"claim SHA line {line_number}: {exc}")
    return paths, errors


def audit_release_closure(
    entries: dict[str, PublicationEntry],
    root: Path,
) -> list[str]:
    """Check static Python, TeX, and claim-SHA closure of the isolated release."""

    included = {name for name, entry in entries.items() if entry.release_included}
    errors: list[str] = []
    for name in sorted(included):
        dependencies: set[str] = set()
        if name.endswith(".py"):
            dependencies.update(_python_dependencies(name, root))
        if name.endswith(".tex"):
            dependencies.update(_tex_dependencies(name, root))
        for dependency in sorted(dependencies - included):
            if (root / dependency).is_file():
                errors.append(f"release dependency missing: {name} -> {dependency}")
    claim_paths, claim_errors = _claim_sha_paths(root)
    errors.extend(claim_errors)
    for name in sorted(claim_paths - included):
        errors.append(f"claim SHA path is not release included: {name}")
    return errors


def audit_ownership(manifest_path: Path, root: Path) -> list[str]:
    """Return fail-closed schema, coverage, filesystem, and release errors."""

    entries = load_publication_manifest(manifest_path)
    actual = repository_paths(root)
    listed = set(entries)
    errors = [f"unclassified repository path: {name}" for name in sorted(actual - listed)]
    errors.extend(
        f"manifest path does not exist: {name}" for name in sorted(listed - actual)
    )
    for name in sorted(actual & listed):
        target = root / PurePosixPath(name)
        if target.is_symlink():
            errors.append(f"repository path must not be a symlink: {name}")
        elif not os.path.lexists(target):
            errors.append(f"repository path is missing from worktree: {name}")
    errors.extend(audit_release_policy(entries))
    errors.extend(audit_release_closure(entries, root))
    return errors


def release_paths(manifest_path: Path) -> set[str]:
    """Return the explicitly reviewed Paper A release file set."""

    return {
        name
        for name, entry in load_publication_manifest(manifest_path).items()
        if entry.release_included
    }


def _format_errors(errors: Iterable[str]) -> str:
    return "\n".join(errors)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Issue 128 root (defaults to the parent of scripts/)",
    )
    args = parser.parse_args(argv)

    try:
        errors = audit_ownership(args.manifest, args.root)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        print(f"publication ownership audit failed: {exc}")
        return 1
    if errors:
        print(_format_errors(errors))
        return 1
    entries = load_publication_manifest(args.manifest)
    included = sum(entry.release_included for entry in entries.values())
    print(
        "publication ownership audit passed: "
        f"{len(entries)} files, {included} Paper A release files"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
