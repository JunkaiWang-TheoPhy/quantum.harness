from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.audit_publication_scope import (
    ALLOWED_LABELS,
    audit_ownership,
    audit_release_closure,
    audit_release_policy,
    load_ownership,
    load_publication_manifest,
    release_paths,
    repository_paths,
)

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "artifacts" / "publication" / "paper-a-file-ownership.json"


def _entry(
    ownership: object = "paper-a",
    release_included: object = False,
) -> dict[str, object]:
    return {
        "ownership": ownership,
        "release_included": release_included,
    }


def _write_manifest(path: Path, data: dict[str, object]) -> None:
    path.write_text(json.dumps(data))


def _init_repository(path: Path) -> None:
    path.mkdir()
    subprocess.run(["git", "init", "-q", str(path)], check=True)


def test_ownership_rejects_unknown_label(tmp_path: Path) -> None:
    path = tmp_path / "ownership.json"
    _write_manifest(path, {"README.md": _entry("misc")})
    with pytest.raises(ValueError, match="unknown ownership label"):
        load_ownership(path)


def test_ownership_rejects_duplicate_release_roles(tmp_path: Path) -> None:
    path = tmp_path / "ownership.json"
    _write_manifest(path, {"README.md": _entry(["paper-a", "obsolete"])})
    with pytest.raises(ValueError, match="exactly one ownership label"):
        load_ownership(path)


def test_ownership_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    path = tmp_path / "ownership.json"
    path.write_text(
        '{"README.md":{"ownership":"paper-a",'
        '"release_included":false,"release_included":true}}'
    )
    with pytest.raises(ValueError, match="duplicate JSON key"):
        load_ownership(path)


@pytest.mark.parametrize(
    "entry",
    [
        {"ownership": "paper-a"},
        {"ownership": "paper-a", "release_included": False, "extra": 1},
        {"ownership": "paper-a", "release_included": 0},
        "paper-a",
    ],
)
def test_ownership_rejects_non_strict_entry_schema(
    tmp_path: Path,
    entry: object,
) -> None:
    path = tmp_path / "ownership.json"
    _write_manifest(path, {"README.md": entry})
    with pytest.raises(ValueError, match="entry must contain exactly|must be a boolean"):
        load_publication_manifest(path)


@pytest.mark.parametrize("name", ["/README.md", "../README.md", "docs//README.md"])
def test_ownership_rejects_noncanonical_paths(tmp_path: Path, name: str) -> None:
    path = tmp_path / "ownership.json"
    _write_manifest(path, {name: _entry()})
    with pytest.raises(ValueError, match="canonical and relative"):
        load_ownership(path)


def test_git_discovery_is_nul_safe_for_unicode_and_newlines(tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    _init_repository(repository)
    root = repository / "scope"
    root.mkdir()
    unusual = "unicode-Δ-and-newline\nfile.txt"
    (root / unusual).write_text("payload")
    subprocess.run(
        ["git", "-C", str(repository), "add", "--", f"scope/{unusual}"],
        check=True,
    )
    assert unusual in repository_paths(root)


def test_repository_root_dot_is_parsed_without_prefix_loss(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repository = tmp_path / "repository"
    _init_repository(repository)
    (repository / "root-file.txt").write_text("payload")
    monkeypatch.chdir(repository)
    assert "root-file.txt" in repository_paths(Path("."))


def test_audit_rejects_tracked_deletion(tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    _init_repository(repository)
    root = repository / "scope"
    root.mkdir()
    tracked = root / "tracked.txt"
    tracked.write_text("payload")
    subprocess.run(
        ["git", "-C", str(repository), "add", "--", "scope/tracked.txt"],
        check=True,
    )
    tracked.unlink()
    manifest = tmp_path / "ownership.json"
    _write_manifest(manifest, {"tracked.txt": _entry()})
    errors = audit_ownership(manifest, root)
    assert "repository path is missing from worktree: tracked.txt" in errors


@pytest.mark.parametrize("broken", [False, True])
def test_audit_rejects_symlinks_including_broken(
    tmp_path: Path,
    broken: bool,
) -> None:
    repository = tmp_path / "repository"
    _init_repository(repository)
    root = repository / "scope"
    root.mkdir()
    target = root / "target.txt"
    if not broken:
        target.write_text("payload")
    os.symlink("target.txt", root / "link.txt")
    manifest = tmp_path / "ownership.json"
    data = {"link.txt": _entry()}
    if not broken:
        data["target.txt"] = _entry()
    _write_manifest(manifest, data)
    errors = audit_ownership(manifest, root)
    assert "repository path must not be a symlink: link.txt" in errors


def test_reviewed_manifest_has_complete_repository_coverage() -> None:
    assert audit_ownership(MANIFEST, ROOT) == []


def test_reviewed_manifest_uses_only_declared_labels() -> None:
    ownership = load_ownership(MANIFEST)
    assert ownership
    assert set(ownership.values()) <= ALLOWED_LABELS


def test_every_release_bit_matches_the_complete_policy() -> None:
    entries = load_publication_manifest(MANIFEST)
    assert audit_release_policy(entries) == []


def test_release_contains_all_formal_paper_a_families() -> None:
    included = release_paths(MANIFEST)
    ownership = load_ownership(MANIFEST)
    required_exact = {
        "certificates/issue128-d4-groups.json",
        "certificates/issue128-d5-groups.json.gz",
        "certificates/issue128-d5-integrated-certificate.json",
        "scripts/reference_verify.py",
        "scripts/verify.py",
        "tests/test_certificate_mutations.py",
        "tests/test_reference_verify.py",
        "docs/manuscript/main.tex",
    }
    assert required_exact <= included
    formal_manuscript = {
        name
        for name in ownership
        if name.startswith("docs/manuscript/")
        and not name.endswith("/.gitignore")
    }
    assert formal_manuscript <= included


def test_release_excludes_every_forbidden_family() -> None:
    included = release_paths(MANIFEST)
    forbidden_tokens = ("fivefold", "d6", "d8", "e9", "processed")
    exceptions = {
        "docs/manuscript/figures/fivefold_gap.pdf",
        "src/trottercert/d6_physical_channels.py",
    }
    assert {
        name
        for name in included
        if any(token in name.lower() for token in forbidden_tokens)
    } == exceptions
    ownership = load_ownership(MANIFEST)
    assert not {
        name for name in included if ownership[name] in {"paper-b", "exploratory"}
    }


def test_shared_ownership_does_not_imply_release_inclusion() -> None:
    entries = load_publication_manifest(MANIFEST)
    shared = {name for name, entry in entries.items() if entry.ownership == "shared"}
    shared_included = {
        name for name in shared if entries[name].release_included
    }
    assert shared_included
    assert shared_included < shared
    assert {
        "src/trottercert/dual_word_manifest.py",
        "src/trottercert/effective_processor.py",
        "src/trottercert/spectral_processing.py",
    }.isdisjoint(shared_included)


def test_fivefold_figure_is_negative_limitations_evidence() -> None:
    entries = load_publication_manifest(MANIFEST)
    figure = entries["docs/manuscript/figures/fivefold_gap.pdf"]
    assert figure.ownership == "paper-a"
    assert figure.release_included


def test_curated_release_decisions_match_reviewed_exceptions() -> None:
    entries = load_publication_manifest(MANIFEST)
    expected = {
        "artifacts/d5-integrated/SHA256SUMS": ("paper-a", False),
        "artifacts/publication/paper-a-claim-SHA256SUMS": ("paper-a", True),
        "certificates/issue128-certificate.json": ("paper-a", True),
        "src/trottercert/d6_physical_channels.py": ("shared", True),
        "scripts/package_delivery.py": ("obsolete", False),
        "tests/test_delivery_package.py": ("obsolete", False),
    }
    assert {
        name: (entries[name].ownership, entries[name].release_included)
        for name in expected
    } == expected


def test_release_has_static_python_tex_and_claim_sha_closure() -> None:
    entries = load_publication_manifest(MANIFEST)
    assert audit_release_closure(entries, ROOT) == []


def _run_release_command(
    release: Path,
    arguments: list[str],
    *,
    timeout: int,
) -> None:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(release / "src")
    completed = subprocess.run(
        [sys.executable, *arguments],
        cwd=release,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    assert completed.returncode == 0, (
        f"isolated release command failed: {' '.join(arguments)}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )


@pytest.mark.slow
def test_isolated_release_smoke(tmp_path: Path) -> None:
    release = tmp_path / "paper-a-release"
    for relative in sorted(release_paths(MANIFEST)):
        source = ROOT / relative
        destination = release / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    copied_entries = load_publication_manifest(
        release / "artifacts/publication/paper-a-file-ownership.json"
    )
    assert audit_release_closure(copied_entries, release) == []
    certificate = "certificates/issue128-d5-integrated-certificate.json"
    _run_release_command(
        release,
        ["scripts/reference_verify.py", certificate],
        timeout=120,
    )
    _run_release_command(
        release,
        ["scripts/verify.py", "--deep", certificate],
        timeout=900,
    )
    _run_release_command(
        release,
        ["docs/manuscript/scripts/validate_claims.py"],
        timeout=120,
    )
    _run_release_command(
        release,
        [
            "-m",
            "pytest",
            "-q",
            "-p",
            "no:cacheprovider",
            "tests/test_certificate_mutations.py",
            "tests/test_reference_verify.py",
        ],
        timeout=300,
    )


def test_required_dual_track_families_are_explicitly_classified() -> None:
    ownership = load_ownership(MANIFEST)
    labels_by_path = {
        "scripts/certify_dual_e9_pairing.py": "paper-b",
        "tests/test_dual_e9_pairing.py": "paper-b",
        "artifacts/processed-kernel-audit/s10-s11-local-audit.json": "exploratory",
        "scripts/build_d6_direct.py": "exploratory",
        "scripts/build_d8_monolithic.py": "exploratory",
        "docs/superpowers/plans/2026-08-01-paper-a-quantum-implementation.md": "paper-a",
        "docs/superpowers/plans/2026-08-01-paper-b-prxq-gates-implementation.md": "paper-b",
    }
    assert {name: ownership.get(name) for name in labels_by_path} == labels_by_path
