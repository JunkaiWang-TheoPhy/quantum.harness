from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.audit_paper_a_claims import audit_claims

CERTIFICATE = "certificates/issue128-d5-integrated-certificate.json"
REFERENCE = "scripts/reference_verify.py"
TRANSCRIPT = "artifacts/d5-integrated/verification-transcript.txt"
MUTATIONS = "tests/test_certificate_mutations.py"


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _record(value: object, location: str) -> dict[str, object]:
    return {
        "value": value,
        "source": CERTIFICATE,
        "checker": REFERENCE,
        "manuscript_locations": [location],
    }


def _write_fixture(root: Path) -> Path:
    for directory in (
        "certificates",
        "scripts",
        "tests",
        "docs/manuscript/sections",
        "artifacts/d5-integrated",
        "artifacts/publication",
    ):
        (root / directory).mkdir(parents=True)

    certificate = {
        "benchmark": {
            "model": "periodic_square_spin_half_isotropic_heisenberg",
            "normalization": "(XX+YY+ZZ)/4",
            "length": 12,
            "time": [1, 1],
            "tolerance": [1, 1_000_000],
        },
        "published_baseline": {"steps": 393, "group_exponentials": 11_791},
        "candidate": {
            "steps": 95,
            "group_exponentials": 2_851,
            "d4_certificate": {"term_count": 75_324, "group_count": 7_576},
            "d5_certificate": {"term_count": 605_832, "group_count": 123_106},
        },
        "claimed_resources": {
            "published_group_exponentials": 11_791,
            "candidate_group_exponentials": 2_851,
        },
        "claims": {"exact_improvement_ratio": [11_791, 2_851]},
    }
    (root / CERTIFICATE).write_text(json.dumps(certificate), encoding="utf-8")
    (root / REFERENCE).write_text("# standard-library checker\n", encoding="utf-8")
    (root / MUTATIONS).write_text("# mutation checker\n", encoding="utf-8")
    (root / "README.md").write_text("not evidence\n", encoding="utf-8")
    (root / TRANSCRIPT).write_text(
        "NORMAL TEST SUITE\nObserved:\n"
        "  146 passed, 12 deselected in 1.0s\n"
        "FOCUSED ADVERSARIAL SUITES\nObserved:\n"
        "  18 passed in 1.0s\n",
        encoding="utf-8",
    )

    sections = {
        "abstract.tex": (
            "periodic $12\\times12$ spin-$1/2$ isotropic Heisenberg model "
            "at $T=1$ and tolerance $10^{-6}$ accepts 95 steps and rejects 94. "
            "The exact result is 11791/2851=4.135741844966678. "
            "D4 and D5 sidecars contain respectively 75,324 and 605,832 terms, "
            "partitioned into 7,576 and 123,106 groups.\n"
        ),
        "problem.tex": r"$h_{uv}=\frac{X_uX_v+Y_uY_v+Z_uZ_v}{4}$" + "\n",
        "introduction.tex": "a pinned 393-step instantiation\n",
        "results.tex": (
            "Pinned published-theorem control & 393 & 11,791 \\\\\n"
            "Proof-carrying model-aware bound & 95 & 2,851 \\\\\n"
        ),
        "reproducibility.tex": (
            "Verifier falsifiability is checked by field-mutation tests; "
            "each mutation must fail.\n"
        ),
        "verification.tex": (
            "Discovery jobs are untrusted producers. A third program uses "
            "only the Python standard library.\n"
        ),
        "limitations.tex": (
            "\\subsection{The fivefold boundary is open}\n"
            "This is arithmetic, not a 78-step error certificate.\n"
            "$r=78$ is therefore not certified.\n"
        ),
    }
    for name, text in sections.items():
        (root / "docs/manuscript/sections" / name).write_text(text, encoding="utf-8")

    claims = {
        "benchmark_model": _record(
            "periodic_square_spin_half_isotropic_heisenberg",
            "docs/manuscript/sections/abstract.tex",
        ),
        "benchmark_normalization": _record(
            "(XX+YY+ZZ)/4", "docs/manuscript/sections/problem.tex"
        ),
        "lattice_size": _record([12, 12], "docs/manuscript/sections/abstract.tex"),
        "simulation_time": _record([1, 1], "docs/manuscript/sections/abstract.tex"),
        "tolerance": _record([1, 1_000_000], "docs/manuscript/sections/abstract.tex"),
        "accepted_steps": _record(95, "docs/manuscript/sections/abstract.tex"),
        "adjacent_rejected_steps": _record(
            94, "docs/manuscript/sections/abstract.tex"
        ),
        "published_steps": _record(393, "docs/manuscript/sections/introduction.tex"),
        "published_groups": _record(11_791, "docs/manuscript/sections/results.tex"),
        "candidate_groups": _record(2_851, "docs/manuscript/sections/results.tex"),
        "published_ratio": _record(
            [11_791, 2_851], "docs/manuscript/sections/abstract.tex"
        ),
        "d4_term_count": _record(75_324, "docs/manuscript/sections/abstract.tex"),
        "d4_group_count": _record(7_576, "docs/manuscript/sections/abstract.tex"),
        "d5_term_count": _record(605_832, "docs/manuscript/sections/abstract.tex"),
        "d5_group_count": _record(123_106, "docs/manuscript/sections/abstract.tex"),
        "test_evidence": {
            "value": {
                "normal_passed": 146,
                "normal_deselected": 12,
                "mutation_passed": 18,
            },
            "source": TRANSCRIPT,
            "checker": MUTATIONS,
            "manuscript_locations": [
                "docs/manuscript/sections/reproducibility.tex"
            ],
        },
        "trusted_computing_base": {
            "value": {
                "discovery_trusted": False,
                "reference_verifier_imports_trottercert": False,
            },
            "source": "docs/manuscript/sections/verification.tex",
            "checker": REFERENCE,
            "manuscript_locations": ["docs/manuscript/sections/verification.tex"],
        },
    }
    matrix = {
        "schema_version": 1,
        "sha256_manifest": "artifacts/publication/paper-a-claim-SHA256SUMS",
        "excluded_pending_claims": [
            {
                "name": "strengthened_control_ratio",
                "reason": "Excluded until independent evidence exists.",
            },
            {
                "name": "fivefold_global_certificate",
                "reason": "Open and not certified by this release.",
            },
        ],
        "claims": claims,
    }
    matrix_path = root / "matrix.json"
    matrix_path.write_text(json.dumps(matrix), encoding="utf-8")
    _rewrite_manifest(root)
    return matrix_path


def _rewrite_manifest(root: Path) -> None:
    relative_paths = [
        CERTIFICATE,
        REFERENCE,
        TRANSCRIPT,
        MUTATIONS,
        *(f"docs/manuscript/sections/{name}" for name in (
            "abstract.tex",
            "problem.tex",
            "introduction.tex",
            "results.tex",
            "reproducibility.tex",
            "verification.tex",
            "limitations.tex",
        )),
    ]
    manifest = root / "artifacts/publication/paper-a-claim-SHA256SUMS"
    manifest.write_text(
        "\n".join(f"{_digest(root / path)}  {path}" for path in relative_paths) + "\n",
        encoding="utf-8",
    )


def _load_matrix(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_matrix(path: Path, matrix: dict[str, object]) -> None:
    path.write_text(json.dumps(matrix), encoding="utf-8")


def test_claim_audit_accepts_bound_fixture(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    assert audit_claims(matrix_path, tmp_path) == []


def test_claim_audit_rejects_wrong_published_ratio(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["claims"]["published_ratio"]["value"] = [1, 2]
    _write_matrix(matrix_path, matrix)
    errors = audit_claims(matrix_path, tmp_path)
    assert "published ratio must equal 11791/2851" in errors


def test_claim_audit_rejects_missing_required_claim(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    del matrix["claims"]["accepted_steps"]
    _write_matrix(matrix_path, matrix)
    assert any("missing required claims accepted_steps" in error for error in audit_claims(matrix_path, tmp_path))


def test_claim_audit_rejects_readme_source_replacement(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["claims"]["published_ratio"]["source"] = "README.md"
    _write_matrix(matrix_path, matrix)
    assert any("source must be" in error for error in audit_claims(matrix_path, tmp_path))


def test_claim_audit_rejects_transcript_count_drift(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["claims"]["test_evidence"]["value"]["mutation_passed"] = 19
    _write_matrix(matrix_path, matrix)
    assert "test_evidence does not match the verification transcript" in audit_claims(matrix_path, tmp_path)


def test_claim_audit_rejects_noncanonical_claim_path(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["claims"]["published_ratio"]["source"] = (
        "certificates/../certificates/issue128-d5-integrated-certificate.json"
    )
    _write_matrix(matrix_path, matrix)
    errors = audit_claims(matrix_path, tmp_path)
    assert any("canonical repository-relative path" in error for error in errors)


def test_claim_audit_rejects_matrix_outside_root(tmp_path: Path) -> None:
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "matrix.json"
    outside.write_text("{}", encoding="utf-8")
    assert audit_claims(outside, root) == ["matrix_path must be inside the Paper A root"]


def test_claim_audit_rejects_integer_substring_marker(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    abstract = tmp_path / "docs/manuscript/sections/abstract.tex"
    abstract.write_text(
        abstract.read_text(encoding="utf-8").replace("accepts 95 steps", "accepts 195 steps"),
        encoding="utf-8",
    )
    _rewrite_manifest(tmp_path)
    errors = audit_claims(matrix_path, tmp_path)
    assert "accepted_steps: stable manuscript marker is absent" in errors


def test_claim_audit_rejects_adjacent_step_digit_suffix(tmp_path: Path) -> None:
    for replacement in ("rejects 949", "rejects 940"):
        case_root = tmp_path / replacement.replace(" ", "-")
        matrix_path = _write_fixture(case_root)
        abstract = case_root / "docs/manuscript/sections/abstract.tex"
        abstract.write_text(
            abstract.read_text(encoding="utf-8").replace("rejects 94", replacement),
            encoding="utf-8",
        )
        _rewrite_manifest(case_root)
        errors = audit_claims(matrix_path, case_root)
        assert "adjacent_rejected_steps: stable manuscript marker is absent" in errors


def test_excluded_strengthened_number_is_forbidden_in_manuscript(
    tmp_path: Path,
) -> None:
    for token in ("10,591", "10591/2851"):
        case_root = tmp_path / token.replace(",", "-").replace("/", "-")
        matrix_path = _write_fixture(case_root)
        abstract = case_root / "docs/manuscript/sections/abstract.tex"
        abstract.write_text(
            abstract.read_text(encoding="utf-8") + f"Excluded control {token}.\n",
            encoding="utf-8",
        )
        _rewrite_manifest(case_root)
        errors = audit_claims(matrix_path, case_root)
        assert any("excluded strengthened-control numeric token" in error for error in errors)


def test_positive_fivefold_claim_is_forbidden_after_hash_refresh(
    tmp_path: Path,
) -> None:
    matrix_path = _write_fixture(tmp_path)
    limitations = tmp_path / "docs/manuscript/sections/limitations.tex"
    limitations.write_text(
        limitations.read_text(encoding="utf-8")
        + "The fivefold target is achieved at r=78.\n",
        encoding="utf-8",
    )
    _rewrite_manifest(tmp_path)
    errors = audit_claims(matrix_path, tmp_path)
    assert any("forbidden positive fivefold claim" in error for error in errors)


def test_fivefold_open_marker_is_required(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    limitations = tmp_path / "docs/manuscript/sections/limitations.tex"
    limitations.write_text("Fivefold discussion removed.\n", encoding="utf-8")
    _rewrite_manifest(tmp_path)
    assert "fivefold open/not-certified manuscript marker is absent" in audit_claims(
        matrix_path, tmp_path
    )


def test_claim_audit_requires_dedicated_manifest_path(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["sha256_manifest"] = "artifacts/d5-integrated/SHA256SUMS"
    (tmp_path / "artifacts/d5-integrated/SHA256SUMS").write_text("", encoding="utf-8")
    _write_matrix(matrix_path, matrix)
    errors = audit_claims(matrix_path, tmp_path)
    assert any("sha256_manifest must be" in error for error in errors)


def test_claim_manifest_rejects_nonclaim_entry(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    manifest = tmp_path / "artifacts/publication/paper-a-claim-SHA256SUMS"
    with manifest.open("a", encoding="utf-8") as handle:
        handle.write(f"{_digest(tmp_path / 'README.md')}  README.md\n")
    errors = audit_claims(matrix_path, tmp_path)
    assert any("claim manifest has non-claim entries README.md" in error for error in errors)


def test_claim_manifest_rejects_missing_claim_entry(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    manifest = tmp_path / "artifacts/publication/paper-a-claim-SHA256SUMS"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    manifest.write_text(
        "\n".join(line for line in lines if not line.endswith("verification.tex"))
        + "\n",
        encoding="utf-8",
    )
    errors = audit_claims(matrix_path, tmp_path)
    assert any("claim manifest missing exact entries" in error for error in errors)


def test_claim_audit_rejects_hash_mismatch(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    (tmp_path / REFERENCE).write_text("# drifted checker\n", encoding="utf-8")
    assert any("SHA-256 mismatch" in error for error in audit_claims(matrix_path, tmp_path))


def test_claim_audit_rejects_pending_claim_value_or_source(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["excluded_pending_claims"][0]["value"] = [1, 2]
    _write_matrix(matrix_path, matrix)
    errors = audit_claims(matrix_path, tmp_path)
    assert any("must contain exactly name and reason" in error for error in errors)


def test_claim_audit_rejects_unknown_top_level_field(tmp_path: Path) -> None:
    matrix_path = _write_fixture(tmp_path)
    matrix = _load_matrix(matrix_path)
    matrix["notes"] = "not schema-bound"
    _write_matrix(matrix_path, matrix)
    assert any("unknown top-level fields notes" in error for error in audit_claims(matrix_path, tmp_path))


def test_repository_claim_matrix_passes() -> None:
    root = Path(__file__).resolve().parents[1]
    matrix = root / "artifacts/publication/paper-a-claim-matrix.json"
    assert audit_claims(matrix, root) == []
