"""Build the evidence-backed Paper A deferred-norming ablation table."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import pairwise
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
if __package__ in {None, ""}:
    sys.path.insert(0, str(ROOT))

from scripts.reference_verify import verify as reference_verify
from trottercert.refined_error import (
    RefinedFourthOrderConstants,
    build_refined_fourth_order_constants,
    certified_d4_cell_coefficients,
    evaluate_refined_fourth_order_bound,
)
from trottercert.resources import (
    fourth_order_four_matching_resources,
    required_steps,
)

STAGE_NAMES = (
    "published",
    "complete_formula",
    "pauli_aggregation",
    "translation_aggregation",
    "d4_anticommuting",
    "finite_ledger",
)
MONOTONIC_CHAIN = (
    "pauli_aggregation",
    "translation_aggregation",
    "d4_anticommuting",
    "finite_ledger",
)
CERTIFICATE_RELATIVE_PATH = Path("certificates/issue128-d5-integrated-certificate.json")
COMPLETE_FORMULA_MISSING_REASON = (
    "no frozen pre-Pauli complete-formula coefficient artifact or "
    "independent recomputation function exists"
)
SCIENTIFIC_SOURCE_NAMES = (
    "refined_error",
    "resources",
    "intervals",
    "local_commutators",
    "rigorous_fourth",
    "higher_order",
    "algebra",
)
STATIC_SOURCE_PATHS = {
    "generator": Path("scripts/build_paper_a_ablation.py"),
    "reference_verifier": Path("scripts/reference_verify.py"),
    "refined_error": Path("src/trottercert/refined_error.py"),
    "resources": Path("src/trottercert/resources.py"),
    "intervals": Path("src/trottercert/intervals.py"),
    "local_commutators": Path("src/trottercert/local_commutators.py"),
    "rigorous_fourth": Path("src/trottercert/rigorous_fourth.py"),
    "higher_order": Path("src/trottercert/higher_order.py"),
    "algebra": Path("src/trottercert/algebra.py"),
}
RECOMPUTE_FUNCTION_ALLOWLIST = {
    "published": (
        "trottercert.resources.required_steps",
        "trottercert.resources.fourth_order_four_matching_resources",
    ),
    "complete_formula": (),
    "pauli_aggregation": (
        "trottercert.refined_error.build_refined_fourth_order_constants",
        "scripts.build_paper_a_ablation._minimum_steps",
        "trottercert.refined_error.evaluate_refined_fourth_order_bound",
        "trottercert.resources.fourth_order_four_matching_resources",
    ),
    "translation_aggregation": (
        "trottercert.refined_error.build_refined_fourth_order_constants",
        "trottercert.refined_error.certified_d4_cell_coefficients",
        "scripts.build_paper_a_ablation._minimum_steps",
        "trottercert.refined_error.evaluate_refined_fourth_order_bound",
        "trottercert.resources.fourth_order_four_matching_resources",
    ),
    "d4_anticommuting": (
        "scripts.reference_verify.verify",
        "trottercert.refined_error.build_refined_fourth_order_constants",
        "scripts.build_paper_a_ablation._minimum_steps",
        "trottercert.refined_error.evaluate_refined_fourth_order_bound",
        "trottercert.resources.fourth_order_four_matching_resources",
    ),
    "finite_ledger": (
        "scripts.reference_verify.verify",
        "trottercert.refined_error.build_refined_fourth_order_constants",
        "scripts.build_paper_a_ablation._minimum_steps",
        "trottercert.refined_error.evaluate_refined_fourth_order_bound",
        "trottercert.resources.fourth_order_four_matching_resources",
    ),
}
_STAGE_FIELDS = {
    "stage",
    "status",
    "bound",
    "steps",
    "groups",
    "leading_d4_site_bound",
    "input_sources",
    "recompute_functions",
    "missing_evidence",
}


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _fraction(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or isinstance(value[0], bool)
        or not isinstance(value[0], int)
        or isinstance(value[1], bool)
        or not isinstance(value[1], int)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be an exact rational pair")
    return Fraction(value[0], value[1])


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_json(payload: object) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def compute_payload_digest(payload: object) -> str:
    """Hash canonical payload bytes with the digest field omitted."""

    if not isinstance(payload, dict):
        raise TypeError("ablation payload must be an object")
    digest_input = dict(payload)
    digest_input.pop("payload_sha256", None)
    return sha256(canonical_json(digest_input).encode()).hexdigest()


def _seal_payload(payload: dict[str, object]) -> dict[str, object]:
    sealed = dict(payload)
    sealed["payload_sha256"] = compute_payload_digest(sealed)
    return sealed


def _source(root: Path, path: Path) -> dict[str, str]:
    resolved = path.resolve()
    return {
        "path": resolved.relative_to(root.resolve()).as_posix(),
        "sha256": _digest(resolved),
    }


def _missing_stage(stage: str, reason: str) -> dict[str, object]:
    return {
        "stage": stage,
        "status": "missing_evidence",
        "bound": None,
        "steps": None,
        "groups": None,
        "leading_d4_site_bound": None,
        "input_sources": [],
        "recompute_functions": [],
        "missing_evidence": reason,
    }


def _computed_stage(
    *,
    stage: str,
    bound: Fraction,
    steps: int,
    n_sites: int,
    leading_d4_site_bound: Fraction | None,
    input_sources: list[str],
    recompute_functions: list[str],
) -> dict[str, object]:
    groups = fourth_order_four_matching_resources(
        n_sites,
        steps,
        stage_count=31,
    ).group_exponentials
    return {
        "stage": stage,
        "status": "computed",
        "bound": _pair(bound),
        "steps": steps,
        "groups": groups,
        "leading_d4_site_bound": (
            None if leading_d4_site_bound is None else _pair(leading_d4_site_bound)
        ),
        "input_sources": input_sources,
        "recompute_functions": recompute_functions,
        "missing_evidence": None,
    }


def _minimum_steps(
    *,
    constants: RefinedFourthOrderConstants,
    n_sites: int,
    tolerance: Fraction,
    upper: int,
    d4_site: Fraction | None,
    d5_site: Fraction | None,
) -> int:
    for steps in range(2, upper + 1):
        try:
            bound = evaluate_refined_fourth_order_bound(
                constants,
                n_sites,
                steps,
                d4_site_override=d4_site,
                d5_site_override=d5_site,
            ).global_error_bound
        except ValueError as exc:
            if "convergence" not in str(exc):
                raise
            continue
        if bound <= tolerance:
            return steps
    raise ValueError("no accepted step count exists within the published baseline")


def _refined_stage(
    *,
    stage: str,
    constants: RefinedFourthOrderConstants,
    n_sites: int,
    tolerance: Fraction,
    upper: int,
    reference_steps: int,
    d4_site: Fraction | None,
    d5_site: Fraction | None,
    input_sources: list[str],
    recompute_functions: list[str],
) -> dict[str, object]:
    steps = _minimum_steps(
        constants=constants,
        n_sites=n_sites,
        tolerance=tolerance,
        upper=upper,
        d4_site=d4_site,
        d5_site=d5_site,
    )
    comparison = evaluate_refined_fourth_order_bound(
        constants,
        n_sites,
        reference_steps,
        d4_site_override=d4_site,
        d5_site_override=d5_site,
    )
    leading = constants.d4_site if d4_site is None else d4_site
    return _computed_stage(
        stage=stage,
        bound=comparison.global_error_bound,
        steps=steps,
        n_sites=n_sites,
        leading_d4_site_bound=leading,
        input_sources=input_sources,
        recompute_functions=recompute_functions,
    )


def _recompute_ablation(certificate: Path, root: Path) -> dict[str, object]:
    """Recompute the unique unsealed payload from authoritative inputs."""

    certificate = certificate.resolve()
    root = root.resolve()
    if certificate != root / CERTIFICATE_RELATIVE_PATH:
        raise ValueError("Paper A ablation requires the frozen certificate path")
    raw = json.loads(certificate.read_text())
    if not isinstance(raw, dict) or raw.get("schema_version") != 3:
        raise ValueError("Paper A ablation requires a schema-v3 certificate")
    reference = reference_verify(certificate)
    if reference.get("valid") is not True:
        raise ValueError("reference verifier did not accept the certificate")

    benchmark = raw["benchmark"]
    published = raw["published_baseline"]
    candidate = raw["candidate"]
    resources = raw["claimed_resources"]
    n_sites = int(benchmark["length"]) ** 2
    tolerance = _fraction(benchmark["tolerance"], "tolerance")
    reference_steps = int(candidate["steps"])
    if reference_steps != 95:
        raise ValueError("frozen Paper A candidate must use 95 steps")
    if int(resources["candidate_group_exponentials"]) != 2_851:
        raise ValueError("frozen Paper A candidate must use 2,851 groups")

    d4_path = certificate.parent / candidate["d4_certificate"]["path"]
    d5_path = certificate.parent / candidate["d5_certificate"]["path"]
    sources = {
        "certificate": _source(root, certificate),
        "d4_sidecar": _source(root, d4_path),
        "d5_sidecar": _source(root, d5_path),
        **{
            name: _source(root, root / relative_path)
            for name, relative_path in STATIC_SOURCE_PATHS.items()
        },
    }

    published_density = _fraction(
        published["site_density_upper"],
        "published site density",
    )
    published_steps = required_steps(
        published_density * n_sites,
        tolerance,
        4,
    )
    if published_steps != int(published["steps"]):
        raise ValueError("published step count regeneration mismatch")
    published_bound = published_density * n_sites / reference_steps**4
    stages: list[dict[str, object]] = [
        _computed_stage(
            stage="published",
            bound=published_bound,
            steps=published_steps,
            n_sites=n_sites,
            leading_d4_site_bound=None,
            input_sources=["certificate", "generator", "resources"],
            recompute_functions=list(RECOMPUTE_FUNCTION_ALLOWLIST["published"]),
        ),
        _missing_stage(
            "complete_formula",
            COMPLETE_FORMULA_MISSING_REASON,
        ),
    ]

    constants = build_refined_fourth_order_constants(
        decimal_digits=int(candidate["coefficient_interval_decimal_digits"]),
        quantization_digits=int(candidate["e5_quantization_digits"]),
    )
    upper = published_steps
    stages.append(
        _refined_stage(
            stage="pauli_aggregation",
            constants=constants,
            n_sites=n_sites,
            tolerance=tolerance,
            upper=upper,
            reference_steps=reference_steps,
            d4_site=None,
            d5_site=None,
            input_sources=[
                "certificate",
                "generator",
                *SCIENTIFIC_SOURCE_NAMES,
            ],
            recompute_functions=list(RECOMPUTE_FUNCTION_ALLOWLIST["pauli_aggregation"]),
        )
    )

    translated = certified_d4_cell_coefficients(
        constants.stages,
        quantization_digits=int(candidate["e5_quantization_digits"]),
    )
    translation_d4_site = (
        sum(
            (coefficient.abs_upper() for coefficient in translated.values()),
            Fraction(),
        )
        / 4
    )
    stages.append(
        _refined_stage(
            stage="translation_aggregation",
            constants=constants,
            n_sites=n_sites,
            tolerance=tolerance,
            upper=upper,
            reference_steps=reference_steps,
            d4_site=translation_d4_site,
            d5_site=None,
            input_sources=[
                "certificate",
                "generator",
                *SCIENTIFIC_SOURCE_NAMES,
            ],
            recompute_functions=list(
                RECOMPUTE_FUNCTION_ALLOWLIST["translation_aggregation"]
            ),
        )
    )

    d4_site = (
        _fraction(
            candidate["d4_certificate"]["cell_norm_upper"],
            "D4 cell norm",
        )
        / 4
    )
    stages.append(
        _refined_stage(
            stage="d4_anticommuting",
            constants=constants,
            n_sites=n_sites,
            tolerance=tolerance,
            upper=upper,
            reference_steps=reference_steps,
            d4_site=d4_site,
            d5_site=None,
            input_sources=[
                "certificate",
                "d4_sidecar",
                "generator",
                "reference_verifier",
                *SCIENTIFIC_SOURCE_NAMES,
            ],
            recompute_functions=list(RECOMPUTE_FUNCTION_ALLOWLIST["d4_anticommuting"]),
        )
    )

    d5_site = _fraction(
        candidate["d5_certificate"]["site_norm_upper"],
        "D5 site norm",
    )
    final = _refined_stage(
        stage="finite_ledger",
        constants=constants,
        n_sites=n_sites,
        tolerance=tolerance,
        upper=upper,
        reference_steps=reference_steps,
        d4_site=d4_site,
        d5_site=d5_site,
        input_sources=[
            "certificate",
            "d4_sidecar",
            "d5_sidecar",
            "generator",
            "reference_verifier",
            *SCIENTIFIC_SOURCE_NAMES,
        ],
        recompute_functions=list(RECOMPUTE_FUNCTION_ALLOWLIST["finite_ledger"]),
    )
    if final["bound"] != candidate["global_error_upper"]:
        raise ValueError("finite-ledger bound regeneration mismatch")
    stages.append(final)

    payload = {
        "schema_version": 1,
        "bound_semantics": "global_operator_error_upper_at_fixed_reference_steps",
        "reference_steps": reference_steps,
        "tolerance": benchmark["tolerance"],
        "monotonic_chain": list(MONOTONIC_CHAIN),
        "sources": sources,
        "stages": stages,
    }
    return payload


def build_ablation(certificate: Path) -> dict[str, object]:
    """Recompute and seal the unique payload from the frozen inputs."""

    return _seal_payload(_recompute_ablation(certificate, ROOT))


def _authoritative_source_paths(root: Path) -> dict[str, Path]:
    certificate = (root / CERTIFICATE_RELATIVE_PATH).resolve()
    raw = json.loads(certificate.read_text())
    candidate = raw["candidate"]
    paths = {
        "certificate": certificate,
        "d4_sidecar": certificate.parent / str(candidate["d4_certificate"]["path"]),
        "d5_sidecar": certificate.parent / str(candidate["d5_certificate"]["path"]),
        **{
            name: root / relative_path
            for name, relative_path in STATIC_SOURCE_PATHS.items()
        },
    }
    resolved: dict[str, Path] = {}
    for name, path in paths.items():
        candidate_path = path.resolve()
        candidate_path.relative_to(root.resolve())
        if not candidate_path.is_file():
            raise ValueError(f"authoritative source is missing: {name}")
        resolved[name] = candidate_path
    return resolved


def _authoritative_fingerprint(root: Path) -> tuple[tuple[str, str], ...]:
    return tuple(
        (name, _digest(path))
        for name, path in sorted(_authoritative_source_paths(root).items())
    )


@lru_cache(maxsize=4)
def _authoritative_payload_bytes(
    root_text: str,
    fingerprint: tuple[tuple[str, str], ...],
) -> bytes:
    """Rebuild once per exact authoritative source fingerprint."""

    if not fingerprint:
        raise ValueError("authoritative source fingerprint is empty")
    root = Path(root_text)
    payload = _seal_payload(_recompute_ablation(root / CERTIFICATE_RELATIVE_PATH, root))
    return canonical_json(payload).encode()


def _validate_stage_shells(payload: dict[str, object]) -> list[str]:
    errors: list[str] = []
    stages = payload.get("stages")
    if not isinstance(stages, list) or any(not isinstance(row, dict) for row in stages):
        return ["ablation stages must be a list of objects"]
    if [row.get("stage") for row in stages] != list(STAGE_NAMES):
        return ["ablation stage order mismatch"]

    for row in stages:
        stage = str(row["stage"])
        if set(row) != _STAGE_FIELDS:
            errors.append(f"{stage}: invalid stage schema")
            continue
        expected_functions = list(RECOMPUTE_FUNCTION_ALLOWLIST[stage])
        if row["recompute_functions"] != expected_functions:
            errors.append(f"{stage}: recompute-function allowlist mismatch")

    complete = stages[1]
    if complete != _missing_stage("complete_formula", COMPLETE_FORMULA_MISSING_REASON):
        errors.append(
            "complete_formula must remain missing_evidence until an "
            "authoritative artifact and recomputation function exist"
        )
    return errors


def validate_ablation(
    payload: object,
    root: Path = ROOT,
    serialized: bytes | None = None,
) -> list[str]:
    """Fail closed against a canonical rebuild from authoritative inputs."""

    if not isinstance(payload, dict):
        return ["ablation payload must be an object"]
    root = root.resolve()
    errors: list[str] = []
    expected_top = {
        "schema_version",
        "bound_semantics",
        "reference_steps",
        "tolerance",
        "monotonic_chain",
        "payload_sha256",
        "sources",
        "stages",
    }
    if set(payload) != expected_top:
        return ["ablation payload has an unexpected top-level schema"]

    if payload["payload_sha256"] != compute_payload_digest(payload):
        errors.append("payload digest mismatch")
    errors.extend(_validate_stage_shells(payload))

    try:
        authoritative_paths = _authoritative_source_paths(root)
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        errors.append(f"cannot resolve authoritative sources: {exc}")
        return errors

    sources = payload.get("sources")
    if not isinstance(sources, dict) or set(sources) != set(authoritative_paths):
        errors.append("ablation source set mismatch")
    else:
        for name, path in authoritative_paths.items():
            expected_source = _source(root, path)
            if sources.get(name) != expected_source:
                errors.append(f"authoritative source binding mismatch for {name}")

    try:
        fingerprint = _authoritative_fingerprint(root)
        expected_bytes = _authoritative_payload_bytes(str(root), fingerprint)
        expected_payload = json.loads(expected_bytes)
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        errors.append(f"authoritative canonical rebuild failed: {exc}")
        return errors

    if payload != expected_payload:
        errors.append("payload fields differ from authoritative canonical rebuild")
    actual_bytes = (
        canonical_json(payload).encode() if serialized is None else serialized
    )
    if actual_bytes != expected_bytes:
        errors.append("payload bytes differ from authoritative canonical rebuild")

    certificate = json.loads(authoritative_paths["certificate"].read_text())
    final = expected_payload["stages"][-1]
    candidate = certificate["candidate"]
    resources = certificate["claimed_resources"]
    if (
        final["bound"] != candidate["global_error_upper"]
        or final["steps"] != candidate["steps"]
        or final["groups"] != resources["candidate_group_exponentials"]
    ):
        errors.append("authoritative final row does not exactly match certificate")
    supplied_stages = payload["stages"]
    if isinstance(supplied_stages, list) and len(supplied_stages) == len(STAGE_NAMES):
        supplied_final = supplied_stages[-1]
        if not isinstance(supplied_final, dict) or (
            supplied_final.get("bound") != candidate["global_error_upper"]
            or supplied_final.get("steps") != candidate["steps"]
            or supplied_final.get("groups") != resources["candidate_group_exponentials"]
        ):
            errors.append("supplied final row does not exactly match certificate")

    chain = [
        next(row for row in expected_payload["stages"] if row["stage"] == name)
        for name in MONOTONIC_CHAIN
    ]
    bounds = [_fraction(row["bound"], f"{row['stage']} bound") for row in chain]
    steps = [int(row["steps"]) for row in chain]
    if any(right > left for left, right in pairwise(bounds)):
        errors.append("authoritative fixed-step bounds are not monotone")
    if any(right > left for left, right in pairwise(steps)):
        errors.append("authoritative certified steps are not monotone")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    payload = build_ablation(arguments.certificate)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(canonical_json(payload))
    print(f"output={arguments.output}")
    for row in payload["stages"]:
        print(
            f"stage={row['stage']} status={row['status']} "
            f"steps={row['steps']} groups={row['groups']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
