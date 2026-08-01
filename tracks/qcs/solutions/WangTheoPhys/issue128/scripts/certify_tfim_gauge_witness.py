"""Build or verify the source-closed TFIM PF4 gauge witness artifact."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import sys
import tempfile
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path

ISSUE_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ISSUE_ROOT / "src"
for import_root in (ISSUE_ROOT, SOURCE_ROOT):
    if str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))

from scripts import reference_tfim_gauge_witness as reference
from trottercert.cubic_field import Cubic
from trottercert.tfim_gauge_witness import (
    CenteredPolynomialWitness,
    PauliOrbitPolynomialCertificate,
    TFIMGaugeWitnessResult,
    build_tfim_centered_polynomial_witnesses,
)

DEFAULT_OUTPUT = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/tfim-gauge-witness.json"
)
SOURCE_ROOTS = (
    "src/trottercert/tfim_gauge_witness.py",
    "scripts/certify_tfim_gauge_witness.py",
    "scripts/reference_tfim_gauge_witness.py",
    "pyproject.toml",
    "requirements-reproducibility.txt",
)


canonical_bytes = reference.canonical_bytes
payload_digest = reference.payload_digest


def _checked_path(relative: str) -> Path:
    if not relative or Path(relative).is_absolute() or Path(relative).as_posix() != relative:
        raise ValueError(f"source path is not canonical relative: {relative!r}")
    candidate = ISSUE_ROOT / relative
    if candidate.is_symlink():
        raise ValueError(f"source path must not be a symlink: {relative}")
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(ISSUE_ROOT.resolve(strict=True))
    except ValueError as error:
        raise ValueError(f"source path escapes issue root: {relative}") from error
    if not resolved.is_file():
        raise ValueError(f"source path is not a file: {relative}")
    return resolved


def _module_sources(module: str) -> set[str]:
    """Return a local module plus every executed parent package initializer."""

    if not module or module.startswith("."):
        return set()
    for base in (ISSUE_ROOT / "src", ISSUE_ROOT):
        parts = module.split(".")
        module_file = base.joinpath(*parts).with_suffix(".py")
        package_file = base.joinpath(*parts, "__init__.py")
        for candidate in (module_file, package_file):
            if candidate.is_file():
                sources = {candidate.relative_to(ISSUE_ROOT).as_posix()}
                for depth in range(1, len(parts)):
                    initializer = base.joinpath(*parts[:depth], "__init__.py")
                    if initializer.is_file():
                        sources.add(initializer.relative_to(ISSUE_ROOT).as_posix())
                return sources
    return set()


def _relative_import_module(source: str, node: ast.ImportFrom) -> str:
    source_path = Path(source)
    if source_path.parts[:1] != ("src",):
        raise ValueError(f"relative import outside src package: {source}")
    module_parts = list(source_path.with_suffix("").parts[1:-1])
    ascend = node.level - 1
    if ascend > len(module_parts):
        raise ValueError(f"relative import escapes package: {source}")
    if ascend:
        module_parts = module_parts[:-ascend]
    if node.module:
        module_parts.extend(node.module.split("."))
    return ".".join(module_parts)


def _local_imports(source: str) -> set[str]:
    path = _checked_path(source)
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=source)
    imports: set[str] = set()
    for node in ast.walk(tree):
        modules: list[tuple[str, tuple[str, ...]]] = []
        if isinstance(node, ast.Import):
            modules.extend((alias.name, ()) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            modules.append(
                (
                    _relative_import_module(source, node)
                    if node.level
                    else (node.module or ""),
                    tuple(alias.name for alias in node.names if alias.name != "*"),
                )
            )
        for module, imported_names in modules:
            local = _module_sources(module)
            if local:
                imports.update(local)
                # ``from package import submodule`` executes the imported
                # submodule even though the AST's ``module`` is only the
                # package name.  Include it when it resolves locally; plain
                # imported attributes correctly resolve to no extra file.
                for imported_name in imported_names:
                    imports.update(_module_sources(f"{module}.{imported_name}"))
            elif module.split(".", 1)[0] in {"trottercert", "scripts"}:
                raise ValueError(f"unresolved local import {module!r} from {source}")
    return imports


def source_closure() -> tuple[str, ...]:
    pending = list(SOURCE_ROOTS)
    closure: set[str] = set()
    while pending:
        relative = Path(pending.pop()).as_posix()
        _checked_path(relative)
        if relative in closure:
            continue
        closure.add(relative)
        if relative.endswith(".py"):
            pending.extend(sorted(_local_imports(relative) - closure))
    result = tuple(sorted(closure))
    if result != reference.SOURCE_PATHS:
        raise ValueError(
            "recursive implementation source closure differs from independent allowlist"
        )
    return result


def implementation_sources() -> dict[str, str]:
    return {
        relative: hashlib.sha256(_checked_path(relative).read_bytes()).hexdigest()
        for relative in source_closure()
    }


def _pair(value: Fraction | int) -> list[int]:
    exact = Fraction(value)
    return [exact.numerator, exact.denominator]


def _cubic(value: Cubic) -> list[list[int]]:
    return [_pair(value.a0), _pair(value.a1), _pair(value.a2)]


def _rational_orbit_terms(
    values: tuple[tuple[tuple[int, int], Fraction], ...]
) -> list[dict[str, object]]:
    return [
        {"monomial": [monomial[0], monomial[1]], "coefficient": _pair(coefficient)}
        for monomial, coefficient in values
    ]


def _cubic_orbit_terms(
    values: tuple[tuple[tuple[int, int], Cubic], ...]
) -> list[dict[str, object]]:
    return [
        {"monomial": [monomial[0], monomial[1]], "coefficient": _cubic(coefficient)}
        for monomial, coefficient in values
    ]


def _orbit_payload(
    certificate: PauliOrbitPolynomialCertificate,
) -> dict[str, object]:
    if not certificate.matches_closed_formula:
        raise ArithmeticError("Pauli orbit polynomial certificate is not closed")
    return {
        "tau_h2": _rational_orbit_terms(certificate.tau_h2_coefficients),
        "tau_h4": _rational_orbit_terms(certificate.tau_h4_coefficients),
        "tau_h_e5": _cubic_orbit_terms(certificate.tau_h_e5_coefficients),
        "tau_h3_e5": _cubic_orbit_terms(
            certificate.tau_h3_e5_coefficients
        ),
    }


def _witness_payload(witness: CenteredPolynomialWitness) -> dict[str, object]:
    return {
        "power": witness.power,
        "identity_coefficient": _pair(witness.identity_coefficient),
        "hamiltonian_coefficient": _pair(witness.hamiltonian_coefficient),
        "tau_w": _pair(witness.tau_w),
        "tau_wh": _pair(witness.tau_wh),
        "commutator_is_zero": witness.commutator_is_zero,
        "tau_w_e5": _cubic(witness.tau_w_e5),
    }


def _instance_payload(result: TFIMGaugeWitnessResult) -> dict[str, object]:
    if result.field != 1 or result.coupling != 1:
        raise ValueError("artifact checked instances are frozen at h=j=1")
    if result.first_nonzero_power != 3:
        raise ArithmeticError("checked TFIM instance lacks the frozen p=3 witness")
    return {
        "length": result.length,
        "regime": result.orbit_formula.regime,
        "parameters": {"h": [1, 1], "j": [1, 1]},
        "tau_h2": _pair(result.tau_h2),
        "tau_h4": _pair(result.tau_h4),
        "tau_h_e5": _cubic(result.tau_h_e5),
        "tau_h3_e5": _cubic(result.tau_h3_e5),
        "orbit_polynomials": _orbit_payload(result.orbit_polynomial_certificate),
        "witnesses": [_witness_payload(witness) for witness in result.witnesses],
        "first_nonzero_power": result.first_nonzero_power,
    }


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _mapping_binding() -> dict[str, object]:
    path = _checked_path(reference.MAPPING_PATH)
    raw = path.read_bytes()
    wrapper = json.loads(raw, object_pairs_hook=_unique_object)
    if not isinstance(wrapper, Mapping) or raw != canonical_bytes(wrapper):
        raise ValueError("PF4 mapping artifact is not canonical JSON")
    mapping = wrapper.get("mapping")
    if not isinstance(mapping, Mapping):
        raise ValueError("PF4 mapping payload is absent")
    mapping_digest = mapping.get("payload_sha256")
    if mapping_digest != reference.MAPPING_PAYLOAD_SHA256:
        raise ValueError("PF4 mapping payload differs from frozen authority")
    return {
        "path": reference.MAPPING_PATH,
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "mapping_payload_sha256": mapping_digest,
        "physical_bridge": "E5=L5",
    }


def _source_digest(sources: Mapping[str, str]) -> str:
    return hashlib.sha256(canonical_bytes(sources)).hexdigest()


def build_payload() -> dict[str, object]:
    instances = []
    for length in reference.CHECKED_LENGTHS:
        print(f"building exact centered-polynomial TFIM witness at L={length}", flush=True)
        result = build_tfim_centered_polynomial_witnesses(
            length, Fraction(1), Fraction(1)
        )
        instances.append(_instance_payload(result))
    sources = implementation_sources()
    payload: dict[str, object] = {
        "schema_version": reference.SCHEMA_VERSION,
        "kind": reference.ARTIFACT_KIND,
        "model": reference.EXPECTED_MODEL,
        "pf4_mapping_binding": _mapping_binding(),
        "gamma": [[37, 900000], [313, 14400000], [29, 1800000]],
        "formulas": reference.EXPECTED_FORMULAS,
        "search": reference.EXPECTED_SEARCH,
        "checked_lengths": list(reference.CHECKED_LENGTHS),
        "checked_instances": instances,
        "claim": reference.EXPECTED_CLAIM,
        "implementation_sources": sources,
        "implementation_sources_digest": _source_digest(sources),
        "payload_sha256": "",
        "hpc_authorized": False,
    }
    payload["payload_sha256"] = payload_digest(payload)
    reference.verify_payload(payload, ISSUE_ROOT)
    return payload


def load_payload(path: Path) -> dict[str, object]:
    return reference.load_payload(path)


def verify_payload(payload: object) -> None:
    reference.verify_payload(payload, ISSUE_ROOT)
    if isinstance(payload, Mapping) and payload.get("implementation_sources") != (
        implementation_sources()
    ):
        raise ValueError("primary implementation source closure mismatch")


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build", action="store_true")
    modes.add_argument("--verify", type=Path, metavar="ARTIFACT")
    modes.add_argument("--verify-full", type=Path, metavar="ARTIFACT")
    parser.add_argument("--output", type=Path, metavar="PATH")
    arguments = parser.parse_args(argv)
    try:
        if arguments.build:
            if arguments.output is None:
                parser.error("--build requires --output PATH")
            _atomic_write(arguments.output, canonical_bytes(build_payload()))
            print(f"wrote canonical TFIM gauge witness: {arguments.output}")
        elif arguments.verify is not None:
            if arguments.output is not None:
                parser.error("--output is valid only with --build")
            verify_payload(load_payload(arguments.verify))
            print("TFIM gauge witness is canonical and source-closed")
        else:
            if arguments.output is not None:
                parser.error("--output is valid only with --build")
            submitted = load_payload(arguments.verify_full)
            verify_payload(submitted)
            if submitted != build_payload():
                raise ValueError("TFIM gauge witness differs from exact full rebuild")
            print("TFIM gauge witness fully recomputed and verified")
    except (ArithmeticError, OSError, TypeError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
