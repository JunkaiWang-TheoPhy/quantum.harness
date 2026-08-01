"""Build and verify the source-closed exact PF4 BCH mapping artifact."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

from trottercert.pf4_bch_mapping import parse_payload, to_payload

ISSUE_ROOT = Path(__file__).resolve().parents[1]
WRAPPER_SCHEMA_VERSION = 2
WRAPPER_KIND = "pf4_bch_mapping_source_closed_artifact"
REFERENCE_ALGORITHM = "recursive_bch_hall_lyndon_v1"
SOURCE_ROOTS = (
    "src/trottercert/pf4_bch_mapping.py",
    "src/trottercert/cubic_field.py",
    "src/trottercert/intervals.py",
    "scripts/derive_pf4_bch_mapping.py",
    "scripts/reference_pf4_bch_mapping.py",
    "pyproject.toml",
    "requirements-reproducibility.txt",
)
WRAPPER_FIELDS = {
    "schema_version",
    "kind",
    "mapping",
    "implementation_sources",
    "reference_algorithm",
    "payload_sha256",
}


def canonical_bytes(payload: object) -> bytes:
    """Return the one accepted compact JSON encoding."""

    return (
        json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def payload_digest(payload: Mapping[str, object]) -> str:
    unsigned = {key: value for key, value in payload.items() if key != "payload_sha256"}
    return hashlib.sha256(canonical_bytes(unsigned)).hexdigest()


def _checked_path(relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ValueError(f"source path is not a relative path: {relative!r}")
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


def _module_source(module: str) -> str | None:
    if not module or module.startswith("."):
        return None
    parts = module.split(".")
    bases = (ISSUE_ROOT / "src", ISSUE_ROOT)
    for base in bases:
        module_file = base.joinpath(*parts).with_suffix(".py")
        package_file = base.joinpath(*parts, "__init__.py")
        for candidate in (module_file, package_file):
            if candidate.is_file():
                return candidate.relative_to(ISSUE_ROOT).as_posix()
    return None


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
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=source)
    except (SyntaxError, UnicodeDecodeError) as error:
        raise ValueError(f"cannot parse source dependency {source}") from error
    imports: set[str] = set()
    for node in ast.walk(tree):
        modules: list[str] = []
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            module = (
                _relative_import_module(source, node)
                if node.level
                else (node.module or "")
            )
            modules.append(module)
        for module in modules:
            local = _module_source(module)
            if local is not None:
                imports.add(local)
            elif module.split(".", 1)[0] in {"trottercert", "scripts"}:
                raise ValueError(
                    f"unresolved local import {module!r} from source {source}"
                )
    return imports


def source_closure() -> tuple[str, ...]:
    """Return normalized roots plus all recursively imported local sources."""

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
    return tuple(sorted(closure))


def implementation_sources() -> dict[str, str]:
    return {
        relative: hashlib.sha256(_checked_path(relative).read_bytes()).hexdigest()
        for relative in source_closure()
    }


def build_payload() -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": WRAPPER_SCHEMA_VERSION,
        "kind": WRAPPER_KIND,
        "mapping": to_payload(),
        "implementation_sources": implementation_sources(),
        "reference_algorithm": REFERENCE_ALGORITHM,
        "payload_sha256": "",
    }
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def _reject_json_numeric_aliases(value: object, path: str = "artifact") -> None:
    if isinstance(value, float):
        raise TypeError(f"{path}: floating-point values are forbidden")
    if isinstance(value, list):
        for index, item in enumerate(value):
            _reject_json_numeric_aliases(item, f"{path}[{index}]")
    elif isinstance(value, dict):
        for key, item in value.items():
            _reject_json_numeric_aliases(item, f"{path}.{key}")


def verify_payload(payload: object) -> None:
    if not isinstance(payload, Mapping):
        raise TypeError("PF4 wrapper artifact must be an object")
    _reject_json_numeric_aliases(payload)
    if set(payload) != WRAPPER_FIELDS:
        raise ValueError("PF4 wrapper field set mismatch")
    if payload.get("schema_version") != WRAPPER_SCHEMA_VERSION:
        raise ValueError("PF4 wrapper schema version mismatch")
    if payload.get("kind") != WRAPPER_KIND:
        raise ValueError("PF4 wrapper kind mismatch")
    if payload.get("reference_algorithm") != REFERENCE_ALGORITHM:
        raise ValueError("PF4 reference algorithm mismatch")
    digest = payload.get("payload_sha256")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
        or digest != payload_digest(payload)
    ):
        raise ValueError("PF4 wrapper payload digest mismatch")
    mapping = payload.get("mapping")
    if not isinstance(mapping, Mapping):
        raise TypeError("PF4 mathematical mapping is missing")
    parse_payload(mapping)
    sources = payload.get("implementation_sources")
    if not isinstance(sources, Mapping) or any(
        not isinstance(path, str)
        or not isinstance(source_digest, str)
        or len(source_digest) != 64
        or any(character not in "0123456789abcdef" for character in source_digest)
        for path, source_digest in sources.items()
    ):
        raise ValueError("PF4 implementation source manifest is malformed")
    if dict(sources) != implementation_sources():
        raise ValueError("PF4 implementation source closure mismatch")
    if dict(payload) != build_payload():
        raise ValueError("PF4 wrapper differs from authoritative rebuild")


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def load_payload(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    try:
        payload = json.loads(raw, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON: {error}") from error
    if not isinstance(payload, dict):
        raise TypeError("PF4 wrapper artifact must be an object")
    if raw != canonical_bytes(payload):
        raise ValueError("PF4 wrapper bytes are not canonical JSON")
    return payload


def _build(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(build_payload()))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--build", type=Path, metavar="OUTPUT")
    modes.add_argument("--verify", type=Path, metavar="ARTIFACT")
    arguments = parser.parse_args(argv)
    try:
        if arguments.build is not None:
            _build(arguments.build)
            print(f"wrote canonical PF4 mapping artifact: {arguments.build}")
        else:
            verify_payload(load_payload(arguments.verify))
            print("PF4 mapping artifact is canonical and source-closed")
    except (ArithmeticError, OSError, TypeError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
