"""Independent standard-library verifier for the TFIM PF4 gauge witness."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path

ISSUE_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 1
ARTIFACT_KIND = "tfim_pf4_centered_polynomial_gauge_witness"
MAPPING_PATH = "docs/experiments/processor-obstruction/pf4-bch-mapping.json"
MAPPING_PAYLOAD_SHA256 = "6e3fb46aaa1bbfc4691fa27e03644870e4409e3d1b789aa5df16731a41eaf54d"
CHECKED_LENGTHS = (4, 6, 8, 10)
SEARCH_ORDER = (2, 3, 4)
GAMMA = (
    Fraction(37, 900000),
    Fraction(313, 14400000),
    Fraction(29, 1800000),
)
SOURCE_PATHS = (
    "pyproject.toml",
    "requirements-reproducibility.txt",
    "scripts/__init__.py",
    "scripts/certify_tfim_gauge_witness.py",
    "scripts/reference_tfim_gauge_witness.py",
    "src/trottercert/__init__.py",
    "src/trottercert/algebra.py",
    "src/trottercert/cubic_field.py",
    "src/trottercert/intervals.py",
    "src/trottercert/pf4_bch_mapping.py",
    "src/trottercert/tfim_gauge_witness.py",
)

TOP_FIELDS = {
    "schema_version",
    "kind",
    "model",
    "pf4_mapping_binding",
    "gamma",
    "formulas",
    "search",
    "checked_lengths",
    "checked_instances",
    "claim",
    "implementation_sources",
    "implementation_sources_digest",
    "payload_sha256",
    "hpc_authorized",
}
MODEL_FIELDS = {
    "hamiltonian",
    "boundary",
    "normalized_trace",
    "checked_parameters",
}
BINDING_FIELDS = {
    "path",
    "file_sha256",
    "mapping_payload_sha256",
    "physical_bridge",
}
SEARCH_FIELDS = {
    "power_order",
    "witness_expression",
    "centering_conditions",
    "commutation_reason",
    "selection_rule",
    "normalization_status",
}
FORMULA_FIELDS = {"definitions", "l4_wraparound", "generic_even_l_ge_6"}
REGIME_FORMULA_FIELDS = {
    "domain",
    "tau_h2",
    "tau_h4",
    "tau_h_e5_over_gamma",
    "tau_h3_e5_over_gamma",
    "p3_hamiltonian_coefficient",
    "tau_w3_e5_over_gamma",
    "nonzero_condition",
}
INSTANCE_FIELDS = {
    "length",
    "regime",
    "parameters",
    "tau_h2",
    "tau_h4",
    "tau_h_e5",
    "tau_h3_e5",
    "orbit_polynomials",
    "witnesses",
    "first_nonzero_power",
}
ORBIT_FIELDS = {"tau_h2", "tau_h4", "tau_h_e5", "tau_h3_e5"}
WITNESS_FIELDS = {
    "power",
    "identity_coefficient",
    "hamiltonian_coefficient",
    "tau_w",
    "tau_wh",
    "commutator_is_zero",
    "tau_w_e5",
}
CLAIM_FIELDS = {
    "qualitative_status",
    "quotient",
    "scope",
    "proof_rule",
    "universal_even_l_ge_6_status",
    "finite_step_status",
    "quantitative_distance_status",
}

EXPECTED_FORMULAS = {
    "definitions": "x=h^2; y=j^2; gamma is bound PF4 cubic coefficient",
    "l4_wraparound": {
        "domain": "L=4 periodic",
        "tau_h2": "4*(x+y)",
        "tau_h4": "40*x^2+64*x*y+64*y^2",
        "tau_h_e5_over_gamma": "512*x*y*(x+8*y/3)",
        "tau_h3_e5_over_gamma": (
            "5120*x^3*y+(40960/3)*x^2*y^2+(65536/3)*x*y^3"
        ),
        "p3_hamiltonian_coefficient": "tau_h4/tau_h2",
        "tau_w3_e5_over_gamma": (
            "(1024/3)*h^4*j^4*(16*j^2-9*h^2)/(x+y)"
        ),
        "nonzero_condition": "h*j!=0 and 9*h^2!=16*j^2",
    },
    "generic_even_l_ge_6": {
        "domain": "periodic even L>=6",
        "tau_h2": "L*(x+y)",
        "tau_h4": "L*((3*L-2)*(x^2+y^2)+(6*L-8)*x*y)",
        "tau_h_e5_over_gamma": "128*L*x*y*(x+8*y/3)",
        "tau_h3_e5_over_gamma": (
            "128*L*x*y*((3*L-2)*x^2+((33*L-34)/3)*x*y+"
            "(8/3)*(3*L-2)*y^2)"
        ),
        "p3_hamiltonian_coefficient": "tau_h4/tau_h2",
        "tau_w3_e5_over_gamma": "(2560/3)*L*h^4*j^6/(x+y)",
        "nonzero_condition": "h*j!=0",
    },
}

EXPECTED_MODEL = {
    "hamiltonian": "H=h*sum_i X_i+j*sum_i Z_i Z_(i+1)",
    "boundary": "periodic_even_chain",
    "normalized_trace": "tau(A)=Tr(A)/2^L",
    "checked_parameters": {"h": [1, 1], "j": [1, 1]},
}
EXPECTED_SEARCH = {
    "power_order": [2, 3, 4],
    "witness_expression": "W_p=H^p-a_p*I-b_p*H",
    "centering_conditions": ["tau(W_p)=0", "tau(W_p*H)=0"],
    "commutation_reason": "W_p is an exact polynomial in H",
    "selection_rule": "first p in frozen order with tau(W_p*E5)!=0",
    "normalization_status": "qualitative_unnormalized_dual_pairing_only",
}
EXPECTED_CLAIM = {
    "qualitative_status": "certified_on_checked_instances",
    "quotient": "image(i ad_H)+span(I,H)",
    "scope": "periodic TFIM h=j=1 at L=4,6,8,10",
    "proof_rule": (
        "[W,H]=0 and tau(W)=tau(WH)=0 and tau(W E5)!=0 imply "
        "E5 not in image(i ad_H)+span(I,H)"
    ),
    "universal_even_l_ge_6_status": (
        "analytic_formula_recorded_induction_not_machine_proved"
    ),
    "finite_step_status": "not_claimed",
    "quantitative_distance_status": "not_claimed",
}


def canonical_bytes(payload: object) -> bytes:
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
        raise TypeError("gauge witness artifact must be an object")
    if raw != canonical_bytes(payload):
        raise ValueError("gauge witness artifact is not canonical JSON")
    return payload


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be an object")
    return value


def _require_fields(
    value: object,
    fields: set[str],
    label: str,
) -> Mapping[str, object]:
    mapping = _mapping(value, label)
    if set(mapping) != fields:
        raise ValueError(f"{label} field set mismatch")
    return mapping


def _reject_numeric_aliases(value: object, path: str = "artifact") -> None:
    boolean_field = path.endswith(
        (".hpc_authorized", ".commutator_is_zero")
    )
    if isinstance(value, float) or (
        isinstance(value, bool) and not boolean_field
    ):
        raise TypeError(f"{path}: forbidden JSON numeric alias")
    if isinstance(value, list):
        for index, item in enumerate(value):
            _reject_numeric_aliases(item, f"{path}[{index}]")
    elif isinstance(value, Mapping):
        for key, item in value.items():
            _reject_numeric_aliases(item, f"{path}.{key}")


def _pair(value: object, label: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
        or value[1] <= 0
    ):
        raise ValueError(f"{label} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if value != [result.numerator, result.denominator]:
        raise ValueError(f"{label} must be a canonical rational pair")
    return result


def _cubic(value: object, label: str) -> tuple[Fraction, Fraction, Fraction]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{label} must have three cubic coordinates")
    return tuple(_pair(item, label) for item in value)  # type: ignore[return-value]


def _cubic_scale(
    value: tuple[Fraction, Fraction, Fraction], scalar: Fraction
) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(coordinate * scalar for coordinate in value)  # type: ignore[return-value]


def _pair_json(value: Fraction | int) -> list[int]:
    exact = Fraction(value)
    return [exact.numerator, exact.denominator]


def _cubic_json(
    value: tuple[Fraction, Fraction, Fraction]
) -> list[list[int]]:
    return [_pair_json(coordinate) for coordinate in value]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _checked_path(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute() or Path(relative).as_posix() != relative:
        raise ValueError(f"source path is not canonical relative: {relative!r}")
    candidate = root / relative
    if candidate.is_symlink():
        raise ValueError(f"source path must not be a symlink: {relative}")
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(root.resolve(strict=True))
    except ValueError as error:
        raise ValueError(f"source path escapes issue root: {relative}") from error
    if not resolved.is_file():
        raise ValueError(f"source path is not a file: {relative}")
    return resolved


def _source_hashes(root: Path) -> dict[str, str]:
    return {relative: _sha256(_checked_path(root, relative)) for relative in SOURCE_PATHS}


def _canonical_digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _load_mapping_binding(root: Path, binding: Mapping[str, object]) -> None:
    if binding.get("path") != MAPPING_PATH:
        raise ValueError("PF4 mapping path mismatch")
    mapping_path = _checked_path(root, MAPPING_PATH)
    raw = mapping_path.read_bytes()
    try:
        wrapper = json.loads(raw, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid PF4 mapping JSON: {error}") from error
    if not isinstance(wrapper, dict) or raw != canonical_bytes(wrapper):
        raise ValueError("bound PF4 mapping artifact is not canonical")
    if binding.get("file_sha256") != hashlib.sha256(raw).hexdigest():
        raise ValueError("PF4 mapping file hash mismatch")
    if set(wrapper) != {
        "schema_version",
        "kind",
        "mapping",
        "implementation_sources",
        "reference_algorithm",
        "payload_sha256",
    }:
        raise ValueError("bound PF4 mapping wrapper field mismatch")
    if wrapper.get("schema_version") != 2 or wrapper.get("kind") != (
        "pf4_bch_mapping_source_closed_artifact"
    ):
        raise ValueError("bound PF4 mapping wrapper identity mismatch")
    if wrapper.get("payload_sha256") != payload_digest(wrapper):
        raise ValueError("bound PF4 mapping wrapper digest mismatch")
    mapping = _mapping(wrapper.get("mapping"), "bound PF4 mapping payload")
    observed_mapping_digest = mapping.get("payload_sha256")
    if observed_mapping_digest != MAPPING_PAYLOAD_SHA256:
        raise ValueError("authoritative PF4 mapping payload hash mismatch")
    if binding.get("mapping_payload_sha256") != observed_mapping_digest:
        raise ValueError("PF4 mapping payload binding mismatch")
    if mapping.get("payload_sha256") != payload_digest(mapping):
        raise ValueError("bound PF4 mathematical mapping digest mismatch")
    bridge = _mapping(mapping.get("physical_bridge"), "PF4 physical bridge")
    if bridge.get("relation") != "E5=L5" or binding.get("physical_bridge") != "E5=L5":
        raise ValueError("PF4 physical mapping bridge mismatch")
    if _cubic(mapping.get("optional_gamma"), "PF4 gamma") != GAMMA:
        raise ValueError("PF4 mapping gamma mismatch")


def _rational_terms(value: object, label: str) -> dict[tuple[int, int], Fraction]:
    if not isinstance(value, list):
        raise TypeError(f"{label} orbit terms must be a list")
    result: dict[tuple[int, int], Fraction] = {}
    previous: tuple[int, int] | None = None
    for record in value:
        fields = _require_fields(record, {"monomial", "coefficient"}, label)
        monomial = fields.get("monomial")
        if (
            not isinstance(monomial, list)
            or len(monomial) != 2
            or any(isinstance(item, bool) or not isinstance(item, int) or item < 0 for item in monomial)
        ):
            raise ValueError(f"{label} monomial is malformed")
        key = (monomial[0], monomial[1])
        if previous is not None and key <= previous:
            raise ValueError(f"{label} orbit monomials are duplicated or unsorted")
        previous = key
        result[key] = _pair(fields.get("coefficient"), label)
    return result


def _cubic_terms(
    value: object, label: str
) -> dict[tuple[int, int], tuple[Fraction, Fraction, Fraction]]:
    if not isinstance(value, list):
        raise TypeError(f"{label} orbit terms must be a list")
    result: dict[tuple[int, int], tuple[Fraction, Fraction, Fraction]] = {}
    previous: tuple[int, int] | None = None
    for record in value:
        fields = _require_fields(record, {"monomial", "coefficient"}, label)
        monomial = fields.get("monomial")
        if (
            not isinstance(monomial, list)
            or len(monomial) != 2
            or any(isinstance(item, bool) or not isinstance(item, int) or item < 0 for item in monomial)
        ):
            raise ValueError(f"{label} monomial is malformed")
        key = (monomial[0], monomial[1])
        if previous is not None and key <= previous:
            raise ValueError(f"{label} orbit monomials are duplicated or unsorted")
        previous = key
        result[key] = _cubic(fields.get("coefficient"), label)
    return result


def _expected_instance(length: int) -> dict[str, object]:
    tau_h2 = Fraction(2 * length)
    gamma = GAMMA
    tau_h_e5_over_gamma = Fraction(1408 * length, 3)
    expected_h2 = {(0, 2): Fraction(length), (2, 0): Fraction(length)}
    expected_h_e5 = {
        (2, 4): _cubic_scale(gamma, Fraction(1024 * length, 3)),
        (4, 2): _cubic_scale(gamma, Fraction(128 * length)),
    }
    if length == 4:
        regime = "l4_wraparound"
        tau_h4 = Fraction(168)
        tau_h3_e5_over_gamma = Fraction(121856, 3)
        tau_w3_e5_over_gamma = Fraction(3584, 3)
        expected_h4 = {
            (0, 4): Fraction(64),
            (2, 2): Fraction(64),
            (4, 0): Fraction(40),
        }
        expected_h3_e5 = {
            (2, 6): _cubic_scale(gamma, Fraction(65536, 3)),
            (4, 4): _cubic_scale(gamma, Fraction(40960, 3)),
            (6, 2): _cubic_scale(gamma, Fraction(5120)),
        }
    else:
        regime = "generic_even_l_ge_6"
        tau_h4 = Fraction(12 * length * (length - 1))
        tau_h3_e5_over_gamma = Fraction(256 * length * (33 * length - 28), 3)
        tau_w3_e5_over_gamma = Fraction(1280 * length, 3)
        expected_h4 = {
            (0, 4): Fraction(length * (3 * length - 2)),
            (2, 2): Fraction(length * (6 * length - 8)),
            (4, 0): Fraction(length * (3 * length - 2)),
        }
        expected_h3_e5 = {
            (2, 6): _cubic_scale(
                gamma, Fraction(1024 * length * (3 * length - 2), 3)
            ),
            (4, 4): _cubic_scale(
                gamma, Fraction(128 * length * (33 * length - 34), 3)
            ),
            (6, 2): _cubic_scale(
                gamma, Fraction(128 * length * (3 * length - 2))
            ),
        }
    b3 = tau_h4 / tau_h2
    zero = (Fraction(), Fraction(), Fraction())
    tau_h_e5 = _cubic_scale(gamma, tau_h_e5_over_gamma)
    tau_h3_e5 = _cubic_scale(gamma, tau_h3_e5_over_gamma)
    tau_w3_e5 = _cubic_scale(gamma, tau_w3_e5_over_gamma)
    witnesses = (
        (2, tau_h2, Fraction(), zero),
        (3, Fraction(), b3, tau_w3_e5),
        (4, tau_h4, Fraction(), zero),
    )
    return {
        "regime": regime,
        "tau_h2": tau_h2,
        "tau_h4": tau_h4,
        "tau_h_e5": tau_h_e5,
        "tau_h3_e5": tau_h3_e5,
        "orbit_h2": expected_h2,
        "orbit_h4": expected_h4,
        "orbit_h_e5": expected_h_e5,
        "orbit_h3_e5": expected_h3_e5,
        "witnesses": witnesses,
    }


def _verify_instance(instance: object, expected_length: int) -> None:
    record = _require_fields(instance, INSTANCE_FIELDS, "checked instance")
    if record.get("length") != expected_length:
        raise ValueError("checked length/order mismatch")
    expected = _expected_instance(expected_length)
    if record.get("regime") != expected["regime"]:
        raise ValueError("checked orbit regime mismatch")
    if record.get("parameters") != {"h": [1, 1], "j": [1, 1]}:
        raise ValueError("checked witness parameters mismatch")
    if _pair(record.get("tau_h2"), "tau_h2") != expected["tau_h2"]:
        raise ValueError("tau_h2 formula mismatch")
    if _pair(record.get("tau_h4"), "tau_h4") != expected["tau_h4"]:
        raise ValueError("tau_h4 formula mismatch")
    if _cubic(record.get("tau_h_e5"), "tau_h_e5") != expected["tau_h_e5"]:
        raise ValueError("tau(H E5) pairing mismatch")
    if _cubic(record.get("tau_h3_e5"), "tau_h3_e5") != expected["tau_h3_e5"]:
        raise ValueError("tau(H^3 E5) pairing mismatch")

    orbit = _require_fields(record.get("orbit_polynomials"), ORBIT_FIELDS, "orbit polynomials")
    if _rational_terms(orbit.get("tau_h2"), "tau_h2 orbit") != expected["orbit_h2"]:
        raise ValueError("tau_h2 orbit polynomial mismatch")
    if _rational_terms(orbit.get("tau_h4"), "tau_h4 orbit") != expected["orbit_h4"]:
        raise ValueError("tau_h4 orbit polynomial mismatch")
    if _cubic_terms(orbit.get("tau_h_e5"), "tau_h_e5 orbit") != expected["orbit_h_e5"]:
        raise ValueError("tau_h_e5 orbit polynomial mismatch")
    if _cubic_terms(orbit.get("tau_h3_e5"), "tau_h3_e5 orbit") != expected["orbit_h3_e5"]:
        raise ValueError("tau_h3_e5 orbit polynomial mismatch")

    witnesses = record.get("witnesses")
    if not isinstance(witnesses, list) or len(witnesses) != 3:
        raise ValueError("witness search record must contain p=2,3,4")
    for submitted, expected_values in zip(witnesses, expected["witnesses"]):
        witness = _require_fields(submitted, WITNESS_FIELDS, "witness")
        power, identity, h_coefficient, pairing = expected_values
        if witness.get("power") != power:
            raise ValueError("witness search order mismatch")
        if _pair(witness.get("identity_coefficient"), "witness identity") != identity:
            raise ValueError("witness identity coefficient mismatch")
        if _pair(witness.get("hamiltonian_coefficient"), "witness Hamiltonian") != h_coefficient:
            raise ValueError("witness Hamiltonian coefficient mismatch")
        if _pair(witness.get("tau_w"), "witness tau_w") != 0:
            raise ValueError("witness is not orthogonal to identity")
        if _pair(witness.get("tau_wh"), "witness tau_wh") != 0:
            raise ValueError("witness is not orthogonal to Hamiltonian")
        if witness.get("commutator_is_zero") is not True:
            raise ValueError("witness commutator field is not certified")
        if _cubic(witness.get("tau_w_e5"), "witness pairing") != pairing:
            raise ValueError("witness pairing mismatch")
    if record.get("first_nonzero_power") != 3:
        raise ValueError("first nonzero witness power mismatch")


def verify_payload(payload: object, root: Path = ISSUE_ROOT) -> None:
    artifact = _require_fields(payload, TOP_FIELDS, "gauge witness artifact")
    _reject_numeric_aliases(artifact)
    if artifact.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("gauge witness schema version mismatch")
    if artifact.get("kind") != ARTIFACT_KIND:
        raise ValueError("gauge witness artifact kind mismatch")
    if artifact.get("payload_sha256") != payload_digest(artifact):
        raise ValueError("gauge witness payload digest mismatch")
    if _require_fields(artifact.get("model"), MODEL_FIELDS, "model") != EXPECTED_MODEL:
        raise ValueError("model field mismatch")
    binding = _require_fields(
        artifact.get("pf4_mapping_binding"), BINDING_FIELDS, "PF4 mapping binding"
    )
    _load_mapping_binding(root, binding)
    if _cubic(artifact.get("gamma"), "gamma") != GAMMA:
        raise ValueError("gauge witness gamma mismatch")

    formulas = _require_fields(artifact.get("formulas"), FORMULA_FIELDS, "formulas")
    _require_fields(formulas.get("l4_wraparound"), REGIME_FORMULA_FIELDS, "L4 formula")
    _require_fields(formulas.get("generic_even_l_ge_6"), REGIME_FORMULA_FIELDS, "generic formula")
    if dict(formulas) != EXPECTED_FORMULAS:
        raise ValueError("closed formulas mismatch")
    search = _require_fields(artifact.get("search"), SEARCH_FIELDS, "search")
    if dict(search) != EXPECTED_SEARCH:
        raise ValueError("search schema/order mismatch")
    if artifact.get("checked_lengths") != list(CHECKED_LENGTHS):
        raise ValueError("checked length list is incomplete, duplicated, or unordered")
    instances = artifact.get("checked_instances")
    if not isinstance(instances, list) or len(instances) != len(CHECKED_LENGTHS):
        raise ValueError("checked instance set mismatch")
    for instance, length in zip(instances, CHECKED_LENGTHS):
        _verify_instance(instance, length)

    claim = _require_fields(artifact.get("claim"), CLAIM_FIELDS, "claim")
    if claim.get("universal_even_l_ge_6_status") != EXPECTED_CLAIM[
        "universal_even_l_ge_6_status"
    ]:
        raise ValueError("universal claim exceeds machine proof")
    if claim.get("finite_step_status") != "not_claimed":
        raise ValueError("finite-step claim exceeds artifact")
    if claim.get("quantitative_distance_status") != "not_claimed":
        raise ValueError("quantitative distance claim exceeds artifact")
    if dict(claim) != EXPECTED_CLAIM:
        raise ValueError("qualitative claim field mismatch")
    if artifact.get("hpc_authorized") is not False:
        raise ValueError("gauge witness artifact does not authorize HPC")

    sources = artifact.get("implementation_sources")
    if not isinstance(sources, Mapping) or tuple(sources) != SOURCE_PATHS:
        raise ValueError("implementation source closure/path order mismatch")
    expected_sources = _source_hashes(root)
    if dict(sources) != expected_sources:
        raise ValueError("implementation source hash mismatch")
    if artifact.get("implementation_sources_digest") != _canonical_digest(sources):
        raise ValueError("implementation source digest mismatch")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", required=True, type=Path, metavar="ARTIFACT")
    arguments = parser.parse_args(argv)
    try:
        verify_payload(load_payload(arguments.verify), ISSUE_ROOT)
    except (ArithmeticError, OSError, TypeError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 1
    print("TFIM PF4 gauge witness is canonical and independently verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
