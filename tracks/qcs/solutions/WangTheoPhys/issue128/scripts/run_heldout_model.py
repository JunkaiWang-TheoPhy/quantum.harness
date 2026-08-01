#!/usr/bin/env python3
"""Freeze an honest held-out TFIM evaluation against predeclared rules."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import subprocess
from decimal import ROUND_CEILING, Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path

import mpmath
import numpy as np
import scipy
from scipy.linalg import expm

from trottercert.algebra import PauliString, PauliSum, to_dense
from trottercert.hamiltonian import SymplecticPauli, tfim_terms
from trottercert.higher_order import (
    _multinomial,
    fourth_order_suzuki_stages,
    nested_commutator,
    weak_compositions,
)
from trottercert.intervals import RationalInterval
from trottercert.rigorous_fourth import fourth_order_suzuki_interval_stages

PAPER_ROOT = Path(__file__).resolve().parents[1]
FROZEN_LENGTHS = (4, 6, 8)
FROZEN_STEPS = (8, 12, 16, 24, 32, 48)
RULE_SCHEMA_PATHS = (
    "artifacts/publication/paper-a-file-ownership.json",
    "artifacts/publication/paper-a-claim-matrix.json",
)
SCIENTIFIC_DEPENDENCY_PATHS = (
    "scripts/run_heldout_model.py",
    "src/trottercert/algebra.py",
    "src/trottercert/hamiltonian.py",
    "src/trottercert/higher_order.py",
    "src/trottercert/intervals.py",
    "src/trottercert/lattice.py",
    "src/trottercert/local_commutators.py",
    "src/trottercert/rigorous_fourth.py",
)
MODEL_NORMALIZATION = "H=-J*sum_i X_i X_{i+1}-h*sum_i Z_i"
BOUND_SIGNIFICANT_DIGITS = 16
ACTUAL_DISPLAY_SIGNIFICANT_DIGITS = 9
THEOREM_CENTER = 6
CANDIDATE_COMPILER_SCOPE = "Heisenberg_four_matching_only"
CANDIDATE_STATUS = "unsupported"
CANDIDATE_STATUS_REASON = "candidate_compiler_model_family_unsupported"
ACTUAL_ERROR_STATUS = "dense_float64_diagnostic_not_a_certificate"
PUBLISHED_BOUND_STATUS = "certified"
_TOP_LEVEL_FIELDS = {
    "schema_version",
    "kind",
    "rules_frozen_at_commit",
    "frozen_rules",
    "field",
    "coupling",
    "periodic",
    "normalization",
    "requested_lengths",
    "requested_steps",
    "candidate_compiler_scope",
    "candidate_grouped_transfer_claimed",
    "actual_error_status",
    "certificates",
    "rows",
    "provenance",
    "provenance_digest",
}


def _canonical_json_bytes(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode(
        "ascii"
    )


def _canonical_digest(payload: object) -> str:
    return hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()


def _fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _strict_fraction_pair(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or type(value[0]) is not int
        or type(value[1]) is not int
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a strict rational pair")
    result = Fraction(value[0], value[1])
    if value != _fraction_pair(result):
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _outward_decimal(value: Fraction) -> str:
    if value < 0:
        raise ValueError("upper-bound serialization requires a nonnegative value")
    with localcontext() as context:
        context.prec = BOUND_SIGNIFICANT_DIGITS
        context.rounding = ROUND_CEILING
        rounded = Decimal(value.numerator) / Decimal(value.denominator)
        return format(rounded, f".{BOUND_SIGNIFICANT_DIGITS - 1}E").lower()


def _diagnostic_decimal(value: float) -> str:
    return format(value, f".{ACTUAL_DISPLAY_SIGNIFICANT_DIGITS - 1}e")


def _git_output(*arguments: str) -> bytes:
    try:
        result = subprocess.run(
            ["git", *arguments],
            cwd=PAPER_ROOT,
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("cannot resolve held-out freeze commit") from error
    return result.stdout


def _frozen_rule_record(rules_frozen_at_commit: str) -> dict[str, object]:
    if not rules_frozen_at_commit:
        raise ValueError("rules_frozen_at_commit is required")
    if not re.fullmatch(r"[0-9a-f]{40}", rules_frozen_at_commit):
        raise ValueError("rules_frozen_at_commit must be a full lowercase commit")
    resolved = (
        _git_output("rev-parse", f"{rules_frozen_at_commit}^{{commit}}")
        .decode("ascii")
        .strip()
    )
    if resolved != rules_frozen_at_commit:
        raise ValueError("rules_frozen_at_commit does not resolve exactly")
    repository_root = Path(
        _git_output("rev-parse", "--show-toplevel").decode().strip()
    ).resolve()

    inputs: dict[str, dict[str, str]] = {}
    for relative in RULE_SCHEMA_PATHS:
        repository_relative = (
            (PAPER_ROOT / relative).resolve().relative_to(repository_root).as_posix()
        )
        frozen = _git_output("show", f"{rules_frozen_at_commit}:{repository_relative}")
        inputs[relative] = {
            "sha256": hashlib.sha256(frozen).hexdigest(),
            "binding": "byte_exact_content_at_rules_frozen_commit",
        }
    unsigned: dict[str, object] = {
        "rules_frozen_at_commit": rules_frozen_at_commit,
        "inputs": inputs,
    }
    return {**unsigned, "rules_digest": _canonical_digest(unsigned)}


def _current_rule_bytes(relative: str) -> bytes:
    return (PAPER_ROOT / relative).read_bytes()


def _assert_current_rules_match(rules_frozen_at_commit: str) -> None:
    repository_root = Path(
        _git_output("rev-parse", "--show-toplevel").decode().strip()
    ).resolve()
    for relative in RULE_SCHEMA_PATHS:
        current_path = PAPER_ROOT / relative
        repository_relative = (
            current_path.resolve().relative_to(repository_root).as_posix()
        )
        frozen = _git_output("show", f"{rules_frozen_at_commit}:{repository_relative}")
        if _current_rule_bytes(relative) != frozen:
            raise ValueError(f"frozen rule mismatch for {relative}")


def _scientific_dependency_hashes() -> dict[str, str]:
    return {
        relative: hashlib.sha256((PAPER_ROOT / relative).read_bytes()).hexdigest()
        for relative in SCIENTIFIC_DEPENDENCY_PATHS
    }


def _runtime_environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "mpmath": mpmath.__version__,
    }


def _pauli_from_masks(key: SymplecticPauli) -> PauliString:
    x_mask, z_mask = key
    sites = x_mask | z_mask
    operators: dict[int, str] = {}
    while sites:
        bit = sites & -sites
        site = bit.bit_length() - 1
        has_x = bool(x_mask & bit)
        has_z = bool(z_mask & bit)
        operators[site] = "Y" if has_x and has_z else ("X" if has_x else "Z")
        sites ^= bit
    return PauliString(operators)


def _tfim_fragments(
    length: int,
    field: Fraction,
    coupling: Fraction,
    periodic: bool,
) -> tuple[PauliSum, PauliSum]:
    terms = tfim_terms(length, field, coupling, periodic)
    coupling_fragment = PauliSum.zero()
    field_fragment = PauliSum.zero()
    for key, coefficient in terms.items():
        term = PauliSum.term(_pauli_from_masks(key), coefficient)
        if key[0]:
            coupling_fragment += term
        else:
            field_fragment += term
    return coupling_fragment, field_fragment


@lru_cache(maxsize=32)
def _published_triangle_constant(
    length: int,
    field: Fraction,
    coupling: Fraction,
    periodic: bool,
    *,
    center: int = THEOREM_CENTER,
    decimal_digits: int = 18,
) -> tuple[Fraction, RationalInterval, int, int]:
    """Direct finite-chain published S4 bound with exact Pauli norms."""

    fragments = _tfim_fragments(length, field, coupling, periodic)
    stages_left, root = fourth_order_suzuki_interval_stages(
        len(fragments), decimal_digits=decimal_digits
    )
    stages = tuple(reversed(stages_left))
    if not 1 <= center <= len(stages):
        raise ValueError("theorem center is outside the merged stage sequence")
    order = 4
    theorem_terms = 0
    weights: dict[tuple[int, ...], Fraction] = {}

    def collect(
        j: int,
        indices: tuple[int, ...],
        composition: tuple[int, ...],
    ) -> None:
        nonlocal theorem_terms
        theorem_terms += 1
        outer: list[int] = []
        scalar = RationalInterval.point(_multinomial(order, composition))
        for stage_index, power in zip(indices, composition):
            stage = stages[stage_index - 1]
            scalar *= stage.coefficient**power
            outer.extend([stage.fragment_index] * power)
        for base_index in range(1, j):
            base_stage = stages[base_index - 1]
            key = tuple(outer) + (base_stage.fragment_index,)
            weights[key] = weights.get(key, Fraction()) + (
                scalar.abs_upper() * base_stage.coefficient.abs_upper()
            )

    for j in range(2, center + 1):
        indices = tuple(range(center, j - 1, -1))
        for composition in weak_compositions(order, len(indices)):
            if composition[-1]:
                collect(j, indices, composition)
    for j in range(center + 1, len(stages) + 1):
        indices = tuple(range(center + 1, j + 1))
        for composition in weak_compositions(order, len(indices)):
            if composition[-1]:
                collect(j, indices, composition)

    cache: dict[tuple[int, ...], PauliSum] = {}
    total = sum(
        (
            weight * nested_commutator(fragments, key, cache).exact_axis_l1()
            for key, weight in weights.items()
        ),
        Fraction(),
    )
    return total / factorial(order + 1), root, theorem_terms, len(weights)


def _certificate_record(
    length: int,
    field: Fraction,
    coupling: Fraction,
    periodic: bool,
) -> dict[str, object]:
    constant, root, theorem_terms, keys = _published_triangle_constant(
        length, field, coupling, periodic
    )
    terms = tfim_terms(length, field, coupling, periodic)
    unsigned: dict[str, object] = {
        "length": length,
        "field": _fraction_pair(field),
        "coupling": _fraction_pair(coupling),
        "periodic": periodic,
        "boundary": "periodic" if periodic else "open",
        "normalization": MODEL_NORMALIZATION,
        "canonical_terms": [
            [x_mask, z_mask, *_fraction_pair(coefficient)]
            for (x_mask, z_mask), coefficient in terms.items()
        ],
        "published_theorem": "high_order_commutator_triangle_bound",
        "norm_method": "direct_finite_chain_exact_Pauli_axis_l1",
        "theorem_center": THEOREM_CENTER,
        "coefficient_interval_decimal_digits": 18,
        "root_interval": {
            "lower": _fraction_pair(root.lower),
            "upper": _fraction_pair(root.upper),
        },
        "constant_upper_exact": _fraction_pair(constant),
        "constant_upper": _outward_decimal(constant),
        "theorem_terms": theorem_terms,
        "expanded_commutator_keys": keys,
    }
    return {**unsigned, "identity_digest": _canonical_digest(unsigned)}


def _dense_diagnostics(
    length: int,
    field: Fraction,
    coupling: Fraction,
    periodic: bool,
    steps: tuple[int, ...],
) -> dict[int, float]:
    fragments = _tfim_fragments(length, field, coupling, periodic)
    fragment_matrices = tuple(to_dense(fragment, length) for fragment in fragments)
    hamiltonian = sum(fragment_matrices, np.zeros_like(fragment_matrices[0]))
    exact = expm(-1j * hamiltonian)
    result: dict[int, float] = {}
    for step_count in steps:
        one_step = np.eye(1 << length, dtype=np.complex128)
        exponential_cache: dict[tuple[int, str], np.ndarray] = {}
        for stage in fourth_order_suzuki_stages(len(fragment_matrices)):
            key = stage.fragment_index, str(stage.coefficient)
            exponential = exponential_cache.get(key)
            if exponential is None:
                exponential = expm(
                    -1j
                    * float(stage.coefficient)
                    * fragment_matrices[stage.fragment_index]
                    / step_count
                )
                exponential_cache[key] = exponential
            one_step = one_step @ exponential
        approximate = np.linalg.matrix_power(one_step, step_count)
        error = float(np.linalg.norm(exact - approximate, ord=2))
        if not math.isfinite(error):
            raise ArithmeticError("TFIM dense diagnostic is not finite")
        result[step_count] = error
    return result


def _heldout_row(
    *,
    certificate: dict[str, object],
    field: Fraction,
    coupling: Fraction,
    periodic: bool,
    step_count: int,
    actual: float,
) -> dict[str, object]:
    constant = Fraction(*certificate["constant_upper_exact"])
    theorem_upper = constant / step_count**4
    published_upper = min(Fraction(2), theorem_upper)
    if actual > float(published_upper):
        raise ArithmeticError("published TFIM bound does not dominate diagnostic")
    return {
        "model": "transverse_field_ising_chain",
        "length": certificate["length"],
        "steps": step_count,
        "field": _fraction_pair(field),
        "coupling": _fraction_pair(coupling),
        "periodic": periodic,
        "boundary": "periodic" if periodic else "open",
        "normalization": MODEL_NORMALIZATION,
        "status": CANDIDATE_STATUS,
        "status_reason": CANDIDATE_STATUS_REASON,
        "status_explanation": (
            "the candidate certificate compiler is restricted to "
            "the Heisenberg four-matching model family"
        ),
        "candidate_certificate_upper": None,
        "actual_operator_error": _diagnostic_decimal(actual),
        "actual_operator_error_status": ACTUAL_ERROR_STATUS,
        "published_upper": _outward_decimal(published_upper),
        "published_upper_exact": _fraction_pair(published_upper),
        "published_bound_status": PUBLISHED_BOUND_STATUS,
        "published_bound_method": (
            "direct_finite_chain_published_triangle_exact_Pauli_l1_"
            "capped_by_unitary_distance"
        ),
        "published_certificate_digest": certificate["identity_digest"],
    }


def run_heldout(
    *,
    field: Fraction,
    coupling: Fraction,
    rules_frozen_at_commit: str,
    periodic: bool = False,
    lengths: tuple[int, ...] = FROZEN_LENGTHS,
    steps: tuple[int, ...] = FROZEN_STEPS,
) -> dict[str, object]:
    """Run the preregistered grid without promoting diagnostics to proofs."""

    frozen_rules = _frozen_rule_record(rules_frozen_at_commit)
    _assert_current_rules_match(rules_frozen_at_commit)
    if not isinstance(field, Fraction) or not isinstance(coupling, Fraction):
        raise TypeError("held-out field and coupling must be exact Fractions")
    if not isinstance(periodic, bool):
        raise TypeError("held-out periodic flag must be bool")
    normalized_lengths = tuple(sorted(set(lengths)))
    normalized_steps = tuple(sorted(set(steps)))
    if not normalized_lengths or not normalized_steps:
        raise ValueError("held-out lengths and steps must be nonempty")
    if any(
        isinstance(value, bool) or not isinstance(value, int)
        for value in normalized_lengths
    ):
        raise TypeError("held-out lengths must be integers")
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value < 1
        for value in normalized_steps
    ):
        raise ValueError("held-out steps must be positive integers")

    certificates: list[dict[str, object]] = []
    rows: list[dict[str, object]] = []
    for length in normalized_lengths:
        certificate = _certificate_record(length, field, coupling, periodic)
        certificates.append(certificate)
        diagnostics = _dense_diagnostics(
            length, field, coupling, periodic, normalized_steps
        )
        for step_count in normalized_steps:
            rows.append(
                _heldout_row(
                    certificate=certificate,
                    field=field,
                    coupling=coupling,
                    periodic=periodic,
                    step_count=step_count,
                    actual=diagnostics[step_count],
                )
            )

    provenance = {
        "scientific_dependency_sha256": _scientific_dependency_hashes(),
        "runtime_environment": _runtime_environment(),
    }
    payload = {
        "schema_version": 1,
        "kind": "paper_a_heldout_tfim_evaluation",
        "rules_frozen_at_commit": rules_frozen_at_commit,
        "frozen_rules": frozen_rules,
        "field": _fraction_pair(field),
        "coupling": _fraction_pair(coupling),
        "periodic": periodic,
        "normalization": MODEL_NORMALIZATION,
        "requested_lengths": list(normalized_lengths),
        "requested_steps": list(normalized_steps),
        "candidate_compiler_scope": CANDIDATE_COMPILER_SCOPE,
        "candidate_grouped_transfer_claimed": False,
        "actual_error_status": ACTUAL_ERROR_STATUS,
        "certificates": certificates,
        "rows": rows,
        "provenance": provenance,
        "provenance_digest": _canonical_digest(provenance),
    }
    verify_heldout_payload(payload)
    return payload


def verify_heldout_payload(payload: dict[str, object]) -> bool:
    """Fail closed on rule, provenance, certificate, and row mutations."""

    if set(payload) != _TOP_LEVEL_FIELDS:
        raise ValueError("held-out top-level schema mismatch")
    if type(payload.get("schema_version")) is not int or payload["schema_version"] != 1:
        raise ValueError("unsupported held-out schema")
    if payload.get("kind") != "paper_a_heldout_tfim_evaluation":
        raise ValueError("held-out kind mismatch")
    if payload.get("normalization") != MODEL_NORMALIZATION:
        raise ValueError("held-out normalization mismatch")
    if payload.get("candidate_compiler_scope") != CANDIDATE_COMPILER_SCOPE:
        raise ValueError("held-out candidate compiler scope mismatch")
    if payload.get("candidate_grouped_transfer_claimed") is not False:
        raise ValueError("held-out grouped transfer claim must remain false")
    if payload.get("actual_error_status") != ACTUAL_ERROR_STATUS:
        raise ValueError("held-out actual error status mismatch")
    commit = payload.get("rules_frozen_at_commit")
    if not isinstance(commit, str):
        raise TypeError("held-out freeze commit is malformed")
    if payload.get("frozen_rules") != _frozen_rule_record(commit):
        raise ValueError("held-out frozen-rule binding mismatch")
    provenance = payload.get("provenance")
    if not isinstance(provenance, dict):
        raise TypeError("held-out provenance is malformed")
    if set(provenance) != {
        "scientific_dependency_sha256",
        "runtime_environment",
    }:
        raise ValueError("held-out provenance schema mismatch")
    if provenance.get("scientific_dependency_sha256") != (
        _scientific_dependency_hashes()
    ):
        raise ValueError("held-out scientific dependency mismatch")
    if provenance.get("runtime_environment") != _runtime_environment():
        raise ValueError("held-out runtime environment mismatch")
    if payload.get("provenance_digest") != _canonical_digest(provenance):
        raise ValueError("held-out provenance digest mismatch")

    field = _strict_fraction_pair(payload["field"], "held-out field")
    coupling = _strict_fraction_pair(payload["coupling"], "held-out coupling")
    periodic = payload.get("periodic")
    if not isinstance(periodic, bool):
        raise TypeError("held-out periodic flag is malformed")
    lengths = tuple(payload.get("requested_lengths", ()))
    steps = tuple(payload.get("requested_steps", ()))
    if lengths != tuple(sorted(set(lengths))) or any(
        isinstance(length, bool) or not isinstance(length, int) or length < 2
        for length in lengths
    ):
        raise ValueError("held-out lengths are not canonical")
    if steps != tuple(sorted(set(steps))) or any(
        isinstance(step, bool) or not isinstance(step, int) or step < 1
        for step in steps
    ):
        raise ValueError("held-out steps are not canonical")
    certificates = payload.get("certificates")
    rows = payload.get("rows")
    if not isinstance(certificates, list) or not isinstance(rows, list):
        raise TypeError("held-out certificates and rows must be lists")
    if len(rows) != len(lengths) * len(steps):
        raise ValueError("held-out grid is incomplete")

    rebuilt_certificates = [
        _certificate_record(length, field, coupling, periodic) for length in lengths
    ]
    if _canonical_json_bytes(certificates) != _canonical_json_bytes(
        rebuilt_certificates
    ):
        raise ValueError("held-out certificate reconstruction mismatch")
    by_digest: dict[str, dict[str, object]] = {}
    for record in rebuilt_certificates:
        if not isinstance(record, dict):
            raise TypeError("held-out certificate is malformed")
        digest = record.get("identity_digest")
        unsigned = {
            key: value for key, value in record.items() if key != "identity_digest"
        }
        if not isinstance(digest, str) or digest != _canonical_digest(unsigned):
            raise ValueError("held-out certificate digest mismatch")
        by_digest[digest] = record

    expected_grid = sorted((length, step) for length in lengths for step in steps)
    diagnostics = {
        length: _dense_diagnostics(length, field, coupling, periodic, steps)
        for length in lengths
    }
    observed_grid: list[tuple[int, int]] = []
    for row in rows:
        if not isinstance(row, dict):
            raise TypeError("held-out row is malformed")
        digest = row.get("published_certificate_digest")
        if not isinstance(digest, str) or digest not in by_digest:
            raise ValueError("held-out row has an unknown certificate")
        record = by_digest[digest]
        length_value = row.get("length")
        step_value = row.get("steps")
        if type(length_value) is not int or type(step_value) is not int:
            raise TypeError("held-out row length and steps must be strict integers")
        length = length_value
        step = step_value
        if length != record["length"]:
            raise ValueError("held-out row references the wrong length certificate")
        actual_value = row.get("actual_operator_error")
        if not isinstance(actual_value, str):
            raise TypeError("held-out actual diagnostic must be a decimal string")
        recomputed_actual = diagnostics[length][step]
        if actual_value != _diagnostic_decimal(recomputed_actual):
            raise ValueError("held-out actual diagnostic mismatch")
        expected_row = _heldout_row(
            certificate=record,
            field=field,
            coupling=coupling,
            periodic=periodic,
            step_count=step,
            actual=recomputed_actual,
        )
        comparable_row = dict(row)
        comparable_expected = dict(expected_row)
        comparable_row.pop("actual_operator_error", None)
        comparable_expected.pop("actual_operator_error", None)
        if _canonical_json_bytes(comparable_row) != _canonical_json_bytes(
            comparable_expected
        ):
            raise ValueError("held-out row canonical reconstruction mismatch")
        observed_grid.append((length, step))
    if observed_grid != expected_grid:
        raise ValueError("held-out rows are not canonical and complete")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules-frozen-at-commit", required=True)
    parser.add_argument("--field", type=Fraction, default=Fraction(1))
    parser.add_argument("--coupling", type=Fraction, default=Fraction(1))
    parser.add_argument("--periodic", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=PAPER_ROOT / "benchmarks/paper-a/heldout-tfim.json",
    )
    arguments = parser.parse_args()
    payload = run_heldout(
        field=arguments.field,
        coupling=arguments.coupling,
        rules_frozen_at_commit=arguments.rules_frozen_at_commit,
        periodic=arguments.periodic,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_bytes(_canonical_json_bytes(payload))
    print(f"output={arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
