#!/usr/bin/env python3
"""Certify the extensive quadratic commutant witness on periodic tori."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

from trottercert.commutant_witness import (
    QuadraticWitnessMoments,
    density_aliases_on_torus,
    heisenberg_symplectic_terms,
    lift_cubic_density,
    pair_lifted_cubic_density,
    quadratic_witness_moments,
    square_real_pauli_terms,
)
from trottercert.cubic_field import Cubic, fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import exact_log_e5_density, exact_matching_density
from trottercert.lattice import SquareLattice
from trottercert.local_commutators import CoordinateRegistry, SymplecticPauli

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/extensive-commutant-witness.json"
)
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/cubic_local.py",
    ISSUE_ROOT / "src/trottercert/commutant_witness.py",
)


@dataclass(frozen=True, slots=True)
class ExtensiveSizeRecord:
    length: int
    n_sites: int
    cells: int
    matching_reconstruction_terms: int
    moments: QuadraticWitnessMoments


def _rational_json(value: Fraction) -> list[int]:
    value = Fraction(value)
    return [value.numerator, value.denominator]


def _cubic_json(value: Cubic) -> list[list[int]]:
    return [_rational_json(value.a0), _rational_json(value.a1), _rational_json(value.a2)]


def _parse_rational(value: object, field: str) -> Fraction:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(not isinstance(entry, int) or isinstance(entry, bool) for entry in value)
        or value[1] <= 0
    ):
        raise ValueError(f"{field} must be a canonical rational pair")
    result = Fraction(value[0], value[1])
    if _rational_json(result) != value:
        raise ValueError(f"{field} must be a canonical rational pair")
    return result


def _parse_cubic(value: object, field: str) -> Cubic:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field} must contain three exact coordinates")
    return Cubic(*(_parse_rational(entry, field) for entry in value))


def _moments_json(moments: QuadraticWitnessMoments) -> dict[str, object]:
    return {
        "tau_h2": _rational_json(moments.tau_h2),
        "tau_h3": _rational_json(moments.tau_h3),
        "tau_h4": _rational_json(moments.tau_h4),
        "identity_coefficient": _rational_json(moments.identity_coefficient),
        "hamiltonian_coefficient": _rational_json(moments.h_coefficient),
        "tau_w2": _rational_json(moments.tau_w2),
        "tau_h_e5": _cubic_json(moments.tau_h_e5),
        "tau_h2_e5": _cubic_json(moments.tau_h2_e5),
        "tau_w_e5": _cubic_json(moments.tau_w_e5),
        "squared_normalized_pairing": _cubic_json(
            moments.squared_normalized_pairing
        ),
    }


def _record_json(record: ExtensiveSizeRecord) -> dict[str, object]:
    return {
        "length": record.length,
        "n_sites": record.n_sites,
        "cells": record.cells,
        "matching_reconstruction_terms": record.matching_reconstruction_terms,
        "moments": _moments_json(record.moments),
        "dual_remainder_threshold_squared": {
            str(steps): _cubic_json(
                record.moments.squared_normalized_pairing / steps**8
            )
            for steps in (78, 95)
        },
    }


def _canonical_digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(ISSUE_ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in SOURCE_PATHS
    }


def build_extensive_payload(
    records: Sequence[ExtensiveSizeRecord],
    *,
    rejected_alias_lengths: Sequence[int],
    implementation_sources: Mapping[str, str],
) -> dict[str, Any]:
    if len(records) < 2:
        raise ValueError("at least two alias-free sizes are required")
    ordered = tuple(sorted(records, key=lambda record: record.length))
    if len({record.length for record in ordered}) != len(ordered):
        raise ValueError("accepted lengths must be unique")

    per_cell_h = ordered[0].moments.tau_h_e5 / ordered[0].cells
    per_cell_h2 = ordered[0].moments.tau_h2_e5 / ordered[0].cells
    per_cell_w = ordered[0].moments.tau_w_e5 / ordered[0].cells
    for record in ordered:
        if record.moments.tau_h_e5 / record.cells != per_cell_h:
            raise ValueError("stable tau(H E5) per-cell relation failed")
        if record.moments.tau_h2_e5 / record.cells != per_cell_h2:
            raise ValueError("stable tau(H^2 E5) per-cell relation failed")
        if record.moments.tau_w_e5 / record.cells != per_cell_w:
            raise ValueError("stable tau(W E5) per-cell relation failed")

    sources = dict(sorted(implementation_sources.items()))
    payload: dict[str, Any] = {
        "schema_version": 1,
        "kind": "issue128_extensive_commutant_witness",
        "model_family": "periodic_even_square_spin_half_heisenberg",
        "formula": "five_copy_fourth_order_suzuki",
        "accepted_lengths": [record.length for record in ordered],
        "rejected_alias_lengths": sorted(set(rejected_alias_lengths)),
        "records": [_record_json(record) for record in ordered],
        "stable_relations": {
            "scope": "exactly_verified_on_all_accepted_sizes",
            "tau_h2": "3 N / 8",
            "tau_h3": "-3 N / 16",
            "tau_h4": "(27 N^2 - 21 N) / 64",
            "tau_w2": "9 N (2 N - 3) / 64",
            "hamiltonian_coefficient": [1, 2],
            "tau_h_e5_per_cell": _cubic_json(per_cell_h),
            "tau_h2_e5_per_cell": _cubic_json(per_cell_h2),
            "tau_w_e5_per_cell": _cubic_json(per_cell_w),
        },
        "finite_step_gate": {
            "criterion": "tau(W R_r)^2 < tau(W E5)^2 / r^8",
            "sufficient_hs_criterion": (
                "tau(R_r^2) < squared_normalized_pairing / r^8"
            ),
            "candidate_steps": [78, 95],
        },
        "claim": {
            "leading_order_status": "no_go",
            "finite_step_status": "inconclusive",
            "reason": "local-log branch and compatible all-order remainder missing",
        },
        "implementation_sources": sources,
        "implementation_sources_digest": _canonical_digest(sources),
        "hpc_authorized": False,
    }
    verify_extensive_payload(payload)
    return payload


def _verify_record(record: Mapping[str, object]) -> tuple[int, Cubic, Cubic, Cubic]:
    length = record.get("length")
    if not isinstance(length, int) or isinstance(length, bool) or length < 6 or length % 2:
        raise ValueError("accepted length must be an even integer at least six")
    n_sites = length * length
    if record.get("n_sites") != n_sites or record.get("cells") != n_sites // 4:
        raise ValueError("length and volume metadata mismatch")
    if record.get("matching_reconstruction_terms") != 6 * n_sites:
        raise ValueError("Hamiltonian reconstruction term count mismatch")
    moments = record.get("moments")
    if not isinstance(moments, Mapping):
        raise ValueError("moment record is missing")
    tau_h2 = _parse_rational(moments.get("tau_h2"), "moment tau_h2")
    tau_h3 = _parse_rational(moments.get("tau_h3"), "moment tau_h3")
    tau_h4 = _parse_rational(moments.get("tau_h4"), "moment tau_h4")
    identity = _parse_rational(
        moments.get("identity_coefficient"), "moment identity coefficient"
    )
    h_coefficient = _parse_rational(
        moments.get("hamiltonian_coefficient"), "moment Hamiltonian coefficient"
    )
    tau_w2 = _parse_rational(moments.get("tau_w2"), "moment tau_w2")
    expected_h2 = Fraction(3 * n_sites, 8)
    expected_h3 = Fraction(-3 * n_sites, 16)
    expected_h4 = Fraction(27 * n_sites**2 - 21 * n_sites, 64)
    expected_w2 = Fraction(9 * n_sites * (2 * n_sites - 3), 64)
    if (tau_h2, tau_h3, tau_h4, tau_w2) != (
        expected_h2,
        expected_h3,
        expected_h4,
        expected_w2,
    ):
        raise ValueError("moment family relation mismatch")
    if identity != -tau_h2 or h_coefficient != Fraction(1, 2):
        raise ValueError("moment witness coefficient mismatch")

    tau_h_e5 = _parse_cubic(moments.get("tau_h_e5"), "moment tau_h_e5")
    tau_h2_e5 = _parse_cubic(moments.get("tau_h2_e5"), "moment tau_h2_e5")
    tau_w_e5 = _parse_cubic(moments.get("tau_w_e5"), "moment tau_w_e5")
    if tau_w_e5 != tau_h2_e5 + h_coefficient * tau_h_e5:
        raise ValueError("moment witness pairing mismatch")
    rho = _parse_cubic(
        moments.get("squared_normalized_pairing"),
        "moment squared normalized pairing",
    )
    if rho != tau_w_e5**2 / tau_w2:
        raise ValueError("moment normalized pairing mismatch")
    thresholds = record.get("dual_remainder_threshold_squared")
    if not isinstance(thresholds, Mapping):
        raise ValueError("finite-step threshold record is missing")
    for steps in (78, 95):
        observed = _parse_cubic(thresholds.get(str(steps)), "finite-step threshold")
        if observed != rho / steps**8:
            raise ValueError("finite-step threshold mismatch")
    return n_sites // 4, tau_h_e5, tau_h2_e5, tau_w_e5


def verify_extensive_payload(payload: Mapping[str, object]) -> None:
    if payload.get("kind") != "issue128_extensive_commutant_witness":
        raise ValueError("extensive witness kind mismatch")
    records = payload.get("records")
    if not isinstance(records, list) or len(records) < 2:
        raise ValueError("at least two extensive records are required")
    accepted = payload.get("accepted_lengths")
    record_lengths = [record.get("length") for record in records if isinstance(record, Mapping)]
    if accepted != record_lengths:
        raise ValueError("accepted length list mismatch")
    rejected = payload.get("rejected_alias_lengths")
    if not isinstance(rejected, list) or any(
        not isinstance(length, int) or length % 2 or length < 4 for length in rejected
    ):
        raise ValueError("rejected alias length list is invalid")
    if set(rejected) & set(record_lengths):
        raise ValueError("a length cannot be both accepted and rejected")

    verified = [_verify_record(record) for record in records]
    cells0, h0, h20, w0 = verified[0]
    expected_stable = {
        "scope": "exactly_verified_on_all_accepted_sizes",
        "tau_h2": "3 N / 8",
        "tau_h3": "-3 N / 16",
        "tau_h4": "(27 N^2 - 21 N) / 64",
        "tau_w2": "9 N (2 N - 3) / 64",
        "hamiltonian_coefficient": [1, 2],
        "tau_h_e5_per_cell": _cubic_json(h0 / cells0),
        "tau_h2_e5_per_cell": _cubic_json(h20 / cells0),
        "tau_w_e5_per_cell": _cubic_json(w0 / cells0),
    }
    for cells, h_pairing, h2_pairing, w_pairing in verified[1:]:
        if (h_pairing / cells, h2_pairing / cells, w_pairing / cells) != (
            h0 / cells0,
            h20 / cells0,
            w0 / cells0,
        ):
            raise ValueError("stable per-cell pairing relation mismatch")
    if payload.get("stable_relations") != expected_stable:
        raise ValueError("stable relation payload mismatch")

    sources = payload.get("implementation_sources")
    if not isinstance(sources, Mapping) or not sources:
        raise ValueError("implementation source manifest is missing")
    if any(
        not isinstance(path, str)
        or not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
        for path, digest in sources.items()
    ):
        raise ValueError("implementation source manifest is invalid")
    if payload.get("implementation_sources_digest") != _canonical_digest(sources):
        raise ValueError("implementation source digest mismatch")
    claim = payload.get("claim")
    if not isinstance(claim, Mapping) or claim.get("leading_order_status") != "no_go":
        raise ValueError("leading-order claim mismatch")
    if claim.get("finite_step_status") != "inconclusive":
        raise ValueError("finite-step claim exceeds the certificate")
    if payload.get("hpc_authorized") is not False:
        raise ValueError("extensive certificate does not authorize HPC")


def compute_size_record(
    registry: CoordinateRegistry,
    e5: Mapping[SymplecticPauli, Cubic],
    length: int,
) -> ExtensiveSizeRecord:
    if density_aliases_on_torus(registry, e5, length):
        raise ValueError(f"E5 density aliases on the periodic L={length} torus")
    lattice = SquareLattice(length)
    hamiltonian = heisenberg_symplectic_terms(lattice)
    reconstructed = {}
    for color in range(4):
        matching_registry, density = exact_matching_density(color)
        lifted = lift_cubic_density(matching_registry, density, lattice)
        for pauli, coefficient in lifted.items():
            if coefficient.a1 or coefficient.a2:
                raise ArithmeticError("matching density is not rational")
            reconstructed[pauli] = reconstructed.get(pauli, Fraction()) + coefficient.a0
    if reconstructed != hamiltonian:
        raise ArithmeticError(f"matching reconstruction failed for L={length}")
    squared = square_real_pauli_terms(hamiltonian)
    moments = quadratic_witness_moments(
        hamiltonian,
        squared,
        pair_lifted_cubic_density(registry, e5, hamiltonian, lattice),
        pair_lifted_cubic_density(registry, e5, squared, lattice),
    )
    return ExtensiveSizeRecord(
        length=length,
        n_sites=lattice.n_sites,
        cells=lattice.n_sites // 4,
        matching_reconstruction_terms=len(reconstructed),
        moments=moments,
    )


def compute_payload(lengths: Sequence[int]) -> dict[str, Any]:
    requested = tuple(sorted(set(lengths)))
    if len(requested) < 2:
        raise ValueError("request at least two even torus lengths")
    print("generating exact fifth-degree logarithm density once", flush=True)
    registry, e5 = exact_log_e5_density(fourth_order_suzuki_cubic_stages(4))
    records: list[ExtensiveSizeRecord] = []
    rejected: list[int] = []
    for length in requested:
        if density_aliases_on_torus(registry, e5, length):
            print(f"rejecting L={length}: periodic coordinate alias", flush=True)
            rejected.append(length)
            continue
        print(f"computing exact extensive witness record for L={length}", flush=True)
        records.append(compute_size_record(registry, e5, length))
    return build_extensive_payload(
        records,
        rejected_alias_lengths=rejected,
        implementation_sources=_source_hashes(),
    )


def _encoded(payload: Mapping[str, object]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", type=int, nargs="+", default=(4, 6, 8, 10, 12))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--verify-full", action="store_true")
    args = parser.parse_args()

    if args.verify and args.verify_full:
        raise SystemExit("choose at most one verification mode")
    if args.verify:
        try:
            payload = json.loads(args.output.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"cannot read extensive witness: {args.output}") from exc
        verify_extensive_payload(payload)
        if payload.get("implementation_sources") != _source_hashes():
            raise SystemExit("extensive witness source hashes do not match")
        print("extensive commutant witness is internally valid", flush=True)
        return

    expected = compute_payload(args.lengths)
    if args.verify_full:
        try:
            submitted = json.loads(args.output.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"cannot read extensive witness: {args.output}") from exc
        verify_extensive_payload(submitted)
        if submitted != expected:
            raise SystemExit("extensive commutant witness artifact mismatch")
        print("extensive commutant witness fully recomputed and valid", flush=True)
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(_encoded(expected))
    print(f"wrote extensive commutant witness: {args.output}", flush=True)


if __name__ == "__main__":
    main()
