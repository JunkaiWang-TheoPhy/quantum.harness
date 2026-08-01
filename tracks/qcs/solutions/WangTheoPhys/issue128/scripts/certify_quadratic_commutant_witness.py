#!/usr/bin/env python3
"""Generate or fully recompute the exact quadratic commutant certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from trottercert.commutant_witness import (
    build_quadratic_witness_payload,
    heisenberg_symplectic_terms,
    lift_cubic_density,
    pair_lifted_cubic_density,
    rational_pauli_pairing,
    square_real_pauli_terms,
    verify_quadratic_witness_payload,
)
from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.cubic_local import exact_log_e5_density, exact_matching_density
from trottercert.lattice import SquareLattice

ISSUE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ISSUE_ROOT
    / "docs/experiments/processor-obstruction/quadratic-commutant-witness.json"
)
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/cubic_local.py",
    ISSUE_ROOT / "src/trottercert/commutant_witness.py",
)


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(ISSUE_ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in SOURCE_PATHS
    }


def compute_payload(length: int = 12) -> dict[str, Any]:
    lattice = SquareLattice(length)
    print(f"building exact H and H^2 for periodic {length}x{length}", flush=True)
    hamiltonian = heisenberg_symplectic_terms(lattice)
    squared = square_real_pauli_terms(hamiltonian)
    tau_h2 = squared.get((0, 0))
    if tau_h2 is None:
        raise ArithmeticError("H^2 has no identity coefficient")
    tau_h3 = rational_pauli_pairing(hamiltonian, squared)

    reconstructed = {}
    for color in range(4):
        registry, density = exact_matching_density(color)
        for pauli, coefficient in lift_cubic_density(
            registry, density, lattice
        ).items():
            if coefficient.a1 or coefficient.a2:
                raise ArithmeticError("matching density is not rational")
            reconstructed[pauli] = reconstructed.get(pauli, 0) + coefficient.a0
    if reconstructed != hamiltonian:
        raise ArithmeticError("matching densities do not reconstruct H")

    print("generating exact fifth-degree logarithm density", flush=True)
    registry, e5 = exact_log_e5_density(fourth_order_suzuki_cubic_stages(4))
    print(
        f"pairing {len(e5)} density terms through {(length // 2) ** 2} translations",
        flush=True,
    )
    tau_h_e5 = pair_lifted_cubic_density(
        registry, e5, hamiltonian, lattice
    )
    tau_h2_e5 = pair_lifted_cubic_density(registry, e5, squared, lattice)
    payload = build_quadratic_witness_payload(
        length=length,
        tau_h2=tau_h2,
        tau_h3=tau_h3,
        tau_h_e5=tau_h_e5,
        tau_h2_e5=tau_h2_e5,
        matching_reconstruction_terms=len(reconstructed),
    )
    payload["implementation_sources"] = _source_hashes()
    verify_quadratic_witness_payload(payload)
    return payload


def _encoded(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=12)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()

    expected = compute_payload(args.length)
    if args.verify:
        try:
            submitted = json.loads(args.output.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise SystemExit(f"cannot read witness artifact: {args.output}") from exc
        verify_quadratic_witness_payload(submitted)
        if submitted != expected:
            raise SystemExit("quadratic commutant witness artifact mismatch")
        print(
            "quadratic commutant witness valid: "
            f"tau_w_e5={submitted['pairings']['tau_w_e5']}",
            flush=True,
        )
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(_encoded(expected))
    print(f"wrote quadratic commutant witness: {args.output}", flush=True)


if __name__ == "__main__":
    main()
