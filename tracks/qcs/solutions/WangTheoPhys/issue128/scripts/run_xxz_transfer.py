"""Run the preregistered exact XXZ transfer-capability sweep."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections.abc import Callable
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from trottercert.hamiltonian import xxz_bond

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"
REFERENCE_VERIFIER = ROOT / "scripts" / "reference_verify.py"
PREREGISTERED_DELTAS = (
    Fraction(0),
    Fraction(1, 4),
    Fraction(1, 2),
    Fraction(1),
    Fraction(3, 2),
    Fraction(2),
    Fraction(4),
)
STATUSES = {"certified", "unsupported", "inconclusive"}
CERTIFICATION_ONLY_FIELDS = (
    "d4_term_count",
    "d4_group_count",
    "d4_norm_bound",
    "d5_term_count",
    "d5_group_count",
    "d5_norm_bound",
    "accepted_steps",
    "source_certificate",
)
TIMING_FIELDS = ("builder_seconds", "verifier_seconds")
_METRIC_FIELDS = {
    "status",
    "reason",
    *CERTIFICATION_ONLY_FIELDS,
    *TIMING_FIELDS,
}
TransferEvaluator = Callable[[Fraction], dict[str, object]]


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_json(payload: object) -> str:
    """Serialize an artifact deterministically with a trailing newline."""

    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def _empty_metrics(status: str, reason: str) -> dict[str, object]:
    return {
        "status": status,
        "reason": reason,
        "d4_term_count": None,
        "d4_group_count": None,
        "d4_norm_bound": None,
        "d5_term_count": None,
        "d5_group_count": None,
        "d5_norm_bound": None,
        "accepted_steps": None,
        "builder_seconds": None,
        "verifier_seconds": None,
        "source_certificate": None,
    }


def evaluate_delta(
    delta: Fraction,
    *,
    timeout_seconds: float = 60.0,
) -> dict[str, object]:
    """Evaluate one delta without extending the frozen compiler's scope."""

    if delta != 1:
        return _empty_metrics(
            "unsupported",
            "the frozen D4/D5 compiler is isotropic and has not been generalized to XXZ",
        )

    try:
        certificate = json.loads(CERTIFICATE.read_text())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return _empty_metrics("inconclusive", f"cannot load isotropic anchor: {exc}")

    candidate = certificate.get("candidate")
    if not isinstance(candidate, dict):
        return _empty_metrics("inconclusive", "isotropic anchor schema is malformed")
    d4 = candidate.get("d4_certificate")
    d5 = candidate.get("d5_certificate")
    if not isinstance(d4, dict) or not isinstance(d5, dict):
        return _empty_metrics("inconclusive", "isotropic D4/D5 metadata is absent")

    try:
        completed = subprocess.run(
            [sys.executable, str(REFERENCE_VERIFIER), str(CERTIFICATE)],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired:
        return _empty_metrics(
            "inconclusive",
            f"reference verification exceeded {timeout_seconds:g} seconds",
        )
    if completed.returncode != 0:
        return _empty_metrics(
            "inconclusive",
            "reference verifier rejected the isotropic anchor",
        )

    result = _empty_metrics(
        "certified",
        "delta=1 exactly reduces to the frozen isotropic Heisenberg certificate",
    )
    result.update(
        {
            "d4_term_count": d4.get("term_count"),
            "d4_group_count": d4.get("group_count"),
            "d4_norm_bound": d4.get("cell_norm_upper"),
            "d5_term_count": d5.get("term_count"),
            "d5_group_count": d5.get("group_count"),
            "d5_norm_bound": d5.get("site_norm_upper"),
            "accepted_steps": candidate.get("steps"),
            "source_certificate": CERTIFICATE.relative_to(ROOT).as_posix(),
        }
    )
    return result


def _validated_metrics(
    delta: Fraction,
    metrics: object,
) -> dict[str, object]:
    if not isinstance(metrics, dict) or set(metrics) != _METRIC_FIELDS:
        raise ValueError("XXZ evaluator must return the exact metric schema")
    status = metrics["status"]
    if status not in STATUSES:
        raise ValueError(f"unknown XXZ status {status!r}")
    if status == "certified" and delta != 1:
        raise ValueError("the current compiler can certify only delta=1")
    if any(metrics[field] is not None for field in TIMING_FIELDS):
        raise ValueError("execution timing must not enter the frozen scientific payload")
    if status != "certified":
        for field in CERTIFICATION_ONLY_FIELDS:
            if metrics[field] is not None:
                raise ValueError(f"uncertified XXZ rows must set {field} to null")
    else:
        for field in CERTIFICATION_ONLY_FIELDS:
            if metrics[field] is None:
                raise ValueError(f"certified XXZ rows require {field}")
    return metrics


def run_xxz_transfer(
    deltas: tuple[Fraction, ...] = PREREGISTERED_DELTAS,
    *,
    evaluator: TransferEvaluator | None = None,
) -> dict[str, object]:
    """Return one explicit result row for every requested exact anisotropy."""

    if not deltas or any(not isinstance(delta, Fraction) for delta in deltas):
        raise TypeError("deltas must be a nonempty tuple of Fractions")
    if len(set(deltas)) != len(deltas):
        raise ValueError("deltas must be unique")
    selected = evaluate_delta if evaluator is None else evaluator
    rows: list[dict[str, object]] = []
    for delta in sorted(deltas):
        metrics = _validated_metrics(delta, selected(delta))
        coefficients = xxz_bond(delta)
        rows.append(
            {
                "delta": _pair(delta),
                "bond_coefficients": {
                    name: _pair(value) for name, value in coefficients.items()
                },
                **metrics,
            }
        )
    return {
        "schema_version": 1,
        "model": "square_lattice_xxz",
        "normalization": "(XX+YY+delta*ZZ)/4",
        "sources": {
            "certificate": {
                "path": CERTIFICATE.relative_to(ROOT).as_posix(),
                "sha256": _digest(CERTIFICATE),
            },
            "reference_verifier": {
                "path": REFERENCE_VERIFIER.relative_to(ROOT).as_posix(),
                "sha256": _digest(REFERENCE_VERIFIER),
            },
            "generator": {
                "path": Path(__file__).resolve().relative_to(ROOT).as_posix(),
                "sha256": _digest(Path(__file__).resolve()),
            },
        },
        "execution_timing_included": False,
        "preregistered_deltas": [_pair(delta) for delta in PREREGISTERED_DELTAS],
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    arguments = parser.parse_args()
    if arguments.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    payload = run_xxz_transfer(
        evaluator=lambda delta: evaluate_delta(
            delta,
            timeout_seconds=arguments.timeout_seconds,
        )
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(canonical_json(payload))
    counts = {
        status: sum(row["status"] == status for row in payload["rows"])
        for status in sorted(STATUSES)
    }
    print(f"output={arguments.output}")
    print(" ".join(f"{status}={count}" for status, count in counts.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
