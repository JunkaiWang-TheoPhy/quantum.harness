#!/usr/bin/env python3
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

from trottercert.refined_error import (
    build_refined_fourth_order_constants,
    evaluate_refined_fourth_order_bound,
)
from trottercert.support_groups import decode_d5_gzip, verify_d5_payload


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "certificates" / "issue128-certificate.json"
D5 = ROOT / "certificates" / "issue128-d5-groups.json.gz"
OUTPUT = ROOT / "certificates" / "issue128-d5-integrated-certificate.json"


def _pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def build() -> dict[str, object]:
    data = json.loads(SOURCE.read_text())
    raw_d5 = D5.read_bytes()
    d5_partition = verify_d5_payload(decode_d5_gzip(raw_d5))
    d5_site = d5_partition.bound / 4
    d4_site = Fraction(*data["candidate"]["d4_certificate"]["cell_norm_upper"]) / 4
    constants = build_refined_fourth_order_constants(
        decimal_digits=int(data["candidate"]["coefficient_interval_decimal_digits"]),
        quantization_digits=int(data["candidate"]["e5_quantization_digits"]),
    )
    if d5_site > constants.d5_site:
        raise ValueError("grouped D5 sidecar does not tighten the generic D5 bound")

    n_sites = int(data["benchmark"]["length"]) ** 2
    tolerance = Fraction(*data["benchmark"]["tolerance"])

    def evaluate(steps: int):
        return evaluate_refined_fourth_order_bound(
            constants,
            n_sites,
            steps,
            d4_site_override=d4_site,
            d5_site_override=d5_site,
        )

    def meets_tolerance(steps: int) -> bool:
        try:
            return evaluate(steps).global_error_bound <= tolerance
        except ValueError as exc:
            if "convergence" not in str(exc):
                raise
            return False

    candidate_steps = next(
        steps
        for steps in range(2, int(data["published_baseline"]["steps"]) + 1)
        if meets_tolerance(steps)
    )
    candidate = evaluate(candidate_steps)
    previous = evaluate(candidate_steps - 1)
    candidate_groups = 30 * candidate_steps + 1
    published_groups = int(data["published_baseline"]["group_exponentials"])
    max_group_size = max(
        (len(group.term_indices) for group in d5_partition.groups), default=0
    )

    record = data["candidate"]
    record["proof_method"] = (
        "local_log_E5_grouped_D4_D5_plus_E7_majorant"
        "_plus_exact_generator_tail"
    )
    record["d5_certificate"] = {
        "path": D5.name,
        "sha256": hashlib.sha256(raw_d5).hexdigest(),
        "coefficient_interval_decimal_digits": 24,
        "max_group_size": max_group_size,
        "term_count": len(d5_partition.paulis),
        "group_count": len(d5_partition.groups),
        "site_norm_upper": _pair(d5_site),
    }
    record["steps"] = candidate_steps
    record["group_exponentials"] = candidate_groups
    record["contributions"] = {
        "degree4": _pair(candidate.degree_four_contribution),
        "degree5": _pair(candidate.degree_five_contribution),
        "degree6": _pair(candidate.degree_six_contribution),
        "degree7": _pair(candidate.degree_seven_contribution),
        "tail": _pair(candidate.tail_contribution),
    }
    record["global_error_upper"] = _pair(candidate.global_error_bound)
    record["previous_step_error_upper"] = _pair(previous.global_error_bound)

    baseline_bonds = published_groups * n_sites // 2
    candidate_bonds = candidate_groups * n_sites // 2
    data["claimed_resources"] = {
        "published_steps": int(data["published_baseline"]["steps"]),
        "candidate_steps": candidate_steps,
        "published_group_exponentials": published_groups,
        "candidate_group_exponentials": candidate_groups,
        "published_bond_propagators": baseline_bonds,
        "candidate_bond_propagators": candidate_bonds,
        "published_cnot_upper": 3 * baseline_bonds,
        "candidate_cnot_upper": 3 * candidate_bonds,
    }
    ratio = Fraction(published_groups, candidate_groups)
    data["claims"] = {
        "global_twofold_target_met": ratio >= 2,
        "global_fourfold_target_met": ratio >= 4,
        "exact_improvement_ratio": _pair(ratio),
    }
    OUTPUT.write_text(json.dumps(data, indent=2) + "\n")
    return {
        "output": str(OUTPUT),
        "candidate_steps": candidate_steps,
        "previous_step_error_upper": str(previous.global_error_bound),
        "candidate_error_upper": str(candidate.global_error_bound),
        "candidate_group_exponentials": candidate_groups,
        "exact_improvement_ratio": str(ratio),
        "improvement": float(ratio),
    }


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
