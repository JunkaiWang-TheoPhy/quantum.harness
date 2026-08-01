#!/usr/bin/env python3
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from scripts.audit_paper_a_claims import audit_claims

with (ROOT / "certificates/issue128-d5-integrated-certificate.json").open() as handle:
    certificate = json.load(handle)
with (ROOT / "certificates/issue128-certificate.json").open() as handle:
    frozen_predecessor = json.load(handle)

baseline = certificate["published_baseline"]
candidate = certificate["candidate"]
resources = certificate["claimed_resources"]

headline_errors: list[str] = []


def require(condition: bool, message: str) -> None:
    if not condition:
        headline_errors.append(message)


require(baseline["steps"] == 393, "published baseline must use 393 steps")
require(candidate["steps"] == 95, "candidate must use 95 steps")
require(
    resources["published_group_exponentials"] == 11791,
    "published resource count must equal 11791 groups",
)
require(
    resources["candidate_group_exponentials"] == 2851,
    "candidate resource count must equal 2851 groups",
)
require(
    resources["published_bond_propagators"] == 848952,
    "published bond-propagator count is inconsistent",
)
require(
    resources["candidate_bond_propagators"] == 205272,
    "candidate bond-propagator count is inconsistent",
)
require(
    resources["published_cnot_upper"] == 2546856,
    "published CNOT upper bound is inconsistent",
)
require(
    resources["candidate_cnot_upper"] == 615816,
    "candidate CNOT upper bound is inconsistent",
)
require(Fraction(11791, 2851) > 4, "published improvement must exceed fourfold")
require(
    certificate["claims"]["exact_improvement_ratio"] == [11791, 2851],
    "exact improvement ratio must equal 11791/2851",
)

tolerance = Fraction(*certificate["benchmark"]["tolerance"])
accepted = Fraction(*candidate["global_error_upper"])
rejected = Fraction(*candidate["previous_step_error_upper"])
require(accepted < tolerance, "95-step bound must be below tolerance")
require(rejected > tolerance, "94-step bound must exceed tolerance")
require(candidate["d4_certificate"]["term_count"] == 75324, "D4 term count drift")
require(candidate["d4_certificate"]["group_count"] == 7576, "D4 group count drift")
require(candidate["d4_certificate"]["max_group_size"] == 10, "D4 max group-size drift")
require(candidate["d5_certificate"]["term_count"] == 605832, "D5 term count drift")
require(candidate["d5_certificate"]["group_count"] == 123106, "D5 group count drift")
require(candidate["d5_certificate"]["max_group_size"] == 10, "D5 max group-size drift")

# The separately hash-bound 97-step certificate remains an auditable
# predecessor; it is not the numerical source for the current manuscript.
require(frozen_predecessor["candidate"]["steps"] == 97, "predecessor step count drift")
require(
    frozen_predecessor["claimed_resources"]["candidate_group_exponentials"] == 2911,
    "predecessor resource count drift",
)

headline_errors.extend(
    audit_claims(
        ROOT / "artifacts/publication/paper-a-claim-matrix.json",
        ROOT,
    )
)
if headline_errors:
    raise SystemExit("\n".join(headline_errors))

print("headline_claims=PASS")
