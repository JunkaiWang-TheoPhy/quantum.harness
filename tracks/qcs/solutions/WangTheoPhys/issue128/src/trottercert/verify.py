from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from .anticommuting import certify_anticommuting_partition
from .baseline import pauli_l1_second_order_constant
from .d6_physical_channels import build_grouped_d6_bound
from .exact_series_certificate import (
    ExactDegreeVerification,
    read_portable_canonical_gzip,
    verify_exact_degree_payload,
)
from .hamiltonian import four_matching_fragments
from .intervals import RationalInterval, cube_root_four_interval
from .lattice import SquareLattice
from .local_commutators import SymplecticPauli
from .resources import (
    fourth_order_four_matching_resources,
    four_matching_resources,
    required_steps,
    three_l_path_resources,
)
from .su2clusters import three_l_path_fragments
from .support_groups import decode_d5_gzip, verify_d5_payload


EXPECTED_NORMALIZATION = "(XX+YY+ZZ)/4"


def _integer(
    value: object,
    label: str,
    *,
    minimum: int | None = None,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{label} must be at least {minimum}")
    return value


@dataclass(frozen=True)
class D4SidecarVerification:
    site_bound: Fraction
    coefficients: dict[SymplecticPauli, RationalInterval]
    term_count: int
    group_count: int
    max_group_size: int


@dataclass(frozen=True)
class D5SidecarVerification:
    site_bound: Fraction
    coefficients: dict[SymplecticPauli, RationalInterval]
    term_count: int
    group_count: int
    max_group_size: int


@dataclass(frozen=True)
class D6SidecarVerification:
    site_bound: Fraction
    artifact: ExactDegreeVerification
    group_count: int = 0
    max_group_size: int = 1


def _fraction(pair: object) -> Fraction:
    if not isinstance(pair, list) or len(pair) != 2:
        raise ValueError("fraction must be [numerator, denominator]")
    numerator, denominator = pair
    if (
        isinstance(numerator, bool)
        or not isinstance(numerator, int)
        or isinstance(denominator, bool)
        or not isinstance(denominator, int)
    ):
        raise ValueError("fraction numerator and denominator must be integers")
    if denominator <= 0:
        raise ValueError("fraction denominator must be positive")
    return Fraction(numerator, denominator)


def _verify_d4_sidecar(
    certificate_path: Path,
    candidate: dict[str, object],
) -> D4SidecarVerification:
    metadata = candidate["d4_certificate"]
    root = certificate_path.resolve().parent
    sidecar_path = (root / str(metadata["path"])).resolve()
    if sidecar_path.parent != root:
        raise ValueError("D4 sidecar path escapes certificate directory")
    raw = sidecar_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != str(metadata["sha256"]):
        raise ValueError("D4 sidecar digest mismatch")

    payload = json.loads(raw)
    if not isinstance(payload, dict):
        raise ValueError("D4 sidecar root must be an object")
    canonical_payload = (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    if raw != canonical_payload:
        raise ValueError("D4 sidecar is not canonical JSON")
    if _integer(payload["schema_version"], "D4 sidecar schema") != 1:
        raise ValueError("unsupported D4 sidecar schema")
    coefficient_denominator = _integer(
        payload["coefficient_denominator"],
        "D4 coefficient denominator",
        minimum=1,
    )
    sqrt_denominator = _integer(
        payload["sqrt_denominator"],
        "D4 square-root denominator",
        minimum=1,
    )

    rows = payload["terms"]
    if not isinstance(rows, list):
        raise ValueError("D4 terms must be a list")
    pauli_rows: list[tuple[SymplecticPauli, RationalInterval]] = []
    for position, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != 4:
            raise ValueError(f"D4 term {position} is malformed")
        pauli = (
            _integer(row[0], "D4 x mask", minimum=0),
            _integer(row[1], "D4 z mask", minimum=0),
        )
        coefficient = RationalInterval(
            Fraction(
                _integer(row[2], "D4 coefficient lower endpoint"),
                coefficient_denominator,
            ),
            Fraction(
                _integer(row[3], "D4 coefficient upper endpoint"),
                coefficient_denominator,
            ),
        )
        pauli_rows.append((pauli, coefficient))
    paulis = tuple(pauli for pauli, _ in pauli_rows)
    if len(paulis) != len(set(paulis)):
        raise ValueError("D4 sidecar coefficient terms are duplicated")
    if paulis != tuple(sorted(paulis)):
        raise ValueError("D4 sidecar coefficient terms are not canonical")
    coefficients = dict(pauli_rows)

    submitted_groups = payload["groups"]
    if not isinstance(submitted_groups, list):
        raise ValueError("D4 groups must be a list")
    groups: list[tuple[SymplecticPauli, ...]] = []
    submitted_bounds: list[Fraction] = []
    for position, row in enumerate(submitted_groups):
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError(f"D4 group {position} is malformed")
        if not isinstance(row[0], list) or not row[0]:
            raise ValueError("D4 group indices must be a nonempty list")
        indices = tuple(
            _integer(index, "D4 group index", minimum=0)
            for index in row[0]
        )
        if any(index >= len(paulis) for index in indices):
            raise ValueError("D4 partition coverage index is out of range")
        groups.append(tuple(paulis[index] for index in indices))
        submitted_bounds.append(
            Fraction(
                _integer(
                    row[1],
                    "D4 group bound numerator",
                    minimum=0,
                ),
                sqrt_denominator,
            )
        )
    regenerated = certify_anticommuting_partition(
        coefficients,
        tuple(groups),
    )
    regenerated_bounds = tuple(
        group.bound for group in regenerated.groups
    )
    if tuple(submitted_bounds) != regenerated_bounds:
        raise ValueError("D4 group bound mismatch")

    max_group_size = max((len(group) for group in groups), default=0)
    if len(paulis) != _integer(
        metadata["term_count"],
        "D4 term count",
        minimum=0,
    ):
        raise ValueError("D4 term count mismatch")
    if len(groups) != _integer(
        metadata["group_count"],
        "D4 group count",
        minimum=0,
    ):
        raise ValueError("D4 group count mismatch")
    if max_group_size != _integer(
        metadata["max_group_size"],
        "D4 maximum group size",
        minimum=0,
    ):
        raise ValueError("D4 maximum group size mismatch")
    cell_bound = _fraction(payload["cell_bound"])
    if cell_bound != regenerated.bound:
        raise ValueError("D4 cell bound mismatch")
    if cell_bound != _fraction(metadata["cell_norm_upper"]):
        raise ValueError("D4 main-certificate bound mismatch")
    return D4SidecarVerification(
        site_bound=cell_bound / 4,
        coefficients=coefficients,
        term_count=len(paulis),
        group_count=len(groups),
        max_group_size=max_group_size,
    )


def _verify_d5_sidecar(
    certificate_path: Path,
    candidate: dict[str, object],
) -> D5SidecarVerification:
    metadata = candidate["d5_certificate"]
    root = certificate_path.resolve().parent
    sidecar_path = (root / str(metadata["path"])).resolve()
    if sidecar_path.parent != root:
        raise ValueError("D5 sidecar path escapes certificate directory")
    raw = sidecar_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != str(metadata["sha256"]):
        raise ValueError("D5 sidecar digest mismatch")
    payload = decode_d5_gzip(raw)
    regenerated = verify_d5_payload(payload)
    max_group_size = max(
        (len(group.term_indices) for group in regenerated.groups),
        default=0,
    )
    if len(regenerated.paulis) != _integer(
        metadata["term_count"],
        "D5 term count",
        minimum=0,
    ):
        raise ValueError("D5 term count mismatch")
    if len(regenerated.groups) != _integer(
        metadata["group_count"],
        "D5 group count",
        minimum=0,
    ):
        raise ValueError("D5 group count mismatch")
    if max_group_size != _integer(
        metadata["max_group_size"],
        "D5 maximum group size",
        minimum=0,
    ):
        raise ValueError("D5 maximum group size mismatch")
    site_bound = regenerated.bound / 4
    if site_bound != _fraction(metadata["site_norm_upper"]):
        raise ValueError("D5 main-certificate bound mismatch")
    return D5SidecarVerification(
        site_bound=site_bound,
        coefficients=dict(zip(regenerated.paulis, regenerated.coefficients)),
        term_count=len(regenerated.paulis),
        group_count=len(regenerated.groups),
        max_group_size=max_group_size,
    )


def _verify_d6_sidecar(
    certificate_path: Path,
    candidate: dict[str, object],
) -> D6SidecarVerification:
    metadata = candidate["d6_certificate"]
    root = certificate_path.resolve().parent
    sidecar_path = (root / str(metadata["path"])).resolve()
    if sidecar_path.parent != root:
        raise ValueError("D6 sidecar path escapes certificate directory")
    raw = sidecar_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != str(metadata["sha256"]):
        raise ValueError("D6 sidecar digest mismatch")
    payload = read_portable_canonical_gzip(sidecar_path)
    verified = verify_exact_degree_payload(
        payload,
        expected_degree=6,
        expected_source_commit=str(metadata["source_commit"]),
    )
    if verified.term_count != _integer(
        metadata["term_count"],
        "D6 term count",
        minimum=0,
    ):
        raise ValueError("D6 term count mismatch")
    if _integer(
        payload["coefficient_interval_decimal_digits"],
        "D6 payload coefficient precision",
        minimum=1,
    ) != _integer(
        metadata["coefficient_interval_decimal_digits"],
        "D6 metadata coefficient precision",
        minimum=1,
    ):
        raise ValueError("D6 coefficient interval precision mismatch")
    if verified.parent_sha256 != str(metadata["parent_sha256"]):
        raise ValueError("D6 parent digest mismatch")
    if "parent_path" in metadata:
        parent_path = (root / str(metadata["parent_path"])).resolve()
        if parent_path.parent != root:
            raise ValueError("D6 parent path escapes certificate directory")
        if hashlib.sha256(parent_path.read_bytes()).hexdigest() != verified.parent_sha256:
            raise ValueError("D6 parent artifact digest mismatch")
    if "groups_path" not in metadata:
        if verified.cell_l1_upper != _fraction(metadata["cell_norm_upper"]):
            raise ValueError("D6 cell bound mismatch")
        if verified.site_l1_upper != _fraction(metadata["site_norm_upper"]):
            raise ValueError("D6 main-certificate bound mismatch")
        return D6SidecarVerification(
            site_bound=verified.site_l1_upper,
            artifact=verified,
        )

    comparison_path = (root / str(metadata["comparison_path"])).resolve()
    if comparison_path.parent != root:
        raise ValueError("D6 comparison path escapes certificate directory")
    comparison_raw = comparison_path.read_bytes()
    if hashlib.sha256(comparison_raw).hexdigest() != str(
        metadata["comparison_sha256"]
    ):
        raise ValueError("D6 comparison audit digest mismatch")
    try:
        comparison = json.loads(comparison_raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("D6 comparison audit is not JSON") from error
    canonical_comparison = (
        json.dumps(comparison, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    if comparison_raw != canonical_comparison:
        raise ValueError("D6 comparison audit is not canonical JSON")
    if (
        comparison.get("schema_version") != 1
        or comparison.get("kind")
        != "issue128_exact_degree_payload_comparison"
        or comparison.get("status") != "complete"
        or comparison.get("degree") != 6
        or comparison.get("exact_coefficient_map_equal") is not True
    ):
        raise ValueError("D6 comparison audit status mismatch")
    if comparison.get("term_count") != verified.term_count:
        raise ValueError("D6 comparison audit term count mismatch")
    if comparison.get("cell_pauli_l1_upper") != payload.get(
        "cell_pauli_l1_upper"
    ) or comparison.get("site_pauli_l1_upper") != payload.get(
        "site_pauli_l1_upper"
    ):
        raise ValueError("D6 comparison audit norm mismatch")
    comparison_lanes = (comparison.get("left"), comparison.get("right"))
    if any(not isinstance(lane, dict) for lane in comparison_lanes):
        raise ValueError("D6 comparison audit lane metadata is malformed")
    if any(
        lane.get("source_commit") != verified.source_commit
        for lane in comparison_lanes
    ):
        raise ValueError("D6 comparison audit source commit mismatch")
    if hashlib.sha256(raw).hexdigest() not in {
        str(lane.get("sha256")) for lane in comparison_lanes
    }:
        raise ValueError("D6 sidecar is absent from comparison audit")
    if verified.parent_sha256 not in {
        str(lane.get("parent_sha256")) for lane in comparison_lanes
    }:
        raise ValueError("D6 parent is absent from comparison audit")

    groups_path = (root / str(metadata["groups_path"])).resolve()
    if groups_path.parent != root:
        raise ValueError("D6 groups path escapes certificate directory")
    groups_raw = groups_path.read_bytes()
    if hashlib.sha256(groups_raw).hexdigest() != str(metadata["groups_sha256"]):
        raise ValueError("D6 groups sidecar digest mismatch")
    try:
        groups_payload = json.loads(groups_raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("D6 groups sidecar is not JSON") from error
    canonical_groups = (
        json.dumps(groups_payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    if groups_raw != canonical_groups:
        raise ValueError("D6 groups sidecar is not canonical JSON")
    if groups_payload.get("schema_version") != 1:
        raise ValueError("unsupported grouped D6 sidecar schema")
    if groups_payload.get("kind") != "issue128_d6_physical_channel_groups":
        raise ValueError("unexpected grouped D6 sidecar kind")
    if groups_payload.get("source_payload_sha256") != hashlib.sha256(raw).hexdigest():
        raise ValueError("grouped D6 source payload digest mismatch")
    if groups_payload.get("source_commit") != verified.source_commit:
        raise ValueError("grouped D6 source commit mismatch")
    if groups_payload.get("degree") != 6:
        raise ValueError("grouped D6 degree mismatch")
    digits = _integer(
        payload["coefficient_interval_decimal_digits"],
        "D6 coefficient precision",
        minimum=1,
    )
    if groups_payload.get("coefficient_interval_decimal_digits") != digits:
        raise ValueError("grouped D6 coefficient precision mismatch")
    candidate_cap = groups_payload.get("candidate_cap")
    if (
        isinstance(candidate_cap, bool)
        or not isinstance(candidate_cap, int)
        or candidate_cap < 1
    ):
        raise ValueError("grouped D6 candidate cap is invalid")

    regenerated = build_grouped_d6_bound(
        verified.terms,
        digits,
        candidate_cap=candidate_cap,
    )
    submitted_groups = groups_payload.get("groups")
    if not isinstance(submitted_groups, list):
        raise ValueError("grouped D6 groups must be a list")
    decoded_groups: list[tuple[int, ...]] = []
    submitted_bounds: list[Fraction] = []
    for group in submitted_groups:
        if not isinstance(group, dict):
            raise ValueError("grouped D6 group is malformed")
        indices = group.get("term_indices")
        if not isinstance(indices, list) or any(
            isinstance(index, bool) or not isinstance(index, int)
            for index in indices
        ):
            raise ValueError("grouped D6 term indices are malformed")
        decoded_groups.append(tuple(indices))
        submitted_bounds.append(_fraction(group["bound"]))
    if tuple(decoded_groups) != regenerated.groups:
        raise ValueError("grouped D6 partition differs from exact regeneration")

    regenerated_bounds = regenerated.group_bounds
    if tuple(submitted_bounds) != regenerated_bounds:
        raise ValueError("grouped D6 group bound mismatch")
    if groups_payload.get("term_count") != regenerated.term_count:
        raise ValueError("grouped D6 term count mismatch")
    if groups_payload.get("group_count") != len(regenerated.groups):
        raise ValueError("grouped D6 group count mismatch")
    if groups_payload.get("max_group_size") != regenerated.max_group_size:
        raise ValueError("grouped D6 maximum group size mismatch")
    if _fraction(groups_payload["cell_pauli_l1_upper"]) != verified.cell_l1_upper:
        raise ValueError("grouped D6 l1 cell baseline mismatch")
    if _fraction(groups_payload["site_pauli_l1_upper"]) != verified.site_l1_upper:
        raise ValueError("grouped D6 l1 site baseline mismatch")
    if _fraction(groups_payload["grouped_cell_bound"]) != regenerated.grouped_cell_bound:
        raise ValueError("grouped D6 cell bound mismatch")
    if _fraction(groups_payload["grouped_site_bound"]) != regenerated.grouped_site_bound:
        raise ValueError("grouped D6 site bound mismatch")
    if _fraction(metadata["l1_site_norm_upper"]) != verified.site_l1_upper:
        raise ValueError("D6 main-certificate l1 baseline mismatch")
    if _fraction(metadata["cell_norm_upper"]) != regenerated.grouped_cell_bound:
        raise ValueError("D6 grouped main-certificate cell bound mismatch")
    if _fraction(metadata["site_norm_upper"]) != regenerated.grouped_site_bound:
        raise ValueError("D6 grouped main-certificate site bound mismatch")
    if _integer(
        metadata["group_count"],
        "D6 group count",
        minimum=0,
    ) != len(regenerated.groups):
        raise ValueError("D6 main-certificate group count mismatch")
    if _integer(
        metadata["max_group_size"],
        "D6 maximum group size",
        minimum=0,
    ) != regenerated.max_group_size:
        raise ValueError("D6 main-certificate maximum group size mismatch")
    return D6SidecarVerification(
        site_bound=regenerated.grouped_site_bound,
        artifact=verified,
        group_count=len(regenerated.groups),
        max_group_size=regenerated.max_group_size,
    )


def _verify_v1(data: dict[str, object]) -> dict[str, object]:
    benchmark = data["benchmark"]
    if benchmark["normalization"] != EXPECTED_NORMALIZATION:
        raise ValueError("Hamiltonian normalization mismatch")
    length = _integer(benchmark["length"], "benchmark length", minimum=1)
    if length % 6:
        raise ValueError("benchmark length must be divisible by six")
    tolerance = _fraction(benchmark["tolerance"])
    reference_length = int(data["proof"]["reference_length"])
    if reference_length != 6:
        raise ValueError("the verifier currently pins reference length L=6")

    lattice = SquareLattice(reference_length)
    baseline_ref = pauli_l1_second_order_constant(
        four_matching_fragments(lattice)
    )
    candidate_ref = pauli_l1_second_order_constant(
        three_l_path_fragments(lattice)
    )
    baseline_density = baseline_ref / lattice.n_sites
    candidate_density = candidate_ref / lattice.n_sites
    if baseline_density != _fraction(data["proof"]["baseline_density"]):
        raise ValueError("baseline density summary is inconsistent")
    if candidate_density != _fraction(data["proof"]["candidate_density"]):
        raise ValueError("candidate density summary is inconsistent")

    n_sites = length * length
    baseline_constant = baseline_density * n_sites
    candidate_constant = candidate_density * n_sites
    baseline_steps = required_steps(baseline_constant, tolerance, 2)
    candidate_steps = required_steps(candidate_constant, tolerance, 2)
    baseline = four_matching_resources(n_sites, baseline_steps)
    candidate = three_l_path_resources(n_sites, candidate_steps)

    claimed = data["claimed_resources"]
    recomputed = {
        "baseline_steps": baseline.steps,
        "candidate_steps": candidate.steps,
        "baseline_local_propagators": baseline.local_propagators,
        "candidate_local_propagators": candidate.local_propagators,
        "baseline_cnot_upper": baseline.cnot_upper,
        "candidate_cnot_upper": candidate.cnot_upper,
    }
    if claimed != recomputed:
        raise ValueError("claimed resources do not match recomputation")
    local_ratio = Fraction(
        candidate.local_propagators,
        baseline.local_propagators,
    )
    return {
        "valid": True,
        "baseline_density": str(baseline_density),
        "candidate_density": str(candidate_density),
        "local_propagator_ratio": str(local_ratio),
        "local_propagator_improvement": float(1 / local_ratio),
        "compiled_cnot_improvement": (
            baseline.cnot_upper / candidate.cnot_upper
        ),
        "twofold_local_target_met": local_ratio <= Fraction(1, 2),
        "twofold_compiled_cnot_target_met": (
            candidate.cnot_upper * 2 <= baseline.cnot_upper
        ),
    }


def _verify_v2(
    data: dict[str, object],
    *,
    deep: bool,
) -> dict[str, object]:
    benchmark = data["benchmark"]
    if benchmark["normalization"] != EXPECTED_NORMALIZATION:
        raise ValueError("Hamiltonian normalization mismatch")
    length = _integer(benchmark["length"], "benchmark length", minimum=1)
    if length % 2:
        raise ValueError("four-matching benchmark requires even length")
    n_sites = length * length
    time = _fraction(benchmark["time"])
    tolerance = _fraction(benchmark["tolerance"])
    if time != 1:
        raise ValueError("schema v2 currently pins T=1")

    proof = data["proof"]
    if proof["formula"] != "five_copy_suzuki_fourth_order":
        raise ValueError("unexpected candidate formula")
    if int(proof["stage_count"]) != 31 or int(proof["center"]) != 17:
        raise ValueError("candidate formula structure mismatch")
    decimal_digits = int(proof["coefficient_interval_decimal_digits"])
    root = cube_root_four_interval(decimal_digits)
    claimed_root = proof["cube_root_four_interval"]
    if (
        root.lower != _fraction(claimed_root["lower"])
        or root.upper != _fraction(claimed_root["upper"])
        or not (root.lower**3 <= 4 <= root.upper**3)
    ):
        raise ValueError("cube-root enclosure is invalid")

    baseline_density = _fraction(proof["baseline_site_density"])
    if baseline_density <= 0 or baseline_density > Fraction(21, 8):
        raise ValueError("strengthened baseline density is outside the pinned range")
    candidate_density = _fraction(proof["candidate_site_density_upper"])
    if candidate_density <= 0:
        raise ValueError("candidate density must be positive")

    deep_verified = False
    if deep:
        from .baseline import strang_commutator_operators
        from .clusters import phase_partitioned_collatz_certificate
        from .rigorous_fourth import fourth_order_rational_pair_certificate

        reference_length = int(proof["baseline_reference_length"])
        if reference_length != 6:
            raise ValueError("deep baseline currently pins L=6")
        reference_lattice = SquareLattice(reference_length)
        baseline_total = Fraction()
        component_bounds: list[list[list[int]]] = []
        for block in strang_commutator_operators(
            four_matching_fragments(reference_lattice)
        ):
            repeated_fragment = phase_partitioned_collatz_certificate(
                block.repeated_fragment,
                reference_lattice,
                iterations=35,
            ).global_bound
            repeated_tail = phase_partitioned_collatz_certificate(
                block.repeated_tail,
                reference_lattice,
                iterations=35,
            ).global_bound
            baseline_total += repeated_fragment / 24 + repeated_tail / 12
            component_bounds.append(
                [
                    [repeated_fragment.numerator, repeated_fragment.denominator],
                    [repeated_tail.numerator, repeated_tail.denominator],
                ]
            )
        regenerated_baseline_density = (
            baseline_total / reference_lattice.n_sites
        )
        if regenerated_baseline_density != baseline_density:
            raise ValueError("deep regenerated baseline differs from certificate")
        if component_bounds != proof["baseline_component_bounds"]:
            raise ValueError("deep baseline component bounds differ")

        regenerated = fourth_order_rational_pair_certificate(
            center=17,
            decimal_digits=decimal_digits,
        )
        if regenerated.site_density_upper != candidate_density:
            raise ValueError("deep regenerated density differs from certificate")
        statistics = proof["statistics"]
        if (
            regenerated.theorem_terms != int(statistics["theorem_terms"])
            or regenerated.paired_terms != int(statistics["paired_terms"])
            or regenerated.singleton_terms != int(statistics["singleton_terms"])
        ):
            raise ValueError("deep proof statistics differ from certificate")
        deep_verified = True

    baseline_steps = required_steps(
        baseline_density * n_sites,
        tolerance,
        2,
        time,
    )
    candidate_steps = required_steps(
        candidate_density * n_sites,
        tolerance,
        4,
        time,
    )
    baseline = four_matching_resources(n_sites, baseline_steps)
    candidate = fourth_order_four_matching_resources(
        n_sites,
        candidate_steps,
        stage_count=31,
    )
    claimed = data["claimed_resources"]
    recomputed = {
        "baseline_steps": baseline.steps,
        "candidate_steps": candidate.steps,
        "baseline_group_exponentials": baseline.group_exponentials,
        "candidate_group_exponentials": candidate.group_exponentials,
        "baseline_bond_propagators": baseline.local_propagators,
        "candidate_bond_propagators": candidate.local_propagators,
        "baseline_cnot_upper": baseline.cnot_upper,
        "candidate_cnot_upper": candidate.cnot_upper,
    }
    if claimed != recomputed:
        raise ValueError("claimed resources do not match recomputation")

    ratio = Fraction(candidate.cnot_upper, baseline.cnot_upper)
    twofold = candidate.cnot_upper * 2 <= baseline.cnot_upper
    if bool(data["claims"]["compiled_cnot_improvement_exceeds_two"]) != twofold:
        raise ValueError("twofold claim is inconsistent")
    return {
        "valid": True,
        "verification_level": "deep" if deep_verified else "fast",
        "deep_proof_regenerated": deep_verified,
        "baseline_site_density": str(baseline_density),
        "candidate_site_density_upper": str(candidate_density),
        "baseline_steps": baseline.steps,
        "candidate_steps": candidate.steps,
        "compiled_cnot_ratio": str(ratio),
        "compiled_cnot_improvement": float(1 / ratio),
        "twofold_compiled_cnot_target_met": twofold,
    }


def _verify_v3(
    data: dict[str, object],
    *,
    deep: bool,
    certificate_path: Path,
) -> dict[str, object]:
    benchmark = data["benchmark"]
    if benchmark["normalization"] != EXPECTED_NORMALIZATION:
        raise ValueError("Hamiltonian normalization mismatch")
    length = _integer(benchmark["length"], "benchmark length", minimum=1)
    if length % 2:
        raise ValueError("four-matching benchmark requires even length")
    n_sites = length * length
    if _fraction(benchmark["time"]) != 1:
        raise ValueError("schema v3 pins T=1")
    tolerance = _fraction(benchmark["tolerance"])

    published = data["published_baseline"]
    if (
        published["formula"] != "five_copy_suzuki_fourth_order"
        or _integer(published["formula_order"], "published formula order") != 4
        or _integer(published["stage_count"], "published stage count") != 31
        or _integer(published["theorem_center"], "published theorem center") != 20
    ):
        raise ValueError("published baseline structure mismatch")
    published_density = _fraction(published["site_density_upper"])
    published_steps = required_steps(
        published_density * n_sites,
        tolerance,
        4,
    )
    if published_steps != _integer(
        published["steps"],
        "published step count",
        minimum=1,
    ):
        raise ValueError("published baseline step count mismatch")
    published_groups = 30 * published_steps + 1
    if published_groups != _integer(
        published["group_exponentials"],
        "published group count",
        minimum=1,
    ):
        raise ValueError("published baseline group count mismatch")

    candidate = data["candidate"]
    legacy_method = (
        "local_log_E5_grouped_D4_plus_E7_majorant"
        "_plus_exact_generator_tail"
    )
    d5_method = (
        "local_log_E5_grouped_D4_D5_plus_E7_majorant"
        "_plus_exact_generator_tail"
    )
    d6_method = (
        "local_log_E5_grouped_D4_D5_exact_D6_plus_E7_majorant"
        "_plus_exact_generator_tail"
    )
    d6_grouped_method = (
        "local_log_E5_grouped_D4_D5_exact_grouped_D6_plus_E7_majorant"
        "_plus_exact_generator_tail"
    )
    if candidate["formula"] != "five_copy_suzuki_fourth_order" or candidate[
        "proof_method"
    ] not in {legacy_method, d5_method, d6_method, d6_grouped_method}:
        raise ValueError("candidate structure mismatch")
    has_d5 = "d5_certificate" in candidate
    has_d6 = "d6_certificate" in candidate
    if has_d5 != (
        candidate["proof_method"]
        in {d5_method, d6_method, d6_grouped_method}
    ):
        raise ValueError("candidate D5 proof method and sidecar disagree")
    if has_d6 != (
        candidate["proof_method"] in {d6_method, d6_grouped_method}
    ):
        raise ValueError("candidate D6 proof method and sidecar disagree")
    d4_verification = _verify_d4_sidecar(
        certificate_path,
        candidate,
    )
    d4_metadata = candidate["d4_certificate"]
    if d4_verification.term_count != _integer(
        d4_metadata["term_count"],
        "D4 term count",
        minimum=0,
    ):
        raise ValueError("D4 term count mismatch")
    if d4_verification.group_count != _integer(
        d4_metadata["group_count"],
        "D4 group count",
        minimum=0,
    ):
        raise ValueError("D4 group count mismatch")
    if (
        d4_verification.max_group_size
        > _integer(
            d4_metadata["max_group_size"],
            "D4 maximum group size",
            minimum=0,
        )
    ):
        raise ValueError("D4 maximum group size mismatch")
    d5_verification = (
        _verify_d5_sidecar(certificate_path, candidate) if has_d5 else None
    )
    d6_verification = (
        _verify_d6_sidecar(certificate_path, candidate) if has_d6 else None
    )

    # Fast verification must reconstruct the complete finite-step ledger rather
    # than merely checking that submitted D4--D7 and tail entries add up.  The
    # expensive coefficient-by-coefficient D4/D5/D6 replay remains a deep-mode
    # responsibility, but every accepted mode independently derives the
    # formula constants, degree contributions, tail, and adjacent-step value.
    from .refined_error import (
        build_refined_fourth_order_constants,
        evaluate_refined_fourth_order_bound,
    )

    constants = build_refined_fourth_order_constants(
        decimal_digits=_integer(
            candidate["coefficient_interval_decimal_digits"],
            "candidate coefficient precision",
            minimum=1,
        ),
        quantization_digits=_integer(
            candidate["e5_quantization_digits"],
            "candidate E5 quantization precision",
            minimum=1,
        ),
    )
    if constants.e5_site_l1 != _fraction(candidate["e5_site_l1_upper"]):
        raise ValueError("E5 site bound regeneration mismatch")
    if constants.e7_site_majorant != _fraction(candidate["e7_site_majorant"]):
        raise ValueError("E7 site majorant regeneration mismatch")

    candidate_steps = _integer(
        candidate["steps"],
        "candidate step count",
        minimum=2,
    )
    candidate_error = _fraction(candidate["global_error_upper"])
    previous_error = _fraction(candidate["previous_step_error_upper"])
    if candidate_error > tolerance:
        raise ValueError("candidate does not meet the requested tolerance")
    if previous_error <= tolerance:
        raise ValueError("candidate step count is not minimal for this certificate")
    candidate_groups = 30 * candidate_steps + 1
    if candidate_groups != _integer(
        candidate["group_exponentials"],
        "candidate group count",
        minimum=1,
    ):
        raise ValueError("candidate group count mismatch")
    contributions = candidate["contributions"]
    contribution_sum = sum(
        (_fraction(contributions[key]) for key in ("degree4", "degree5", "degree6", "degree7", "tail")),
        Fraction(),
    )
    if contribution_sum != candidate_error:
        raise ValueError("candidate contribution sum mismatch")
    expected_degree_four = (
        Fraction(n_sites)
        * d4_verification.site_bound
        / (5 * candidate_steps**4)
    )
    if _fraction(contributions["degree4"]) != expected_degree_four:
        raise ValueError("candidate grouped D4 contribution mismatch")
    if d5_verification is not None:
        expected_degree_five = (
            Fraction(n_sites)
            * d5_verification.site_bound
            / (6 * candidate_steps**5)
        )
        if _fraction(contributions["degree5"]) != expected_degree_five:
            raise ValueError("candidate grouped D5 contribution mismatch")
    if d6_verification is not None:
        expected_degree_six = (
            Fraction(n_sites)
            * d6_verification.site_bound
            / (7 * candidate_steps**6)
        )
        if _fraction(contributions["degree6"]) != expected_degree_six:
            raise ValueError("candidate exact D6 contribution mismatch")

    rebuilt = evaluate_refined_fourth_order_bound(
        constants,
        n_sites,
        candidate_steps,
        d4_site_override=d4_verification.site_bound,
        d5_site_override=(
            d5_verification.site_bound
            if d5_verification is not None
            else None
        ),
        d6_site_override=(
            d6_verification.site_bound
            if d6_verification is not None
            else None
        ),
    )
    rebuilt_previous = evaluate_refined_fourth_order_bound(
        constants,
        n_sites,
        candidate_steps - 1,
        d4_site_override=d4_verification.site_bound,
        d5_site_override=(
            d5_verification.site_bound
            if d5_verification is not None
            else None
        ),
        d6_site_override=(
            d6_verification.site_bound
            if d6_verification is not None
            else None
        ),
    )
    regenerated_contributions = {
        "degree4": rebuilt.degree_four_contribution,
        "degree5": rebuilt.degree_five_contribution,
        "degree6": rebuilt.degree_six_contribution,
        "degree7": rebuilt.degree_seven_contribution,
        "tail": rebuilt.tail_contribution,
    }
    for name, regenerated in regenerated_contributions.items():
        if _fraction(contributions[name]) != regenerated:
            raise ValueError(f"candidate {name} contribution regeneration mismatch")
    if rebuilt.global_error_bound != candidate_error:
        raise ValueError("candidate bound regeneration mismatch")
    if rebuilt_previous.global_error_bound != previous_error:
        raise ValueError("candidate adjacent-step regeneration mismatch")

    claimed = data["claimed_resources"]
    baseline_bonds = published_groups * n_sites // 2
    candidate_bonds = candidate_groups * n_sites // 2
    recomputed = {
        "published_steps": published_steps,
        "candidate_steps": candidate_steps,
        "published_group_exponentials": published_groups,
        "candidate_group_exponentials": candidate_groups,
        "published_bond_propagators": baseline_bonds,
        "candidate_bond_propagators": candidate_bonds,
        "published_cnot_upper": 3 * baseline_bonds,
        "candidate_cnot_upper": 3 * candidate_bonds,
    }
    if claimed != recomputed:
        raise ValueError("schema v3 resource summary mismatch")

    deep_verified = False
    baseline_centers_scanned = 0
    if deep:
        from .refined_error import (
            certified_d4_cell_coefficients,
        )
        from .rigorous_fourth import (
            fourth_order_published_triangle_center_scan,
        )

        published_scan = fourth_order_published_triangle_center_scan(
            decimal_digits=_integer(
                published["coefficient_interval_decimal_digits"],
                "published coefficient precision",
                minimum=1,
            ),
        )
        baseline_centers_scanned = len(published_scan)
        published_rebuilt = min(
            published_scan,
            key=lambda certificate: certificate.site_density_upper,
        )
        if published_rebuilt.center != _integer(
            published["theorem_center"],
            "published theorem center",
        ):
            raise ValueError("deep published baseline center scan mismatch")
        if published_rebuilt.site_density_upper != published_density:
            raise ValueError("deep published baseline regeneration mismatch")
        if (
            published_rebuilt.theorem_terms
            != _integer(
                published["theorem_terms"],
                "published theorem term count",
                minimum=0,
            )
            or published_rebuilt.expanded_commutator_keys
            != _integer(
                published["expanded_commutator_keys"],
                "published expanded commutator count",
                minimum=0,
            )
        ):
            raise ValueError("deep published baseline statistics mismatch")

        regenerated_d4 = certified_d4_cell_coefficients(
            constants.stages,
            quantization_digits=_integer(
                candidate["e5_quantization_digits"],
                "candidate E5 quantization precision",
                minimum=1,
            ),
        )
        if regenerated_d4 != d4_verification.coefficients:
            raise ValueError("deep D4 coefficient regeneration mismatch")
        if d5_verification is not None:
            from .cubic_field import fourth_order_suzuki_cubic_stages
            from .cubic_local import exact_d5_density, exact_log_e5_density

            exact_stages = fourth_order_suzuki_cubic_stages(4)
            registry, exact_e5 = exact_log_e5_density(exact_stages)
            exact_d5 = exact_d5_density(registry, exact_e5)
            d5_root = cube_root_four_interval(
                _integer(
                    candidate["d5_certificate"][
                        "coefficient_interval_decimal_digits"
                    ],
                    "D5 coefficient precision",
                    minimum=1,
                )
            )
            regenerated_d5 = {
                pauli: coefficient.enclose(d5_root)
                for pauli, coefficient in exact_d5.items()
            }
            if regenerated_d5 != d5_verification.coefficients:
                raise ValueError("deep D5 coefficient regeneration mismatch")
        if d6_verification is not None:
            from .cubic_field import fourth_order_suzuki_cubic_stages
            from .cubic_local import exact_right_generator_local_series
            from .hpc_artifacts import (
                coordinate_decode_terms,
                coordinate_encode_terms,
            )

            exact_stages = fourth_order_suzuki_cubic_stages(4)
            registry, exact_series = exact_right_generator_local_series(
                exact_stages,
                6,
            )
            regenerated_d6 = coordinate_decode_terms(
                coordinate_encode_terms(registry, exact_series[6])
            )
            if regenerated_d6 != d6_verification.artifact.terms:
                raise ValueError("deep D6 coefficient regeneration mismatch")
        deep_verified = True

    ratio = Fraction(published_groups, candidate_groups)
    claimed_ratio = _fraction(data["claims"]["exact_improvement_ratio"])
    if ratio != claimed_ratio:
        raise ValueError("claimed improvement ratio mismatch")
    if bool(data["claims"]["global_twofold_target_met"]) != (ratio >= 2):
        raise ValueError("global twofold claim mismatch")
    if bool(data["claims"]["global_fourfold_target_met"]) != (ratio >= 4):
        raise ValueError("global fourfold claim mismatch")
    result = {
        "valid": True,
        "verification_level": "deep" if deep_verified else "fast",
        "deep_proof_regenerated": deep_verified,
        "finite_step_bound_recomputed": True,
        "baseline_centers_scanned": baseline_centers_scanned,
        "published_steps": published_steps,
        "candidate_steps": candidate_steps,
        "published_group_exponentials": published_groups,
        "candidate_group_exponentials": candidate_groups,
        "candidate_error_upper": str(candidate_error),
        "previous_step_error_upper": str(previous_error),
        "d4_cell_norm_upper": str(
            4 * d4_verification.site_bound
        ),
        "d4_term_count": d4_verification.term_count,
        "d4_group_count": d4_verification.group_count,
        "d4_max_group_size": d4_verification.max_group_size,
        "exact_improvement_ratio": str(ratio),
        "improvement": float(ratio),
        "global_twofold_target_met": ratio >= 2,
        "global_fourfold_target_met": ratio >= 4,
    }
    if d5_verification is not None:
        result.update(
            {
                "d5_site_norm_upper": str(d5_verification.site_bound),
                "d5_term_count": d5_verification.term_count,
                "d5_group_count": d5_verification.group_count,
                "d5_max_group_size": d5_verification.max_group_size,
            }
        )
    if d6_verification is not None:
        result.update(
            {
                "d6_site_norm_upper": str(d6_verification.site_bound),
                "d6_term_count": d6_verification.artifact.term_count,
                "d6_parent_sha256": d6_verification.artifact.parent_sha256,
            }
        )
    return result


def verify_certificate(
    path: str | Path,
    *,
    deep: bool = False,
) -> dict[str, object]:
    certificate_path = Path(path)
    data = json.loads(certificate_path.read_text())
    schema = data.get("schema_version")
    if isinstance(schema, bool) or not isinstance(schema, int):
        raise ValueError("certificate schema version must be an integer")
    if schema == 1:
        if deep:
            raise ValueError("deep verification is available only for schema v2")
        return _verify_v1(data)
    if schema == 2:
        return _verify_v2(data, deep=deep)
    if schema == 3:
        return _verify_v3(
            data,
            deep=deep,
            certificate_path=certificate_path,
        )
    raise ValueError("unsupported certificate schema")
