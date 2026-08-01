import copy
from dataclasses import replace
from fractions import Fraction

import pytest

from scripts.certify_extensive_commutant_witness import (
    ExtensiveSizeRecord,
    build_extensive_payload,
    verify_extensive_payload,
)
from trottercert.commutant_witness import QuadraticWitnessMoments
from trottercert.cubic_field import Cubic


def _record(length: int) -> ExtensiveSizeRecord:
    n_sites = length * length
    cells = n_sites // 4
    tau_h2 = Fraction(3 * n_sites, 8)
    tau_h3 = Fraction(-3 * n_sites, 16)
    tau_h4 = Fraction(27 * n_sites**2 - 21 * n_sites, 64)
    tau_w2 = Fraction(9 * n_sites * (2 * n_sites - 3), 64)
    tau_h_e5 = cells * Cubic(1, 2, 3)
    tau_h2_e5 = cells * Cubic(4, 5, 6)
    tau_w_e5 = tau_h2_e5 + tau_h_e5 / 2
    return ExtensiveSizeRecord(
        length=length,
        n_sites=n_sites,
        cells=cells,
        matching_reconstruction_terms=6 * n_sites,
        moments=QuadraticWitnessMoments(
            tau_h2=tau_h2,
            tau_h3=tau_h3,
            tau_h4=tau_h4,
            identity_coefficient=-tau_h2,
            h_coefficient=Fraction(1, 2),
            tau_w2=tau_w2,
            tau_h_e5=tau_h_e5,
            tau_h2_e5=tau_h2_e5,
            tau_w_e5=tau_w_e5,
            squared_normalized_pairing=tau_w_e5**2 / tau_w2,
        ),
    )


def test_extensive_payload_certifies_stable_per_cell_pairing() -> None:
    payload = build_extensive_payload(
        (_record(6), _record(8)),
        rejected_alias_lengths=(4,),
        implementation_sources={"source.py": "a" * 64},
    )
    verify_extensive_payload(payload)
    assert payload["accepted_lengths"] == [6, 8]
    assert payload["rejected_alias_lengths"] == [4]
    assert payload["stable_relations"]["tau_w_e5_per_cell"] == [
        [9, 2],
        [6, 1],
        [15, 2],
    ]
    assert payload["claim"]["leading_order_status"] == "no_go"
    assert payload["claim"]["finite_step_status"] == "inconclusive"
    assert payload["hpc_authorized"] is False


@pytest.mark.parametrize(
    ("path", "value", "message"),
    (
        (("records", 0, "moments", "tau_h2"), [1, 1], "moment"),
        (("accepted_lengths", 0), 10, "length"),
        (("stable_relations", "tau_w_e5_per_cell", 0), [0, 1], "stable"),
        (("claim", "finite_step_status"), "proved", "finite-step"),
        (("implementation_sources", "source.py"), "b" * 64, "source"),
    ),
)
def test_extensive_payload_rejects_mutations(path, value, message) -> None:
    payload = build_extensive_payload(
        (_record(6), _record(8)),
        rejected_alias_lengths=(4,),
        implementation_sources={"source.py": "a" * 64},
    )
    forged = copy.deepcopy(payload)
    target = forged
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError, match=message):
        verify_extensive_payload(forged)


def test_extensive_payload_rejects_inconsistent_size_family() -> None:
    inconsistent = _record(8)
    inconsistent = replace(
        inconsistent,
        moments=replace(
            inconsistent.moments,
            tau_w_e5=inconsistent.moments.tau_w_e5 + Cubic.one(),
        ),
    )
    with pytest.raises(ValueError, match="stable"):
        build_extensive_payload(
            (_record(6), inconsistent),
            rejected_alias_lengths=(4,),
            implementation_sources={"source.py": "a" * 64},
        )
