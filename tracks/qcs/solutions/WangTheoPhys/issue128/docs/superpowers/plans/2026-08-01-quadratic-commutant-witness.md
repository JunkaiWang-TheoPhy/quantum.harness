# Quadratic Commutant Witness Implementation Plan

**Goal:** Certify an independent polynomial commutant witness for the 12-by-12
Issue-128 leading defect and promote only the justified calibrated leading-order
claim.

- [ ] Add failing tests for exact unit-cell lifting, Hamiltonian reconstruction,
  rational `H^2` moments, cubic-field serialization, and fail-closed payload
  mutation checks.
- [ ] Implement the sparse exact torus and moment kernel in
  `src/trottercert/commutant_witness.py`.
- [ ] Add the deterministic generation/full-verification CLI and generate the
  12-by-12 witness artifact.
- [ ] Update the gauge-aware audit to bind the witness artifact, report
  `leading_order_no_go`, and retain `finite_step_status = inconclusive`.
- [ ] Run focused tests, the explicit full certificate recomputation, and the
  complete non-slow Issue-128 suite.
- [ ] Mark this plan complete and commit only the scoped files.
