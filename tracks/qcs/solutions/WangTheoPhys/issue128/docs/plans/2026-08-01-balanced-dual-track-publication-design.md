# Issue 128 balanced dual-track publication design

**Date:** 2026-08-01
**Status:** Approved
**Program:** Quantum-ready Paper A plus gated PRX Quantum Paper B sprint

## 1. Decision

Adopt a balanced dual-track program:

- allocate approximately 55% of near-term effort to Paper A and 45% to Paper B;
- produce a Quantum-ready Paper A in 10--14 weeks;
- evaluate Paper B against explicit theorem and finite-step gates over 4--5 months;
- increase Paper B to approximately 65% of effort only after its family-theorem
  and finite-step gates pass;
- never delay Paper A for D8, a fivefold resource claim, or an unfinished Paper B
  theorem.

The current certified result remains the principal protected asset:

\[
393\to95\text{ steps},\qquad
11791\to2851\text{ merged groups},\qquad
\frac{11791}{2851}=4.135741844966678\ldots.
\]

Relative to a strengthened approximately 353-step control, the same candidate
corresponds to an approximately 3.71-fold reduction.  These two ratios must
always be labeled by their distinct denominators.

## 2. Program architecture

The two papers share exact-algebra and provenance infrastructure but have
separate scientific claims, certificates, manuscripts, and completion gates.

```text
Exact word and Lie algebra
        |
        +-- Paper A: proof-carrying upper bounds
        |      +-- D4/D5/D6 data
        |      +-- finite-step error ledger
        |      +-- independent verifier
        |
        +-- Paper B: spectral lower bounds and obstructions
               +-- primal processors
               +-- commutant witnesses
               +-- E5/E7/E9 dual contractions
               +-- E11-and-higher analytic envelope
```

Shared components may include exact word arithmetic, cubic-field arithmetic,
Pauli canonicalization, rational intervals, manifest formats, hashing, and test
conventions.  Neither paper may cite an exploratory artifact as a frozen result
solely because it uses shared code.

Before moving files, every existing file must be classified as `paper-a`,
`paper-b`, `shared`, `exploratory`, or `obsolete`.  The current dirty worktree
must not be reorganized with bulk moves until that ownership manifest has been
reviewed.

## 3. Paper A: proof-carrying finite-step simulation

### 3.1 Target and title

Primary target: **Quantum**.  PR Research or PRA are transfer options, not the
design target.

Working title:

> Proof-Carrying Finite-Step Product-Formula Simulation of Local Quantum
> Many-Body Systems

### 3.2 Central claim

The paper introduces a proof-carrying compiler that turns a supported local
Pauli Hamiltonian and product formula into a finite-step error and resource
certificate that an independent implementation can reject or accept.

The 12 by 12 periodic Heisenberg calculation is the flagship instance.  The
paper must not claim that 95 is the physically necessary step count, that the
4.1357-fold ratio applies against every modern simulation method, or that an
unfrozen 94-step candidate is established.

### 3.3 Main theorems

#### A1. Certificate soundness

For every supported input, acceptance must imply the stated global error bound:

\[
\operatorname{Verify}(H,S_r,C)=\operatorname{accept}
\Longrightarrow
\|e^{-iHT}-S_r(T)\|_\infty\le\epsilon_C.
\]

The proof must cover formal-log coefficients, homogeneous Lie defects, Pauli
evaluation, exact aggregation, grouped norm certificates, the finite-step
ledger, the rational tail, and resource arithmetic.

#### A2. Deferred-norming dominance

Prove the validity and monotonicity of:

1. combining identical Pauli strings before absolute values;
2. aggregating translation orbits before norm inequalities;
3. using certified pairwise-anticommuting groups;
4. excluding discovery heuristics from the trusted computing base.

#### A3. Finite-lattice transfer

State explicit conditions for local support, periodic wraparound, orbit
multiplicity, aliasing, and the transfer from unit-cell densities to finite
lattices.  Include a deliberately aliased small case that the implementation
rejects or treats separately.

#### A4. Fail-closed verification

The verifier must reject drift in coefficients, coverage, group metadata,
integer types, tails, adjacent-step values, resource counts, ratios, hashes,
and anticommutation claims.

### 3.4 Evidence matrix

Use a deliberately small but meaningful generalization matrix:

| Role | Model or system | Purpose |
|---|---|---|
| Flagship | 12 by 12 isotropic Heisenberg | Large strict resource certificate |
| Parameter family | XXZ at a preregistered set of anisotropies | Test dependence on SU(2) symmetry |
| Held-out model | Small TFIM or generic XYZ | Test transfer after freezing rules |
| Exact calibration | 2 by 3, 2 by 4, and 3 by 3 | Actual-error versus certificate gap |

Report certified steps, merged groups, builder and verifier time, peak memory,
certificate size, baseline definitions, actual/certified error gap when exact
diagonalization is feasible, and a full deferred-norming ablation.

### 3.5 Figures

1. Compiler and independent-verifier data flow.
2. Ablation from the published bound to the final finite-step certificate.
3. Actual error, new certificate, and published bound versus step count.
4. Model transfer, certificate size, and verification cost.

### 3.6 Paper A completion gate

Paper A is complete only when all of the following hold:

- a unique 95-step authoritative release exists;
- A1--A4 are in reviewable form;
- a clean clone reproduces primary, deep, and reference verification;
- nondegenerate exact-diagonalization checks are archived;
- the XXZ family and one held-out model are reported;
- a mutation corpus covers all advertised failure modes;
- every manuscript number is generated or validated from frozen data;
- source, certificate, data, and reproduction transcripts have a DOI;
- the manuscript has been submitted to Quantum.

Paper A does not wait for D8, fivefold performance, or Paper B finite-step
closure.

## 4. Paper B: spectral obstructions to low-cost correction

### 4.1 Target and title

Primary target after all gates pass: **PRX Quantum**.

Working title:

> Spectral Obstructions to Additive-Cost Corrected Product Formulas

### 4.2 Gauge and central question

Define the allowed leading-order gauge

\[
\mathcal G_H=operatorname{Im}(i\,\operatorname{ad}_H)
+\operatorname{span}\{I,H\}.
\]

The three directions represent basis rotation, global phase, and time or energy
rescaling.  The paper studies when a defect is nonzero in
\(\mathcal A/\mathcal G_H\), and what circuit or repetition cost is required to
change that spectral class.

### 4.3 Main theorems

#### B1. Gauge-aware correctability duality

Give a finite-dimensional quotient or dual characterization of

\[
\inf_{Q,a,b}\|E+i[Q,H]+aI+bH\|_\infty
\]

using trace-norm-bounded commutant witnesses satisfying

\[
[W,H]=0,\qquad \operatorname{Tr}W=0,\qquad \operatorname{Tr}(WH)=0.
\]

The statement must handle degenerate energy blocks.

#### B2. Telescoping-cost trilemma

Formalize that conjugation processors telescope but do not change the kernel
spectrum, while spectrum-changing nonconjugate corrections generally repeat at
each step.  Derive a quantitative tradeoff between spectral residual and
failure to commute or telescope.

#### B3. PF2 universal obstruction or complete exception classification

For genuinely noncommuting two-fragment Hamiltonians, prove a universal
leading spectral obstruction, or characterize every exception explicitly.

#### B4. PF4 symmetry-family theorem

Derive and verify the exact trace quadratic form, and prove a strict obstruction
for an appropriate symmetry-protected family rather than only at the isotropic
Heisenberg point.

#### B5. TFIM scaling optimality

Prove for a TFIM parameter family that pure symplectic endpoint correctors
cannot improve the structural \(\alpha^2\) residual scaling.  Coefficients and
finite-size counting identities must be exact and independently checkable.

#### B6. Finite-step spectral lower bound

Combine exact E5, E7, and E9 witness contractions with an independently
verifiable E11-and-higher envelope.  For the effective local-log generator,
obtained by dividing the one-step logarithm by \(h\), promotion requires the
strictly positive signed margin

\[
|h^4q_5+h^6q_7+h^8q_9|-|R_{\ge11}^{\rm dual}|>0.
\]

The un-divided one-step logarithm has the corresponding odd powers; the two
normalizations must not be mixed in either the certificate or the paper.

An exact E9 result alone never promotes the finite-step claim.

### 4.4 Sharpness

Include both:

- at least one family with a certified obstruction;
- at least one exactly or rigorously correctable positive example.

The intended positive seed is the tuned one-qubit example, with a local
extensive version added if it remains within scope.  The theory should act as a
classification rather than an indiscriminate no-go statement.

### 4.5 Figures

1. Corrector taxonomy and telescoping costs.
2. Quotient geometry with primal processors and dual witnesses.
3. TFIM scaling floor and finite-step crossover.
4. Correctability-cost frontier with positive and negative examples.

### 4.6 Paper B promotion gates

Paper B becomes a PRX Quantum submission only if it has:

- B1 and B2;
- at least one general family theorem from B3 or B4;
- B5;
- B6 with a positive finite-step margin;
- a positive correctable example;
- an exact verifier;
- fair comparison with corrected and processed product formulas.

If B1--B5 hold but B6 remains inconclusive, target Quantum or PR Research and
do not delay Paper A.  If only a Heisenberg leading-order witness survives, use
it as a strong Paper A limitation or a narrower follow-up rather than claiming
a flagship theory.

## 5. Twelve-week execution design

### Weeks 1--2: baseline and theorem skeletons

- classify all current files by paper ownership;
- create a Paper A claim-to-artifact matrix;
- reproduce the 95-step bundle in a clean clone;
- freeze the Paper B gauge and claim boundary;
- draft A1's theorem dependency graph;
- complete the finite-dimensional form of B1 on exact small examples.

### Weeks 3--5: Quantum evidence

- complete A2 and A3;
- add exhaustive wraparound and aliasing tests;
- produce nondegenerate exact-diagonalization calibration data;
- derive PF2 and PF4 exact identities independently;
- run the preregistered XXZ family;
- freeze and run the held-out model;
- derive and verify the TFIM counting formulas.

### Weeks 6--8: ablation and finite-step obstruction

- produce a computed Paper A ablation table;
- launch exact E9 only after all manifest and regression gates pass;
- recover 64 of 64 shards, independently rerun selected shards, and verify the
  reduced E9 artifact in a clean checkout;
- combine E9 with the centered dual-log E11-and-higher envelope;
- record a positive, negative, or interval-crossing signed margin without
  overriding an inconclusive result.

### Weeks 9--12: submission and sprint decision

- freeze Paper A prose, figures, data, and release bundle;
- complete external clean-clone reproduction;
- perform a structured mock review;
- mint the Paper A DOI and submit to Quantum;
- if B6 passes, raise Paper B effort to approximately 65% and begin the PRX
  Quantum sprint;
- otherwise keep Paper B at leading-order scope and redirect effort to the
  strongest valid journal target.

## 6. HPC policy

E9 production is permitted only after:

- an immutable source commit exists;
- 64 manifests and their index verify locally;
- E5 and E7 manifest paths reproduce frozen pairings;
- the reachability prune agrees exactly with unpruned regression fixtures;
- load balancing covers all inner patterns in every shard;
- Slurm scripts pass static checks;
- a new explicit run root and recovery plan are recorded.

The run must preserve the manifest index, all worker outputs, scheduler
accounting, the reduced artifact, hashes, and at least two independent shard
reruns.  E9 computation is not authorization to start D8.  D8 remains paused
unless a separate integer-step or tail-replacement decision gate is approved.

## 7. Failure handling and scope control

- If D6 does not change an integer step, retain it as scaling infrastructure
  and do not delay Paper A.
- If E9 production fails operationally, Paper A continues.
- If the finite-step margin crosses zero, keep the claim inconclusive.
- Permit only one preregistered tightening of an identifiable dominant tail
  slack; do not enumerate higher degrees without a feasibility calculation.
- If XXZ or the held-out model shows weak improvement, publish the boundary
  rather than deleting the result.
- If a verifier mutation succeeds, stop release, repair the schema or trusted
  path, regenerate manifests, and rerun all affected evidence.
- Never use a test count, green manifest, or hash as evidence for a claim whose
  relevant semantics that check does not cover.

## 8. Success map

| Valid completed package | Primary target |
|---|---|
| Paper A completion gate | Quantum |
| Paper A without convincing transfer evidence | PR Research or PRA transfer |
| B1--B5 without finite-step closure | Quantum or PR Research |
| B1--B6 plus sharpness and exact verification | PRX Quantum |
| Above plus a general gate-count, range, or depth lower bound | Consider PRL |

The program is successful only if it protects a real Quantum submission while
testing the higher-risk PRX Quantum thesis against explicit gates.  Moving a
resource ratio by a small amount is not a substitute for either completion
definition.
