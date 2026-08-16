# Wander — Issue #276: Mechanism-Dependent Geometric ETH

![Quantum geometry over an exactly degenerate manifold](https://raw.githubusercontent.com/JunkaiWang-TheoPhy/quantum.harness/664506886413b68cdbc5a362a6075c78fd5d2c46/docs/showcase/ranger-archive/assets/missions-v2/07-quantum-geometry-curved-state-space-v2.png)

> **EN** The spectrum is silent; the geometry is not.
>
> **中文** 能谱沉默，几何仍然说话。

## Team

| Field | Value |
|---|---|
| **Team name** | Wander |
| **Members** | Chenxi Wan, Yedi Shen, Junkai Wang |
| **Contact email** | WangTheoPhys@outlook.com |

## Challenge and Answer

This submission addresses [QuantumBFS/quantum.harness#276](https://github.com/QuantumBFS/quantum.harness/issues/276): **What probes quantum chaos inside an exactly degenerate eigenspace?**

The executable answer is the **non-Abelian quantum geometry of the protected projector over coupling space**. If (P(\lambda)) is the degenerate projector and (X_a=(1-P)\partial_aP), then

$$\mathcal Q_{ab}=X_a^\dagger X_b,\qquad g_{ab}=\frac{\mathcal Q_{ab}+\mathcal Q_{ba}}{2},\qquad F_{ab}=i(\mathcal Q_{ab}-\mathcal Q_{ba}).$$

Ordinary level statistics is undefined inside a flat multiplet, but (g), (F), gauge-invariant four-channel contractions, Chern numbers, and Wilson holonomy remain nontrivial.

## Independent Protection Mechanisms

The release no longer rests on one fractional-quantum-Hall example.

| Mechanism | Protected manifold | Exact response structure | Role |
|---|---|---|---|
| Kapit–Mueller/Laughlin parent (H=B^\dagger B) | Frustration-free, gapped zero modes | One-sided response (X=-H_\perp^+B^\dagger\delta BP) | Local topological benchmark, Jacobi curvature, fixed-Chern holonomy |
| Generic cubic (\mathcal N=2) SYK (H=\{Q,Q^\dagger\}) | Charge-resolved harmonic cohomology | Orthogonal exact/coexact response (X=X_-\oplus X_+) | Independent supersymmetric mechanism and sealed covariance test |
| Continuum Moore–Read parent | Three-body clustering zero modes | Fixed-four-quasihole ranks (D=42,120) | Positive complete-covariance statistic at (N=4,6) |
| Lattice supersymmetry on (C_6^{\sqcup m}) | Independence-complex cohomology | Local exact/coexact response, rank (2^m) | Independent local cohomological transfer test |
| Periodic X-cube code | Commuting stabilizer constraints | Exact zero or scalar curvature under two deformations | Counterexample: moving projectors need not be geometrically chaotic |

Nilpotency gives (X_-^\dagger X_+=0) exactly. The two-sided Hodge decomposition is therefore a physical response-complex structure, not a basis rewrite of the Laughlin calculation.

## Registered Scientific Test

For eight coupling tangents, the code whitens the channel covariance and evaluates the gauge-invariant tensor

$$\mathcal T_{abcd}=\frac1D\operatorname{Tr}(\widehat X_a^\dagger\widehat X_b\widehat X_c^\dagger\widehat X_d).$$

The (N=8,10,12) pilot uses complete disorder realizations as the uncertainty unit. A held-out (N=14) calculation freezes two prediction models before any four-channel outcome is opened:

1. a collapsed separable covariance Gaussian null;
2. an exact/coexact Hodge-resolved separable covariance Gaussian null.

Safe covariates and predictions are SHA-256 sealed. Unsealing is an explicit, separate command and is never launched by the scheduler dependency chain. The seal was independently checked before opening, and the frozen result is:

| $N=14$ sparse sector | Physical median (95% bootstrap) | Collapsed null (97.5% prediction) | Hodge null (97.5% prediction) |
|---|---:|---:|---:|
| Adjacent | 0.301529 [0.291527, 0.312061] | [0.111789, 0.111852] | [0.112344, 0.112513] |
| Central | 0.374993 [0.368980, 0.380473] | [0.111338, 0.111353] | [0.111333, 0.111348] |

Both nulls miss both primary sectors, selecting `cohomological_non_gaussian_class` for the SYK branch.

## Cross-Mechanism Result

The v12 comparison applies one projector-response convention and a covariance-complete three-pairing U-statistic to clustering, local-cohomology, and stabilizer mechanisms. Moore–Read gives positive fixed-direction estimates (0.1908758\pm0.0277690) at (N=4,D=42) and (0.1936048\pm0.0248286) at (N=6,D=120), with one-sided (p=3.13\times10^{-12}) and (3.15\times10^{-15}). The same registered direction does not pass for lattice supersymmetry at (m=1,2,3). X-cube coefficient changes give (g=F=0), while a local isospectral orbit gives scalar (F=-I/2) and exactly zero connected curvature variance.

The supported v12 branch is `domain_limited_geometric_eth`: the positive Moore–Read result is conditional on exchangeability of 24 local tangent panels at one base Hamiltonian. It is not an independent Hamiltonian ensemble or an asymptotic law.

## Claim Boundary

This submission establishes multiple independent exact-degeneracy mechanisms and finite-size geometric response memory. It also establishes an exact negative control: projector motion alone is insufficient for chaos. It does **not** claim asymptotic Geometric ETH, a universal cross-mechanism law, conventional energy-resolved ETH, real-time chaos, or a thermodynamic-limit theorem. The Moore–Read inference is conditional on local-panel exchangeability, and the SYK null rejection does not by itself prove intrinsic non-Gaussianity after matching every entrywise nonseparable covariance.

## Deliverables

- [Research and reproduction guide](research/README.md)
- [Task-level guide](research/01_task_folder/task_05/README.md)
- [Letter source](research/overleaf_sync/cohomological_geometric_eth/main.tex)
- [Supplement source](research/overleaf_sync/cohomological_geometric_eth/supplement.tex)
- [Compiled Letter](research/01_task_folder/task_05/script/output/response_complex_memory_v7.pdf)
- [Compiled Supplemental Material](research/01_task_folder/task_05/script/output/response_complex_memory_supplement_v7.pdf)
- [Mechanism-dependent main article](research/01_task_folder/task_05/script/output/mechanism_dependent_geometric_eth_v12.pdf)
- [Mechanism-dependent Supplemental Material](research/01_task_folder/task_05/script/output/mechanism_dependent_geometric_eth_supplement_v12.pdf)
- [v12 cross-mechanism evidence](research/01_task_folder/task_05/script/output/cross_mechanism_geometric_eth_v12.json)
- [Exact-data result report](research/01_task_folder/task_05/script/output/susy_hodge_geometric_eth_report_v7.md)
- [Scientific ceiling and literature memo](research/docs/2026-08-01-scientific-ceiling-strategy.md)
- [Public v12 research branch](https://github.com/JunkaiWang-TheoPhy/Chaos-of-Quantum-Geometry/tree/codex/geometric-eth-mechanisms-v12)

## Verification

```bash
cd tracks/mps/solutions/Wander-276
bash verify.sh
```

The verification path checks the exact Hodge identity, analytic controls, pilot-grid completeness, outcome-blind sealing, Moore–Read ranks and gaps, lattice-SUSY cohomology, X-cube identities, complete-covariance inference, figure provenance, and both audited paper packages.

## Next Scientific Gate

The finite-size result is designed to expose the next theorem rather than hide it. The two high-ceiling targets are:

- an asymptotic concentration/scaling law for the covariance-controlled response tensor;
- a spatially local nilpotent-supercharge model with a stable zero-mode count, open gap, and moving cohomological projector.

These are journal-ceiling gates, not assumptions built into the current claim.
