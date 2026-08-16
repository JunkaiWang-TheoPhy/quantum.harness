<!-- wander-artwork -->
![Wander — Issue #276](https://raw.githubusercontent.com/JunkaiWang-TheoPhy/quantum.harness/664506886413b68cdbc5a362a6075c78fd5d2c46/docs/showcase/ranger-archive/assets/missions-v2/07-quantum-geometry-curved-state-space-v2.png)

> **EN** The spectrum is silent; the geometry is not.
>
> **中文** 能谱沉默，几何仍然说话。
<!-- /wander-artwork -->

## Team

| Field | Value |
|---|---|
| **Team name** | Wander |
| **Members** | Chenxi Wan, Yedi Shen, Junkai Wang |
| **Contact email** | WangTheoPhys@outlook.com |

## Challenge

Addresses #276: **What probes quantum chaos inside an exactly degenerate eigenspace?**

The executable answer is the non-Abelian quantum geometry of the protected projector over coupling space. For a degenerate projector $P(\lambda)$, the response $X_a=(1-P)\partial_aP$ generates

$$\mathcal Q_{ab}=X_a^\dagger X_b,\qquad g_{ab}=\tfrac12(\mathcal Q_{ab}+\mathcal Q_{ba}),\qquad F_{ab}=i(\mathcal Q_{ab}-\mathcal Q_{ba}).$$

Exact degeneracy silences internal level statistics, while the quantum metric, Berry curvature, gauge-invariant four-channel tensors, Chern numbers, and Wilson holonomy remain defined.

## Results

The submission now covers five protection mechanisms or controls rather than one FQH model.

| Mechanism | Exact-degeneracy origin | Result |
|---|---|---|
| Kapit–Mueller/Laughlin | Frustration-free two-body parent | Structured topological benchmark, fixed-Chern holonomy, and one-sided response |
| Generic cubic $\mathcal N=2$ SYK | Charge-resolved harmonic cohomology | Sealed $N=14$ rejection of both registered separable covariance nulls |
| Continuum Moore–Read | Three-body clustering constraints | Exact ranks $D=42,120$, open gaps, and positive complete-covariance statistic at $N=4,6$ |
| Lattice SUSY on $C_6^{\sqcup m}$ | Independence-complex cohomology | Exact ranks $2,4,8$ and nonzero response, but no registered positive transfer at $m=1,2,3$ |
| Periodic X-cube | Commuting stabilizer constraints | Exact nonchaotic controls: $F=0$ or scalar $F=-I/2$, both with zero connected curvature variance |

### Sealed cohomological test

The complete $N=8,10,12$ cubic-SYK pilot and separately sealed held-out $N=14$ sparse pair use complete disorder realizations as the uncertainty unit. The prediction seal was checked before explicit unsealing.

| $N=14$ sector | Physical median (95% bootstrap) | Collapsed null (97.5% prediction) | Hodge null (97.5% prediction) |
|---|---:|---:|---:|
| Adjacent | 0.301529 [0.291527, 0.312061] | [0.111789, 0.111852] | [0.112344, 0.112513] |
| Central | 0.374993 [0.368980, 0.380473] | [0.111338, 0.111353] | [0.111333, 0.111348] |

This selects `cohomological_non_gaussian_class` relative to the frozen collapsed and Hodge-separable covariance models.

### Mechanism-dependent v12 test

For Moore–Read, the complete three-pairing U-statistic over all 276 unordered pairs among 24 local tangent panels is $0.1908758\pm0.0277690$ at $(N,D)=(4,42)$ and $0.1936048\pm0.0248286$ at $(6,120)$, with one-sided $p=3.13\times10^{-12}$ and $3.15\times10^{-15}$. The same fixed direction does not pass for lattice SUSY at $m=1,2,3$.

The combined branch is `domain_limited_geometric_eth`. The positive Moore–Read inference is conditional on exchangeability of local tangent panels at a fixed base Hamiltonian. The X-cube calculation proves that even a moving projector can have completely nonchaotic scalar geometry.

## Claim Boundary

- Established: multiple independent exact-degeneracy mechanisms; exact Hodge, clustering, cohomology, and stabilizer structures; finite-size geometric response memory; and exact counterexamples to “moving projector implies chaos.”
- Not established: asymptotic Geometric ETH, an independent Moore–Read Hamiltonian/disorder ensemble, a universal cross-mechanism limiting law, thermalization, real-time chaos, or a thermodynamic-limit theorem.
- The SYK result rejects the frozen separable covariance nulls, not every possible nonseparable entrywise Gaussian covariance model.

## Deliverables

- [Harness solution overview](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/README.md)
- [Research and reproduction guide](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/research/README.md)
- [Mechanism-dependent main article](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/research/01_task_folder/task_05/script/output/mechanism_dependent_geometric_eth_v12.pdf)
- [Mechanism-dependent Supplemental Material](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/research/01_task_folder/task_05/script/output/mechanism_dependent_geometric_eth_supplement_v12.pdf)
- [Machine-readable v12 inference](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/research/01_task_folder/task_05/script/output/cross_mechanism_geometric_eth_v12.json)
- [Cubic-SYK Letter](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/research/01_task_folder/task_05/script/output/response_complex_memory_v7.pdf)
- [Cubic-SYK Supplemental Material](https://github.com/JunkaiWang-TheoPhy/quantum.harness/blob/codex/issue-276-quantum-geometry/tracks/mps/solutions/Wander-276/research/01_task_folder/task_05/script/output/response_complex_memory_supplement_v7.pdf)
- [Public v12 research branch](https://github.com/JunkaiWang-TheoPhy/Chaos-of-Quantum-Geometry/tree/codex/geometric-eth-mechanisms-v12)

## Verification

```bash
cd tracks/mps/solutions/Wander-276
bash verify.sh
```

The verifier checks the existing outcome-blind v7 chain and the new Moore–Read, lattice-SUSY, X-cube, complete-covariance, artifact-hash, and manuscript gates.

@OkongOyangO, please review the exact-degeneracy interpretation, mechanism comparison, and claim boundary as Wander's Issue #276 submission.
