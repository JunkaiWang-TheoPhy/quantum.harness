# PF4--TFIM exact BCH and endpoint-conjugation ledger

## Certified statement

For the two-fragment five-copy fourth-order Suzuki formula, let

\[
\log S_4(s)=s(A+B)+s^5E_5+O(s^7),
\]

with exact stage weights in
`Q(alpha)`, `alpha^3=4`, at the positive real embedding.  Define

\[
C=[A,[A,B]],\qquad D=[B,[B,A]].
\]

The exact cyclic free-trace identity is

\[
\operatorname{Tr}((A+B)E_5)
=\gamma\operatorname{Tr}\!\left(C^2-4CD+\frac83D^2\right),
\]

where

\[
\gamma=
\frac{37}{900000}
+\frac{313}{14400000}\alpha
+\frac{29}{1800000}\alpha^2.
\]

All three rational coordinates are positive and `alpha>0`, hence
`gamma>0` without floating-point inference.

## Cyclic free-word proof

`src/trottercert/pf4_bch_mapping.py` expands the exact Suzuki logarithm through
degree five.  Degrees two, three, and four vanish exactly.  After multiplying
the degree-five term by `A+B`, both sides of the claimed identity reduce to the
following ten cyclic degree-six classes.  Each displayed entry is the exact
multiple of `gamma`.

| Cyclic word | Coefficient divided by `gamma` |
|---|---:|
| `000011` | `2` |
| `000101` | `-8` |
| `000111` | `-8` |
| `001001` | `6` |
| `001011` | `12` |
| `001101` | `12` |
| `001111` | `16/3` |
| `010101` | `-16` |
| `010111` | `-64/3` |
| `011011` | `16` |

Equality of this full mapping is the authoritative proof.  Seeded exact
rational symmetric matrices in dimensions two and three independently verify
the evaluated trace identity.  Tests mutate each of the three quadratic
coefficients and require exact class inequality.

## Physical-unitary convention

Substituting the anti-Hermitian generators `-iA` and `-iB` into a homogeneous
degree-five Lie polynomial contributes `(-i)^5=-i`.  Multiplication of the
principal logarithm by `i/s` therefore produces the same Hermitian `E5`
coefficient in the effective Hamiltonian.  No sign or factor is inferred from
a numerical logarithm.

## Periodic TFIM specialization

For every even periodic length `L >= 4` and exact rational couplings,

\[
A=h\sum_iX_i,\qquad B=j\sum_i Z_iZ_{i+1}.
\]

Exact Pauli orthogonality and local multiplicities give

\[
\frac{\operatorname{Tr}(C^2)}d=128Lh^4j^2,
\quad
\frac{\operatorname{Tr}(CD)}d=0,
\quad
\frac{\operatorname{Tr}(D^2)}d=128Lh^2j^4.
\]

Consequently

\[
\frac{\operatorname{Tr}(HE_5)}d
=128L\gamma h^2j^2\left(h^2+\frac83j^2\right).
\]

This element is strictly positive exactly when `h*j != 0`.  If either
coupling is zero, the two fragments commute and the pairing vanishes.

## Endpoint-conjugation lower bound

For every Hermitian endpoint generator `Q`, cyclicity gives

\[
\tau\!\left(H\,i[Q,H]\right)=0,
\]

where `tau=Tr/d`.  For `X=E5+i[Q,H]`, normalized trace duality implies

\[
|\tau(H E_5)|=|\tau(HX)|
\le \tau(|H|)\,\|X\|_\infty
\le L(|h|+|j|)\,\|X\|_\infty.
\]

The artifact therefore stores the exact algebraic coefficient

\[
\inf_Q\|E_5+i[Q,H]\|_\infty
\ge
\frac{128\gamma |h^2j^2(h^2+8j^2/3)|}{|h|+|j|},
\]

with zero defined at the zero Hamiltonian.  This is a leading-order
coefficient; the physical small-step contribution carries the corresponding
fourth power of the step size in the effective generator.

## Claim boundary

The schema-v2 certificate proves:

- the exact Suzuki PF4 cyclic BCH/free-trace identity;
- the all-even-length periodic TFIM trace-moment formulas;
- a strict leading-order obstruction to pure endpoint conjugation when both
  couplings are nonzero;
- an exact algebraic operator-norm lower-bound coefficient at fixed time and
  fixed Hamiltonian normalization.

It does not prove:

- an obstruction after quotienting arbitrary `a I + b H` directions;
- a finite-step TFIM no-go;
- a total-time eigenphase lower bound modulo `2 pi`;
- any part of the independent Heisenberg E9 finite-step certificate.

The old research-plan form with coefficients `1/2, 14/3, 4/3` was never
connected to the actual Suzuki BCH defect.  It is superseded and retained only
in a negative mutation test that proves the canonical verifier rejects it.

## Machine evidence

```text
artifact:
docs/experiments/processor-obstruction/tfim-family-obstruction.json

file sha256:
a1f628104f816ea2cd35c7846b64c1205df76fae73eecf2b82cc69a1e9f789e4

payload sha256:
ffd5bf98bf40138f221dde86b990eb43e97e4e101332bd898465bcbae53c09d3

PF4 identity digest:
bb77fd1dfc2f1289e3695dcfc1791b57b3b11d2a1366a7c1b28eeaa1b1a7975c
```

The artifact binds these implementation sources by SHA-256:

- `scripts/certify_tfim_obstruction.py`
- `src/trottercert/algebra.py`
- `src/trottercert/cubic_field.py`
- `src/trottercert/cubic_local.py`
- `src/trottercert/intervals.py`
- `src/trottercert/pf4_bch_mapping.py`
- `src/trottercert/tfim_obstruction.py`
- `src/trottercert/trace_obstruction.py`

From `tracks/qcs/solutions/WangTheoPhys/issue128`, verify with

```bash
PYTHONPATH=src:. python scripts/certify_tfim_obstruction.py \
  --verify docs/experiments/processor-obstruction/tfim-family-obstruction.json

PYTHONPATH=src:. python -m pytest -q \
  tests/test_pf4_bch_mapping.py \
  tests/test_tfim_obstruction.py
```

The observed focused gate is `42 passed`.  Ruff reports no findings for the
mapping, trace, TFIM, certificate, and focused-test files.
