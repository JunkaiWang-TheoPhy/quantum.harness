# Processed spectral discovery screen

## Decision

The deterministic four-fragment dense diagnostic passes for all three local
candidates.  It reproduces the defining separation: the raw kernel has lower
operator order, while its eigenphases converge at effective order six.

| kernel | repeated operator slope | repeated phase slope |
|---|---:|---:|
| `s8` | `2.00151769` | `5.99184313` |
| `s10` | `2.00155973` | `5.99712310` |
| `s11` | `4.00258240` | `5.95571683` |

The screen uses four seeded noncommuting `8 x 8` Hermitian fragments, each
normalized to spectral norm `0.27`, and repetitions
`r = 2,3,4,5,6,8,10,12`.  Phase-slope fitting excludes values within
`64 * eps(binary64)` of the eigensolver floor; all raw errors remain in the
JSON artifact.

## Reproduction

From the issue-128 directory:

```bash
PYTHONPATH=src python -u -m scripts.screen_processed_spectra \
  --output docs/experiments/processed-spectral-screen/random-four-fragment.json
python -m pytest -q tests/test_processed_spectral_screen.py
```

## Claim boundary

This is a non-rigorous discovery diagnostic.  It checks the Eq. (1.8) stage
order, phase matching, binary64 floor handling, and qualitative candidate
behavior.  It is not a locality bound, does not scale the dense calculation
to the periodic `12 x 12` benchmark, and cannot be substituted for the
processed finite-step spectral certificate or the exact physical E7/D6
sidecars.
