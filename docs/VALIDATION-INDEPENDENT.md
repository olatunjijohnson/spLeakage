# Real-data validation against true independent error

Does spLeakage predict the **real** optimism gap — not just simulated optimism?
Using three **real spatial fields** (meuse log-zinc, Paraná rainfall, ca20 calcium;
`data-raw/validate_independent.R`), each replicate trains on part of the data and
evaluates on a **genuinely independent held-out set** of real observations. Two
deployment scenarios span the spectrum:

- **extrapolation** — hold out a spatially coherent region (predict into unsampled
  space, the canonical leakage scenario);
- **interpolation** — hold out a random interspersed subset (deployment ≈ CV).

Because the fields are real (not Gaussian simulations), this is an out-of-
(simulation)-distribution test of the whole approach. 240 replicates.

## Results

**The problem is real and large.** Under genuine extrapolation deployment, naive
random 10-fold CV **understates the true independent error by a median of 36%**
(`E_cv / E_indep = 0.64`). A practitioner reporting random-CV accuracy on real
clustered data is materially overconfident.

**The SLI tracks real optimism and separates the scenarios.**

| | extrapolation | interpolation |
|---|---|---|
| mean **real** optimism | **+0.35** | −0.04 |
| mean **SLI_rho** | **+0.28** | −0.01 |

The SLI (given the true deployment target) correctly distinguishes "deploying to a
new region → 35% optimistic" from "interpolation deployment → not optimistic", and
its values mirror the real optimism. Across all 240 replicates,
`cor(SLI_rho, real optimism) = +0.52` — strong for real, noisy, n≈150 datasets
(cf. r ≈ 0.97 in the clean GP simulation).

**The de-leaked estimate corrects a meaningful fraction — and a target-aware control
corrects most of it.** Under extrapolation (`E_cv / E_indep = 0.62`, naive understates
true error by 38%):

| estimate | `E / E_indep` | gap closed |
|---|---|---|
| naive random CV | 0.62 | — |
| de-leak (fixed block) | 0.80 | 58% |
| **de-leak (target-aware, `target=`)** | **0.89** | **71%** |

Passing the deployment target to `deleak_estimate(target=)` selects the control
scheme whose `SLI_rho` against the target is closest to zero (NNDM matching as a
selection rule), recovering **89%** of the true independent error — without ever
seeing the independent data.

## Remaining limitations (and what they teach)

- Even the target-aware de-leak does not fully close the gap (71%, not 100%) for the
  hardest extrapolation: a finite sample cannot fully diagnose performance arbitrarily
  far beyond its footprint — this is exactly where **area-of-applicability flagging**
  (abstain rather than extrapolate) is the honest answer.
- The earlier finding that *fixed* block control is target-agnostic motivated the now-
  implemented `target=` argument (turn #23); it is what lifts gap closure from 58% to
  71%.

## Takeaway for the paper

On real fields, naive random CV understates true extrapolation error by ~36%, and the
target-aware SLI predicts this (separating interpolation from extrapolation deployment
with `cor = 0.52`). The de-leaked estimate corrects ~40% of the gap. This is the
real-data counterpart to the GP theory (`docs/THEORY-RESULTS.md`) and the NNDM
benchmark, and it honestly bounds the method's current reach: **target-matched
control** is the next improvement needed to fully recover true extrapolation accuracy.
