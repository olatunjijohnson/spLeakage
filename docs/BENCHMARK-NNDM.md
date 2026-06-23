# Benchmark — signed SLI vs NNDM's unsigned W statistic

The most likely reviewer objection is "why not just use NNDM/kNNDM's Wasserstein `W`
statistic?" This benchmark answers it. Ground truth is the **exact GP optimism**
(excess explained variance, `docs/THEORY.md` eq. 2), computed from covariance algebra
across 288 configs spanning sampling design × split (random/block) × target
(grid/interpolation) × range × signal. **64% of configs are pessimistic
(optimism < 0)** — the common Wadoux regime — so this is not a contrived setup.

Script: `data-raw/benchmark_nndm.R`.

## Results

| Task | signed `SLI_rho` | NNDM `W` (unsigned) |
|---|---|---|
| **Predict signed optimism** (optimistic vs pessimistic) | **r = +0.977** | r = −0.090 |
| **Predict magnitude `|optimism|`** | **r = 0.955** | r = 0.672 |
| **Get the direction right** (sign match) | **97%** | silent by construction |

## The point

- **NNDM's `W` cannot sign the bias.** It is an *unsigned* mismatch magnitude, so it
  is essentially uncorrelated with *signed* optimism (r = −0.09). It conflates the
  Milà case (CV closer than deployment → optimistic) with the Wadoux case (CV farther
  than deployment → pessimistic) — exactly the distinction the whole framework exists
  to make. A practitioner with a high `W` does not know whether to trust their numbers
  more or less.
- **The signed SLI recovers the direction** (97% sign agreement) and so predicts
  signed optimism almost perfectly (r = 0.98).
- **Even on magnitude, SLI beats `W`** (0.955 vs 0.672), because `W` is pure geometry
  while `SLI_rho` carries the signal proportion `w` that scales optimism
  (`optimism ≈ V w² · (…)`, eq. 6). `W` ignores how much spatial signal there is.

## Takeaway for the paper

NNDM/kNNDM solve a *fold-construction* problem (make `W` small). `spLeakage` solves a
*diagnosis* problem: given any split, **how much and in which direction is the
reported accuracy biased.** `W` is necessary machinery but insufficient as a
diagnostic — it is direction-blind and signal-blind. This is the empirical
justification for the signed, dependence-aware SLI, and a direct, quantitative
response to the "isn't this just NNDM?" critique.
