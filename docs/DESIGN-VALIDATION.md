# Optimal validation design — inverting the SLI

The diagnostic question is "is my validation biased?"; the **design** question is
"where should I put a budget of independent validation points to estimate true map
accuracy as precisely as possible?" `design_validation()` answers the second by
**inverting the SLI**: rather than detecting a CV-vs-deployment mismatch, it *places*
validation points to match deployment and minimise the variance of the accuracy
estimate. This unifies leakage diagnosis with spatial survey design (and the author's
NTD survey-design work).

## Method

Given training locations, a deployment target, and a budget `n_v`:

1. **Stratify by distance to training** (strata defined on the *deployment*
   distance-to-training distribution, so the strata represent where the map will be
   used).
2. **Per-stratum error model** from the GP theory: `L(d) = V[1 − w² ρ(d)²]`
   (`docs/THEORY.md`), giving the pointwise MSE — large in under-sampled (far) strata,
   small near training.
3. **Allocate the budget** across strata:
   - `proportional` — match deployment (unbiased with a simple mean);
   - `optimal` (Neyman) — allocate ∝ `W_h · S_h` (stratum weight × error SD) to
     **minimise the variance** of the estimated mean accuracy, with **inclusion
     weights** (`W_h / n_h`) returned so the weighted-mean estimate stays unbiased.
4. **Spatially spread** the points within each stratum (k-means centroids).

## Validation (`data-raw/validate_design.R`)

120 replicates; estimate the true map RMSE (known over a dense grid) from a budget of
30 validation points placed randomly vs by the optimal design:

| placement | bias | SD of the accuracy estimate |
|---|---|---|
| random | −0.013 | 0.102 |
| **optimal design** | −0.016 | **0.056** |

**The optimal design estimates true map accuracy with 45% lower variance for the same
budget** (theory predicts ≥24% from the within-stratum term alone; the empirical gain
is larger because stratification also removes between-stratum variance). Bias stays
small. With 30 points instead of, say, 100, you get the same precision — a direct,
quantified survey-saving.

## Why it matters

- A genuinely **novel** contribution: optimal experimental design for *validation*
  (not prediction), which the spatial-CV literature does not address.
- It closes the loop opened by turn ⑳ (true accuracy needs independent points) and
  the meta-audit (the field is overconfident): not only *detect* and *correct*
  optimism, but *design the cheapest data collection to measure true accuracy*.
- Directly applicable to disease-mapping / NTD survey design: "given a budget, where
  to put validation surveys to pin down map accuracy" — the original `spSurveyAI`
  idea, now grounded in the leakage theory.
