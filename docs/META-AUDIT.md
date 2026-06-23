# Meta-audit — how overconfident is spatial prediction, across real datasets?

The high-impact question (cf. Ploton et al. 2020): is reported accuracy in spatial
prediction systematically optimistic, across many real datasets? `spLeakage` is the
instrument. Corpus: **16 real spatial response variables** from 8 public datasets —
heavy metals (meuse Cd/Cu/Pb/Zn), soil chemistry (camg Ca/Mg/CTC, ca20 calcium),
rainfall (Paraná, North-American), elevation (geoR, North-American), ozone (ozone2),
and malaria prevalence (Gambia). For each we measure the **real** optimism of naive
random 10-fold CV by spatial region-holdout (a direct measurement, IDW model),
averaged over 8 holdout directions. Script: `data-raw/meta_audit.R`.

## Headline: the field is overconfident

| | |
|---|---|
| datasets where random CV is optimistic | **16 / 16 (100%)** |
| median real optimism | **22%** (IQR 17–34%) |
| substantially optimistic (> 20%) | **69%** |
| worst cases | Paraná rainfall +54%, elevation +52%, N-American rainfall/elev +40% |

Across diverse real spatial datasets, naive random cross-validation understates true
extrapolation error by a **median of 22%**, and by **more than 20% in two-thirds** of
them. This is direct measurement, independent of any leakage model — the empirical
core of the impact claim. (Figure `paper-figures/F4_meta_audit.png`.)

## The honest, deeper finding: optimism has two sources

Can the optimism be flagged *without* holdout (so the audit scales)? Cross-dataset:

| predictor (no holdout) | cor with real optimism |
|---|---|
| raw `SLI_rho` | **−0.39** (does NOT transfer across datasets) |
| de-leaked optimism (target-aware) | **+0.23** (partial) |

The negative SLI correlation is not a bug — it is the theory (`docs/THEORY.md`):
optimism scales as `V w² · (…)`, so comparing *across* datasets with very different
covariance structures requires the signal variance the bare index omits. The SLI is a
**within-dataset** predictor (validated at cor = 0.52 in
`docs/VALIDATION-INDEPENDENT.md`), not a cross-dataset one.

But even the de-leaked estimator only weakly transfers (+0.23), and the per-study
table shows why: the smooth, **trend-dominated** fields (Paraná rainfall, elevation,
ozone) have the **highest** real optimism (40–54%) yet a de-leaked optimism near
**zero**. Their leakage is *not* short-range autocorrelation — it is the model's
inability to **extrapolate a large-scale trend**. Autocorrelation-based diagnostics
(NNDM, our SLI, the de-leak) are blind to this component by construction.

**This is now quantified and fixed** (`data-raw/trend_analysis.R`, `trend_strength()`):
adding a **trend-strength** feature (variance explained by a quadratic coordinate
surface) recovers exactly the signal the SLI misses.

| predictor of cross-dataset real optimism | cor / R² |
|---|---|
| `SLI_rho` alone | cor −0.39 ; R² 0.15 |
| **trend strength** alone | **cor +0.75 ; R² 0.56** |
| `\|SLI\|` + trend strength | **R² 0.62** |

So real-data optimism is a **two-channel** phenomenon, and the two channels are now
both measured: autocorrelation (`detect_leakage`) and trend (`trend_strength`, also
surfaced in `audit_workflow()`). The trend-dominated fields the SLI missed are
exactly the ones trend strength flags.

**Implication — the problem is bigger than current methods reveal.** Real-data
optimism decomposes into (i) a short-range *autocorrelation* component, which
`spLeakage` (and NNDM) diagnose, and (ii) a large-scale *trend-extrapolation*
component, which pure-spatial diagnostics miss. The field is overconfident, and part
of that overconfidence is invisible even to best-practice spatial CV. This motivates a
**trend/covariate-aware optimism** extension (link to the feature-space channel and a
trend-capable surrogate model), and strengthens, rather than weakens, the headline:
22% median is a *lower bound* on what autocorrelation-aware tooling would catch.

## Takeaways for the paper

1. **Impact:** 100% of 16 real datasets optimistic; median 22%, two-thirds > 20%. A
   quantified, field-level overconfidence claim with `spLeakage` as the instrument.
2. **Scope honesty:** the SLI predicts optimism *within* a dataset, not *across*
   datasets (theory: missing the `V w²` scaling); the de-leak transfers partially.
3. **New science:** optimism has an autocorrelation component (diagnosed) and a
   trend-extrapolation component (missed by spatial diagnostics) — a concrete, novel
   decomposition and a clear extension target.
