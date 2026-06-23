# spLeakage — Paper Plan

Living outline for the methods/software paper. Keep in sync with `docs/VISION.md`.

## Working title
*Detecting and quantifying spatial information leakage in predictive modelling:
a diagnostic framework and the `spLeakage` R package.*

## Core claim / gap
Existing spatial-CV tools **generate corrected folds**; none **diagnose an existing
split**, **quantify the optimism** it causes, or **reconcile** the Milà-vs-Wadoux
debate. `spLeakage` does all three — held together by one thesis:

> **Leakage and optimism are only well-posed relative to a declared
> `design × estimand × target`.** Make them well-posed, then diagnose (C1),
> quantify (C2), and prescribe (C3).

This reframing is the paper's real intellectual contribution: it dissolves the
Milà-vs-Wadoux controversy (each is right in a different cell of the space) and
removes the circular dependency a reviewer would otherwise attack. The load-bearing
*methodological* novelty is the optimism emulator (C2b). See `VISION.md` §3 and the
reviewer-objection register §8.

## Contributions (see VISION §3)
1. **SLI** — Spatial Leakage Index: a post-hoc, [0,1] leakage score for *any* split,
   from the divergence between observed test→train and target prediction→sample
   NN-distance ECDFs.
2. **Optimism estimator** — empirical (refit) + a cheap **simulation-calibrated
   emulator** mapping data/split geometry to expected accuracy inflation.
3. **Design-aware recommendation** — a decision framework stating which validation
   is correct (and when spatial CV is *not*).

## Planned figures/experiments
- F1. Conceptual: random split leaks; ECDF mismatch picture.
- F2. Simulation: SLI & emulator vs *true* optimism across autocorrelation range ×
  sampling clustering × prediction target (GRF study, shipped in `data-raw/`).
- F3. Design × target decision matrix with empirical support (shows random CV best
  for probability samples; spatial needed for clustered).
- F4–F6. Real case studies (target ≥1 disease-mapping/NTD dataset) re-analysed:
  "reported accuracy was X% optimistic."

## Results obtained so far (real numbers, to be scaled for the paper)

**Emulator calibration & validation** (`data-raw/simulate_optimism.R`, ~2000 denoised
rows from 1000 configs × 6 realizations; Matérn × smoothness × signal × {random,
clustered, preferential} × n × {gaussian, poisson, binomial} × {idw, rf, gam}):
- Held-out (by config) **R² = 0.76**, **in-AOA 97%**, **90%-interval coverage ≈ 94%**
  (gradient-boosted point + config-grouped split-conformal intervals). Figure F1.
- **`cor(SLI_rho, optimism) = 0.67`** — direct evidence for the C1→C2 mechanism (F2).
- **Design effect (reconciles Milà vs Wadoux), F3:** the *ordering* is robust —
  clustered/preferential sampling makes random CV more optimistic than a
  probability-like random sample (mean optimism random −0.06 < clustered −0.02). For
  these well-covered grid targets random CV is fine/slightly pessimistic (Wadoux),
  and clustering shifts it toward optimism (Milà). Absolute level depends on target
  geometry/model; the relative design effect is the paper-relevant claim.
- The two optimism routes agree on a held-in example (emulator vs empirical refit).

**Nigeria malaria case study** (`docs/CASE-STUDY-NIGERIA.md`; MAP Pf prevalence,
n = 66 geolocated):
- `SLI_rho = +0.212` on a naive random 10-fold split; variogram `signal_prop ≈ 0`.
- **Leakage attributed to co-located repeat surveys**: deduplicating collapses
  `SLI_rho` to **+0.006**; the grouped-leakage channel reports **21.2%** of test
  points leaking via a shared location — *numerically equal to* `c_obs = 0.212`,
  the two diagnostics agreeing — and `group_kfold()` drives it to **0%**.
- Empirical optimism a modest **+1.8%** (co-located records differ by year; weak
  signal) — the C1-vs-C2 distinction in the wild.
- The emulator **correctly refused** the out-of-distribution query (n below training
  envelope), demonstrating the AOA guard on real data.

**Multi-channel leakage implemented:** geographic (SLI), grouped/duplicated-location,
feature-space (covariate AOA), and temporal (lookahead) — each with a paired fix
(`group_kfold`, `temporal_kfold`, design-matched control).

This already supplies the spine of F2 (emulator validation), F3 (design split), and
F4 (real case study). For the paper: scale configs/reps, add covariate-trend +
spatial-GLMM realism and normalised conformal, and add 2–3 more real datasets.

## Target venues (in order)
Methods in Ecology and Evolution · Geoscientific Model Development · Journal of
Statistical Software · Ecography.

## Datasets (see `docs/DATA.md`)
- **Anchor (author's domain):** Malaria Atlas Project prevalence via `malariaAtlas`
  (open) — Binomial, clustered/preferential, wall-to-wall.
- **Continuous:** CAST `cookfarm` (interpolation) and `splotdata` (wall-to-wall).
- **Binary/SDM:** blockCV Australia presence/absence.
- **Benchmark/simulation:** NNDMpaper harness (Zenodo `10.5281/zenodo.6366985`) +
  our own generator (`docs/METHOD-EMULATOR.md`).

## To decide
- A real, shareable **probability-sample** dataset for the design-based /
  negative-optimism demonstration (else simulation-only).
- Whether to include Ploton AGB data (licence/size).

## Theory (extension ①, DONE — `docs/THEORY.md`, `docs/THEORY-RESULTS.md`)
Optimism is formalised as **excess explained variance** under a GP (a real estimand),
with a single-nearest-neighbour closed form `optimism ≈ V w² (E_t[ρ²(g)] − E_p[ρ²(f)])`
and a **Wasserstein-1 bound** `|optimism| ≤ ||L'|| · W₁(P_cv, P_dep)`. Numerically
validated (192 GRF configs, exact covariance algebra):
- the closed form predicts exact optimism with **r = 0.975**;
- **`SLI_rho` is a near-sufficient statistic for GP optimism: r = 0.976** (~95% of
  variance) — and the unsquared form beats the single-NN squared form (0.841),
  vindicating the package;
- the Wasserstein bound holds in **100%** for the single-NN optimism (multi-NN
  amplifies the exact by ~1.6×).
This gives the SLI a formal estimand, a theory-grounded (emulator-free) optimism
estimate, and a bridge to distribution-shift generalization theory — the upgrade that
opens a top-tier (stats/ML-theory) venue.

## Benchmark vs NNDM (DONE — `docs/BENCHMARK-NNDM.md`)
Answers "isn't this just NNDM's W?". Against exact GP optimism (288 configs, 64%
pessimistic): signed `SLI_rho` predicts *signed* optimism at **r = 0.98**; NNDM's
*unsigned* `W` at **r = −0.09** (it cannot sign the bias — conflates Milà optimism
with Wadoux pessimism). SLI also beats `W` on magnitude (0.96 vs 0.67, since `W`
ignores the signal proportion). 97% sign agreement. This is the quantitative
justification for the signed, dependence-aware index.

## Data-driven recommendation (extension ③, DONE — `R/rank-schemes.R`)
`rank_cv_schemes()` turns C3 from an expert table into a computation: it builds each
candidate scheme, measures its `SLI_rho` against the declared target, and ranks by
deployment match. On clustered data it correctly shows random CV optimistic, buffered
CV *over-corrected into pessimism*, block CV best — nuance a static table cannot give.

## Real-data validation vs true independent error (turn ⑳, DONE — `docs/VALIDATION-INDEPENDENT.md`)
Three real fields (meuse, Paraná, ca20), 240 replicates, train/hold-out with a
genuinely independent test set across interpolation vs extrapolation deployments:
- **naive random CV understates true extrapolation error by a median 36%** on real
  data — the problem is real and large;
- the target-aware **SLI separates the scenarios as real optimism does** (SLI +0.28
  extrapolation vs −0.01 interpolation; real optimism +0.35 vs −0.04;
  `cor(SLI, real optimism) = 0.52`);
- the de-leaked estimate closes **58%** of the gap with a fixed block control, and
  **71%** with the **target-aware** control (`deleak_estimate(target=)`, turn #23),
  recovering 0.89x of true independent error without seeing it.
- Remaining 29% is the hardest extrapolation, where AOA abstention is the honest
  answer.

## Meta-audit (turn ④, DONE — `docs/META-AUDIT.md`, figure F4)
The impact result. 16 real spatial response variables (metals, soil chemistry,
rainfall, elevation, ozone, malaria); real optimism measured by spatial holdout:
- **100% of datasets optimistic; median 22% (IQR 17–34%); 69% > 20%.** Field-level
  overconfidence, quantified, with spLeakage as the instrument.
- **New decomposition:** optimism has a short-range *autocorrelation* component
  (diagnosed by SLI/de-leak) and a large-scale *trend-extrapolation* component
  (missed by all spatial-CV diagnostics; smooth fields like rainfall/elevation show
  40–54% real optimism but ~0 de-leak signal). Cross-dataset: raw SLI does NOT
  transfer (−0.39, per the V·w² theory); de-leak transfers partially (+0.23). → 22%
  is a *lower bound*; motivates trend/covariate-aware optimism (task #24).

## Optimal validation design (extension ⑥, DONE — `docs/DESIGN-VALIDATION.md`)
Inverts the SLI: `design_validation()` places a budget of independent validation
points (distance-stratified Neyman allocation, per-stratum error from the GP theory
`L(d)`, Horvitz-Thompson inclusion weights) to estimate true map accuracy with
minimum variance. Validated (120 reps): **45% lower SD than random placement for the
same budget**, small bias. Unifies leakage diagnosis with spatial survey design (the
`spSurveyAI` idea, grounded in the leakage theory) — a novel "design for validation,
not prediction" contribution.

## Method specs (locked)
- Signed SLI: `docs/METHOD-SLI.md`. Optimism emulator: `docs/METHOD-EMULATOR.md`.
- Theory: `docs/THEORY.md` (+ results `docs/THEORY-RESULTS.md`).

## Reference list
Maintained in `docs/VISION.md` §9 — verify all before submission.
