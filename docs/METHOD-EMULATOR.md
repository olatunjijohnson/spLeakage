# Locked spec — the optimism emulator & its simulation generator

Status: **locked** (v1). The load-bearing methodological novelty (C2b). Answers
objections C2-① … C2-④ in `VISION.md` §8. Built in phase P3; this spec is fixed now
so P1/P2 outputs (SLI, geometry features) are emulator-ready.

## 0. What the emulator is

A pre-trained, compact regression model

```
emulator : features(data, split, declaration)  →  expected optimism  (+ interval)
```

shipped with the package, so a user gets an instant optimism estimate **without
refitting their model**. It is calibrated once, offline, on a large simulation study
where the *true* optimism is known by construction.

## 1. The estimand: optimism defined against ground truth, NOT against NNDM

This is the fix for the fatal objection C2-①. Because the data are simulated, we
know the true field, so we can compute the **true generalisation error on the
declared target** directly — we never treat NNDM as "truth."

For one simulated replicate, a fitted model, a chosen metric `E` (RMSE, MAE, log-
score, Brier/AUC, …), and a declared `(design, estimand, target)`:

```
E_cv    = error the user's CV scheme reports        (e.g. random 10-fold)
E_true  = error on fresh target locations drawn from the declared target geometry,
          scored against the known truth
optimism_abs = E_true − E_cv
optimism_rel = (E_true − E_cv) / E_true        (scale-free; the modelled response)
```

`optimism_rel > 0` ⇒ the user's scheme was optimistic; **`< 0` ⇒ pessimistic**
(legitimately occurs for probability samples / design-based estimand — the Wadoux
case). The emulator predicts `optimism_rel`, which can be negative. This is what
makes "optimism" well-posed and design-relative rather than NNDM-relative.

## 2. Simulation generator — the factor grid

A realism ladder, not a single toy. Factors crossed (fractional/space-filling
design over the grid; not full factorial):

| # | factor | levels / range | why it matters |
|---|--------|----------------|----------------|
| 1 | autocorrelation range `φ/L` | 0.02, 0.05, 0.1, 0.2, 0.5 (rel. to domain `L`) | sets leakage length scale |
| 2 | smoothness (Matérn `ν`) | 0.5, 1.5, 2.5 | field roughness |
| 3 | spatial-signal proportion `σ²/(σ²+τ²)` | 0.1 → 0.9 | **no signal ⇒ no optimism** |
| 4 | sampling design | uniform-random, mild cluster, strong cluster, **preferential** (intensity ∝ field; Diggle), transect, gridded | design drives both leakage *and* range-estimation error |
| 5 | sample size `n` | 100, 300, 1000 | small-`n` instability |
| 6 | response type | Gaussian, **Poisson (counts)**, **Binomial (prevalence)** via spatial GLMM | disease-mapping realism (objection C2-②) |
| 7 | covariates | none, smooth trend, covariate-with-own-spatial-structure | separates trend vs autocorrelation; confounding |
| 8 | prediction target | within-sample interpolation, wall-to-wall, new-region extrapolation | the deployment geometry |
| 9 | user CV scheme (evaluated) | random k-fold, spatial block, buffered LOO, NNDM/kNNDM | the thing whose optimism we predict |
| 10 | model class fitted | kriging/GP, Random Forest, GAM/GLM | optimism is model-dependent (objection C2-⑤) |

Replicates per cell: enough for stable Monte-Carlo truth (target ≥ 50; tune).
Everything seeded and scripted in `data-raw/` (reproducible).

Key design choices that answer reviewers:
- **Preferential sampling is in the grid** (factor 4) because that is exactly where
  range estimation breaks and leakage is worst (objections C2-②/④).
- **Non-Gaussian responses via spatial GLMM** (factor 6) so the calibration universe
  covers counts/prevalence, not just Gaussian fields (objection C2-②).
- **Model class is a factor** (factor 10) so the emulator predicts *model-conditional*
  optimism (objection C2-⑤).

## 3. Features (emulator inputs) — estimated, not oracle

Computed from the sample/split exactly as a real user would (so the emulator learns
to tolerate realistic input noise — objection C2-④):

- `φ̂` (estimated range), `ν̂`, estimated spatial-signal proportion (from variogram).
- **Signed `SLI_ρ` and `SLI_d`**, `δ` (directionality), `W` — the C1 outputs.
- Sampling-clustering index (Clark–Evans / nearest-neighbour index; Ripley-`L`
  summary), preferential-sampling diagnostic.
- `n`, point density, domain shape summary.
- One-hot: response type, target type, CV scheme, model class, declared design.

Crucially the features are **the same ones P1/P2 already compute** — so the emulator
consumes the package's own diagnostics. That is why this spec is locked before P1.

## 4. The emulator model

- **Form:** gradient-boosted trees or a GP regressor on the feature→`optimism_rel`
  table; quantile loss (or full predictive distribution) to get the **interval**.
- **Uncertainty:** combine emulator predictive uncertainty with **propagated input
  uncertainty** (resample `φ̂` etc. from their estimation uncertainty, push through
  the emulator). Report `optimism_rel` with a calibrated interval; check interval
  coverage on held-out cells.
- **Compactness:** ship a small fitted object in `data/` (built by `data-raw/`),
  loadable without heavy dependencies at predict time.

## 5. The emulator's own area-of-applicability (objection C2-③)

The emulator is itself a spatial-ish predictor and must not commit the sin the
package polices. So we ship, alongside it:
- the **training feature envelope** + a Meyer–Pebesma-style dissimilarity (AOA)
  check;
- at predict time, compute the query's dissimilarity; if **outside the area of
  applicability**, the package **refuses to give a point estimate** and falls back to
  the empirical/refit route (C2a) with a clear message.

This self-AOA guard is a deliberate, reviewer-facing feature ("quis custodiet").

## 6. Validation plan (paper — Figure 2 and friends)

1. **Within-simulation:** held-out *cells* of the factor grid — emulator-predicted
   vs true `optimism_rel` (calibration plot, interval coverage).
2. **Out-of-distribution (the answer to "GRF is a toy"):** evaluate on the **real
   datasets** in `DATA.md`, where `E_true` is approximated by the strongest available
   leakage-controlled estimate; report degradation honestly.
3. **Ablations:** drop signal (factor 3 → 0) ⇒ optimism → 0 (sanity); show
   `optimism_rel` rises with `SLI_ρ × signal` (the C1↔C2 mechanistic link).
4. **Negative-optimism demonstration:** probability-sample + design-based estimand
   cells produce `optimism_rel < 0` and the design-aware recommendation (C3) avoids
   it — the head-to-head that reconciles Milà vs Wadoux.

## 7. Dependencies / build

Simulation in `data-raw/` may use `gstat`/`fields`/`RandomFields`-style GRF
simulation, `sf`/`terra`, a spatial-GLMM tool for non-Gaussian responses, and
`ranger`/GP libs for model classes — all **build-time only**. The shipped artefacts
(fitted emulator + AOA reference) live in `data/`; runtime depends only on a light
predictor. Keep the heavy simulation deps out of `Imports`.
