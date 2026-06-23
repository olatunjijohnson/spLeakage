# CLAUDE.md — spLeakage

Guidance for working on the `spLeakage` R package. Read `docs/VISION.md` first for
the full methodology, feature catalogue, and roadmap.

## What this package is

A **diagnostic** package for *spatial information leakage*. Given spatial data and a
train/test split (or an existing modelling workflow), it detects, quantifies, and
explains leakage, estimates the optimism it causes, and recommends the appropriate
validation strategy for the user's sampling design and prediction target.

**It is NOT a fold generator.** That space is mature (CAST, blockCV, spatialsample,
mlr3spatiotempcv). We *audit and quantify*, and *delegate* fold generation to those
packages. Keep this boundary sharp — it is the whole reason the package exists.

## Dual end goal

Every change serves **both** a strong methods paper **and** a mature R package.
Before adding a feature, ask: does it strengthen the paper's argument, the
software's usefulness, or (ideally) both? See `docs/PAPER.md`.

## The central thesis (everything serves this)

**Leakage and optimism are only well-posed relative to a declared
`design × estimand × target`.** The same split can be leaking, correct, or
pessimistic depending on those three declarations. The package's job is to make
leakage well-posed given the user's declarations, then diagnose, quantify, and
prescribe. This reframing *is* the reconciliation of the Milà-vs-Wadoux debate and
removes the circular dependency between contributions. See `docs/VISION.md` §3.

Consequences to hold in mind while coding:
- **Optimism is defined relative to the design-matched correct validation, not
  NNDM.** It can legitimately be *negative* (pessimism) for probability samples.
- **Design-basedness is elicited, never inferred from geometry.** Point patterns
  can't distinguish a designed cluster sample from a convenience sample. Geometry
  only flags risk.
- **Estimand first:** population-mean map accuracy (design-based) ≠ conditional
  predictive skill (model-based). Never recommend across estimands.
- **Geographic distance ≠ statistical dependence:** offer covariance-aware variants;
  plain distance is only the cheap default.
- **The emulator must pass its own area-of-applicability check** — don't commit the
  sin the package polices.

The load-bearing novelty is the optimism emulator (C2b) + this reconciliation
framing. `docs/VISION.md` §8 is the reviewer-objection register — treat its
hardening moves as design constraints, not optional polish.

## Conventions

- **Language/tooling:** R (≥ 4.1). `roxygen2` for docs, `testthat` (3e) for tests,
  `devtools`/`usethis` workflow, `pkgdown` site, `lintr`+`styler` (tidyverse style).
- **Dependencies:** prefer light, well-maintained deps. Geometry: `sf`, `terra`,
  `nabor`/`FNN`/`RANN`. Variograms: `gstat` (or own light impl) — keep `gstat` in
  Suggests if possible. Adaptors (`rsample`, `mlr3`, `caret`, `blockCV`) go in
  **Suggests**, guarded by `rlang::check_installed()` / `requireNamespace()`.
- **Naming:** snake_case functions and args; verbs for actions
  (`detect_leakage`, `estimate_optimism`, `recommend_validation`,
  `audit_workflow`, `report_leakage`). Constructors return classed S3 objects with
  `print`/`summary`/`plot`/`autoplot` methods.
- **Output objects:** return rich S3 objects, not bare values; provide tidy
  (`broom`-style) accessors. Side-effect-free core; plotting/reporting separate.
- **CRS/distance:** never silently assume planar coords. Detect CRS; use geodesic
  distance for geographic CRS, projected otherwise; warn on missing CRS.
- **Reproducibility:** seed all stochastic routines; expose `seed`/RNG args.

## Repo layout

```
spLeakage/
  R/             # package source
  tests/testthat # tests (write alongside each feature)
  man/           # roxygen-generated (do not hand-edit)
  vignettes/     # long-form usage + the case studies
  data-raw/      # scripts to build bundled/simulated datasets
  docs/          # VISION.md, PAPER.md (planning; not the pkgdown site)
  DESCRIPTION, NAMESPACE
```

## Dev workflow

- `devtools::load_all()` to iterate; `devtools::document()` after roxygen changes
  (never hand-edit `NAMESPACE`/`man/`); `devtools::test()`; `devtools::check()`
  must pass clean before any milestone.
- Add a test with every exported function. Simulation-based tests use fixed seeds.
- Commits: only when the user asks. This folder is not yet a git repo on its own.

## Research extensions (post-P4, in turns; roadmap in docs/PAPER.md)

- **① Theory DONE** (`docs/THEORY.md`, `docs/THEORY-RESULTS.md`): optimism formalised
  as excess explained variance under a GP; single-NN closed form + Wasserstein bound.
  Validated (192 GRF configs, exact covariance algebra): closed form predicts exact
  optimism r=0.975; **SLI_rho is a near-sufficient statistic for GP optimism
  (r=0.976)**; unsquared SLI beats the single-NN squared form (vindicates the
  package); bound holds 100% for single-NN, multi-NN amplifies ~1.6x.
- **② De-leaked estimator DONE** (`R/deleak.R`): `deleak_estimate()` -- corrected
  accuracy + fold-bootstrap CI; corrects a *reported* metric value without refitting
  (meta-audit enabler). Wired into `report_leakage(deleak=)`.
- Remaining turns: ③ data-driven recommend, ④ meta-audit, ⑤a/b extraction+ST leakage,
  ⑥ optimal validation design, proper scoring rules, conditional conformal, NNDM
  benchmark, real independent-validation, richer emulator generators, multi-NN SLI.

## Current status

Phase **P4 multi-channel mostly done** (R CMD check 0/0/0; ~67 tests; 15 exported
functions). Channels: geographic (`detect_leakage`), grouped/duplicated-location
(`detect_group_leakage`/`group_kfold`), feature-space/covariate-AOA
(`detect_feature_leakage`), temporal lookahead (`detect_temporal_leakage`/
`temporal_kfold`). `.parse_split` accepts fold vectors, fold lists,
`list(test=,train=)`, pre-built fold lists, and tidymodels `rsample` `rset` objects.
CRAN polish: examples on key fns, `NEWS.md`, `cran-comments.md`. Remaining: optional
preprocessing/static-code leakage (stretch), more real datasets, paper write-up.

### Earlier: emulator
Phase **P3 full emulator study done** (R CMD check clean). The optimism
emulator (C2b) is built, validated, and shipped: `R/emulator.R`
(`predict_optimism()`, shared `.emulator_features()`, category-conditional
[model x response-type] instance-based emulator with predictive intervals + a
Meyer-Pebesma **area-of-applicability guard**), calibrated by
`data-raw/simulate_optimism.R`. Artifact in `R/sysdata.rda`; training table in
`data-raw/`.

The study (≈900 denoised rows from 450 configs x 4 realizations, parallelised):
Matern fields x smoothness {0.5,1.5,2.5} x signal x {random, clustered (varied
tightness), preferential} sampling x n {150,300,600} x {gaussian, poisson, binomial}
responses x {idw, rf(ranger), gam(mgcv)} learners x {grid, interpolation} targets.
**Learner: gradient-boosted trees (xgboost) for the point estimate +
**normalized** split-conformal intervals (config-grouped; interval width scales with
a local difficulty sigma(x) = mean residual of the k_sigma feature-space neighbours,
so intervals adapt -- cor(width, |error|) = 0.52, was 0); AOA stays instance-based.** Held-out (by config)
validation (paper-scale, 1000 configs x 6 reps): **R2 ~= 0.76, in-AOA 97%,
90%-interval coverage ~94%**; `cor(SLI_rho, optimism) ~= 0.67`. Figures in
`paper-figures/`; site config `_pkgdown.yml` (builds to `pkgdown/`). The shipped xgboost model is serialised as raw
bytes in `R/sysdata.rda` (**85 KB** total; ranger QRF was 5 MB -> over the CRAN
limit, hence xgboost). `xgboost` is a Suggests, guarded at predict time. Validation
metrics are stored as an attribute on the emulator and shown by `print()`.

Four correctness facts established during P3 (keep these -- they were real bugs):
1. Ground-truth optimism must score the target against **noisy** observations (same
   nugget floor as E_cv), else a noise asymmetry fakes pessimism.
2. Single-realization labels are too noisy (held-out R2 0.19); **averaging reps per
   config** to denoise the label raised R2 to ~0.5.
3. Sampling locations from the candidate **grid floors the NN distance**, so
   clustered samples looked regular (Clark-Evans >= 1) -> jitter with continuous
   noise; and the clustered design needs **varied tightness** to cover realistic
   `nn_index` down to ~0.5.
4. `sli_d` (= A/phi) is numerically unstable (phi can be tiny) and uninformative ->
   dropped from the emulator features; `sli_rho` is the SLI feature.

**Real-data case study done** (`data-raw/case_study_nigeria.R`,
`docs/CASE-STUDY-NIGERIA.md`): Nigeria MAP Pf prevalence (n=66 geolocated; DHS coords
withheld). Exercised the geographic/geodesic path; the package **attributed** the
leakage to co-located repeat surveys (SLI_rho +0.21 with duplicates -> +0.006
deduplicated; audit flagged 8 duplicated coords), and the emulator correctly refused
the out-of-distribution query (n below training envelope) and deferred to the
empirical route. This motivated the now-implemented grouped/duplicate-location
channel (C4): `detect_group_leakage()` + `group_kfold()` (`R/group-leakage.R`), wired
into `audit_workflow()`. On the Nigeria data it recovers the same 21.2% co-location
leakage the SLI flagged (numerically c_obs); `group_kfold()` fixes it to 0%. This
opens phase P4 (multi-channel leakage); remaining channels: feature-space (AOA-style),
temporal/spatiotemporal, preprocessing.

Remaining for the paper: scale configs/reps, add covariate trends + spatial-GLMM
realism, locally-adaptive (normalised) conformal to close the last coverage gap
(85% -> 90%), and formal out-of-distribution validation on the real case-study
datasets. Feature contract + AOA architecture are fixed. The Milà-vs-Wadoux split
reproduces in the data (clustered -> +optimism, random -> ~0).

### Prior phases
Phase **P2 core done** (R CMD check clean: 0/0/0). On top of the polished
P1 diagnostic (geometry, dependence+anisotropy+param covariance, `detect_leakage`
with signed `SLI_d`/`SLI_rho`, crossing decomposition, per-point/fold leakage,
Monte-Carlo `n_boot` uncertainty), P2 added:
- `estimate_optimism()` empirical route (`R/optimism.R`) -- optimism = user-CV error
  minus **design-matched controlled** error (block / buffered-block), NOT NNDM;
  can be negative; model-agnostic via `predict_fun` (IDW default). C2-(1) fix baked in.
- CV-scheme infra (`R/cv-schemes.R`): `spatial_block_cv()`, buffered folds, CV runner,
  IDW learner, metrics.
- `recommend_validation()` (`R/recommend.R`): estimand-first design x target decision
  engine; design is **elicited, never inferred**; Clark-Evans NN index is a risk flag
  only. C3 hardening baked in.
- `audit_workflow()` + `report_leakage()` scorecard (`R/audit-report.R`).
Tests: `test-optimism.R`, `test-recommend-audit.R`.

Notable design corrections made during P2 (keep in mind): buffered CV must run on
spatially *contiguous* folds (buffering a dispersed random fold empties training);
buffers are capped at 25% of domain diameter (range >= domain => no independent data).

Next: a first **vignette** + the real disease-mapping case study (`malariaAtlas`,
see `docs/DATA.md`) to hit the P2 milestone, optional `gstat` backend and ECDF
threshold sweep, then **P3** (simulation study + optimism emulator, the load-bearing
novelty -- see `docs/METHOD-EMULATOR.md`). Don't start the emulator before the
simulation harness and `optimism_rel` ground-truth are in place.

## Working style here

- This is a research project: think before coding, document decisions in `docs/`.
- When unsure about a methodological choice, surface it in `docs/VISION.md` §8
  (open questions) rather than silently picking.
- Verify references against publisher pages before citing in the paper.
