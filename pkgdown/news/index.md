# Changelog

## spLeakage 0.1.0

First release: a diagnostic toolkit for spatial information leakage.

### Diagnose

- [`detect_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_leakage.md)
  — signed Spatial Leakage Index (`SLI_rho` dependence form and `SLI_d`
  distance form), crossing decomposition, per-point/per-fold leakage
  map, geometric anisotropy, Monte-Carlo uncertainty (`n_boot`), and an
  optional multi-neighbour retained-correlation form (`k`; `k = 1`
  default is recommended).
- [`estimate_dependence()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_dependence.md)
  — variogram, correlation function, practical range, optional
  anisotropy.

### Quantify optimism

- [`estimate_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/estimate_optimism.md)
  — empirical route: optimism relative to a design-matched
  leakage-controlled scheme (can be negative for probability samples).
  Supports proper scoring rules (`brier`/`logloss` for binary, `poisson`
  deviance for counts) so optimism is meaningful for disease-prevalence
  and count responses.
- [`predict_optimism()`](https://olatunjijohnson.github.io/spLeakage/reference/predict_optimism.md)
  — fast emulator route (gradient-boosted point estimate +
  **normalized** split-conformal intervals that widen in harder regions)
  with an area-of-applicability guard. Calibrated on a simulation study
  (`data-raw/simulate_optimism.R`).
- [`deleak_estimate()`](https://olatunjijohnson.github.io/spLeakage/reference/deleak_estimate.md)
  — bias-corrected (“de-leaked”) accuracy estimate with a fold bootstrap
  CI; can also correct a published metric value without refitting, and
  (with `target=`) uses a target-matched control for extrapolation
  deployment.

### Recommend & report

- [`recommend_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/recommend_validation.md)
  — estimand-first, design × target decision engine; design is elicited,
  never inferred from geometry.
- [`rank_cv_schemes()`](https://olatunjijohnson.github.io/spLeakage/reference/rank_cv_schemes.md)
  — data-driven: ranks candidate CV schemes by how well their geometry
  matches the declared deployment target (smallest `|SLI_rho|`).
- [`audit_workflow()`](https://olatunjijohnson.github.io/spLeakage/reference/audit_workflow.md),
  [`report_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/report_leakage.md)
  — submission-ready scorecard.

### Multi-channel leakage

- [`detect_group_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_group_leakage.md)
  /
  [`group_kfold()`](https://olatunjijohnson.github.io/spLeakage/reference/group_kfold.md)
  — grouped / duplicated-location leakage.
- [`detect_feature_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_feature_leakage.md)
  — covariate-space (area-of-applicability) leakage.
- [`detect_extraction_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_extraction_leakage.md)
  — covariate-extraction-overlap leakage (focal/buffer/ kernel windows
  that straddle the train/test boundary; common in remote sensing/SDM).
- [`detect_temporal_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_temporal_leakage.md)
  /
  [`temporal_kfold()`](https://olatunjijohnson.github.io/spLeakage/reference/temporal_kfold.md)
  — temporal lookahead leakage.
- [`detect_st_leakage()`](https://olatunjijohnson.github.io/spLeakage/reference/detect_st_leakage.md)
  — joint spatiotemporal leakage in a scaled space-time metric (catches
  leakage that space-only or time-only analysis misses).
- [`trend_strength()`](https://olatunjijohnson.github.io/spLeakage/reference/trend_strength.md)
  — large-scale trend channel (extrapolation-optimism a trend-blind
  model misses; the second optimism source found in the meta-audit).
  Surfaced in
  [`audit_workflow()`](https://olatunjijohnson.github.io/spLeakage/reference/audit_workflow.md).

### Validation design

- [`design_validation()`](https://olatunjijohnson.github.io/spLeakage/reference/design_validation.md)
  — optimal placement of a budget of independent validation points
  (distance-stratified Neyman allocation + inclusion weights) to
  estimate true map accuracy with minimum variance (~45% lower SD than
  random placement in tests).

### Interoperability

- Accepts plain fold vectors, fold lists, `list(test=, train=)`,
  pre-built fold lists, and tidymodels `rsample` `rset` objects.

### Validation

- End-to-end real-data case study on Nigeria malaria prevalence
  (`docs/CASE-STUDY-NIGERIA.md`).
