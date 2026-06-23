# Package index

## Overview

- [`spLeakage`](spLeakage-package.md)
  [`spLeakage-package`](spLeakage-package.md) : spLeakage: Detect and
  Quantify Spatial Information Leakage

## Diagnose leakage

Detect and quantify spatial information leakage in a split.

- [`detect_leakage()`](detect_leakage.md) : Detect and quantify spatial
  leakage in a train/test split
- [`sli()`](sli.md) : Extract the Spatial Leakage Index
- [`plot(`*`<leakage_diagnosis>`*`)`](plot.leakage_diagnosis.md) : Plot
  a leakage diagnosis
- [`estimate_dependence()`](estimate_dependence.md) : Estimate spatial
  dependence (variogram, range, correlation, anisotropy)
- [`prediction_target()`](prediction_target.md) : Declare the prediction
  target (deployment geometry)

## Quantify optimism

How inflated is the reported accuracy?

- [`estimate_optimism()`](estimate_optimism.md) : Estimate the optimism
  (accuracy inflation) of a cross-validation scheme
- [`predict_optimism()`](predict_optimism.md) : Predict optimism from
  the calibrated emulator (fast, no model refitting)
- [`deleak_estimate()`](deleak_estimate.md) : De-leaked (bias-corrected)
  accuracy estimate

## Recommend & report

- [`recommend_validation()`](recommend_validation.md) : Recommend a
  validation strategy for a sampling design and prediction target
- [`rank_cv_schemes()`](rank_cv_schemes.md) : Rank candidate
  cross-validation schemes by deployment match
- [`audit_workflow()`](audit_workflow.md) : Audit a cross-validation
  workflow for spatial leakage
- [`report_leakage()`](report_leakage.md) : Assemble a leakage report
  (scorecard)

## Multi-channel leakage

Grouped/duplicated-location, feature-space, and temporal leakage.

- [`detect_group_leakage()`](detect_group_leakage.md) : Detect grouped /
  duplicated-location leakage in a split
- [`group_kfold()`](group_kfold.md) : Group-aware k-fold assignment (the
  fix for grouped leakage)
- [`detect_feature_leakage()`](detect_feature_leakage.md) : Detect
  feature-space (covariate) leakage
- [`detect_extraction_leakage()`](detect_extraction_leakage.md) : Detect
  covariate-extraction-overlap leakage
- [`detect_temporal_leakage()`](detect_temporal_leakage.md) : Detect
  temporal (lookahead) leakage in a split
- [`temporal_kfold()`](temporal_kfold.md) : Forward-chaining temporal
  cross-validation (the fix for lookahead leakage)
- [`detect_st_leakage()`](detect_st_leakage.md) : Detect joint
  spatiotemporal leakage
- [`trend_strength()`](trend_strength.md) : Trend strength: large-scale
  structure that resists extrapolation

## Validation design

Where to place independent validation points to estimate true accuracy.

- [`design_validation()`](design_validation.md) : Design an optimal
  independent validation sample

## Cross-validation schemes

- [`spatial_block_cv()`](spatial_block_cv.md) : Spatial block
  cross-validation folds
