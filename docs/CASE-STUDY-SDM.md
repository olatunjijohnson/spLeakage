# Ecological case study — species distribution model (presence/absence)

The flagship *ecological* case study for the paper (lead case study, ahead of the
malaria example). Script: `data-raw/case_study_sdm.R`. Data: the Valavi et al. (2019,
*MEE*) south-east Australia presence/absence dataset shipped with **blockCV** — 500
records (243 presence / 257 absence) with four bioclim covariates (bio_4, bio_5,
bio_12, bio_15). This is a canonical SDM benchmark, from an MEE paper, so it is
in-scope and familiar to reviewers.

## Setup

An ecologist builds a wall-to-wall habitat-suitability map over SE Australia and
validates it with a **naive random 10-fold split** (the common default). Deployment
target = 4,000 sampled cells of the bioclim raster (the map's actual footprint).
Response is binary presence/absence, so optimism is measured with the **Brier score**
(a proper scoring rule).

## Results

| Diagnostic | Result | Reading |
|---|---|---|
| `SLI_rho` (random 10-fold) | **+0.103** (90% CI +0.071, +0.142) | optimistic geographic leakage; CI excludes 0 |
| Brier-score optimism vs block control | **+41.6%** | the random split makes the Brier score ~42% too good |
| Feature-space leakage (bioclim) | feature `SLI = +2.69`; test-in-AOA 93%; reach test 1.04 vs **deployment 3.73** | the wall-to-wall map **extrapolates the covariates** far beyond where CV tested |
| Co-location / duplicated coords | **0%** / 0 | the leak is *not* a survey artefact (contrast malaria) |
| `rank_cv_schemes` | **block** best (SLI −0.004); random +0.103 (optimistic); buffered −0.220 (over-corrected to pessimism) | block matches deployment; buffered over-corrects |
| `recommend_validation` (clustered, grid) | NNDM/kNNDM or block; avoid random; NN index 0.44 (clustered) | design-aware advice |
| `audit_workflow` grade | **C** | autocorrelation flagged; trend mild (0.26); CRS undocumented |

## Why this case study earns its place

1. **Canonically ecological and reproducible** — a presence/absence SDM from an MEE
   paper, shipped openly in blockCV; no licensing or redistribution issues.
2. **Exercises a different leakage signature from malaria.** Here the optimism is
   *genuine spatial autocorrelation* (SLI +0.10, +42% Brier optimism) **plus
   covariate-space extrapolation** (the map predicts well outside the bioclim range
   the model was tested on — an area-of-applicability problem). There are **no**
   co-located duplicates. In the malaria case the geographic SLI was instead traced
   to co-located repeat surveys (group channel), with near-zero true autocorrelation.
   Together the two case studies show the multi-channel audit attributing leakage to
   the *right* cause.
3. **Reproduces the paper's headline nuance**: `rank_cv_schemes` shows random CV
   optimistic, buffered CV *over-corrected into pessimism*, and block CV best — the
   nuance a static decision table cannot give.
4. **Ecologically actionable**: the feature-space result is a concrete AOA warning an
   SDM practitioner cares about (don't trust suitability where bioclim is
   extrapolated).

## Note on CRS
The blockCV coordinates are projected (equal-area, metres); distances are therefore
Euclidean in metres and the variogram range `phi ≈ 1.85e6 m`. `audit_workflow()`
flags `CRS documented: FALSE` only because the plain data.frame carries no CRS
attribute; supplying an `sf` object with the projection set clears the flag.
