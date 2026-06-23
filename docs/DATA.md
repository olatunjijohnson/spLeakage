# Datasets — inventory & plan

What we found by checking the cited papers, and the plan for case studies + the
emulator's out-of-distribution validation. Goal: cover the **design × estimand ×
target** space across **response types** (Gaussian / binary / count-prevalence) using
**openly licensed** data we can reach from vignettes without redistributing.

## Finding: the cited papers give us simulation harnesses + a few real datasets, but **no disease-mapping data**

| Source (cited) | What it provides | Usable how |
|---|---|---|
| **Milà et al. 2022 (NNDM)** — Zenodo `10.5281/zenodo.6366985`, GitHub `carlesmila/NNDMpaper` | Landscape data + **simulated** outcomes (random fields; virtual species), CSV outputs of sims 1 & 2 | **Benchmark harness.** Reuse to compare spLeakage apples-to-apples with NNDM; a template for our own simulation generator. Mostly *simulated*, not a shareable real dataset. |
| **CAST package** (`cookfarm`) | Soil-moisture logger data, 42 locations, spatio-temporal, 2007–2013; **continuous** | Real **interpolation / digital-soil-mapping** case study. Ships with CAST (Suggests). |
| **CAST package** (`splotdata`) | sPlotOpen vegetation plots, South America, **species richness** + predictors; **continuous** | Real **wall-to-wall mapping** case study. sPlotOpen is openly licensed. |
| **Valavi et al. 2019 (blockCV)** | Australia presence/absence, 243 pres / 257 abs (500 pts) + raster covariates; **binary** | Real **SDM / classification** case study (Brier/AUC optimism). Ships with blockCV. |
| **Linnenbrink et al. 2024 (kNNDM)** — GMD (code/data policy ⇒ Zenodo) | k-fold extension + study data | Secondary benchmark; pull exact DOI at write-up. |
| **Ploton et al. 2020** (Nat. Commun.) | Tropical forest AGB mapping data | Possible high-profile **wall-to-wall** case (check licence/extent before use). |

## The gap and the fix: a disease-mapping anchor (the author's domain)

None of the cited papers ship disease/NTD data. The author's flagship domain needs a
real **prevalence** example. Best open option (not from the cited papers, but the
right anchor):

- **Malaria Atlas Project** point-prevalence surveys via the **`malariaAtlas`** R
  package — real, point-referenced parasite-rate surveys, **openly licensed**,
  pulled via API (no redistribution needed). Gives a **Binomial/prevalence** response
  with genuinely **clustered/preferential** sampling — exactly the regime where
  leakage and range-estimation problems bite hardest. *This is the proposed flagship
  case study.* (Confirm current API/licence at use time.)
- Alternatives if needed: `PrevMap`/`geostatsp` example prevalence data; DHS cluster
  data (access-gated — avoid for a reproducible vignette).

## Coverage matrix (what each case study exercises)

| Case study | Response | Sampling design | Target / estimand | Source |
|---|---|---|---|---|
| Malaria prevalence (anchor) | Binomial / prevalence | clustered / preferential | wall-to-wall map, model-based | `malariaAtlas` |
| Soil moisture | Gaussian | logger network | interpolation, model-based | CAST `cookfarm` |
| Species richness | Gaussian | plot network | wall-to-wall, model-based | CAST `splotdata` |
| SDM presence/absence | Binary | opportunistic | wall-to-wall, model-based | blockCV (Australia) |
| Simulated benchmark | all | all (incl. probability sample ⇒ design-based) | all | NNDMpaper harness + our generator |

The **probability-sample / design-based** cell — needed to demonstrate *negative*
optimism (the Wadoux case, where spatial CV is over-pessimistic) — is covered by the
**simulation** (we can generate a true probability sample with known inclusion
density); finding a real, shareable probability-sample dataset is an open item.

## Licensing / distribution policy

- Do **not** vendor datasets into the package. Depend on `CAST`, `blockCV`,
  `malariaAtlas` via **Suggests** and load in vignettes/tests guarded by
  `requireNamespace()`.
- Simulated data are generated reproducibly in `data-raw/` (seeded); only compact
  derived artefacts (the fitted emulator + AOA reference) ship in `data/`.

## Status: Nigeria case study DONE
`malariaAtlas` confirmed open and working (`getPR`/`getShp`). Nigeria Pf gives **66
geolocated points** (DHS coords withheld by MAP). Case study built and run -- see
`docs/CASE-STUDY-NIGERIA.md` and `data-raw/case_study_nigeria.R`. Key finding: the
leakage is **co-located repeat surveys**, not autocorrelation (signal_prop ~ 0);
`fillDHSCoordinates()` + `rdhs` would add many more points for a richer run.

## Open data items (track in VISION §8)
- For a stronger headline run: geolocated DHS clusters (rdhs), narrow epoch/age
  filter, and covariates so the model beats IDW.
- Decide whether to use Ploton AGB data (licence/size).
- Find (or accept simulation-only) a real **probability-sample** dataset for the
  design-based / negative-optimism demonstration.

## Sources
- NNDMpaper / Zenodo: https://doi.org/10.5281/zenodo.6366985 ·
  https://github.com/carlesmila/NNDMpaper
- CAST datasets (`cookfarm`, `splotdata`):
  https://cran.r-project.org/web/packages/CAST/
- blockCV (Australia presence/absence): Valavi et al. 2019, MEE; package vignettes.
- kNNDM: https://doi.org/10.5194/gmd-17-5897-2024
- Malaria Atlas Project: https://malariaatlas.org (R package `malariaAtlas`).
