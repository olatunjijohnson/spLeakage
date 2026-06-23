# Case study — spatial leakage in Nigeria malaria prevalence mapping

Real-data demonstration of `spLeakage` on **Malaria Atlas Project** *P. falciparum*
parasite-rate point surveys for Nigeria (open data, pulled via `malariaAtlas`).
Reproducible script: `data-raw/case_study_nigeria.R`. This also serves as an
**out-of-distribution test of the optimism emulator** on real data.

## Data

`malariaAtlas::getPR("Nigeria", "Pf")` returns 2,672 survey records, but MAP cannot
share DHS cluster coordinates, so **66 geolocated points** (58 unique locations)
remain, spanning the country (lon 3.0–13.1, lat 4.6–12.5), prevalence mean 0.41,
surveys 1985–2018. Small, clustered, and pooled across decades — a realistic,
messy compilation. The response is the empirical parasite rate `pr`. This is the
first exercise of the package's **geographic-CRS / geodesic-distance** path.

## What the workflow found

Naive baseline: a random 10-fold CV split (what most papers do).

| Step | Result |
|------|--------|
| **1. `detect_leakage`** | `SLI_rho = +0.212` (grade **D**, optimistic), `c_obs = 0.212` vs `c_pred = 0.000` |
| **2. `estimate_optimism`** (block control) | random-CV RMSE 0.265 vs controlled 0.270 → **+1.8% optimistic** |
| **3. `predict_optimism`** (emulator) | **refused — outside area of applicability** (defers to the empirical route) |
| **3b. attribution** (deduplicate) | SLI_rho `+0.212` (n=66) → **`+0.006` (n=58 deduplicated)** |
| **4. `recommend_validation`** | convenience × grid → **spatial CV (NNDM/block); avoid random**; clustering flag NN index 0.58 |
| **5. `audit_workflow`** | grade D; **8 duplicated coordinates flagged**; CRS documented ✓ |

## The story the package tells

This case is more instructive than a generic "random CV was X% optimistic," because
the tools **diagnose the *type* and *source* of the leakage**:

1. **There is little continuous spatial autocorrelation here.** The fitted variogram
   gives `signal_prop ≈ 0`: pooling 33 years of surveys (varied methods, ages,
   epochs) swamps the spatial signal with heterogeneity. So the leakage is *not*
   range-based.

2. **The leakage is co-location.** With `signal_prop ≈ 0`, the only way `c_obs`
   (mean retained correlation of test→train) reaches 0.212 is exact co-locations:
   ~21% of test points share a site with a training point (repeat surveys at the
   same place, split across folds → distance 0 → ρ = 1). The **degenerate CI
   `[+0.212, +0.212]`** is *correct* — co-location is fixed geometry with no
   parametric uncertainty, unlike an estimated range.

3. **Deduplication proves it.** Averaging prevalence per unique location collapses
   `SLI_rho` from **+0.212 to +0.006** — the leakage was entirely co-located repeat
   surveys. The data-hygiene audit independently flags the **8 duplicated
   coordinates** that cause it.

4. **The consequence is modest (+1.8%)** because co-located records are from
   *different years* (same site, different prevalence), so even an exact-location
   training neighbour predicts the test value imperfectly — and the overall spatial
   signal is weak. This is exactly the **C1-vs-C2 distinction the package is built
   on**: the geometry leaks (SLI), but the *consequence* (optimism) depends on how
   much spatial information actually exists. Reporting only one would mislead.

5. **The emulator behaves honestly.** It **refuses** this query (n = 66 is below its
   training envelope of 150–600, and the real feature profile is unlike the
   simulation), correctly flagging out-of-distribution and deferring to the
   empirical `estimate_optimism()` — the area-of-applicability guard doing its job
   on real data.

## Takeaways

- `spLeakage` ran end-to-end on real, geographic prevalence data and produced a
  coherent, *attributed* diagnosis, not just a number.
- The flagship lesson is **leakage-source attribution**: here the culprit is
  **duplicated-location records**. This motivated the now-implemented C4 channel:
  `detect_group_leakage()` reports **14/66 (21.2%)** of test points leaking via a
  shared location -- numerically identical to the SLI's `c_obs = 0.212`, confirming
  the two diagnostics agree -- and `group_kfold()` drives it to **0%**.
- Honest limitations: n = 66 (DHS coordinates withheld); decades-pooled surveys give
  weak resolvable autocorrelation. A cleaner run would filter to a narrow epoch/age
  and ideally use geolocated DHS clusters (via `malariaAtlas::fillDHSCoordinates()`
  + `rdhs`), and add covariates so the predictive model is stronger than IDW.

## Reproduce

```r
devtools::load_all(); source("data-raw/case_study_nigeria.R")
```
Requires `malariaAtlas` (network) and `sf`.
