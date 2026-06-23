# spLeakage — Vision, Methodology & Roadmap

> **One sentence.** `spLeakage` is a *diagnostic* package: given spatial data and a
> train/test split (or an existing modelling workflow), it **detects, quantifies,
> and explains spatial information leakage**, estimates the resulting optimism in
> reported accuracy, and recommends the *appropriate* validation strategy for the
> user's sampling design and prediction target.

The end goal is **two coupled deliverables**: a strong methods paper and a mature,
well-tested R package. Every feature below is judged against both.

---

## 0. Why this package, and why now

Spatial predictive modelling (species distribution, disease mapping, digital soil
mapping, remote-sensing land cover) routinely reports accuracy from **random**
train/test splits. When observations are spatially autocorrelated, nearby points
land in both train and test sets, so the test set is not independent of the
training set. The reported error is then **optimistically biased**: it measures
near-interpolation, not the out-of-sample prediction the map is actually used for.

This is now a recognised, high-impact problem (Ploton et al. 2020 in *Nature
Communications* showed large-scale ecological maps are far less accurate than
reported once spatial structure is respected). A family of corrected
cross-validation methods exists. **But there is a gap and a live controversy:**

- **The gap.** Existing tools *generate corrected folds* (CAST, blockCV,
  spatialsample, mlr3spatiotempcv). None take an *arbitrary split a user already
  made* and tell them *how much it leaks* and *how inflated their numbers are*.
  There is no leakage **diagnostic / audit / scorecard**.

- **The controversy.** Milà et al. (2022) and Linnenbrink et al. (2024) argue you
  should match the cross-validation geometry to the prediction geometry (NNDM /
  kNNDM). Wadoux et al. (2021) argue spatial CV "has no theoretical underpinning"
  and that for **probability samples** design-based inference (i.e. ordinary
  random validation) is the *correct* unbiased approach, with spatial CV
  introducing *pessimistic* bias. **Both are right, in different regimes.** A
  practitioner has no tool that tells them which regime they are in and therefore
  which validation is correct.

`spLeakage` occupies exactly this gap: not "here are better folds," but
**"here is what is wrong with your evaluation, by how much, and what you should do
instead — given your design and your prediction target."**

---

## 1. Competitive landscape (what already exists — do not rebuild)

| Tool | What it does | Why it is *not* spLeakage |
|------|--------------|---------------------------|
| **CAST** (Meyer, Milà) | NNDM, kNNDM, Area of Applicability, geodist plots, `global_validation` | Generates folds & quantifies AOA; not a leakage audit of a user's own split; takes a strong "spatial CV is right" stance |
| **blockCV** (Valavi et al. 2019) | Spatial/environmental blocking folds, range estimation | Fold generator only |
| **spatialsample** (tidymodels) | `spatial_block_cv`, `spatial_nndm_cv`, clustering folds | Fold generator within tidymodels |
| **mlr3spatiotempcv** | Spatial/temporal resampling for mlr3 | Fold generator within mlr3 |
| **sperrorest** (Brenning) | Spatial error estimation by spatial CV | Runs spatial CV; not a diagnostic of leakage magnitude |

**Design principle:** `spLeakage` *consumes* the outputs of these tools and audits
them; it interoperates rather than competes. We will provide adaptors for
`rsample`/`rset`, `mlr3` `ResamplingSpCV*`, `caret` index lists, `blockCV` folds,
and plain index vectors / `sf` objects.

---

## 2. Methodological foundations to build on

These are the established ideas the package stands on. Each maps to a feature.

1. **NN-distance distribution matching (NNDM / kNNDM).** Milà et al. 2022 (MEE);
   Linnenbrink et al. 2024 (GMD). Core insight: the *test→training* nearest-
   neighbour distance distribution during CV should match the *prediction→sample*
   NN-distance distribution at deployment. Fold quality is measured by the
   Wasserstein W statistic between the two ECDFs. **We reuse this comparison as a
   diagnostic rather than a fold-construction objective.**

2. **Design-based vs model-based inference.** Brus & de Gruijter; Wadoux et al.
   2021; de Bruin et al. 2022. For *probability samples*, design-based validation
   gives unbiased accuracy and random CV is appropriate; spatial CV can be
   *pessimistically* biased. For *non-probability / clustered* samples, this breaks
   down. **This is the backbone of our recommendation engine.**

3. **Buffered / dead-zone LOO CV.** Le Rest et al. 2014; Telford & Birks 2009;
   Ploton et al. 2020. Exclude training points within a buffer of each test point;
   the buffer is naturally set by the autocorrelation range.

4. **Block / clustered CV.** Roberts et al. 2017 (Ecography — the canonical
   reference); Valavi et al. 2019 (blockCV). Block size driven by the empirical
   variogram range.

5. **Variogram-based length scales.** Empirical/fitted variogram → practical range
   = the natural distance over which leakage operates. Anisotropy → directional
   leakage.

6. **Area of Applicability & feature-space dissimilarity.** Meyer & Pebesma 2021
   (MEE). Leakage is not only geographic: test points close to training points *in
   covariate space* also leak. We generalise leakage to feature space.

7. **Residual spatial autocorrelation.** Moran's I / correlograms on model
   residuals as evidence that dependence remains (a leakage symptom).

8. **Effective sample size under autocorrelation.** Clifford–Richardson; Dale &
   Fortin. Autocorrelation inflates the *effective* training–test overlap; useful
   for framing optimism.

---

## 3. Novel contributions (the paper's spine)

### The thesis: leakage is only well-posed relative to *design × estimand × target*

The central claim that unifies the package and the paper:

> **"Spatial leakage" and "optimism" are not properties of a dataset or a split in
> isolation. They are only well-defined *relative to a declared sampling design, a
> declared estimand, and a declared prediction target.* The same split can be
> leaking, correct, or pessimistic depending on these three declarations.**

This reframing is what makes the contribution coherent and defensible, and it
*reconciles* the Milà-vs-Wadoux debate instead of taking a side:

- **Sampling design** (probability / clustered / convenience) decides whether
  random validation is unbiased (Wadoux et al.) or optimistic (Milà et al.).
- **Estimand** — *population-mean map accuracy over a region* (design-based) vs
  *conditional predictive skill at a location* (model-based) — decides what
  "accuracy" even means. These are different questions; conflating them is a common
  error and a frequent reviewer objection. We separate them explicitly.
- **Prediction target** (within-sample interpolation / wall-to-wall map /
  new region) sets the geometry the validation must imitate.

So the package's job is not "detect leakage" in the abstract; it is: **given the
user's declared (design, estimand, target), make leakage well-posed, then diagnose
it, quantify it, and prescribe the matched validation.** The three contributions
below are instantiations of this single thesis, not three independent tools — and
they are deliberately *conditioned on* the declarations, which removes the circular
dependency a reviewer would otherwise attack (see §8). Honest novelty assessment:
the **load-bearing** novelty is the optimism emulator (C2b) plus this reconciliation
framing; C1 and C3 are valuable but must be sold as the workflow/decomposition and
the reconciliation, not as a new statistic or an expert lookup table.

### C1 — A Spatial Leakage Index (SLI) for an *arbitrary* split *(diagnostic)*
A signed, interpretable leakage score for a split the user *already has*, from a
divergence between:
- the **observed** test→train nearest-neighbour ECDF of the user's split, and
- the **target** prediction→sample nearest-neighbour ECDF implied by the declared
  prediction target.

Design notes (each answers a known objection — see §8):
- **Not just the kNNDM W statistic renamed.** The contribution is the *post-hoc
  audit of an arbitrary split* (which fold-generators cannot do), the **per-fold and
  per-point decomposition**, and the **spatial map** of leakage ("where does my
  split leak?"). Sell the workflow, not the metric.
- **Signed / directional**, not a symmetric divergence: it must distinguish
  *optimistic* leakage (test closer than deployment) from *pessimistic* anti-leakage
  (test farther) — the Wadoux failure mode — via stochastic-dominance direction.
- **Covariance-aware variant**, because geographic distance ≠ statistical
  dependence. Offer a dependence-scaled distance (h/range, or the fitted ρ(h)) so
  the index measures *correlation* leakage under anisotropy/non-stationarity, with
  the plain-distance version as the cheap default. Show they agree under
  stationarity and diverge otherwise.
- **Intrinsic normalisation** by the variogram range (not the map extent), so the
  index is comparable across datasets; state and test the invariance claimed.
- **Validated via its monotone link to true optimism (C2)** — C1 is, by design,
  subordinate to C2 rather than a standalone unvalidated heuristic.

### C2 — An optimism-bias estimator *(the number users want)*
Translate leakage into the practitioner's currency: *"your reported RMSE is likely
X% optimistic (90% interval a–b%)."*

**Estimand fix (defuses the single most dangerous objection).** Optimism is defined
**relative to the validation that is *correct for the declared design × estimand ×
target* — not relative to NNDM.** For a probability sample with a population-mean
estimand, the "correct" baseline is design-based/random, so NNDM-vs-random would
register as *pessimism (negative optimism)*, which the tool will report as such.
This conditioning is what makes "optimism" well-posed and stops a Wadoux-camp
reviewer from detonating the paper.

Two routes:
- **(a) Empirical/refit:** re-evaluate the user's model under the *design-matched*
  scheme and report the gap. Model-conditional (optimism depends on the learner) —
  reported as such. Accurate but requires refits.
- **(b) Cheap surrogate emulator (the load-bearing novelty):** a pre-trained
  predictor mapping `(autocorrelation range, sampling clustering, split geometry /
  signed SLI, sample size, signal-to-noise, design, target)` → expected optimism,
  **calibrated once via a large simulation study** shipped with the package; gives
  an instant estimate with *no refitting*. Hardening (see §8):
  - **Richer generators**, not just stationary GRFs: non-Gaussian responses
    (spatial GLMM for counts/prevalence), covariate trends, non-stationarity, and
    *preferential sampling* — then report **out-of-distribution** performance
    against held-out real datasets, not only within-simulation.
  - **The emulator gets its own area-of-applicability check** and refuses/flags
    extrapolated inputs — the package must not commit the sin it polices.
  - **Propagate input-estimation uncertainty** (range etc. are hard to estimate
    from exactly the clustered data that leaks) into the reported interval; show
    interval coverage.

### C3 — A design-aware recommendation engine *(reconciles the controversy)*
Turns the thesis into actionable, *conditional* guidance. **Estimand first, then
design × target** — because design-based and model-based validation answer different
questions:

Step 1 — elicit the **estimand**: population-mean map accuracy over a region
(design-based) *or* conditional predictive skill at locations (model-based)?

Step 2 — given the estimand, recommend over **design × target**:

| Sampling design \ target | Within-sample interpolation | Wall-to-wall mapping | New-region extrapolation |
|---|---|---|---|
| Probability / design-based | Random CV (design-based) | Design-based estimator | Caution; AOA |
| Clustered / preferential | NNDM / buffered | NNDM / kNNDM | Block + AOA |
| Convenience / opportunistic | Buffered LOO | kNNDM | Block + AOA, flag risk |

Design rules (each answers a known objection — see §8):
- **Design-basedness is NOT inferred from geometry.** Whether data are a probability
  sample is a fact about collection (known inclusion probabilities), unrecoverable
  from coordinates — a *designed* cluster sample and a *convenience* sample can have
  identical point patterns. Geometry (Clark–Evans, Ripley's K, inclusion-density)
  only **flags risk**; the design is an **elicited input**. The engine is decision
  *support*, not decision *automation*.
- **Every cell is backed** by a citation or a row of the simulation study — no
  unsupported "expert-system" cells.
- **Conditional, not neutral-while-secretly-partisan:** guidance is phrased "*if*
  your estimand is X under framework Y, *then*…", surfacing the assumption rather
  than hiding a choice of side.
- Uncovered designs (spatiotemporal, citizen-science, informative/preferential
  intensity, multi-source fusion) are **flagged and referred**, not forced into the
  3×3.
- Still states, loudly, **when spatial CV is NOT appropriate** — the nuance no
  current tool surfaces.

### C4 — Multi-channel leakage detection
Beyond geographic distance:
- **Feature-space leakage** (covariate-space NN dissimilarity / AOA-style).
- **Temporal & spatiotemporal leakage** (space–time NN distances).
- **Grouped / duplicated-location leakage** (same site, repeated measures, plots
  split across folds).
- **Preprocessing leakage** (scaling/imputation/feature-selection/target encoding
  fit *before* splitting; spatial covariate extraction with neighbourhoods/buffers
  that straddle the split).

### C5 — Workflow & code auditing
Inspect *objects* (`rsample`, `mlr3` resamplings, `caret` indices, `blockCV`) and,
optionally, *scripts* (static heuristics) to flag leakage-prone patterns:
`sample()`/`initial_split()` on spatial data, pre-split transforms, missing seeds,
undocumented CRS. Produces a reproducibility-style scorecard. (Static code analysis
is a stretch goal; object auditing is core.)

---

## 4. Full feature catalogue (everything we'd like it to do)

Grouped by subsystem. `[Px]` = target phase (see §6). Not a commitment to build
all — a menu to pull from.

### A. Data ingestion & geometry `[P1]`
- Accept `sf`, `data.frame`+coords, `terra`/`SpatRaster` prediction grids,
  `sftime`/space-time.
- CRS handling: detect missing/geographic CRS, warn, compute geodesic vs projected
  distances correctly.
- A **prediction-target specification** API: `target = "grid"` (wall-to-wall),
  `"newdata"`, `"interpolation"`, or an explicit set of prediction locations.

### B. Distance & dependence engine `[P1]`
- Fast NN-distance computation (`nabor`/`FNN`/`RANN` k-d trees) for
  test→train, prediction→sample, train→train.
- Empirical variogram + automatic practical-range estimation; anisotropy detection.
- Clustering / sampling-design diagnostics (NN index / Clark–Evans, Ripley's K/L,
  inclusion-density heuristics).
- Effective sample size estimate.

### C. Core leakage diagnostics `[P1–P2]`
- `detect_leakage(split, data, target)` → leakage report object.
- Spatial Leakage Index (SLI) with fold-level decomposition (C1).
- Distance-threshold view: "% of test points within d km of a training point,"
  sweepable over d; default d = estimated range.
- ECDF overlays (observed vs target NN-distance), the central NNDM-style plot.
- Per-test-point leakage scores → mappable.

### D. Optimism quantification `[P2–P3]`
- `estimate_optimism(...)` empirical route (refit under controlled scheme) (C2a).
- Surrogate emulator route + bundled GRF simulation study + calibration (C2b).
- Uncertainty intervals on the optimism estimate.

### E. Recommendation engine `[P2–P3]`
- `recommend_validation(data, target, design = NULL)` → ranked strategies +
  rationale + *what to avoid* + when spatial CV is inappropriate (C3).
- Auto-generate corrected folds by *delegating* to blockCV/CAST/spatialsample.

### F. Multi-channel leakage `[P3–P4]`
- Feature-space leakage (C4), temporal/spatiotemporal (C4), grouped/duplicate (C4),
  preprocessing/pipeline leakage (C4).

### G. Auditing & reporting `[P2–P4]`
- `audit_workflow(rset | resampling | caret_index | blockCV)` (C5).
- `audit_script(path)` static heuristics (stretch) (C5).
- `report_leakage(...)` → reproducible HTML/PDF (Quarto) scorecard for journal
  submission: leakage grade, optimism estimate, recommendation, figures.
- Residual diagnostics: Moran's I / correlogram on supplied residuals.

### H. Interoperability & ergonomics `[ongoing]`
- Adaptors: tidymodels (`rsample`/`spatialsample`/`workflows`), mlr3, caret.
- S3 print/plot/summary/autoplot methods; tidy `broom`-style output.
- Sensible defaults; everything driven by one estimated length scale by default.

---

## 5. Public API sketch (subject to change)

```r
library(spLeakage)

# 1. Describe the deployment scenario
tgt <- prediction_target(grid = pred_grid)          # or newdata=, or "interpolation"

# 2. Diagnose an existing split
lk  <- detect_leakage(split = my_rsplit, data = obs_sf, target = tgt)
lk                                                   # printed scorecard
plot(lk)                                             # NN-distance ECDF overlay
sli(lk)                                              # the Spatial Leakage Index

# 3. Quantify the damage
opt <- estimate_optimism(lk, model = my_model, metric = "rmse")
opt                                                  # "RMSE ~14% optimistic (8–22%)"

# 4. Get told what to do instead
recommend_validation(obs_sf, target = tgt, design = "clustered")

# 5. Audit a whole workflow / produce a submission-ready report
audit_workflow(my_resampling)
report_leakage(lk, opt, file = "leakage_report.html")
```

---

## 6. Phased roadmap (build order)

Each phase ends in something publishable/usable. We pursue them one at a time.

- **P0 — Scaffolding & design** ✅ *done*: package skeleton, vision doc, paper
  outline, method specs (`METHOD-SLI.md`, `METHOD-EMULATOR.md`), dataset survey.
- **P1 — Geometry & core diagnostic (MVP)** ✅ *done + polished*: distance/variogram
  engine, prediction-target API, `detect_leakage()` + signed `SLI_d`/`SLI_rho` +
  crossing decomposition + per-point/fold map + `sf` adaptor; covariance-aware
  **anisotropy** and **Monte-Carlo SLI uncertainty**. *Milestone met: audits any
  split and reports a (uncertainty-quantified) leakage score.*
- **P2 — Optimism (empirical) + recommendation engine** 🔶 *core done*:
  `estimate_optimism()` (design-matched controlled route, model-agnostic),
  `recommend_validation()` (estimand-first decision engine, elicited design),
  `audit_workflow()`, `report_leakage()` scorecard. ✅ *Milestone met:* end-to-end
  real disease-mapping case study (`data-raw/case_study_nigeria.R`,
  `docs/CASE-STUDY-NIGERIA.md`) on Nigeria MAP prevalence + the usage vignette.
  The case study attributed the leakage to **co-located repeat surveys** (SLI_rho
  +0.21 -> +0.006 on deduplication), exercising the geographic/geodesic path.
- **P3 — Simulation study + surrogate emulator** ✅ *study done*: full harness
  (`data-raw/simulate_optimism.R`) -- Matern fields x smoothness x signal x {random,
  clustered, preferential} x n x {gaussian, poisson, binomial} x {idw, rf, gam} x
  {grid, interpolation}, ~900 denoised rows (450 configs x 4 reps), ground-truth
  `optimism_rel`. `predict_optimism()` emulator: **gradient-boosted (xgboost) point
  estimate + config-grouped split-conformal intervals**, instance-based Meyer-Pebesma
  AOA guard (`R/emulator.R`), shipped as raw bytes in `R/sysdata.rda`. **Paper-scale
  run (1000 configs x 6 reps, ~2000 denoised rows): held-out R2 ~= 0.76, in-AOA 97%,
  interval coverage ~94%; `cor(SLI_rho,optimism) ~= 0.67`.** Design ordering (clustered
  > random optimism) reconciles Milà-vs-Wadoux; figures in `paper-figures/`.
  *Remaining for the paper: scale configs/reps, add covariate-trend + spatial-GLMM
  realism, normalised conformal to close 85%->90% coverage, formal out-of-distribution
  validation on the real case-study datasets.*
- **P4 — Multi-channel leakage + polish** 🔶 *mostly done*: **grouped/duplicated-
  location** (`detect_group_leakage`/`group_kfold`), **feature-space / covariate-AOA**
  (`detect_feature_leakage`), and **temporal lookahead** (`detect_temporal_leakage`/
  `temporal_kfold`) channels implemented + tested; `rsample` adaptor; examples on key
  functions; `NEWS.md`, `cran-comments.md`. R CMD check 0/0/0. *Remaining:
  preprocessing/pipeline (static-code) leakage -- a stretch goal -- additional real
  datasets, and the paper write-up.*

---

## 7. The paper (target & framing)

- **Working title:** *"Detecting and quantifying spatial information leakage in
  predictive modelling: a diagnostic framework and the `spLeakage` R package."*
- **Likely venues:** *Methods in Ecology and Evolution*, *Geoscientific Model
  Development* (software+method, where NNDM/kNNDM landed), *Journal of Statistical
  Software*, or *Ecography*.
- **Narrative:** problem (leakage & inflated maps) → the controversy (Milà vs
  Wadoux) → our reframing as *diagnosis matched to design × target* → SLI + optimism
  emulator (methods) → simulation validation + real case studies (disease mapping,
  ecology, soil) → the package.
- **Evidence plan:** (i) GRF simulation across autocorrelation range × sampling
  clustering × prediction target showing SLI and the emulator track true optimism;
  (ii) 2–3 real datasets re-analysed, reporting how inflated the original numbers
  were; (iii) head-to-head showing spLeakage's design-aware advice avoids *both*
  the optimism of random CV *and* the pessimism Wadoux warns about.

---

## 8. Reviewer objection register (design must answer these)

A pre-emptive Reviewer 2 pass. Severity: **fatal** (can sink the paper) ·
**wounding** (forces major revision) · **minor**. Each row has the objection and the
**hardening move** the design commits to. These are the constraints that steer P1+.

### Cross-cutting (the one to watch)
- **[fatal] Circular entanglement.** C2's "optimism" depends on C3's "correct
  target," which depends on a design fact C3 can't infer. → **Resolve by the §3
  thesis:** define leakage/optimism *relative to declared (design, estimand,
  target)*. The reconciliation becomes the contribution; the circularity disappears
  because the declarations are inputs, not inferences.
- **[wounding] Novelty concentration.** A reviewer may read C1 as incremental and
  C3 as an expert system. → Be explicit that the load-bearing novelty is the
  emulator (C2b) + the reconciliation framing; position C1/C3 accordingly.

### C1 — Spatial Leakage Index
- **[wounding] "Just kNNDM's W renamed."** → Contribution is the post-hoc audit of
  an *arbitrary* split + per-fold/per-point **decomposition** + spatial **map**, not
  the statistic. Sell the workflow.
- **[fatal-ish] "Geographic distance ≠ statistical dependence."** Anisotropy /
  non-stationarity / covariate-mediated dependence break distance-based leakage. →
  Provide a **covariance-aware variant** (h/range or fitted ρ(h)); show agreement
  with plain distance under stationarity, divergence otherwise.
- **[wounding] "Conflates leakage with pessimism."** Symmetric divergence flags both
  optimism and anti-leakage. → Make SLI **signed/directional** (stochastic
  dominance).
- **[wounding] "Normalisation not invariant."** → Normalise by **variogram range**,
  not extent; state and test the invariance claim.
- **[wounding] "Validated against what?"** → C1 validated by its **monotone link to
  true optimism (C2)**; accept C1 is subordinate to C2.

### C2 — Optimism-bias estimator
- **[fatal] "Your optimism assumes NNDM is the truth."** For a probability sample,
  random CV may be *correct* and NNDM *pessimistic*, so the "gap" can be negative. →
  **Define optimism relative to the design-matched correct validation**, not NNDM;
  report pessimism (negative optimism) when that's what occurs.
- **[wounding] "GRF simulations are a toy."** → Enrich generators: spatial-GLMM
  (counts/prevalence), covariate trends, non-stationarity, **preferential
  sampling**; report **out-of-distribution** performance on held-out real data.
- **[wounding] "Does the emulator have its own AOA problem?"** → Ship an
  **area-of-applicability check on the emulator itself**; flag/refuse extrapolation.
- **[wounding] "Inputs are the hardest things to estimate."** Range etc. are badly
  estimated from clustered data. → **Propagate input-estimation uncertainty** into
  the interval; show coverage; sensitivity analysis.
- **[minor] "Optimism is model-dependent."** → Report empirical route as
  model-conditional; emulator predicts *expected* optimism for a stated model class.

### C3 — Design-aware recommendation engine
- **[fatal] "Clustering tests detect pattern, not design."** Designed cluster
  samples and convenience samples can look identical. → **Do not infer
  design-basedness from geometry**; make design an **elicited input**; geometry only
  flags risk. Decision *support*, not automation.
- **[fatal-ish] "Design-based and model-based target different estimands."**
  Population-mean accuracy ≠ conditional predictive skill. → **Elicit the estimand
  first**, then recommend; never recommend across estimands.
- **[wounding] "Matrix is expert opinion dressed as method."** → **Back every cell**
  with a citation or a simulation-study row; no unsupported cells.
- **[wounding] "Claims neutrality but adopts Wadoux/Brus as judge."** → Phrase all
  guidance **conditionally** ("if estimand X under framework Y, then…").
- **[minor] "Real designs don't fit a 3×3."** → **Flag and refer** uncovered designs
  (spatiotemporal, citizen-science, preferential intensity, multi-source).

### Resolved & locked (see dedicated specs)
- **Signed-SLI math** — LOCKED in `docs/METHOD-SLI.md` v1: signed area between
  NN-distance ECDFs (`SLI_d`, range-normalised) + the dependence-form `SLI_ρ ∈
  [−1,1]` (headline), crossing decomposition, per-point/fold leakage map, reduction
  to NNDM `W`, design-relative sign convention.
- **Emulator generator** — LOCKED in `docs/METHOD-EMULATOR.md` v1: ground-truth
  optimism (not NNDM-relative; can be negative), 10-factor generator with
  preferential sampling + non-Gaussian (spatial-GLMM) responses + model-class factor,
  estimated (non-oracle) features = the C1 outputs, self-AOA guard, OOD validation.
- **Datasets** — surveyed in `docs/DATA.md`: cited papers give simulation harnesses
  (NNDMpaper/Zenodo) + real continuous (`cookfarm`, `splotdata`) and binary (blockCV
  Australia) cases; **no disease data in the cited papers** → anchor with Malaria
  Atlas Project via `malariaAtlas` (open). Depend via Suggests, don't vendor.

### Still-open design questions (no committed answer yet)
- How far to push static code analysis (fragile) vs object auditing (robust).
- Default distance metric & CRS policy (geodesic vs projected) — leaning geodesic for
  geographic CRS; confirm in P1.
- A real, shareable **probability-sample** dataset for the design-based /
  negative-optimism demonstration (else simulation-only).
- `malariaAtlas` licence/API still open + specific country/year; Ploton AGB licence.

---

## 9. Key references (verify before citing in the paper)

- Milà, Mateu, Pebesma, Meyer (2022) *Nearest neighbour distance matching LOO CV
  for map validation.* **Methods in Ecology and Evolution** 13:1304–1316.
- Linnenbrink, Milà, Ludwig, Meyer (2024) *kNNDM CV: k-fold nearest-neighbour
  distance matching CV for map accuracy estimation.* **Geoscientific Model
  Development** 17:5897–5912.
- Wadoux, Heuvelink, de Bruin, Brus (2021) *Spatial cross-validation is not the
  right way to evaluate map accuracy.* **Ecological Modelling** 457:109692.
- de Bruin, Brus, Heuvelink, van Ebbenhorst Tengbergen, Wadoux (2022) *Dealing with
  clustered samples for assessing map accuracy by cross-validation.* **Ecological
  Informatics.**
- Ploton et al. (2020) *Spatial validation reveals poor predictive performance of
  large-scale ecological mapping models.* **Nature Communications** 11:4540.
- Roberts et al. (2017) *Cross-validation strategies for data with temporal,
  spatial, hierarchical, or phylogenetic structure.* **Ecography** 40:913–929.
- Meyer & Pebesma (2021) *Predicting into unknown space? Estimating the area of
  applicability of spatial prediction models.* **Methods in Ecology and Evolution.**
- Valavi, Elith, Lahoz-Monfort, Guillera-Arroita (2019) *blockCV.* **MEE.**
- Le Rest, Pinaud, Monestiez, Chadoeuf, Bretagnolle (2014) *Spatial leave-one-out
  cross-validation.* **Global Ecology and Biogeography.**
- Meyer, Milà, Ludwig, Linnenbrink, Schumacher (2024) *The CAST package for
  training and assessment of spatial prediction models in R.* (arXiv:2404.06978.)

*(Authors/years checked against publisher pages June 2026; confirm volume/pages at
write-up time.)*
