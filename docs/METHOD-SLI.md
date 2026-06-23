# Locked spec — the Signed Spatial Leakage Index (SLI)

Status: **locked** (v1). Supersedes the sketch in `VISION.md` §3/C1. Changes go
through a version bump here. Answers objections C1-① … C1-⑤ in `VISION.md` §8.

## 0. Notation

- `S` = sample (observed) locations, `|S| = n`.
- A CV scheme assigns, to each observation `t`, a training set `Tr(t)` used to
  predict `t` when `t` is held out. (LOO: `Tr(t) = S \ {t}` minus any buffer.
  k-fold: `Tr(t) = S \ fold(t)`. Buffered methods enlarge the excluded set.)
- **Prediction target** `P` = the locations the model is actually deployed on,
  derived from the user's *declared* target (see §5). `|P| = m`.
- `d(·,·)` = the chosen ground distance (geodesic for geographic CRS, projected
  otherwise; anisotropic/dependence variants in §6).

Two nearest-neighbour functions, the objects everything is built on:

- **Observed (CV) NN distance:** `gobs(t) = min_{s ∈ Tr(t)} d(t, s)`.
  ECDF `Ĝobs(r) = (1/n) Σ_t 1{gobs(t) ≤ r}`.  *(This is the NNDM `Ĝ` function.)*
- **Target (deployment) NN distance:** `fpred(p) = min_{s ∈ S} d(p, s)`.
  ECDF `Ĝpred(r) = (1/m) Σ_p 1{fpred(p) ≤ r}`.

Intuition: leakage ⇔ test points are *closer* to their training data during CV
than deployment points are to the sample ⇔ `Ĝobs` shifted to **smaller** distances
than `Ĝpred` ⇔ `Ĝobs(r) ≥ Ĝpred(r)`.

## 1. The distance-form index `SLI_d` (cheap, model-free default)

Define the **signed area between the ECDFs**:

```
A = ∫_0^∞ ( Ĝobs(r) − Ĝpred(r) ) dr
```

By the layer-cake identity `E[X] = ∫_0^∞ (1 − F(x)) dx` for non-negative `X`,

```
A = ∫_0^∞ [ (1 − Ĝpred) − (1 − Ĝobs) ] dr
  = mean(fpred over P) − mean(gobs over S)
  = (mean deployment NN distance) − (mean CV NN distance).
```

So **A is literally how much farther deployment reaches than your CV pretended**,
in distance units. `A > 0` ⇒ optimistic leakage; `A < 0` ⇒ pessimism.

Normalise by the **intrinsic dependence length** — the variogram practical range
`φ` (NOT the map extent; this is the invariance fix for objection C1-④):

```
SLI_d = A / φ            (signed; positive = optimistic leakage)
```

Interpretation: the CV-vs-deployment reach gap *as a fraction of the correlation
range*. Beyond `φ` there is no dependence, so a distance gap there causes no
leakage — which is exactly why `φ` is the right scale.

### Crossing decomposition (don't let the integral cancel)
ECDFs can cross (leakage at short range, pessimism at long range). Report the parts:

```
A+ = ∫ max(Ĝobs − Ĝpred, 0) dr        (optimistic mass)
A− = ∫ max(Ĝpred − Ĝobs, 0) dr        (pessimistic mass)
A  = A+ − A−        (signed, = numerator of SLI_d · φ)
W  = A+ + A−        (total = the unsigned NNDM/kNNDM Wasserstein-W statistic)
```

- **`SLI_d` reduces (in magnitude) to the NNDM/kNNDM `W`** when the ECDFs don't
  cross (`A− = 0` ⇒ `|A| = W`). This is the explicit bridge to the established
  literature — and the reason the contribution is the *signing + decomposition +
  post-hoc-on-arbitrary-split use*, not a new statistic (objection C1-①).
- Report a **directionality ratio** `δ = A / W ∈ [−1, 1]`: `δ ≈ +1` pure leakage,
  `δ ≈ −1` pure pessimism, `δ ≈ 0` mixed/crossing → tell the user to look at the
  ECDF plot rather than trust a single number.

## 2. The dependence-form index `SLI_ρ` (headline; covariance-aware)

Distance ≠ statistical dependence (objection C1-②). Replace raw distances by the
**fitted spatial correlation** `ρ(h)` (from the variogram: Matérn/exponential, with
`ρ(0)=1`, `ρ(∞)=0`). Define the *retained spatial information*:

```
c_obs  = (1/n) Σ_t  ρ( gobs(t) )         mean correlation CV keeps to training
c_pred = (1/m) Σ_p  ρ( fpred(p) )        mean correlation available at deployment
SLI_ρ  = c_obs − c_pred ∈ [−1, 1]        (signed; positive = optimistic leakage)
```

Properties:
- **Bounded, dimensionless, signed.** `ρ ∈ [0,1]` ⇒ `SLI_ρ ∈ [−1,1]` with no ad-hoc
  normalisation. `+` = CV retains *more* spatial information than deployment has =
  optimistic; `−` = pessimism.
- **Covariance-aware:** under anisotropy/non-stationarity use a directional/local
  `ρ` (see §6); the index then tracks *dependence*, not geometry.
- **Mechanistic link to optimism (gives C1 its referent, objection C1-⑤).** In a GP,
  the predictive-variance reduction at a location from its neighbours grows with the
  retained neighbour correlation. The excess `c_obs − c_pred` is, to first order,
  the excess variance-reduction CV enjoys over deployment — i.e. the *mechanism* of
  optimism. Heuristic to validate in simulation:
  `optimism ≈ increasing function of ( SLI_ρ × spatial-signal proportion )`.
- **Agreement with `SLI_d`:** under isotropic stationarity `SLI_ρ` is a monotone
  transform of `A` (both driven by the NN-distance shift). We *verify* this in
  simulation and *expect divergence* under anisotropy/non-stationarity — that
  divergence is itself a reported diagnostic.

`SLI_ρ` is the headline index; `SLI_d` is the model-free fallback when no reliable
variogram is available.

## 3. Spatial decomposition (the "where does it leak?" map)

Per-observation leakage contribution:

```
ℓ(t) = ρ( gobs(t) ) − c_pred           (excess retained correlation at t)
SLI_ρ = (1/n) Σ_t ℓ(t)
```

- `ℓ(t) > 0`: point `t` is evaluated more easily than deployment → a leaking point.
- **Map `ℓ(t)`** over space → the leakage hotspot map (a headline figure).
- **Per fold:** `SLI_ρ(fold k) = mean_{t ∈ fold k} ℓ(t)` → which folds leak.

This per-point/per-fold decomposition (which fold-generators cannot produce for an
arbitrary split) is the substance of contribution C1.

## 4. Sign convention & the link to design (this is the crux)

| sign | name | meaning |
|------|------|---------|
| `SLI > 0` | optimistic leakage | CV easier than deployment; reported error too low |
| `SLI ≈ 0` | matched | CV imitates deployment |
| `SLI < 0` | pessimism | CV harder than deployment; reported error too high |

**Whether a non-zero SLI is a *problem* depends on the declared design × estimand ×
target (the §3 thesis).** For a probability sample with a population-mean estimand,
the deployment target `P` is the *population under the inclusion density*, and a
negative SLI versus a spatially-blocked scheme is the *correct* signal that spatial
CV is over-pessimistic (the Wadoux case). SLI is therefore always reported **against
the declared target**, never against NNDM-as-truth. This is what couples C1 to C2/C3
and removes the circularity.

## 5. Constructing the target ECDF `Ĝpred` from the declaration

- **`grid` / wall-to-wall:** `P` = prediction grid cells over the mapping region;
  `fpred(p) = d(p, nearest sample)`. (Distribution is insensitive to grid
  resolution above the sample spacing; sensitive — correctly — to region extent and
  masking.)
- **`newdata`:** `P` = user-supplied locations.
- **`interpolation` (within-sample):** `P` = unsampled locations inside the sampled
  domain; approximate by a fine grid masked to the sampling support (density/hull),
  or by the sample-to-sample LOO NN distances.
- **`design-based` (probability sample, population-mean estimand):** `P` drawn from
  the region weighted by the inclusion density `π(·)`; this is the case where SLI is
  read relative to the design-based correct validation (see §4).

The target is an **explicit, declared input**, so SLI is a property of (split, data,
declaration) — never of the split alone. (This is the answer to "garbage-in":
the assumption is surfaced, not hidden.)

## 6. Anisotropy, non-stationarity, and the distance metric

- **Distance:** geodesic for geographic CRS, projected otherwise; warn on missing
  CRS. Duplicated coordinates ⇒ `gobs = 0 ⇒ ρ = 1` ⇒ flagged as extreme leakage and
  handed to the grouped/duplicate channel (C4).
- **Anisotropy:** fit a directional variogram; use the anisotropic (Mahalanobis-type)
  distance inside `gobs`/`fpred`, or directly a directional `ρ(h, θ)`.
- **Non-stationarity:** allow a locally-varying `ρ` (moving-window or covariate-
  driven range); `ℓ(t)` then uses the local correlation at `t`.

## 7. Uncertainty (objection C2-④ applied to C1's inputs)

`φ` and `ρ` are *estimated*, often from clustered/preferential data where the range
is poorly identified. Therefore:
- Default = **plug-in** SLI with a fitted variogram.
- Optional = **Monte-Carlo propagation:** sample variogram parameters from their
  (profile-likelihood or Bayesian) uncertainty, recompute SLI for each draw → report
  `SLI` with an interval. Coverage of these intervals is checked in the simulation
  study.

## 8. What gets returned (API contract for P1)

`detect_leakage()` returns an S3 object exposing:
`SLI_rho` (point + interval), `SLI_d`, `A+`, `A−`, `W`, `δ` (directionality),
`c_obs`, `c_pred`, the two ECDFs (for `plot`), `ℓ(t)` per point and per fold (for the
map), the declared target, and the fitted `φ`/`ρ`. `print`/`summary`/`plot`/
`autoplot` methods follow.

## 9. Validation plan (paper)

1. **Monotonicity:** `SLI_ρ` increases with true optimism across the simulation grid
   (§ METHOD-EMULATOR). This is C1's validation (objection C1-⑤).
2. **Reduction:** `|A| = W` (NNDM) under non-crossing — numerical check.
3. **`d` vs `ρ` agreement** under stationarity; divergence under anisotropy.
4. **Invariance:** `SLI` stable under rescaling/reprojection once range-normalised.
