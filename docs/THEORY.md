# Theory — optimism as a distribution-shift gap, and why the SLI predicts it

Status: **derivation v1** (extension ① in `docs/PAPER.md`). Goal: give the Spatial
Leakage Index a formal estimand, derive expected optimism under a Gaussian process,
show the SLI is its leading term, and bound optimism by a Wasserstein distance
between two nearest-neighbour-distance laws. This converts the package from a
validated heuristic into a principled method and connects it to distribution-shift
generalization theory.

## 1. Model

Latent zero-mean Gaussian process `Z` with stationary covariance
`Cov(Z(s), Z(s')) = sigma^2 rho(||s - s'||)`, observed with iid nugget noise:

```
Y(s) = mu + Z(s) + eps(s),   eps ~ N(0, tau^2),   independent of Z.
```

Total variance `V = sigma^2 + tau^2`; **signal proportion** `w = sigma^2 / V` (this
is exactly the emulator's `signal` feature). `rho(0)=1`, `rho` decreasing to 0.

## 2. Kriging MSE at a location

Predict the (noisy) observable `Y(s0)` from a training set `T` of `m` observations.
Simple kriging of the latent gives `Zhat(s0) = c0' Sigma^{-1} (y - mu)` with
`c0 = sigma^2 rho(s0, T)` and `Sigma = sigma^2 R + tau^2 I` (R = correlation matrix
of `T`). The mean-squared error of predicting `Y(s0) = Z(s0) + eps0` is

```
MSE(s0 | T) = V - c0' Sigma^{-1} c0 .                                        (1)
```

Call `EV(s0 | T) = c0' Sigma^{-1} c0 >= 0` the **explained variance**: the reduction
in predictive variance contributed by the training data. It is large when `s0` is
strongly correlated with (near) the training points, small when far.

## 3. Optimism is a difference in explained variance

A cross-validation scheme reports, per test point `t`, the error of predicting it
from its fold's training set `Tr(t)`; deployment predicts at target locations `p`
from the full sample `S`. Averaging (1), the `V` cancels and

```
E_cv   = V - mean_t  EV(t | Tr(t))
E_true = V - mean_p  EV(p | S)
optimism_abs = E_true - E_cv = mean_t EV(t | Tr(t)) - mean_p EV(p | S).        (2)
```

**Optimism is exactly the excess explained variance the CV scheme enjoys over
deployment.** This is the formal estimand: an interpretable population quantity, not
a heuristic. (2) holds for *any* GP — no approximation yet.

## 4. The nearest-neighbour closed form (where the SLI appears)

`EV` depends on the whole local training configuration. Take the **single-nearest-
neighbour approximation**: keep only the nearest training point, at distance `d`.
Then `c0 = sigma^2 rho(d)`, `Sigma = V`, and

```
EV(s0) ≈ sigma^4 rho(d)^2 / V = V w^2 rho(d)^2 .                              (3)
```

So the pointwise loss is a function of the nearest-training distance only:

```
L(d) = MSE ≈ V [ 1 - w^2 rho(d)^2 ] .                                         (4)
```

Substituting (3) into (2), with `g = gobs(t)` (test→train NN distance) and
`f = fpred(p)` (deployment→sample NN distance):

```
optimism_abs ≈ V w^2 ( E_t[ rho(g)^2 ] - E_p[ rho(f)^2 ] ).                   (5)
```

### 4.1 This *is* the SLI (and the numerics vindicate the unsquared form)
The package's dependence-form index is `SLI_rho = w (E_t[rho_Z(g)] - E_p[rho_Z(f)])`.
The single-NN derivation suggests the **squared** correlation should appear:

```
SLI_2 := E_t[ rho_Z(g)^2 ] - E_p[ rho_Z(f)^2 ],   optimism_abs ≈ (V w^2) * SLI_2.  (6)
```

**But the numerical validation (§7) overturns this in the realistic regime.** Against
the *exact* (multi-neighbour) GP optimism (2), the package's unsquared `SLI_rho`
correlates **r = 0.976** while the single-NN `SLI_2` correlates only **0.841**. The
squaring is an artefact of the single-nearest-neighbour limit; with realistic dense
training (many correlated neighbours), explained variance saturates and the *linear*
(unsquared) correlation difference is the better sufficient statistic. So the
package's original `SLI_rho` is *vindicated by the exact computation*: it explains
~95% of the variance in true GP optimism. We keep `SLI_rho` as the headline feature;
`SLI_2` is reported as the single-NN-limit diagnostic.

### 4.2 Why the emulator's feature set is exactly right
Equation (6) says optimism is a function of: the correlation function `rho` (hence
**range** and smoothness), the **signal proportion** `w`, and the **NN-distance
distributions** of the split and target (summarised by the **SLI** and the
clustering geometry). These are precisely the emulator's inputs. The emulator was
empirically rediscovering (6); the theory explains *why* those features and predicts
the `V w^2` scaling — a strong, falsifiable narrative.

## 5. Optimism as a distribution-shift gap, and a Wasserstein bound

Write `P_cv` for the law of the test→train NN distance and `P_dep` for the law of the
deployment→sample NN distance. By (4),

```
optimism_abs ≈ E_{d~P_dep}[L(d)] - E_{d~P_cv}[L(d)] = ∫ L d(P_dep - P_cv).      (7)
```

This is a **distribution-shift generalization gap**: the same pointwise loss `L`
evaluated under two different distributions of the governing variable `d`. Because
`L` is bounded-Lipschitz (it is smooth and monotone increasing in `d`, with
`||L'||_inf = V w^2 sup_d |d/dd rho(d)^2|`), Kantorovich–Rubinstein gives

```
|optimism_abs| <= ||L'||_inf * W_1(P_cv, P_dep),                              (8)
```

the **Wasserstein-1** distance between the two NN-distance laws. The package's
distance-form index `SLI_d` is (up to the range normalisation) exactly the signed
`W_1` between these ECDFs, so **(8) bounds optimism by an SLI_d-type quantity** and
the directionality `delta = A/W` tells you whether the bound is tight (no crossing)
or loose (crossing). This is the theoretical bridge: `SLI_rho` gives the
*value* (6); `SLI_d` gives a rigorous *bound* (8).

**Validation caveat (honest):** (8) is derived from the single-NN loss `L`, so it
rigorously bounds the *single-NN* optimism — verified to hold in **100%** of configs
(§7), confirming the Kantorovich–Rubinstein step. The *exact* multi-neighbour
optimism is ~1.5× larger (multi-NN amplification), so the bound as stated holds for
the exact optimism in ~79% of configs; bounding the exact optimism requires inflating
`||L'||` by the empirical multi-NN factor (~1.5–1.6).

This framing also imports the broader distribution-shift / domain-adaptation toolkit
(importance weighting, H-divergence, MMD bounds) — a route to model-agnostic
optimism bounds beyond the GP, and to a much wider (ML-theory) readership.

## 6. Scope, assumptions, and where it breaks (be honest)

- **Single-NN approximation (3):** exact when training is sparse relative to the
  range; for dense training it *underestimates* `EV` (ignores extra neighbours), so
  (6) is a lower bound on explained variance and the optimism approximation degrades.
  The multi-neighbour correction (`task #22`) replaces `rho(d)^2` by the full
  `c0' Sigma^{-1} c0`; the exact (2) always holds.
- **Stationarity & Gaussianity:** (1)–(2) need only second-order structure; the GP
  is used for the MSE form. Non-Gaussian responses (counts/prevalence) enter through
  the link and are handled empirically by the emulator — the theory is the
  Gaussian backbone, not the whole story.
- **Model = optimal kriging:** real models (RF, GAM) are sub-optimal; their EV is
  `<=` kriging's, so (2) is an *upper bound* on their realisable optimism. The
  emulator's model-class factor absorbs this; the theory gives the envelope.
- **Estimand caveat:** "optimism" is defined relative to the declared deployment
  target `P_dep` (the §3 thesis). For a probability sample with a population-mean
  estimand, `P_dep` is the inclusion-weighted population law and the sign can flip
  (the Wadoux case) — consistent with (7).

## 7. Numerical validation (this repo)

`data-raw/validate_theory.R` checks, across GRF configs with **known** `sigma^2,
tau^2, rho`:
1. **Exact identity (2):** measured kriging optimism == `mean_t EV - mean_p EV`
   computed from the covariance algebra (should match to numerical error).
2. **Single-NN approximation (6):** `(V w^2) * SLI_2` vs the exact optimism (2) —
   how good is the closed form, and across what range/density regimes.
3. **`SLI_2` vs `SLI_rho`:** both correlate with optimism; `SLI_2` should track the
   exact optimism more tightly (the theory's prediction).
4. **Wasserstein bound (8):** verify `|optimism| <= ||L'||_inf * W_1` holds and
   report tightness vs `delta`.

See `docs/THEORY-RESULTS.md` for the numbers once run.

## 8. What this buys the paper

- The SLI gets a **formal estimand** (excess explained variance, eq. 2) and a
  **closed-form link to optimism** (eq. 6) — answering the top-stats reviewer's
  "what does it estimate?".
- A **rigorous bound** (eq. 8) connecting the distance-form index to optimism.
- A principled, **GRF-free** optimism estimate (plug `rho`, `w`, and the empirical
  NN-distance ECDFs into (6)) that does not depend on the simulation-trained
  emulator — mitigating the out-of-distribution fragility (`task #21`).
- A bridge to **distribution-shift generalization theory**, broadening the audience
  from spatial statistics to ML.
