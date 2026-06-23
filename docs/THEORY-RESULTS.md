# Theory — numerical validation results

From `data-raw/validate_theory.R` (192 GRF configs: range × signal × n × {random,
clustered}, exponential covariance, V = 1; explained variance computed exactly from
the covariance algebra — no field realisation). Validates `docs/THEORY.md`.

## Headline results

| Claim | Result | Reading |
|---|---|---|
| **Eq. 6 — single-NN form predicts exact GP optimism** | `cor = 0.975` | The closed form `V w² (E[ρ²(g)] − E[ρ²(f)])` tracks exact optimism almost perfectly in *shape* |
| **Multi-NN amplification** | `median(exact / nn_form) = 1.56` | Single-NN *underestimates magnitude* ~1.6× (it ignores extra neighbours), as predicted; one calibration constant recovers it |
| **SLI_rho is a near-sufficient statistic for optimism** | `cor(SLI_rho, exact) = 0.976` | The package's index explains **~95%** of the variance in true GP optimism |
| **Unsquared beats squared (in the realistic regime)** | `cor(SLI_2, exact) = 0.841 < 0.976` | The single-NN limit suggests squaring, but with dense training the *unsquared* `SLI_rho` is better — vindicating the package's original choice |
| **Eq. 8 — Wasserstein bound (single-NN optimism)** | holds in **100%** | Confirms the Kantorovich–Rubinstein step is exact for the quantity it bounds |
| **Eq. 8 — bound vs *exact* (multi-NN) optimism** | holds in **79%**, tightness 0.34 | Multi-NN amplifies optimism beyond the single-NN Lipschitz constant; inflate `||L'||` by ~1.6 to bound the exact |

## What this establishes (for the paper)

1. **The SLI has a formal estimand** — excess explained variance, `mean_t EV(t|Tr) −
   mean_p EV(p|S)` (eq. 2) — not a heuristic. This answers the top-stats reviewer.
2. **The SLI is a near-sufficient statistic for GP optimism** (r ≈ 0.98). This is the
   headline theoretical+empirical result: a cheap, model-free geometric index
   captures essentially all of the true optimism in the Gaussian regime.
3. **A rigorous Wasserstein bound** (eq. 8) holds exactly for the single-NN
   component, with a quantified multi-NN inflation for the exact quantity.
4. **A theory-derived, GRF-free optimism estimate**: plug `ρ`, `w`, and the empirical
   NN-distance ECDFs into eq. (6) × ~1.6 — no simulation-trained emulator needed.
   This is the principled alternative that sidesteps the emulator's out-of-
   distribution fragility (`task #21`).
5. **Honest scope:** the single-NN form is biased low in magnitude (motivating the
   multi-neighbour refinement, `task #22`); the theory is the Gaussian backbone, with
   non-Gaussian responses and sub-optimal learners handled empirically.

## Refinement investigated: multi-neighbour SLI (task #22)

The single-NN form underestimates exact optimism in *magnitude* (~1.3-1.6x), which
motivated a `k`-neighbour version of `SLI_rho` (noisy-OR retained correlation, exposed
as `detect_leakage(k=)`). Honest result across the GRF configs:

| | cor with exact optimism | magnitude ratio exact/SLI |
|---|---|---|
| `k = 1` (single-NN, default) | **0.977** | 1.28 (underestimates) |
| `k = 5` (noisy-OR) | 0.699 | 0.85 (overestimates) |

`k > 1` improves the magnitude but **reduces the correlation** — the noisy-OR
over-counts mutually correlated neighbours and saturates. Since correlation (ranking /
predicting optimism) is what matters, **`k = 1` is kept as the default and is
near-sufficient**; the magnitude gap is better closed by a single calibration constant
(or the exact multi-NN `c0' Sigma^-1 c0`) than by the cheap noisy-OR. A clean negative
result: the original single-NN choice is vindicated again.

## Connection to distribution-shift theory

Optimism = `∫ L d(P_dep − P_cv)` with pointwise loss `L(d) = V[1 − w² ρ(d)²]` and the
shift carried by the NN-distance laws. This is a textbook distribution-shift
generalization gap; eq. (8) is its Kantorovich–Rubinstein bound. The route to
model-agnostic optimism bounds (beyond GP) runs through this framing — importance
weighting / H-divergence / MMD — and broadens the audience from spatial statistics to
ML theory.
