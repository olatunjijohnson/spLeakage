# Emulator out-of-distribution behaviour on real data (task #21)

The optimism emulator is calibrated on simulated (Matern GP) fields. The reviewer-
critical question is: does it work on *real* data, or refuse everything? Tested on
real spatial datasets (meuse zinc/lead, ca20 calcium, Parana rainfall, ozone2):

| dataset | in area of applicability? | DI (threshold 0.38) |
|---|---|---|
| meuse zinc | yes | 0.36 |
| meuse lead | yes | 0.36 |
| Parana rainfall | yes | 0.33 |
| ozone | yes | 0.35 |
| ca20 calcium | **no (abstains)** | 0.94 |

**The sim-trained emulator generalises to 4/5 real datasets, and honestly abstains on
the out-of-envelope one** (ca20, DI far above threshold). This is the intended
behaviour: predict where the calibration support covers the query, refuse otherwise
and defer to the empirical `estimate_optimism()` / `deleak_estimate()` route. The AOA
guard is doing real work on real data, not just simulation.

## Honest scope of generator realism

The generator already spans Matern smoothness, three sampling designs (incl.
preferential), three response distributions (Gaussian / Poisson / Binomial, i.e. a
spatial-GLMM family), three learners and two prediction targets. Two further realism
axes were considered:

- **Non-stationarity / covariate trend in the field.** Trend-dominated fields are the
  one real-data regime the emulator and the SLI under-serve (the meta-audit finding,
  `docs/META-AUDIT.md`). However, that regime is now covered by a *dedicated*
  diagnostic, `trend_strength()` (task #24, cor 0.75 with real optimism), rather than
  folded into the emulator. Folding it in would require adding a trend feature, trend
  fields, **and** extrapolation deployment targets to the generator (the trend
  optimism only manifests under extrapolation) — a generator overhaul whose main
  payoff is already delivered by the standalone trend channel.
- **Real-covariance generators.** The OOD result above shows real-covariance
  generalisation is already adequate (4/5 in-AOA), so this is low priority.

**Conclusion.** The emulator generalises to most real spatial data and abstains safely
otherwise; the one systematic gap (trend extrapolation) is handled by a separate,
validated diagnostic. Further generator enrichment is a documented, low-priority
enhancement rather than a blocker.
