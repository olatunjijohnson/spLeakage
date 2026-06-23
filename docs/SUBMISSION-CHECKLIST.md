# CRAN submission checklist

Status after this session: **`R CMD check --as-cran` passes with 0 errors,
0 warnings, 2 NOTEs** — both benign (see below). Version bumped to **0.1.0**.

## The 2 NOTEs (and what they mean)

1. **CRAN incoming feasibility** — "New submission" (expected for any first
   submission) plus two **404 GitHub URLs**. These resolve once the public repo
   exists at `github.com/olatunjijohnson/spLeakage`. *Action: create the repo (or
   edit `URL`/`BugReports` in DESCRIPTION to the real location).*
2. **HTML manual / "no command 'tidy' found"** — a *local environment* note only
   (HTML Tidy isn't installed in this sandbox). CRAN's build machines have it, so it
   will not appear there. No action needed.

## Remaining steps (maintainer-owned; not done here because they are outward-facing)

- [ ] Create the public GitHub repo so the DESCRIPTION URLs are valid.
- [ ] Run on the winbuilder / macbuilder / R-hub multi-platform services
      (`devtools::check_win_devel()`, `rhub::rhub_check()`).
- [ ] Confirm `cran-comments.md` is accurate, then submit with
      `devtools::submit_cran()` (or the web form). **This is an irreversible
      outward-facing action and is intentionally left to you.**

## Build / check commands used

```r
devtools::document("spLeakage")
tar <- devtools::build("spLeakage", path = "/tmp/spbuild", vignettes = TRUE)
# then, in a shell:
# R CMD check --as-cran /tmp/spbuild/spLeakage_0.1.0.tar.gz
```

## Paper-scale extensions still open (not blockers for CRAN)

Larger simulation grid, covariate-trend + spatial-GLMM realism, normalised conformal,
additional real datasets, and the manuscript — tracked in `docs/VISION.md` and
`docs/PAPER.md`.
