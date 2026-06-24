# Ecological case study: species distribution model (presence/absence)
# ---------------------------------------------------------------------------
# Audits the spatial validation of an SDM built on the Valavi et al. (2019, MEE)
# south-east Australia presence/absence dataset shipped with blockCV: 500 records
# (243 presence / 257 absence) and four bioclim covariates. Demonstrates the
# spLeakage workflow on a canonical ecological mapping problem and contrasts the
# leakage signature with the malaria case study (here: genuine autocorrelation +
# covariate-space extrapolation; there: co-located repeat surveys).
#
# Reproducible; requires Suggests: blockCV, terra.
suppressMessages({library(spLeakage); library(terra)})
set.seed(1)

## ---- data ----------------------------------------------------------------
sp  <- read.csv(system.file("extdata/species.csv", package = "blockCV"))
ed  <- system.file("extdata", package = "blockCV")
r   <- terra::rast(list.files(file.path(ed, "au"), full.names = TRUE))  # bio_4,5,12,15
covn <- names(r)

cov <- terra::extract(r, sp[, c("x", "y")])[, -1]
d   <- cbind(sp, cov)
d   <- d[stats::complete.cases(d), ]
message(sprintf("n = %d points (%d presence); covariates: %s",
                nrow(d), sum(d$occ), paste(covn, collapse = ", ")))

## ---- deployment target: wall-to-wall map over the covariate extent --------
grid_all <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
grid_xy  <- grid_all[sample(nrow(grid_all), min(4000, nrow(grid_all))), c("x", "y")]
tgt      <- prediction_target(grid = as.matrix(grid_xy), type = "grid")

## ---- the naive analysis an ecologist might run: random 10-fold ------------
folds <- sample(rep_len(1:10, nrow(d)))

res <- list()

## 1. Geographic leakage of the random split
res$leak <- detect_leakage(d, split = folds, target = tgt, response = "occ",
                           coords = c("x", "y"), n_boot = 200)

## 2. Optimism in the Brier score (proper scoring rule for presence/absence)
res$opt <- estimate_optimism(d, folds, response = "occ", coords = c("x", "y"),
                             metric = "brier", control = "block")

## 3. Feature-space (covariate) leakage: does the wall-to-wall map extrapolate
##    the bioclim covariates beyond where the model was validated?
res$feat <- detect_feature_leakage(
  d, folds, covariates = covn,
  newdata = grid_all[sample(nrow(grid_all), 3000), covn])

## 4. Design-aware recommendation and data-driven ranking
res$rec  <- recommend_validation(d, estimand = "prediction", design = "clustered",
                                 target = "grid", coords = c("x", "y"))
res$rank <- rank_cv_schemes(d, tgt, response = "occ", coords = c("x", "y"))

## 5. One-screen multi-channel audit
res$audit <- audit_workflow(d, split = folds, target = tgt, response = "occ",
                            coords = c("x", "y"))

print(res$leak); cat("\n"); print(res$opt); cat("\n"); print(res$feat); cat("\n")
print(res$rec);  cat("\n"); print(res$rank); cat("\n"); print(res$audit)

invisible(res)
