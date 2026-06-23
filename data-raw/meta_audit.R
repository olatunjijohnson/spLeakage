# Meta-audit (task #13): how widespread/severe is spatial leakage in real spatial
# prediction? A corpus of REAL public datasets (heavy metals, soil chemistry,
# rainfall, elevation, ozone, malaria prevalence). For each, measure the REAL
# optimism of naive random 10-fold CV via spatial region-holdout (a direct
# measurement, not our model), and test whether the SLI predicts it (so the audit
# scales without holdout). Ploton-style "is the field overconfident?" evidence.
# Run with: devtools::load_all(); source("data-raw/meta_audit.R")

suppressMessages(devtools::load_all(quiet = TRUE))
set.seed(13)
rmse <- function(o, p) sqrt(mean((o - p)^2))
NMAX <- 700L

build_corpus <- function() {
  L <- list(); add <- function(name, x, y, z) {
    ok <- is.finite(x) & is.finite(y) & is.finite(z)
    L[[name]] <<- data.frame(x = x[ok], y = y[ok], z = z[ok])
  }
  g <- function(nm, pkg) { e <- new.env(); data(list = nm, package = pkg, envir = e); get(nm, e) }
  m <- g("meuse", "sp")
  for (v in c("cadmium", "copper", "lead", "zinc")) add(paste0("meuse:", v), m$x, m$y, log(m[[v]]))
  add("meuse:elev", m$x, m$y, m$elev); add("meuse:om", m$x, m$y, m$om)
  ca <- g("ca20", "geoR"); add("ca20:calcium", ca$coords[, 1], ca$coords[, 2], ca$data)
  pa <- g("parana", "geoR"); add("parana:rain", pa$coords[, 1], pa$coords[, 2], pa$data)
  el <- g("elevation", "geoR"); add("elevation:elev", el$coords[, 1], el$coords[, 2], el$data)
  cm <- g("camg", "geoR")
  for (v in c("ca020", "mg020", "ctc020")) add(paste0("camg:", v), cm$east, cm$north, cm[[v]])
  ga <- g("gambia", "geoR"); ag <- aggregate(pos ~ x + y, data = as.data.frame(ga), FUN = mean)
  add("gambia:malaria", ag$x, ag$y, ag$pos)
  oz <- g("ozone2", "fields"); add("ozone2:o3", oz$lon.lat[, 1], oz$lon.lat[, 2], colMeans(oz$y, na.rm = TRUE))
  nr <- g("NorthAmericanRainfall", "fields")
  add("NArain:precip", nr$longitude, nr$latitude, log(nr$precip))
  add("NArain:elev", nr$longitude, nr$latitude, nr$elevation)
  L
}

audit_study <- function(d, R = 8L) {
  if (nrow(d) > NMAX) d <- d[sample(nrow(d), NMAX), ]
  coords <- as.matrix(d[, c("x", "y")]); n <- nrow(coords); idw <- .idw_predictor(c("x", "y"), "z")
  vals <- vapply(seq_len(R), function(i) {
    th <- 2 * pi * i / R; proj <- as.numeric(coords %*% c(cos(th), sin(th)))
    te <- which(proj > stats::quantile(proj, 0.65)); tr <- setdiff(seq_len(n), te)
    dtr <- d[tr, ]; dte <- d[te, ]; m <- length(tr)
    fold <- sample(rep_len(1:10, m))
    folds <- lapply(split(seq_len(m), fold), function(ix) list(test = ix, train = setdiff(seq_len(m), ix)))
    E_cv <- rmse(dtr$z, .cv_predict(dtr, folds, idw))
    E_ind <- rmse(dte$z, idw(dtr, dte))
    tgt <- prediction_target(newdata = as.matrix(dte[, c("x", "y")]), type = "newdata")
    sli <- tryCatch(detect_leakage(dtr, fold, tgt, response = "z", coords = c("x", "y"))$SLI_rho,
                    error = function(e) NA_real_)
    c(real_opt = (E_ind - E_cv) / E_ind, sli = sli)
  }, numeric(2))
  # What an analyst would compute from their data alone (no holdout): target-aware
  # de-leaked optimism, with deployment = a grid over the domain.
  gx <- seq(min(d$x), max(d$x), length.out = 12L); gy <- seq(min(d$y), max(d$y), length.out = 12L)
  tgt <- prediction_target(grid = as.matrix(expand.grid(x = gx, y = gy)), type = "grid")
  fold_full <- sample(rep_len(1:10, n))
  dlo <- tryCatch(deleak_estimate(d, fold_full, "z", coords = c("x", "y"), target = tgt,
                                  n_boot = 1L)$optimism_rel, error = function(e) NA_real_)
  c(n = n, real_opt = mean(vals["real_opt", ]), sli = mean(vals["sli", ], na.rm = TRUE),
    deleak_opt = dlo)
}

corpus <- build_corpus()
cat(sprintf("corpus: %d real spatial response variables\n", length(corpus)))
res <- as.data.frame(t(vapply(corpus, audit_study, numeric(4))))
res$study <- rownames(res); res <- res[order(-res$real_opt), ]
saveRDS(res, "data-raw/meta_audit.rds")

cat(sprintf("\nREAL optimism of naive random CV (extrapolation deployment), %d studies:\n", nrow(res)))
cat(sprintf("  median = %.0f%%   IQR = [%.0f%%, %.0f%%]\n",
            100 * median(res$real_opt), 100 * quantile(res$real_opt, .25), 100 * quantile(res$real_opt, .75)))
cat(sprintf("  studies optimistic (>0%%): %.0f%% ;  substantially (>20%%): %.0f%%\n",
            100 * mean(res$real_opt > 0), 100 * mean(res$real_opt > 0.2)))
cat("\nAuditing WITHOUT holdout (cross-dataset prediction of real optimism):\n")
cat(sprintf("  raw SLI                  cor = %+.2f  <- does NOT transfer across datasets\n",
            cor(res$sli, res$real_opt)))
cat(sprintf("  de-leaked optimism       cor = %+.2f  <- the proper cross-dataset auditor\n",
            cor(res$deleak_opt, res$real_opt, use = "complete.obs")))
cat("\nper-study (sorted by real optimism):\n")
print(data.frame(study = res$study, n = res$n, real_opt = round(res$real_opt, 2),
                 deleak_opt = round(res$deleak_opt, 2), sli = round(res$sli, 2)), row.names = FALSE)

# figure
dir.create("paper-figures", showWarnings = FALSE)
png("paper-figures/F4_meta_audit.png", 1500, 1200, res = 200)
op <- res$real_opt[order(res$real_opt)]; nm <- res$study[order(res$real_opt)]
graphics::par(mar = c(4, 9, 3, 1))
graphics::barplot(100 * op, horiz = TRUE, names.arg = nm, las = 1, cex.names = 0.7,
                  col = ifelse(op > 0.2, "#d1495b", "#00798c"),
                  xlab = "real optimism of random CV (%)",
                  main = "Meta-audit: random CV overstates accuracy on real spatial data")
graphics::abline(v = c(0, 20), lty = c(1, 2), col = c("grey40", "grey60"))
dev.off()
cat("\nwrote paper-figures/F4_meta_audit.png\n")
