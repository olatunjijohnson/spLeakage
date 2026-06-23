# Real-data validation (task #20). REAL spatial fields (meuse zinc, Parana rainfall,
# ca20 calcium). For each replicate we draw a CLUSTERED training set (as real surveys
# cluster) and hold out the rest as a genuine INDEPENDENT test set. The question:
# does spLeakage predict the real CV-vs-independent optimism gap -- and does the
# de-leaked estimate (computed WITHOUT the independent data) recover the true
# independent error? This is out-of-(simulation)-distribution: real fields, not GRFs.
# Run with: devtools::load_all(); source("data-raw/validate_independent.R")

suppressMessages({devtools::load_all(quiet = TRUE)})
set.seed(20)

get_field <- function(name) {
  if (name == "meuse") { data(meuse, package = "sp")
    return(data.frame(x = meuse$x, y = meuse$y, z = log(meuse$zinc))) }
  if (name == "parana") { e <- new.env(); data("parana", package = "geoR", envir = e)
    p <- e$parana; return(data.frame(x = p$coords[, 1], y = p$coords[, 2], z = p$data)) }
  if (name == "ca20") { e <- new.env(); data("ca20", package = "geoR", envir = e)
    p <- e$ca20; return(data.frame(x = p$coords[, 1], y = p$coords[, 2], z = p$data)) }
}
rmse <- function(o, p) sqrt(mean((o - p)^2))

# Two deployment scenarios spanning interpolation -> extrapolation:
#  - "region": hold out a spatially coherent region (far ~35% along a random
#     direction) = genuine extrapolation deployment (leakage should bite).
#  - "interp": hold out a random ~35% interspersed with training = interpolation
#     deployment (deployment ~ CV, little leakage).
make_split <- function(coords, scenario) {
  n <- nrow(coords)
  if (scenario == "region") {
    th <- runif(1, 0, 2 * pi); proj <- as.numeric(coords %*% c(cos(th), sin(th)))
    te <- which(proj > stats::quantile(proj, 0.65)); list(tr = setdiff(seq_len(n), te), te = te)
  } else {
    te <- sample(n, round(0.35 * n)); list(tr = setdiff(seq_len(n), te), te = te)
  }
}

run_dataset <- function(name, B = 40L) {
  d <- get_field(name); coords <- as.matrix(d[, c("x", "y")]); idw <- .idw_predictor(c("x", "y"), "z")
  scen <- rep(c("region", "interp"), each = B)
  out <- lapply(scen, function(scenario) {
    s <- make_split(coords, scenario); dtr <- d[s$tr, ]; dte <- d[s$te, ]; m <- nrow(dtr)
    fold <- sample(rep_len(1:10, m))
    folds <- lapply(split(seq_len(m), fold), function(i) list(test = i, train = setdiff(seq_len(m), i)))
    E_cv <- rmse(dtr$z, .cv_predict(dtr, folds, idw))           # naive reported error
    E_ind <- rmse(dte$z, idw(dtr, dte))                         # TRUE independent error
    tgt <- prediction_target(newdata = as.matrix(dte[, c("x", "y")]), type = "newdata")
    dl  <- deleak_estimate(dtr, fold, "z", coords = c("x", "y"), n_boot = 1L)          # block
    dlt <- deleak_estimate(dtr, fold, "z", coords = c("x", "y"), target = tgt, n_boot = 1L)  # matched
    lk <- detect_leakage(dtr, fold, tgt, response = "z", coords = c("x", "y"))
    data.frame(scenario = scenario, real_opt = (E_ind - E_cv) / E_ind,
               pred_opt = dl$optimism_rel, sli = lk$SLI_rho,
               E_cv = E_cv, E_ind = E_ind, E_dl = dl$deleaked, E_dlt = dlt$deleaked)
  })
  cbind(dataset = name, do.call(rbind, out))
}

res <- do.call(rbind, lapply(c("meuse", "parana", "ca20"), run_dataset))
saveRDS(res, "data-raw/independent_validation.rds")

reg <- res[res$scenario == "region", ]
cat(sprintf("replicates: %d (3 real datasets x 2 deployment scenarios)\n\n", nrow(res)))
cat("Does the SLI predict the REAL optimism gap (across interp + extrapolation)?\n")
cat(sprintf("  cor(SLI_rho, real optimism)        = %+.3f\n", cor(res$sli, res$real_opt)))
cat(sprintf("  cor(de-leaked pred, real optimism) = %+.3f\n", cor(res$pred_opt, res$real_opt)))
cat("\nDoes the SLI separate the two deployment scenarios (as real optimism does)?\n")
cat(sprintf("  mean real optimism : extrapolation %+.2f vs interpolation %+.2f\n",
            mean(reg$real_opt), mean(res$real_opt[res$scenario == "interp"])))
cat(sprintf("  mean SLI_rho       : extrapolation %+.2f vs interpolation %+.2f\n",
            mean(reg$sli), mean(res$sli[res$scenario == "interp"])))
cat("\nExtrapolation deployment -- does de-leak recover the TRUE error? (ratio to E_indep)\n")
cat(sprintf("  naive CV              E_cv  / E_indep = %.2f  (understates by %.0f%%)\n",
            median(reg$E_cv / reg$E_ind), 100 * median(1 - reg$E_cv / reg$E_ind)))
cat(sprintf("  de-leak (block)       E_dl  / E_indep = %.2f  (closes %.0f%% of the gap)\n",
            median(reg$E_dl / reg$E_ind),
            100 * median((reg$E_dl - reg$E_cv) / (reg$E_ind - reg$E_cv))))
cat(sprintf("  de-leak (target-aware) E_dlt / E_indep = %.2f  (closes %.0f%% of the gap)\n",
            median(reg$E_dlt / reg$E_ind),
            100 * median((reg$E_dlt - reg$E_cv) / (reg$E_ind - reg$E_cv))))
