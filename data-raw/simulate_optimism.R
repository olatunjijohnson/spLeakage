# Simulation study calibrating the optimism emulator (contribution C2b).
# See docs/METHOD-EMULATOR.md. Ground-truth optimism is known by construction:
#   optimism_rel = (E_true - E_cv) / E_true
# scored against NOISY target observations (same nugget floor as E_cv).
#
# Rich generators: Matern fields x smoothness x signal x {random, clustered,
# preferential} sampling x sample size x {gaussian, poisson, binomial} response,
# and three learners {idw, rf, gam}. A space-filling random design over the factor
# grid (not full factorial). Parallelised over scenarios.
#
# Run with: devtools::load_all(); source("data-raw/simulate_optimism.R")

library(spLeakage)
stopifnot(requireNamespace("fields"), requireNamespace("ranger"), requireNamespace("mgcv"))

N_CONFIG    <- 1000L         # distinct factor combinations (paper-scale run)
N_REP       <- 6L            # realizations per config, averaged to denoise the label
N_CORES     <- max(1L, parallel::detectCores() - 1L)
GRIDN       <- 26L           # candidate-location grid (GRIDN^2 points)
set.seed(2024)

rmse <- function(o, p) sqrt(mean((o - p)^2, na.rm = TRUE))

# ---- learners: predict_fun(train, test) on columns x, y, z --------------------
fit_idw <- .idw_predictor(c("x", "y"), "z")
fit_rf <- function(train, test) {
  m <- ranger::ranger(z ~ x + y, data = train, num.trees = 200, num.threads = 1)
  predict(m, test)$predictions
}
fit_gam <- function(train, test) {
  k <- min(60L, max(10L, floor(nrow(train) / 5)))
  m <- mgcv::gam(z ~ s(x, y, k = k), data = train)
  as.numeric(mgcv::predict.gam(m, newdata = test))
}
LEARNERS <- list(idw = fit_idw, rf = fit_rf, gam = fit_gam)

# ---- factor combinations (configs) ---------------------------------------------
make_config <- function(id) {
  list(id = id,
       range_rel = sample(c(0.05, 0.10, 0.20, 0.40), 1),
       nu        = sample(c(0.5, 1.5, 2.5), 1),
       signal    = runif(1, 0.3, 0.9),
       design    = sample(c("random", "clustered", "preferential"), 1),
       n_clusters = sample(c(5L, 8L, 12L), 1),
       cl_spread  = sample(c(0.03, 0.05, 0.08), 1),   # tightness -> spans realistic clustering
       n         = sample(c(150L, 300L, 600L), 1),
       rtype     = sample(c("gaussian", "poisson", "binomial"), 1),
       model     = sample(c("idw", "rf", "gam"), 1),
       ttype     = sample(c("grid", "interpolation"), 1))
}

# ---- one realization of a config -> up to 2 rows (random & block split) ---------
sim_scenario <- function(cfg, seed) {
  set.seed(seed)
  range_rel <- cfg$range_rel; nu <- cfg$nu; signal <- cfg$signal; design <- cfg$design
  n <- cfg$n; rtype <- cfg$rtype; model <- cfg$model; ttype <- cfg$ttype

  gs <- seq(0.02, 0.98, length.out = GRIDN); G <- as.matrix(expand.grid(x = gs, y = gs))
  nG <- nrow(G)
  if (ttype == "grid") {
    ts <- seq(0.04, 0.96, length.out = 13L); TG <- as.matrix(expand.grid(x = ts, y = ts))
  } else {
    TG <- cbind(x = runif(150, 0.05, 0.95), y = runif(150, 0.05, 0.95))
  }
  allp <- rbind(G, TG)
  D <- as.matrix(stats::dist(allp))
  C <- fields::Matern(D, range = range_rel, smoothness = nu)
  S <- tryCatch(as.numeric(t(chol(C + diag(1e-7, nrow(C)))) %*% rnorm(nrow(C))),
                error = function(e) NULL)
  if (is.null(S)) return(NULL)
  Sg <- S[seq_len(nG)]; St <- S[(nG + 1L):nrow(allp)]

  # sample n locations from the candidate grid by design
  idx <- switch(design,
    random = sample(nG, n),
    clustered = {
      ce <- G[sample(nG, cfg$n_clusters), , drop = FALSE]
      dmin <- apply(G, 1L, function(p) min(colSums((t(ce) - p)^2)))
      sample(nG, n, prob = exp(-dmin / (2 * cfg$cl_spread^2)) + 1e-6)
    },
    preferential = sample(nG, n, prob = exp(1.5 * as.numeric(scale(Sg)))))
  # Jitter off the candidate grid with continuous noise, otherwise the grid spacing
  # floors the nearest-neighbour distance and clustered samples look regular
  # (Clark-Evans index never drops below ~1).
  spacing <- (0.96 / (GRIDN - 1))
  sxy <- G[idx, , drop = FALSE] + matrix(rnorm(2 * n, 0, 0.4 * spacing), n, 2)
  sxy[] <- pmin(pmax(sxy, 0), 1)
  Ss <- Sg[idx]

  # responses (sample + noisy target observations), scored on the response scale
  nug <- sqrt((1 - signal) / signal)
  mk_resp <- function(s, m) switch(rtype,
    gaussian = s + rnorm(m, 0, nug),
    poisson  = rpois(m, exp(1 + 0.6 * s)),
    binomial = rbinom(m, 1, plogis(0.6 * s)))
  z  <- mk_resp(Ss, n)
  zt <- mk_resp(St, length(St))
  df <- data.frame(x = sxy[, 1], y = sxy[, 2], z = z)
  if (stats::sd(z) < 1e-8) return(NULL)

  pf <- LEARNERS[[model]]
  dep <- tryCatch(estimate_dependence(df, "z", c("x", "y")), error = function(e) NULL)
  if (is.null(dep)) return(NULL)
  tgt <- prediction_target(grid = TG, type = "grid")

  # E_true: fit on full sample, predict at target, score vs noisy target obs
  E_true <- tryCatch(rmse(zt, pf(df, data.frame(x = TG[, 1], y = TG[, 2]))),
                     error = function(e) NA_real_)
  if (!is.finite(E_true) || E_true <= 0) return(NULL)

  out <- list()
  for (sp in c("random", "block")) {
    ids <- if (sp == "random") sample(rep_len(1:10, n)) else spatial_block_cv(df, 10, c("x", "y"))
    folds <- lapply(split(seq_len(n), ids),
                    function(te) list(test = te, train = setdiff(seq_len(n), te)))
    E_cv <- tryCatch(rmse(z, .cv_predict(df, folds, pf)), error = function(e) NA_real_)
    feats <- tryCatch(.emulator_features(df, ids, tgt, "z", c("x", "y"), dependence = dep),
                      error = function(e) NULL)
    if (is.null(feats) || !is.finite(E_cv)) next
    out[[length(out) + 1L]] <- data.frame(
      as.list(feats), model = model, response_type = rtype, split = sp,
      optimism_rel = (E_true - E_cv) / E_true,
      config = cfg$id, design = design, nu = nu, target = ttype,
      stringsAsFactors = FALSE)
  }
  if (length(out)) do.call(rbind, out) else NULL
}

# ---- run (parallel): N_CONFIG configs x N_REP realizations ----------------------
cat(sprintf("running %d configs x %d reps on %d cores...\n", N_CONFIG, N_REP, N_CORES))
configs <- lapply(seq_len(N_CONFIG), make_config)
jobs <- expand.grid(ci = seq_len(N_CONFIG), rep = seq_len(N_REP))
res <- parallel::mclapply(seq_len(nrow(jobs)), function(j)
  sim_scenario(configs[[jobs$ci[j]]], seed = 1000L * jobs$ci[j] + jobs$rep[j]),
  mc.cores = N_CORES)
raw <- do.call(rbind, res[!vapply(res, is.null, logical(1))])

# Average features and label across realizations within (config, split): denoises
# the single-realization Monte-Carlo error in optimism_rel.
key <- interaction(raw$config, raw$split, drop = TRUE)
num <- c(.EMU_FEATURES, "optimism_rel")
agg <- lapply(split(seq_len(nrow(raw)), key), function(ix) {
  r0 <- raw[ix[1], ]
  data.frame(as.list(colMeans(raw[ix, num, drop = FALSE])),
             model = r0$model, response_type = r0$response_type,
             config = r0$config, design = r0$design, nu = r0$nu,
             target = r0$target, n_rep = length(ix), stringsAsFactors = FALSE)
})
tab <- do.call(rbind, agg)
cat(sprintf("denoised rows: %d (from %d raw realizations)\n", nrow(tab), nrow(raw)))
cat("cor(sli_rho, optimism):", round(cor(tab$sli_rho, tab$optimism_rel), 3), "\n")
print(round(tapply(tab$optimism_rel, tab$design, mean), 3))

# ---- held-out validation (by config, within-simulation) ------------------------
set.seed(1)
scen <- unique(tab$config); ho <- sample(scen, floor(0.2 * length(scen)))
train_tab <- tab[!tab$config %in% ho, ]; test_tab <- tab[tab$config %in% ho, ]
em_val <- .fit_optimism_emulator(train_tab)
bv <- .emu_booster(em_val)
pv <- t(vapply(seq_len(nrow(test_tab)), function(i) {
  r <- test_tab[i, ]
  p <- .emu_predict_one(em_val, unlist(r[.EMU_FEATURES]), r$model, r$response_type, booster = bv)
  c(pred = p$optimism_rel, lo = p$optimism_ci[1], hi = p$optimism_ci[2], in_aoa = p$in_aoa)
}, numeric(4)))
val <- data.frame(truth = test_tab$optimism_rel, pv)
ina <- val$in_aoa == 1 & is.finite(val$pred)
validation <- list(
  n_test = nrow(val), aoa_coverage = mean(ina),
  R2 = 1 - sum((val$truth[ina] - val$pred[ina])^2) / sum((val$truth[ina] - mean(val$truth[ina]))^2),
  rmse = rmse(val$truth[ina], val$pred[ina]),
  interval_coverage = mean(val$truth[ina] >= val$lo[ina] & val$truth[ina] <= val$hi[ina]))
cat(sprintf("VALIDATION: in-AOA %.0f%%  R2=%.2f  RMSE=%.3f  90%%-cover=%.0f%%\n",
            100 * validation$aoa_coverage, validation$R2, validation$rmse,
            100 * validation$interval_coverage))

# ---- fit final emulator on all rows and ship -----------------------------------
optimism_emulator <- .fit_optimism_emulator(tab)
attr(optimism_emulator, "validation") <- validation
print(optimism_emulator)

optimism_sim_table <- tab
save(optimism_emulator, file = "R/sysdata.rda", compress = "xz")
saveRDS(optimism_sim_table, file = "data-raw/optimism_sim_table.rds")
cat("saved R/sysdata.rda and data-raw/optimism_sim_table.rds\n")
