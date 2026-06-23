# Numerical validation of design_validation (task #6). Does placing validation points
# by the optimal design estimate true map accuracy more precisely than random
# placement, for the same budget? Truth = RMSE of the model over a dense evaluation
# grid (known). We estimate it from `budget` validation points chosen randomly vs by
# the design, and compare bias and variance of the estimator over replicates.
# Run with: devtools::load_all(); source("data-raw/validate_design.R")

suppressMessages(devtools::load_all(quiet = TRUE))
set.seed(6)
idw <- .idw_predictor(c("x", "y"), "z"); BUD <- 30L; B <- 120L
s <- seq(0.02, 0.98, length.out = 22); EVAL <- as.matrix(expand.grid(x = s, y = s))
field <- function(xy) sin(2 * pi * xy[, 1]) + cos(2 * pi * xy[, 2]) + 0.5 * sin(4 * pi * xy[, 1] * xy[, 2])

one <- function(b) {
  nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, 300, TRUE)
  xy <- ce[cl, ] + cbind(rnorm(300, 0, .05), rnorm(300, 0, .05))
  xy <- cbind(pmin(pmax(xy[, 1], 0), 1), pmin(pmax(xy[, 2], 0), 1))
  d <- data.frame(x = xy[, 1], y = xy[, 2], z = field(xy) + rnorm(300, 0, .15))
  ev <- data.frame(x = EVAL[, 1], y = EVAL[, 2])
  e2 <- (field(EVAL) + rnorm(nrow(EVAL), 0, .15) - idw(d, ev))^2     # squared errors at eval
  true_acc <- sqrt(mean(e2))
  tgt <- prediction_target(grid = EVAL, type = "grid")
  rs <- sample(nrow(EVAL), BUD); est_rand <- sqrt(mean(e2[rs]))
  des <- design_validation(d, tgt, budget = BUD, response = "z", coords = c("x", "y"),
                           allocation = "optimal", seed = b)
  est_opt <- sqrt(stats::weighted.mean(e2[des$index], des$weights))
  c(true = true_acc, rand = est_rand, opt = est_opt)
}
R <- t(vapply(seq_len(B), one, numeric(3)))

err_rand <- R[, "rand"] - R[, "true"]; err_opt <- R[, "opt"] - R[, "true"]
cat(sprintf("replicates: %d, budget: %d validation points, eval grid: %d\n", B, BUD, nrow(EVAL)))
cat("\nEstimating true map accuracy from a few validation points:\n")
cat(sprintf("  random placement : bias %+.4f,  SD %.4f\n", mean(err_rand), sd(err_rand)))
cat(sprintf("  optimal design   : bias %+.4f,  SD %.4f\n", mean(err_opt), sd(err_opt)))
cat(sprintf("\n  optimal design reduces the estimate's SD by %.0f%% (same budget).\n",
            100 * (1 - sd(err_opt) / sd(err_rand))))
