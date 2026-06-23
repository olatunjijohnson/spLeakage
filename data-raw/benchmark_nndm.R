# Benchmark: do the signed SLI / theory predict TRUE optimism better than NNDM's
# (unsigned) Wasserstein W statistic? Ground truth = exact GP optimism (excess
# explained variance, docs/THEORY.md eq. 2), computed from covariance algebra.
# Splits and targets are varied so signed optimism spans pessimism (<0) to optimism
# (>0) -- the regime where W's lack of sign should hurt it.
# Run with: devtools::load_all(); source("data-raw/benchmark_nndm.R")

suppressMessages(devtools::load_all(quiet = TRUE))
set.seed(11)

sample_coords <- function(n, design) {
  if (design == "random") return(cbind(runif(n), runif(n)))
  nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, n, TRUE)
  xy <- ce[cl, ] + cbind(rnorm(n, 0, .05), rnorm(n, 0, .05))
  cbind(pmin(pmax(xy[, 1], 0), 1), pmin(pmax(xy[, 2], 0), 1))
}
EV <- function(Dq, Dtr, w, r) {
  Cq <- w * exp(-Dq / r); S <- w * exp(-Dtr / r); diag(S) <- 1
  rowSums((Cq %*% solve(S)) * Cq)
}
pdist <- function(A, B) sqrt(outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B))

grid <- expand.grid(r = c(.05, .1, .2, .4), w = c(.3, .6, .9), n = 250L,
                    design = c("random", "clustered"),
                    split = c("random", "block"), targ = c("grid", "interp"),
                    rep = 1:3, stringsAsFactors = FALSE)
gs <- seq(.04, .96, length.out = 13); GRIDT <- as.matrix(expand.grid(gs, gs))

rows <- lapply(seq_len(nrow(grid)), function(i) {
  cfg <- grid[i, ]; r <- cfg$r; w <- cfg$w; n <- cfg$n
  S <- sample_coords(n, cfg$design); DSS <- as.matrix(dist(S))
  TG <- if (cfg$targ == "grid") GRIDT else cbind(runif(160, .05, .95), runif(160, .05, .95))
  DTS <- pdist(TG, S)
  ev_p <- EV(DTS, DSS, w, r); f <- apply(DTS, 1, min)
  ids <- if (cfg$split == "random") sample(rep_len(1:10, n)) else spatial_block_cv(S, 10)
  ev_t <- numeric(n); g <- numeric(n)
  for (k in unique(ids)) {
    te <- which(ids == k); tr <- which(ids != k); Dq <- DSS[te, tr, drop = FALSE]
    ev_t[te] <- EV(Dq, DSS[tr, tr, drop = FALSE], w, r); g[te] <- apply(Dq, 1, min)
  }
  ar <- .sli_areas(g, f)
  data.frame(exact = mean(ev_t) - mean(ev_p),                 # signed true optimism
             sli_rho = w * (mean(exp(-g / r)) - mean(exp(-f / r))),
             W = ar$W)                                         # NNDM unsigned statistic
})
res <- do.call(rbind, rows)
saveRDS(res, "data-raw/benchmark_nndm.rds")

cat(sprintf("rows: %d   (%.0f%% pessimistic / optimism<0)\n",
            nrow(res), 100 * mean(res$exact < 0)))
cat("\n--- predicting SIGNED optimism (what determines optimistic vs pessimistic) ---\n")
cat(sprintf("  cor(SLI_rho, optimism) = %+.3f\n", cor(res$sli_rho, res$exact)))
cat(sprintf("  cor(NNDM W, optimism)  = %+.3f   <- unsigned, cannot sign the bias\n",
            cor(res$W, res$exact)))
cat("\n--- predicting the MAGNITUDE |optimism| ---\n")
cat(sprintf("  cor(|SLI_rho|, |optimism|) = %.3f\n", cor(abs(res$sli_rho), abs(res$exact))))
cat(sprintf("  cor(NNDM W,   |optimism|) = %.3f\n", cor(res$W, abs(res$exact))))
cat("\n--- sign agreement (does the diagnostic get optimistic-vs-pessimistic right?) ---\n")
cat(sprintf("  SLI_rho sign matches optimism sign: %.0f%%\n",
            100 * mean(sign(res$sli_rho) == sign(res$exact))))
cat("  (NNDM W has no sign, so it is silent on direction by construction.)\n")
