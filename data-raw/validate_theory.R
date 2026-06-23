# Numerical validation of docs/THEORY.md. Explained variance EV = c0' Sigma^-1 c0 is
# pure covariance algebra (no field realisation needed), so the exact GP optimism
# (eq. 2) is computed directly and compared to the single-NN closed form (eq. 6) and
# the Wasserstein bound (eq. 8). Latent correlation rhoZ(d)=exp(-d/r); V=1.
# Run with: devtools::load_all(); source("data-raw/validate_theory.R")

suppressMessages(devtools::load_all(quiet = TRUE))
set.seed(7)

sample_coords <- function(n, design) {
  if (design == "random") return(cbind(runif(n), runif(n)))
  nc <- 8; ce <- cbind(runif(nc, .1, .9), runif(nc, .1, .9)); cl <- sample(nc, n, TRUE)
  xy <- ce[cl, ] + cbind(rnorm(n, 0, .05), rnorm(n, 0, .05))
  cbind(pmin(pmax(xy[, 1], 0), 1), pmin(pmax(xy[, 2], 0), 1))
}
EV <- function(Dq, Dtr, w, r) {           # explained variance per query row
  Cq <- w * exp(-Dq / r)                   # q x ntr cross-covariance
  S  <- w * exp(-Dtr / r); diag(S) <- 1    # ntr x ntr (= wR + (1-w)I)
  Si <- solve(S)
  rowSums((Cq %*% Si) * Cq)
}
nnd <- function(Dq) apply(Dq, 1, min)

grid <- expand.grid(r = c(.05, .1, .2, .4), w = c(.3, .6, .9),
                    n = c(150L, 300L), design = c("random", "clustered"),
                    rep = 1:4, stringsAsFactors = FALSE)
gs <- seq(.04, .96, length.out = 13); TG <- as.matrix(expand.grid(gs, gs))

rows <- lapply(seq_len(nrow(grid)), function(i) {
  cfg <- grid[i, ]; r <- cfg$r; w <- cfg$w; n <- cfg$n
  S <- sample_coords(n, cfg$design)
  DSS <- as.matrix(dist(S))
  DTS <- as.matrix(proxy_dist <- outer(seq_len(nrow(TG)), seq_len(n),
                  Vectorize(function(a, b) sqrt(sum((TG[a, ] - S[b, ])^2)))))
  # target explained variance + NN distance (deployment)
  ev_p <- EV(DTS, DSS, w, r); f <- nnd(DTS)
  # random 10-fold: per-fold test EV + NN distance
  fold <- sample(rep_len(1:10, n)); ev_t <- numeric(n); g <- numeric(n)
  for (k in 1:10) {
    te <- which(fold == k); tr <- which(fold != k)
    Dq <- DSS[te, tr, drop = FALSE]
    ev_t[te] <- EV(Dq, DSS[tr, tr, drop = FALSE], w, r)
    g[te] <- apply(Dq, 1, min)
  }
  rhoZ2 <- function(d) exp(-2 * d / r)
  data.frame(
    exact   = mean(ev_t) - mean(ev_p),                       # eq. 2 (exact optimism)
    nn_form = w^2 * (mean(rhoZ2(g)) - mean(rhoZ2(f))),        # eq. 6 (single-NN form)
    sli2    = mean(rhoZ2(g)) - mean(rhoZ2(f)),                # squared-correlation SLI
    sli_rho = w * (mean(exp(-g / r)) - mean(exp(-f / r))),    # package-style SLI_rho
    W1      = .sli_areas(g, f)$W,                             # Wasserstein-1 (unsigned)
    Lprime  = w^2 * 2 / r,                                    # ||L'||_inf (exponential)
    r = r, w = w, n = n, design = cfg$design)
})
res <- do.call(rbind, rows)
res$bound <- res$Lprime * res$W1

cat(sprintf("configs: %d\n", nrow(res)))
cat(sprintf("[eq.6] cor(single-NN form, exact optimism) = %.3f\n",
            cor(res$nn_form, res$exact)))
cat(sprintf("[eq.6] median ratio exact/nn_form          = %.2f (multi-NN amplification)\n",
            median(res$exact / res$nn_form)))
cat(sprintf("[feature] cor(SLI_rho, exact) = %.3f   cor(SLI_2, exact) = %.3f\n",
            cor(res$sli_rho, res$exact), cor(res$sli2, res$exact)))
cat(sprintf("[eq.8] bound holds for single-NN optimism in %.0f%% (Kantorovich, must be 100%%)\n",
            100 * mean(abs(res$nn_form) <= res$bound + 1e-9)))
cat(sprintf("[eq.8] bound holds for EXACT (multi-NN) optimism in %.0f%%\n",
            100 * mean(abs(res$exact) <= res$bound)))
cat(sprintf("[eq.8] median tightness |nn_form|/bound    = %.2f\n",
            median(abs(res$nn_form) / res$bound)))
saveRDS(res, "data-raw/theory_validation.rds")
