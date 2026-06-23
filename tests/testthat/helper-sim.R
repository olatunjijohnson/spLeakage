# Deterministic-ish spatially autocorrelated test data on the unit square.
# A smooth latent field (sum of low-frequency sinusoids) + small nugget noise gives
# a well-defined variogram range, without needing a GRF simulator dependency.
sim_spatial <- function(n = 400L, noise = 0.15, seed = 1L) {
  set.seed(seed)
  xy <- cbind(runif(n), runif(n))
  field <- sin(2 * pi * xy[, 1]) + cos(2 * pi * xy[, 2]) +
    0.5 * sin(4 * pi * (xy[, 1] + xy[, 2]))
  z <- field + rnorm(n, 0, noise)
  data.frame(x = xy[, 1], y = xy[, 2], z = z)
}

# Clustered sample of the same field: points gathered around a few centres, leaving
# under-sampled gaps. This is the canonical leakage regime (cf. Ploton et al. 2020):
# random CV neighbours are within-cluster (tiny), but a wall-to-wall grid reaches the
# gaps (large), so random CV is strongly optimistic.
sim_spatial_clustered <- function(n = 500L, n_clusters = 8L, spread = 0.04,
                                  noise = 0.15, seed = 3L) {
  set.seed(seed)
  centers <- cbind(runif(n_clusters, 0.1, 0.9), runif(n_clusters, 0.1, 0.9))
  cl <- sample(n_clusters, n, replace = TRUE)
  xy <- centers[cl, ] + cbind(rnorm(n, 0, spread), rnorm(n, 0, spread))
  xy[, 1] <- pmin(pmax(xy[, 1], 0), 1)
  xy[, 2] <- pmin(pmax(xy[, 2], 0), 1)
  field <- sin(2 * pi * xy[, 1]) + cos(2 * pi * xy[, 2]) +
    0.5 * sin(4 * pi * (xy[, 1] + xy[, 2]))
  z <- field + rnorm(n, 0, noise)
  data.frame(x = xy[, 1], y = xy[, 2], z = z)
}

# Gaussian-random-field sample matching the emulator's training generator
# (exponential covariance), so emulator queries land inside its area of
# applicability. range/signal are within the training grid by default.
sim_gp <- function(n = 300L, range = 0.15, signal = 0.6,
                   design = c("clustered", "random"), seed = 11L) {
  design <- match.arg(design)
  set.seed(seed)
  if (design == "clustered") {
    nc <- 8L; ce <- cbind(runif(nc, 0.1, 0.9), runif(nc, 0.1, 0.9))
    cl <- sample(nc, n, replace = TRUE)
    xy <- ce[cl, ] + cbind(rnorm(n, 0, 0.04), rnorm(n, 0, 0.04))
    xy <- cbind(pmin(pmax(xy[, 1], 0), 1), pmin(pmax(xy[, 2], 0), 1))
  } else {
    xy <- cbind(runif(n), runif(n))
  }
  D <- as.matrix(dist(xy)); L <- chol(exp(-D / range) + diag(1e-8, n))
  f <- as.numeric(t(L) %*% rnorm(n))
  data.frame(x = xy[, 1], y = xy[, 2], z = f + rnorm(n, 0, sqrt((1 - signal) / signal)))
}

# Random k-fold assignment.
random_folds <- function(n, k = 10L, seed = 2L) {
  set.seed(seed)
  sample(rep_len(seq_len(k), n))
}

# Spatially blocked k-fold: assign folds by a coarse grid of blocks.
block_folds <- function(xy, k = 10L) {
  nb <- ceiling(sqrt(k))
  bx <- cut(xy[, 1], breaks = nb, labels = FALSE)
  by <- cut(xy[, 2], breaks = nb, labels = FALSE)
  block <- (bx - 1L) * nb + by
  match(block, sort(unique(block)))
}

# A coarse wall-to-wall prediction grid over the unit square.
grid_target <- function(g = 20L) {
  s <- seq(0.025, 0.975, length.out = g)
  as.matrix(expand.grid(x = s, y = s))
}
