# Build Figure 5 for the MEE paper: the SDM case study (4 panels).
# Reproduces the audit in data-raw/case_study_sdm.R and draws:
#  (a) per-point leakage map, (b) CV vs deployment NN-distance ECDFs,
#  (c) signed SLI by candidate scheme, (d) Brier score random vs block control.
suppressMessages({library(spLeakage); library(terra)})
set.seed(1)

sp  <- read.csv(system.file("extdata/species.csv", package = "blockCV"))
ed  <- system.file("extdata", package = "blockCV")
r   <- terra::rast(list.files(file.path(ed, "au"), full.names = TRUE))
covn <- names(r)
cov <- terra::extract(r, sp[, c("x", "y")])[, -1]
d   <- cbind(sp, cov); d <- d[stats::complete.cases(d), ]

ga  <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
gi  <- ga[sample(nrow(ga), 4000), c("x", "y")]
tgt <- prediction_target(grid = as.matrix(gi), type = "grid")
folds <- sample(rep_len(1:10, nrow(d)))

lk  <- detect_leakage(d, folds, tgt, response = "occ", coords = c("x", "y"), n_boot = 200)
opt <- estimate_optimism(d, folds, response = "occ", coords = c("x", "y"),
                         metric = "brier", control = "block")
rk  <- rank_cv_schemes(d, tgt, response = "occ", coords = c("x", "y"))

teal <- "#00798c"; red <- "#d1495b"; grey <- "grey80"

draw <- function() {
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, 0.8, 0),
            cex.main = 1.15, font.main = 2)
  on.exit(par(op))

  ## (a) leakage map
  lp <- lk$leak_point; co <- lk$coords; keep <- !is.na(lp)
  rng <- max(abs(lp[keep]))
  pal <- grDevices::colorRampPalette(c(teal, "grey92", red))(100)
  idx <- as.integer(cut(lp[keep], breaks = seq(-rng, rng, length.out = 101),
                        include.lowest = TRUE))
  plot(co[keep, 1] / 1e3, co[keep, 2] / 1e3, col = pal[idx], pch = 19, cex = 0.8,
       asp = 1, xlab = "easting (km)", ylab = "northing (km)",
       main = "(a) Where the random split leaks")
  legend("topright", bty = "n", pch = 19, col = c(red, teal),
         legend = c("leaking (test easy)", "not leaking"), cex = 0.85)

  ## (b) NN-distance ECDFs
  e <- lk$ecdf
  plot(e$r / 1e3, e$Gobs, type = "s", col = red, lwd = 2.5, ylim = c(0, 1),
       xlab = "nearest-neighbour distance (km)", ylab = "ECDF",
       main = "(b) CV vs deployment reach")
  lines(e$r / 1e3, e$Gpred, type = "s", col = teal, lwd = 2.5)
  abline(v = lk$phi / 1e3, lty = 3, col = "grey50")
  legend("bottomright", bty = "n", lwd = 2.5, col = c(red, teal, "grey50"),
         lty = c(1, 1, 3), legend = c("CV (test->train)", "deployment (pred->sample)",
         expression(phi)), cex = 0.85)

  ## (c) signed SLI by scheme
  rr <- rk$ranking[match(c("random", "block", "buffered"), rk$ranking$scheme), ]
  cols <- ifelse(rr$sli_rho > 0.02, red, ifelse(rr$sli_rho < -0.02, teal, grey))
  bp <- barplot(rr$sli_rho, names.arg = rr$scheme, col = cols, border = NA,
                ylim = range(c(rr$sli_rho, 0)) * 1.25,
                ylab = expression(SLI[rho] ~ "(signed)"),
                main = "(c) Which scheme matches deployment")
  abline(h = 0, col = "grey40")
  text(bp, rr$sli_rho, sprintf("%+.2f", rr$sli_rho),
       pos = ifelse(rr$sli_rho >= 0, 3, 1), cex = 0.85, xpd = NA)
  text(bp[2], 0, "best\n(matched)", pos = 1, cex = 0.75, col = teal, xpd = NA)

  ## (d) Brier: reported (random) vs honest (block)
  vals <- c(opt$E_cv, opt$E_control)
  bp2 <- barplot(vals, names.arg = c("random CV\n(reported)", "block CV\n(honest)"),
                 col = c(red, teal), border = NA, ylim = c(0, max(vals) * 1.25),
                 ylab = "Brier score (lower = better)",
                 main = "(d) Reported vs honest accuracy")
  text(bp2, vals, sprintf("%.3f", vals), pos = 3, cex = 0.9, xpd = NA)
  arrows(bp2[1], vals[1] + 0.02, bp2[2], vals[2] + 0.02, length = 0, col = "grey40", lty = 2)
  text(mean(bp2), max(vals) * 1.12,
       sprintf("optimism +%.0f%%", 100 * opt$optimism_rel), cex = 0.95, font = 2)
}

for (f in c("paper-mee/figures/F5_sdm.png", "paper-figures/F5_sdm.png")) {
  grDevices::png(f, width = 2000, height = 1750, res = 210)
  draw()
  grDevices::dev.off()
  message("wrote ", f)
}
