# Headline figures for the paper, from the calibrated emulator simulation table.
# Run with: devtools::load_all(); source("data-raw/make_figures.R")

suppressMessages(devtools::load_all(quiet = TRUE))
tab <- readRDS("data-raw/optimism_sim_table.rds")
dir.create("paper-figures", showWarnings = FALSE)

dcol <- c(random = "#00798c", clustered = "#d1495b", preferential = "#edae49")
pt <- dcol[tab$design]

# ---- F1: emulator calibration (held-out configs) -------------------------------
set.seed(1)
ho <- sample(unique(tab$config), floor(0.2 * length(unique(tab$config))))
tr <- tab[!tab$config %in% ho, ]; te <- tab[tab$config %in% ho, ]
em <- .fit_optimism_emulator(tr); bb <- .emu_booster(em)
pred <- vapply(seq_len(nrow(te)), function(i)
  .emu_predict_one(em, unlist(te[i, .EMU_FEATURES]), te$model[i], te$response_type[i])$optimism_rel,
  numeric(1))
ok <- is.finite(pred)
R2 <- 1 - sum((te$optimism_rel[ok] - pred[ok])^2) / sum((te$optimism_rel[ok] - mean(te$optimism_rel[ok]))^2)

png("paper-figures/F1_calibration.png", 1400, 1300, res = 200)
rng <- range(c(te$optimism_rel[ok], pred[ok]))
plot(te$optimism_rel[ok], pred[ok], pch = 19, col = "#00798c88",
     xlab = "true optimism (held-out simulation)", ylab = "emulator prediction",
     main = sprintf("Optimism emulator calibration (held-out configs, R2 = %.2f)", R2),
     xlim = rng, ylim = rng)
abline(0, 1, lwd = 2, col = "grey40"); abline(h = 0, v = 0, lty = 3, col = "grey70")
dev.off()

# ---- F2: SLI_rho vs optimism (the C1 -> C2 mechanism) --------------------------
png("paper-figures/F2_sli_vs_optimism.png", 1500, 1300, res = 200)
plot(tab$sli_rho, tab$optimism_rel, pch = 19, col = paste0(pt, "99"),
     xlab = expression("Spatial Leakage Index  SLI"[rho]),
     ylab = "optimism (relative)",
     main = sprintf("Leakage predicts optimism (r = %.2f)", cor(tab$sli_rho, tab$optimism_rel)))
abline(h = 0, lty = 3, col = "grey60"); abline(lm(optimism_rel ~ sli_rho, tab), lwd = 2)
legend("topleft", names(dcol), col = dcol, pch = 19, bty = "n", title = "sampling")
dev.off()

# ---- F3: optimism by sampling design (Mila vs Wadoux) --------------------------
png("paper-figures/F3_design.png", 1300, 1300, res = 200)
boxplot(optimism_rel ~ factor(design, names(dcol)), data = tab,
        col = dcol, ylab = "optimism (relative)", xlab = "sampling design",
        main = "Random CV optimism by sampling design")
abline(h = 0, lty = 2, lwd = 2, col = "grey30")
dev.off()

cat(sprintf("figures written. calibration R2 = %.2f, cor(SLI,opt) = %.2f\n",
            R2, cor(tab$sli_rho, tab$optimism_rel)))
print(round(tapply(tab$optimism_rel, tab$design, mean), 3))
