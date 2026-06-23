# Real-data case study: spatial leakage in Nigeria Pf malaria prevalence mapping.
# Data: Malaria Atlas Project parasite-rate point surveys (open). This also serves
# as an out-of-distribution test of the optimism emulator on real data.
#
# Run with: devtools::load_all(); source("data-raw/case_study_nigeria.R")

suppressMessages({library(spLeakage); library(malariaAtlas); library(sf)})

# ---- data ----------------------------------------------------------------------
pr <- malariaAtlas::getPR(country = "Nigeria", species = "Pf")
g <- pr[!is.na(pr$latitude) & !is.na(pr$longitude) & !is.na(pr$pr) & pr$examined > 0,
        c("longitude", "latitude", "pr", "examined", "positive", "year_start")]
saveRDS(g, "data-raw/nigeria_pr.rds")
sfd <- sf::st_as_sf(g, coords = c("longitude", "latitude"), crs = 4326)
cat(sprintf("Nigeria Pf prevalence: n = %d points (%d unique), prevalence mean %.2f\n",
            nrow(sfd), nrow(unique(g[, 1:2])), mean(g$pr)))

# ---- prediction target: wall-to-wall grid over Nigeria -------------------------
shp <- tryCatch(malariaAtlas::getShp(country = "Nigeria"), error = function(e) NULL)
bb <- sf::st_bbox(sfd)
grd <- expand.grid(longitude = seq(bb["xmin"], bb["xmax"], length.out = 40),
                   latitude  = seq(bb["ymin"], bb["ymax"], length.out = 40))
gsf <- sf::st_as_sf(grd, coords = c("longitude", "latitude"), crs = 4326)
if (!is.null(shp)) {
  shp <- sf::st_make_valid(sf::st_as_sf(shp))
  gsf <- gsf[lengths(sf::st_intersects(gsf, shp)) > 0, ]
}
tgt <- prediction_target(grid = gsf, type = "grid")
cat(sprintf("prediction grid: %d cells over Nigeria\n", nrow(gsf)))

# ---- the naive analysis: random 10-fold CV ------------------------------------
set.seed(1)
folds <- sample(rep_len(1:10, nrow(sfd)))

cat("\n==== 1. DETECT LEAKAGE (random 10-fold CV) ====\n")
lk <- detect_leakage(sfd, split = folds, target = tgt, response = "pr", n_boot = 500)
print(lk)

cat("\n==== 2. QUANTIFY OPTIMISM (empirical, block control) ====\n")
opt <- estimate_optimism(sfd, split = folds, response = "pr", control = "block")
print(opt)

cat("\n==== 3. EMULATOR (out-of-distribution real-data test) ====\n")
pe <- tryCatch(predict_optimism(sfd, split = folds, target = tgt, response = "pr",
                                model = "idw", response_type = "gaussian"),
               warning = function(w) { message("  emulator: ", conditionMessage(w)); NULL })
if (!is.null(pe)) print(pe)

cat("\n==== 3b. ATTRIBUTING THE LEAKAGE: deduplicate co-located surveys ====\n")
# Nigeria PR pools surveys 1985-2018; the variogram finds almost no continuous
# autocorrelation (signal_prop ~ 0), so SLI is driven by EXACT co-locations (repeat
# surveys at the same site split across train/test -> distance 0 -> rho = 1).
# Averaging prevalence per unique location removes that channel; SLI should collapse.
dd <- aggregate(pr ~ longitude + latitude, data = g, FUN = mean)
sfd_dd <- sf::st_as_sf(dd, coords = c("longitude", "latitude"), crs = 4326)
set.seed(1); folds_dd <- sample(rep_len(1:10, nrow(sfd_dd)))
lk_dd <- detect_leakage(sfd_dd, folds_dd, tgt, response = "pr")
cat(sprintf("  with duplicates : n=%d, SLI_rho=%+.3f\n", nrow(sfd), lk$SLI_rho))
cat(sprintf("  deduplicated    : n=%d, SLI_rho=%+.3f  -> leakage was co-location\n",
            nrow(sfd_dd), lk_dd$SLI_rho))

cat("\n==== 4. RECOMMENDATION (convenience survey, prevalence map) ====\n")
rec <- recommend_validation(sfd, estimand = "prediction", design = "convenience",
                            target = "grid")
print(rec)

cat("\n==== 5. AUDIT + SCORECARD ====\n")
print(audit_workflow(sfd, split = folds, target = tgt, response = "pr"))
cat("\n")
print(report_leakage(lk, optimism = opt, recommendation = rec))
