# Hypothesis test for task #24: is the optimism the SLI misses explained by TREND
# STRENGTH (large-scale structure IDW cannot extrapolate)? Compute a trend R^2 per
# meta-audit study and test whether (SLI + trend) predicts real optimism much better
# than SLI alone. Run with: devtools::load_all(); source("data-raw/trend_analysis.R")

suppressMessages(devtools::load_all(quiet = TRUE))

# rebuild the corpus (same as data-raw/meta_audit.R)
build_corpus <- function() {
  L <- list(); add <- function(name, x, y, z) { ok <- is.finite(x) & is.finite(y) & is.finite(z)
    L[[name]] <<- data.frame(x = x[ok], y = y[ok], z = z[ok]) }
  g <- function(nm, pkg) { e <- new.env(); data(list = nm, package = pkg, envir = e); get(nm, e) }
  m <- g("meuse", "sp")
  for (v in c("cadmium", "copper", "lead", "zinc")) add(paste0("meuse:", v), m$x, m$y, log(m[[v]]))
  add("meuse:elev", m$x, m$y, m$elev); add("meuse:om", m$x, m$y, m$om)
  ca <- g("ca20", "geoR"); add("ca20:calcium", ca$coords[, 1], ca$coords[, 2], ca$data)
  pa <- g("parana", "geoR"); add("parana:rain", pa$coords[, 1], pa$coords[, 2], pa$data)
  el <- g("elevation", "geoR"); add("elevation:elev", el$coords[, 1], el$coords[, 2], el$data)
  cm <- g("camg", "geoR"); for (v in c("ca020", "mg020", "ctc020")) add(paste0("camg:", v), cm$east, cm$north, cm[[v]])
  ga <- g("gambia", "geoR"); ag <- aggregate(pos ~ x + y, data = as.data.frame(ga), FUN = mean)
  add("gambia:malaria", ag$x, ag$y, ag$pos)
  oz <- g("ozone2", "fields"); add("ozone2:o3", oz$lon.lat[, 1], oz$lon.lat[, 2], colMeans(oz$y, na.rm = TRUE))
  nr <- g("NorthAmericanRainfall", "fields")
  add("NArain:precip", nr$longitude, nr$latitude, log(nr$precip)); add("NArain:elev", nr$longitude, nr$latitude, nr$elevation)
  L
}

# trend strength = variance of z explained by a quadratic coordinate surface.
trend_r2 <- function(d) {
  s <- as.data.frame(scale(d[, c("x", "y")])); s$z <- d$z
  summary(stats::lm(z ~ poly(x, 2) + poly(y, 2) + x:y, data = s))$r.squared
}

corpus <- build_corpus()
tr <- vapply(corpus, trend_r2, numeric(1))
res <- readRDS("data-raw/meta_audit.rds")
res$trend <- tr[res$study]

cat("Does trend strength explain the optimism the SLI misses?\n")
cat(sprintf("  cor(SLI,   real optimism) = %+.2f\n", cor(res$sli, res$real_opt)))
cat(sprintf("  cor(trend, real optimism) = %+.2f   <- large-scale trend\n", cor(res$trend, res$real_opt)))
m1 <- summary(lm(real_opt ~ abs(sli), res))$r.squared
m2 <- summary(lm(real_opt ~ trend, res))$r.squared
m3 <- summary(lm(real_opt ~ abs(sli) + trend, res))$r.squared
cat(sprintf("\nPredicting real optimism (R^2):  |SLI| only %.2f | trend only %.2f | BOTH %.2f\n", m1, m2, m3))
cat("\nper-study:\n")
print(data.frame(study = res$study, real_opt = round(res$real_opt, 2),
                 sli = round(res$sli, 2), trend = round(res$trend, 2)), row.names = FALSE)
