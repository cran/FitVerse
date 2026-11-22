## ============================================================================
## FitVerse: Replication Script
## Reproduces all examples from the JSS paper:
##   "FitVerse: Parametric Distribution Fitting and Analysis"
##   Reddy and Khan, Journal of Statistical Software
## ============================================================================
##
## Requirements:
##   install.packages(c("FitVerse", "fitdistrplus", "evd"))
##
## Run from an R session in the directory that contains this script.
## Figures are saved as PDF files in the working directory.
## ============================================================================

library("FitVerse")
library("fitdistrplus")
library("evd")

## ----------------------------------------------------------------------------
## Section 3.2 - Groundbeef: basic three-method fit
## ----------------------------------------------------------------------------

data("groundbeef", package = "fitdistrplus")
fit <- fitverse(
  groundbeef$serving,
  dists   = c("gamma", "lognormal", "weibull", "loglogistic"),
  method  = c("MLE", "MOM", "LMOM"),
  xlab    = "Serving size (g)",
  verbose = FALSE
)

## Ranking table (top results by AIC)
ranking(fit)[, c("Distribution", "Method", "AIC", "Delta_AIC", "AD_p")]

## Summary
s <- summary(fit)
print(s)

## Accessor functions
head(ranking(fit))

## Diagnostic figure
save_plot(fit, file = "fig_groundbeef.pdf", width = 10, height = 5)
message("Saved: fig_groundbeef.pdf")


## ----------------------------------------------------------------------------
## Section 3.3 - Port Pirie sea levels: extreme-value analysis
## ----------------------------------------------------------------------------

data("portpirie", package = "evd")
fit_sea <- fitverse(
  as.numeric(portpirie),
  dists   = c("gev", "gumbel", "normal", "logistic"),
  method  = c("MLE", "LMOM"),
  xlab    = "Annual maximum sea level (m)",
  verbose = FALSE
)

summary(fit_sea)

## Return levels
return_level(fit_sea, return_periods = c(10, 50, 100, 200))

## Bootstrap confidence intervals (B = 999 for paper; use B = 99 for quick run)
set.seed(1)
ci <- bootstrap_ci(fit_sea, B = 999, return_periods = c(10, 50, 100, 200))
ci

## Diagnostic figure
save_plot(fit_sea, file = "fig_portpirie.pdf", width = 10, height = 5)
message("Saved: fig_portpirie.pdf")


## ----------------------------------------------------------------------------
## Section 3.4 - Censored data: fitverse_censored
## ----------------------------------------------------------------------------

set.seed(1)
times   <- rweibull(100, shape = 1.8, scale = 50)
status  <- as.integer(times <= 60)
times   <- pmin(times, 60)
fit_cens <- fitverse_censored(
  times,
  status = status,
  dists  = c("weibull", "lognormal", "gamma"),
  method = "MLE"
)
summary(fit_cens)


## ----------------------------------------------------------------------------
## Section 3.5 - Multi-method comparison: simulated Gamma data
## ----------------------------------------------------------------------------

set.seed(42)
x <- rgamma(150, shape = 2, rate = 0.5)
fit_sim <- fitverse(
  x,
  dists   = c("lognormal", "gamma", "weibull", "loglogistic"),
  method  = c("MLE", "MOM", "LMOM"),
  xlab    = "Simulated Gamma(2, 0.5)",
  verbose = FALSE
)

ranking(fit_sim)[, c("Distribution", "Method", "AIC", "Delta_AIC", "AD_p")]

## Inspect per-method Gamma fits using accessor
gamma_fits <- all_fits(fit_sim)[
  grep("^gamma_", names(all_fits(fit_sim)), value = TRUE)
]
sapply(gamma_fits, function(f) c(AIC = f$aic, AD_p = f$gof$P_value[2]))


## ----------------------------------------------------------------------------
## Section 3.6 - Grouped / interval-censored data: Danish fire insurance
## ----------------------------------------------------------------------------

## Danish fire insurance losses 1980-1990, grouped (McNeil 1997, Table 1)
## Losses in millions of DKK; only claims >= 1 million DKK are recorded.
## The file danish_fire_grouped.csv is bundled with the package.
g_file <- system.file("extdata", "danish_fire_grouped.csv", package = "FitVerse")
if (nchar(g_file) > 0 && file.exists(g_file)) {
  g <- read_grouped_csv(g_file)
  fit_dk <- fitverse_grouped(
    lower  = g$lower,
    upper  = g$upper,
    counts = g$counts,
    dists  = c("pareto", "lomax", "lognormal", "burr", "loglogistic"),
    method = "MLE",
    xlab   = "Loss (millions DKK)"
  )
  print(fit_dk)
  ranking(fit_dk)[, c("Distribution", "Method", "AIC", "Delta_AIC")]
} else {
  message("Note: danish_fire_grouped.csv not found in package extdata; skipping grouped example.")
}


## ----------------------------------------------------------------------------
## Section 4 - Batch fitting: fitverse_batch
## ----------------------------------------------------------------------------

set.seed(1)
df <- data.frame(
  x1 = rgamma(200, shape = 2, scale = 3),
  x2 = rlnorm(200, meanlog = 1, sdlog = 0.5)
)
batch_res <- fitverse_batch(df, method = "MLE")
batch_res$summary[, c("Column", "Distribution", "AIC")]


## ----------------------------------------------------------------------------
## Section 5 - compare_dists utility
## ----------------------------------------------------------------------------

set.seed(1)
x_demo <- rgamma(200, shape = 2, scale = 3)
tbl <- compare_dists(x_demo, method = "MLE")
head(tbl)


## ============================================================================
## End of replication script
## ============================================================================
message("\nAll replication examples completed successfully.")
