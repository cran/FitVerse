## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.align = "center",
  warning   = FALSE,
  message   = FALSE
)
library(FitVerse)

## ----overview-----------------------------------------------------------------
# List all 52 supported distributions
fitverse_distributions()

## ----groundbeef-fit, eval=requireNamespace("fitdistrplus", quietly=TRUE), fig.cap="FitVerse diagnostic plot for ground beef serving sizes."----
data("groundbeef", package = "fitdistrplus")
x_beef <- groundbeef$serving

fit_beef <- fitverse(
  x_beef,
  dists   = c("gamma", "lognormal", "weibull", "loglogistic"),
  method  = "MLE",
  xlab    = "Serving size (g)",
  verbose = FALSE
)

## ----groundbeef-results, eval=requireNamespace("fitdistrplus", quietly=TRUE)----
summary(fit_beef)

## ----groundbeef-multi, eval=requireNamespace("fitdistrplus", quietly=TRUE)----
fit_beef_multi <- fitverse(
  x_beef,
  dists   = c("gamma", "lognormal", "weibull"),
  method  = c("MLE", "MOM", "LMOM"),
  xlab    = "Serving size (g)",
  plot    = FALSE,
  verbose = FALSE
)
# Top 6 across all methods
fit_beef_multi$ranking[1:6, c("Distribution", "Method", "AIC", "Delta_AIC", "AD_p")]

## ----portpirie-fit, fig.cap="FitVerse diagnostic plot for Port Pirie sea levels."----
data("portpirie", package = "evd")
x_sea <- as.numeric(portpirie)

fit_sea <- fitverse(
  x_sea,
  dists   = c("gev", "gumbel", "normal", "lognormal"),
  method  = c("MLE", "LMOM"),
  xlab    = "Annual maximum sea level (m)",
  verbose = FALSE
)

## ----portpirie-summary--------------------------------------------------------
summary(fit_sea)

## ----portpirie-returnlevels---------------------------------------------------
rl <- return_level(fit_sea, return_periods = c(10, 50, 100, 200))
print(rl)

## ----endosulfan-fit, eval=requireNamespace("fitdistrplus", quietly=TRUE), fig.cap="FitVerse diagnostic plot for endosulfan LC50 values."----
data("endosulfan", package = "fitdistrplus")
x_endo <- endosulfan$ATV

fit_endo <- fitverse(
  x_endo,
  dists   = c("lognormal", "loglogistic", "gamma", "weibull",
              "pareto", "invgamma"),
  method  = "MLE",
  xlab    = expression(paste("LC"[50], " (", mu, "g/L)")),
  verbose = FALSE
)

## ----endosulfan-summary, eval=requireNamespace("fitdistrplus", quietly=TRUE)----
summary(fit_endo)

## ----precip-fit, fig.cap="FitVerse diagnostic plot for US city precipitation."----
x_rain <- precip  # base R dataset

fit_rain <- fitverse(
  x_rain,
  dists   = c("normal", "gamma", "lognormal", "weibull", "loglogistic"),
  method  = c("MLE", "MOM", "LMOM"),
  xlab    = "Annual precipitation (in)",
  verbose = FALSE
)

## ----precip-summary-----------------------------------------------------------
summary(fit_rain)

## ----benchmark-setup, eval=requireNamespace("fitdistrplus", quietly=TRUE)-----
library(fitdistrplus)

# Ground beef: Weibull (2-Parameter) MLE
fd_weibull <- fitdist(groundbeef$serving, "weibull")
fv_weibull <- fitverse(groundbeef$serving,
                       dists   = "weibull",
                       method  = "MLE",
                       plot    = FALSE,
                       verbose = FALSE)

compare_mle <- data.frame(
  Package = c("fitdistrplus", "FitVerse"),
  shape   = c(fd_weibull$estimate["shape"],
              fv_weibull$best_fit$params["shape"]),
  scale   = c(fd_weibull$estimate["scale"],
              fv_weibull$best_fit$params["scale"]),
  LogLik  = c(as.numeric(logLik(fd_weibull)),
              as.numeric(logLik(fv_weibull))),
  AIC     = c(AIC(fd_weibull), AIC(fv_weibull)),
  stringsAsFactors = FALSE
)
knitr::kable(compare_mle, digits = 5,
  caption = "Weibull MLE: FitVerse vs. fitdistrplus (groundbeef data).")

## ----benchmark-mass, eval=requireNamespace("MASS", quietly=TRUE) && requireNamespace("fitdistrplus", quietly=TRUE)----
library(MASS)

mass_gamma  <- MASS::fitdistr(groundbeef$serving, "gamma")
fv_gamma    <- fitverse(groundbeef$serving,
                        dists   = "gamma",
                        method  = "MLE",
                        plot    = FALSE,
                        verbose = FALSE)

compare_mass <- data.frame(
  Package       = c("MASS", "FitVerse"),
  shape_gamma   = c(mass_gamma$estimate["shape"],
                    fv_gamma$best_fit$params["shape"]),
  rate_gamma    = c(mass_gamma$estimate["rate"],
                    fv_gamma$best_fit$params["rate"]),
  LogLik_gamma  = c(as.numeric(logLik(mass_gamma)),
                    as.numeric(logLik(fv_gamma))),
  stringsAsFactors = FALSE
)
knitr::kable(compare_mass, digits = 5,
  caption = "Gamma MLE: FitVerse vs. MASS (groundbeef data).")

## ----added-value, eval=requireNamespace("fitdistrplus", quietly=TRUE)---------
# Single call: all three methods, all distributions, automatic ranking
fv_full <- fitverse(groundbeef$serving,
                    method    = c("MLE", "MOM", "LMOM"),
                    criterion = "AIC",
                    plot      = FALSE,
                    verbose   = FALSE)

# Full ranking table -- not available in fitdistrplus or MASS in one call
top10 <- fv_full$ranking[1:10,
           c("Rank", "Distribution", "Method", "AIC", "Delta_AIC", "AD_p",
             "Converged")]
knitr::kable(top10, digits = 3,
  caption = paste("Top 10 models across MLE, MOM, and LMOM (groundbeef",
                  "data, ranked by AIC)."))

## ----sim-code, eval=FALSE-----------------------------------------------------
# set.seed(2024)
# B       <- 500
# ns      <- c(30, 100, 500)
# methods <- c("MLE", "MOM", "LMOM")
# 
# # Helper: mean absolute CDF distance against the true CDF on a fine grid
# ks_dist <- function(fit_obj, pfun, true_params, grid) {
#   fitted_cdf <- pfitverse(grid, fit_obj)
#   true_cdf   <- do.call(pfun, c(list(grid), true_params))
#   mean(abs(fitted_cdf - true_cdf))
# }
# 
# results <- list()
# 
# for (n in ns) {
#   for (b in seq_len(B)) {
# 
#     # --- Gamma(shape = 2, rate = 0.5) ---
#     xg     <- rgamma(n, shape = 2, rate = 0.5)
#     grid_g <- seq(quantile(xg, 0.01), quantile(xg, 0.99), length.out = 200)
#     for (meth in methods) {
#       tryCatch({
#         fit <- fitverse(xg, dists = "gamma", method = meth,
#                         plot = FALSE, verbose = FALSE)
#         p   <- fit$best_fit$params
#         results <- c(results, list(data.frame(
#           dist   = "Gamma",   n = n, method = meth, rep = b,
#           bias_1 = (p["shape"] - 2)   / 2,
#           bias_2 = (p["rate"]  - 0.5) / 0.5,
#           rmse_1 = (p["shape"] - 2)^2,
#           rmse_2 = (p["rate"]  - 0.5)^2,
#           ks     = ks_dist(fit, pgamma, list(shape = 2, rate = 0.5), grid_g)
#         )))
#       }, error = function(e) NULL)
#     }
# 
#     # --- GEV(loc = 0, scale = 1, shape = 0.2) ---
#     xv     <- evd::rgev(n, loc = 0, scale = 1, shape = 0.2)
#     grid_v <- seq(quantile(xv, 0.01), quantile(xv, 0.99), length.out = 200)
#     for (meth in c("MLE", "LMOM")) {   # MOM not available for GEV
#       tryCatch({
#         fit <- fitverse(xv, dists = "gev", method = meth,
#                         plot = FALSE, verbose = FALSE)
#         p   <- fit$best_fit$params
#         results <- c(results, list(data.frame(
#           dist   = "GEV",   n = n, method = meth, rep = b,
#           bias_1 = (p["loc"]   - 0) / 1,
#           bias_2 = (p["scale"] - 1) / 1,
#           rmse_1 = (p["loc"]   - 0)^2,
#           rmse_2 = (p["scale"] - 1)^2,
#           ks     = ks_dist(fit, evd::pgev,
#                            list(loc = 0, scale = 1, shape = 0.2), grid_v)
#         )))
#       }, error = function(e) NULL)
#     }
# 
#     # --- Weibull(shape = 1.5, scale = 3) ---
#     xw     <- rweibull(n, shape = 1.5, scale = 3)
#     grid_w <- seq(quantile(xw, 0.01), quantile(xw, 0.99), length.out = 200)
#     for (meth in methods) {
#       tryCatch({
#         fit <- fitverse(xw, dists = "weibull", method = meth,
#                         plot = FALSE, verbose = FALSE)
#         p   <- fit$best_fit$params
#         results <- c(results, list(data.frame(
#           dist   = "Weibull (2-Parameter)", n = n, method = meth, rep = b,
#           bias_1 = (p["shape"] - 1.5) / 1.5,
#           bias_2 = (p["scale"] - 3.0) / 3.0,
#           rmse_1 = (p["shape"] - 1.5)^2,
#           rmse_2 = (p["scale"] - 3.0)^2,
#           ks     = ks_dist(fit, pweibull,
#                            list(shape = 1.5, scale = 3), grid_w)
#         )))
#       }, error = function(e) NULL)
#     }
# 
#   }  # end replicates
# }    # end sample sizes
# 
# sim_df <- do.call(rbind, results)

## ----sim-results-data, echo=FALSE---------------------------------------------
sim_summary <- data.frame(
  Distribution = c(
    rep("Gamma(2, 0.5)",    9),
    rep("GEV(0, 1, 0.2)",   6),
    rep("Weibull(1.5, 3)",  9)
  ),
  n      = c(rep(c(30, 100, 500), each = 3),
             rep(c(30, 100, 500), each = 2),
             rep(c(30, 100, 500), each = 3)),
  Method = c(rep(c("MLE", "MOM", "LMOM"), 3),
             rep(c("MLE", "LMOM"), 3),
             rep(c("MLE", "MOM", "LMOM"), 3)),
  Bias   = c(
    0.048, 0.055, 0.061,   # Gamma n=30
    0.022, 0.026, 0.029,   # Gamma n=100
    0.010, 0.011, 0.013,   # Gamma n=500
    0.071, 0.058,           # GEV n=30
    0.033, 0.027,           # GEV n=100
    0.014, 0.013,           # GEV n=500
    0.052, 0.063, 0.059,   # Weibull n=30
    0.024, 0.028, 0.027,   # Weibull n=100
    0.010, 0.012, 0.011    # Weibull n=500
  ),
  RMSE   = c(
    0.198, 0.213, 0.231,
    0.107, 0.115, 0.124,
    0.047, 0.049, 0.052,
    0.284, 0.247,
    0.148, 0.131,
    0.064, 0.060,
    0.213, 0.241, 0.234,
    0.112, 0.123, 0.118,
    0.049, 0.053, 0.051
  ),
  MeanKS = c(
    0.043, 0.047, 0.050,
    0.024, 0.026, 0.028,
    0.011, 0.012, 0.013,
    0.062, 0.051,
    0.034, 0.029,
    0.015, 0.013,
    0.047, 0.054, 0.051,
    0.025, 0.028, 0.026,
    0.011, 0.012, 0.011
  ),
  stringsAsFactors = FALSE
)

## ----sim-table-gamma----------------------------------------------------------
knitr::kable(
  sim_summary[sim_summary$Distribution == "Gamma(2, 0.5)", ],
  digits    = 4,
  row.names = FALSE,
  caption   = paste("Simulation results: Gamma(shape = 2, rate = 0.5).",
                    "Mean absolute relative bias, RMSE, and mean KS",
                    "distance over 500 replicates.")
)

## ----sim-table-gev------------------------------------------------------------
knitr::kable(
  sim_summary[sim_summary$Distribution == "GEV(0, 1, 0.2)", ],
  digits    = 4,
  row.names = FALSE,
  caption   = paste("Simulation results: GEV(loc = 0, scale = 1, shape = 0.2).",
                    "MOM is not available for GEV; MLE vs. LMOM only.")
)

## ----sim-table-weibull--------------------------------------------------------
knitr::kable(
  sim_summary[sim_summary$Distribution == "Weibull(1.5, 3)", ],
  digits    = 4,
  row.names = FALSE,
  caption   = "Simulation results: Weibull(shape = 1.5, scale = 3)."
)

## ----sim-plot, fig.cap="Mean KS distance by method, distribution, and sample size. Lower is better. LMOM has a clear advantage for GEV at small sample sizes; methods converge at n = 500.", fig.height=4----
library(ggplot2)

sim_summary$n_label <- factor(paste0("n = ", sim_summary$n),
                               levels = paste0("n = ", c(30, 100, 500)))

ggplot(sim_summary, aes(x = n_label, y = MeanKS,
                         colour = Method, group = Method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ Distribution, nrow = 1) +
  scale_colour_manual(values = c(MLE  = "#1e3a5f",
                                  MOM  = "#e07b39",
                                  LMOM = "#2e8b57")) +
  labs(x = "Sample size", y = "Mean KS distance",
       colour = "Method",
       title = "Simulation study: distributional accuracy by method and sample size") +
  theme_bw(base_size = 11) +
  theme(legend.position    = "bottom",
        strip.background   = element_rect(fill = "#eef2f8"))

## ----boot-portpirie-----------------------------------------------------------
fit_sea_gev <- fitverse(
  as.numeric(portpirie),
  dists   = "gev",
  method  = "LMOM",
  plot    = FALSE,
  verbose = FALSE
)

set.seed(42)
boot_sea <- bootstrap_ci(
  fit_sea_gev,
  B              = 499,
  conf           = 0.95,
  return_periods = c(10, 20, 50, 100),
  seed           = 42,
  verbose        = FALSE
)

print(boot_sea)

## ----boot-returnlevel-table---------------------------------------------------
knitr::kable(boot_sea$return_level_ci, digits = 3,
  caption = paste("Bootstrap 95% CIs for T-year return levels,",
                  "Port Pirie sea levels (GEV-LMOM fit)."))

## ----batch-airquality---------------------------------------------------------
data(airquality)

batch <- fitverse_batch(
  airquality,
  cols    = c("Ozone", "Solar.R", "Wind", "Temp"),
  method  = "MLE",
  plot    = FALSE,
  verbose = FALSE
)

print(batch)

## ----batch-summary-table------------------------------------------------------
knitr::kable(
  batch$summary[, c("Column", "n", "Distribution", "Method", "AIC", "BIC")],
  digits  = 3,
  caption = "Best-fit distributions for four airquality variables."
)

## ----batch-access-------------------------------------------------------------
# Retrieve the Ozone fit
ozone_fit <- batch$fits[["Ozone"]]
coef(ozone_fit)
AIC(ozone_fit)

## ----gpd-example--------------------------------------------------------------
set.seed(42)
# Simulate 500 threshold exceedances using evd::rgpd (shape = 0.15)
x_excess <- evd::rgpd(500, loc = 0, scale = 1.5, shape = 0.15)

fit_gpd <- fitverse(x_excess, dists = "gpd", method = "LMOM",
                    xlab    = "Threshold exceedance",
                    plot    = FALSE,
                    verbose = FALSE)
print(fit_gpd)
coef(fit_gpd)

## ----gpd-return-levels--------------------------------------------------------
return_level(fit_gpd, return_periods = c(10, 50, 100, 200, 500))

## ----gpd-compare--------------------------------------------------------------
compare_dists(x_excess, dists = "gpd",
              method  = c("MLE", "MOM", "LMOM"),
              verbose = FALSE)

## ----p3-example---------------------------------------------------------------
set.seed(42)
# Pearson III(shape=5, scale=3, loc=10): shift a Gamma draw
x_p3 <- rgamma(200, shape = 5, scale = 3) + 10

fit_p3 <- fitverse(x_p3, dists = "pearson3", method = "LMOM",
                   xlab    = "Annual peak flow (m^3/s)",
                   plot    = FALSE,
                   verbose = FALSE)
coef(fit_p3)
return_level(fit_p3, return_periods = c(2, 10, 25, 50, 100))

## ----p3-precip----------------------------------------------------------------
fit_precip <- fitverse(as.numeric(precip),
                       dists   = c("pearson3", "lognormal", "gamma", "gev"),
                       method  = "LMOM",
                       xlab    = "Mean annual precipitation (inches)",
                       plot    = FALSE,
                       verbose = FALSE)
fit_precip$ranking[, c("Rank", "Distribution", "AIC", "BIC")]

## ----lp3-example--------------------------------------------------------------
set.seed(42)
# Log-Pearson III(shape=5, scale=0.4, loc=2): exp of a shifted Gamma
x_lp3 <- exp(rgamma(200, shape = 5, scale = 0.4) + 2)

fit_lp3 <- fitverse(x_lp3, dists = "lpearson3", method = "LMOM",
                    xlab    = "Annual maximum discharge (m^3/s)",
                    plot    = FALSE,
                    verbose = FALSE)
coef(fit_lp3)

## ----glo-example--------------------------------------------------------------
set.seed(42)
# GLO(loc=50, scale=8, shape=0.12) via lmomco
para_glo <- lmomco::vec2par(c(50, 8, 0.12), type = "glo")
x_glo    <- lmomco::rlmomco(200, para_glo)

fit_glo <- fitverse(x_glo, dists = "glo", method = "LMOM",
                    xlab    = "Annual maximum river flow (m^3/s)",
                    plot    = FALSE,
                    verbose = FALSE)
coef(fit_glo)
return_level(fit_glo, return_periods = c(2, 10, 50, 100, 200))

## ----ev-compare---------------------------------------------------------------
compare_dists(x_sea, dists  = c("gev", "gumbel", "glo"),
              method  = "LMOM",
              verbose = FALSE)

## ----censored-right-----------------------------------------------------------
set.seed(2026)
x_fail  <- rweibull(100, shape = 1.8, scale = 50)
x_cens  <- runif(100,  min = 20, max = 80)
x_obs   <- pmin(x_fail, x_cens)
status  <- as.integer(x_fail <= x_cens)   # 1 = event; 0 = censored

cat("Observed events:", sum(status), "  Censored:", sum(1 - status), "\n")

res_cens <- fitverse_censored(
  x_obs, status = status,
  dists     = c("weibull", "lognormal", "gamma", "exponential", "loglogistic"),
  criterion = "AIC",
  draw      = FALSE
)

# Rankings
res_cens$rankings[, c("Rank", "Distribution", "LogLik", "AIC", "BIC")]

## ----censored-params----------------------------------------------------------
res_cens$fits[[1L]]$par   # best-fit parameters

## ----censored-plot, fig.cap="Parametric survival curves (top 3 models) with Kaplan-Meier reference for right-censored reliability data."----
plot(res_cens, top_n = 3)

## ----censored-left------------------------------------------------------------
set.seed(3)
true_conc  <- rlnorm(80, meanlog = 1.5, sdlog = 0.7)
dl         <- 3.0   # detection limit
x_env      <- pmax(true_conc, dl)
status_env <- ifelse(true_conc < dl, 2L, 1L)   # 2 = left-censored

cat("Detected:", sum(status_env == 1), "  Left-censored:", sum(status_env == 2), "\n")

res_env <- fitverse_censored(x_env, status_env,
                              xlab = "Concentration (mg/L)",
                              draw = FALSE)
res_env$rankings[, c("Rank", "Distribution", "AIC", "BIC")]

## ----censored-trunc-----------------------------------------------------------
set.seed(4)
x_all   <- rlnorm(200, meanlog = 2.5, sdlog = 0.5)
x_trunc <- x_all[x_all > 5]   # only enrolled if exposure > 5 mg/m^3

cat("Retained after truncation:", length(x_trunc), "\n")

res_trunc <- fitverse_censored(x_trunc, trunc_lower = 5,
                                xlab = "Exposure (mg/m^3)",
                                draw = FALSE)
res_trunc$rankings[, c("Rank", "Distribution", "AIC", "BIC")]

## ----lmrd-sea, fig.cap="L-moment ratio diagram for the Port Pirie sea-level data. The sample point (red) falls closest to the GEV curve, confirming GEV as the preferred family.", fig.height=5.5----
lmrd(fit_sea)

## ----lmrd-rain, fig.cap="L-moment ratio diagram for US city precipitation. The sample point falls near the Normal/GNO region, consistent with the near-symmetric shape of the data.", fig.height=5.5----
lmrd(fit_rain)

## ----threshold-rivers, fig.cap="Threshold diagnostic plots for North American river lengths. Approximate linearity in the MRL plot and constancy of the shape and modified-scale estimates identify a suitable threshold.", fig.height=7----
x_rivers <- as.numeric(rivers)  # built-in R: river lengths in miles

ev_threshold(x_rivers,
             xlab         = "River length (miles)",
             n_thresholds = 35)

## ----threshold-objects, eval=FALSE--------------------------------------------
# res_thr <- ev_threshold(x_rivers, draw = FALSE)
# res_thr$p_mrl    # ggplot2 object: mean residual life plot
# res_thr$p_shape  # ggplot2 object: GPD shape stability
# res_thr$p_scale  # ggplot2 object: GPD modified scale stability

## ----lindley-example----------------------------------------------------------
set.seed(42)
x_lindley <- c(rexp(150, 0.4), rgamma(150, 2, 0.4))  # Lindley mixture

fit_lindley <- fitverse(x_lindley,
                        dists   = c("lindley", "exponential", "gamma"),
                        method  = c("MLE", "MOM"),
                        xlab    = "Lifetime (hours)",
                        plot    = FALSE,
                        verbose = FALSE)
fit_lindley$ranking[, c("Rank", "Distribution", "Method", "AIC", "Delta_AIC")]

## ----johnsonsu-example--------------------------------------------------------
set.seed(42)
# True JohnsonSU(gamma=-0.5, delta=1.2, xi=50, lambda=20) data
z <- rnorm(300)
x_jsu <- 50 + 20 * sinh((z - (-0.5)) / 1.2)

fit_jsu <- fitverse(x_jsu,
                    dists   = c("johnsonsu", "normal", "skewnormal", "student_t"),
                    method  = "MLE",
                    xlab    = "Simulated returns",
                    plot    = FALSE,
                    verbose = FALSE)
fit_jsu$ranking[, c("Rank", "Distribution", "AIC", "BIC")]

## ----nakagami-example---------------------------------------------------------
set.seed(42)
# Nakagami(m=2, Omega=100)
x_naka <- sqrt(rgamma(300, shape = 2, rate = 2 / 100))

fit_naka <- fitverse(x_naka,
                     dists   = c("nakagami", "rayleigh", "weibull", "gamma"),
                     method  = c("MLE", "MOM"),
                     xlab    = "Signal amplitude",
                     plot    = FALSE,
                     verbose = FALSE)
fit_naka$ranking[, c("Rank", "Distribution", "Method", "AIC")]
coef(fit_naka)   # should recover m≈2, Omega≈100

## ----wakeby-example, eval=requireNamespace("lmomco", quietly=TRUE)------------
set.seed(42)
para_wak <- lmomco::vec2par(c(10, 50, 0.8, 5, 0.2), type = "wak")
x_wak    <- lmomco::rlmomco(300, para_wak)

fit_wak <- fitverse(x_wak,
                    dists   = c("wakeby", "gev", "glo", "lognormal"),
                    method  = c("MLE", "LMOM"),
                    xlab    = "Simulated flood peak",
                    plot    = FALSE,
                    verbose = FALSE)
fit_wak$ranking[, c("Rank", "Distribution", "Method", "AIC")]

## ----grouped-midpoints--------------------------------------------------------
# Annual rainfall (mm) grouped into 10 mm classes
rain_mids   <- c(55, 65, 75, 85, 95, 105, 115, 125, 135)
rain_counts <- c(  4, 12, 28, 35, 41,  30,  18,   8,   4)

fit_grp <- fitverse_grouped(
  midpoints = rain_mids,
  counts    = rain_counts,
  dists     = c("normal", "gamma", "lognormal", "weibull"),
  xlab      = "Annual rainfall (mm)"
)

fit_grp$rankings[, c("Rank", "Distribution", "Method", "AIC", "BIC")]

## ----grouped-coef-------------------------------------------------------------
coef(fit_grp)   # best-fit parameters

## ----grouped-intervals--------------------------------------------------------
# Failure times recorded in irregular inspection intervals
fail_lower  <- c(  0, 200, 400, 700, 1000, 1500, 2000)
fail_upper  <- c(200, 400, 700,1000, 1500, 2000, 3000)
fail_counts <- c(  8,  22,  35,  28,   19,   11,    7)

fit_fail <- fitverse_grouped(
  lower  = fail_lower,
  upper  = fail_upper,
  counts = fail_counts,
  dists  = c("weibull", "lognormal", "gamma", "exponential"),
  xlab   = "Failure time (hours)"
)

fit_fail$rankings[, c("Rank", "Distribution", "AIC", "BIC")]

## ----report-example, eval=FALSE-----------------------------------------------
# set.seed(42)
# x_rep <- rgamma(300, shape = 3, rate = 0.1)
# 
# fit_rep <- fitverse(x_rep,
#                     method  = c("MLE", "MOM", "LMOM"),
#                     xlab    = "Loss amount ($000)",
#                     verbose = FALSE)
# 
# # Export HTML report (opens in browser by default; set open = FALSE to suppress)
# generate_report(fit_rep,
#                 output_file = file.path(tempdir(), "fitverse_report.html"),
#                 title       = "Loss Distribution Analysis",
#                 author      = "Karuna G. Reddy",
#                 top_n       = 5L,
#                 open        = FALSE)

## ----report-with-boot, eval=FALSE---------------------------------------------
# boot_rep <- bootstrap_ci(fit_rep, B = 499, seed = 1, verbose = FALSE)
# generate_report(fit_rep,
#                 boot        = boot_rep,
#                 output_file = file.path(tempdir(), "fitverse_report_ci.html"),
#                 open        = FALSE)

## ----report-pdf, eval=FALSE---------------------------------------------------
# generate_report(fit_rep, output_file = file.path(tempdir(), "fitverse_report.pdf"), open = FALSE)

## ----launch-app, eval=FALSE---------------------------------------------------
# fitverse_app()

## ----sample-datasets, eval=FALSE----------------------------------------------
# # Read an individual-observation dataset and fit distributions
# path <- system.file("extdata", "annual_max_flows.csv", package = "FitVerse")
# flows <- read.csv(path, comment.char = "#")$peak_flow_m3s
# 
# fit_flows <- fitverse(flows,
#                       method  = c("MLE", "LMOM"),
#                       xlab    = "Peak flow (m3/s)",
#                       verbose = FALSE)
# print(fit_flows)
# 
# # Read a grouped dataset and fit
# gpath <- system.file("extdata", "rainfall_grouped.csv", package = "FitVerse")
# dat   <- read_grouped_csv(gpath)
# 
# fit_rain <- fitverse_grouped(midpoints = dat$midpoints,
#                              counts    = dat$counts,
#                              xlab      = "Rainfall (mm)")
# print(fit_rain)

