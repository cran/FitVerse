# ============================================================
# FitVerse: Parametric Distribution Fitting and Analysis
# ============================================================
# File       : s3_methods.R
# Description: S3 generic methods for the 'fitverse' class:
#                print.fitverse   - concise one-line summary
#                summary.fitverse - full diagnostic report
#                coef.fitverse    - extract best-fit parameters
#                logLik.fitverse  - extract log-likelihood
#                AIC.fitverse     - extract AIC
#                BIC.fitverse     - extract BIC
#
#              These make fitverse objects play nicely with
#              standard R workflows and model-comparison tools.
#
# Authors:
#   Karuna G. Reddy <karuna.reddy@auckland.ac.nz>
#   University of Auckland, New Zealand
#
#   M. G. M. Khan <mgm.khan@usp.ac.fj>
#   University of the South Pacific, Fiji
#
# Package : FitVerse
# Version : 1.0-1
# Date    : 2026-06-17
# License : GPL (>= 3)
# ============================================================


# -- print.fitverse ------------------------------------------------------------

#' Print a concise summary of a fitverse object
#'
#' @param x   A \code{fitverse} object.
#' @param ... Ignored.
#' @export
#' @method print fitverse
print.fitverse <- function(x, ...) {
  fit   <- x$best_fit
  entry <- get_dist_entry(fit$dist_key)

  cat("\n")
  cat("-----------------------------------------------------\n")
  cat("  FitVerse - Best-Fit Distribution\n")
  cat("-----------------------------------------------------\n")
  cat(sprintf("  Distribution : %s\n",  entry$display_name))
  cat(sprintf("  Family       : %s\n",  entry$family))
  cat(sprintf("  Method       : %s\n",  fit$method))
  cat(sprintf("  n            : %d\n",  x$n))
  cat("\n")
  cat("  Estimated parameters:\n")
  for (i in seq_along(fit$params)) {
    cat(sprintf("    %-12s = %10.6f\n", names(fit$params)[i], fit$params[i]))
  }
  cat("\n")
  cat(sprintf("  Log-likelihood : %10.4f\n", fit$loglik))
  cat(sprintf("  AIC            : %10.4f\n", fit$aic))
  cat(sprintf("  BIC            : %10.4f\n", fit$bic))

  # GoF summary
  if (!is.null(fit$gof)) {
    cat("\n  Goodness-of-fit:\n")
    for (i in seq_len(nrow(fit$gof))) {
      row <- fit$gof[i, ]
      cat(sprintf("    %-25s  stat = %7.4f  p = %.4f  %s\n",
                  row$Test, row$Statistic, row$P_value, row$Stars))
    }
    cat("  Signif. codes: *** p<0.001  ** p<0.01  * p<0.05  . p<0.10\n")
  }

  cat("-----------------------------------------------------\n")
  cat(sprintf("  (Based on %d models fitted; ranked by %s)\n",
              nrow(x$ranking), x$criterion))
  cat("\n")

  invisible(x)
}


# -- summary.fitverse ----------------------------------------------------------

#' Full diagnostic report for a fitverse object
#'
#' Prints the best-fit details and the full ranking table of all
#' distributions that were fitted.
#'
#' @param object A \code{fitverse} object.
#' @param top_n  Number of rows to show in the ranking table. Default 10.
#' @param ...    Ignored.
#' @export
#' @method summary fitverse
summary.fitverse <- function(object, top_n = 10L, ...) {

  cat("\n")
  cat("===============================================================\n")
  cat("  FitVerse: Distribution Fitting Summary\n")
  cat("===============================================================\n")

  # Use data_summary if present (fitverse >= 1.0), else compute on-the-fly
  if (!is.null(object$data_summary)) {
    print_data_summary(object$data_summary)
  } else {
    x <- object$data
    cat(sprintf("  Sample size  : %d\n", object$n))
    cat(sprintf("  Mean         : %.4f\n", mean(x)))
    cat(sprintf("  Std Dev      : %.4f\n", stats::sd(x)))
    cat(sprintf("  Skewness     : %.4f\n", skewness_stat(x)))
    cat(sprintf("  Kurtosis     : %.4f  (excess)\n", kurtosis_stat(x)))
    cat(sprintf("  Range        : [%.4f, %.4f]\n", min(x), max(x)))
  }
  cat(sprintf("  Methods used : %s\n", paste(object$method, collapse = ", ")))
  cat(sprintf("  Criterion    : %s\n",  object$criterion))

  # Best-fit section
  cat("\n")
  print(object)   # reuse print method for best-fit block

  # Ranking table
  cat("\n  Full ranking table (top", min(top_n, nrow(object$ranking)), "models):\n\n")

  tbl <- head(object$ranking, top_n)
  # Select display columns
  show_cols <- c("Rank", "Distribution", "Method", "AIC", "Delta_AIC",
                 "BIC", "KS_p", "AD_p")
  show_cols <- show_cols[show_cols %in% names(tbl)]
  print(tbl[, show_cols], row.names = FALSE)

  cat("\n")
  cat("===============================================================\n")

  invisible(object)
}


# -- coef.fitverse ------------------------------------------------------------

#' Extract estimated parameters from the best-fit distribution
#'
#' @param object A \code{fitverse} object.
#' @param ...    Ignored.
#' @return Named numeric vector of parameter estimates.
#' @export
#' @method coef fitverse
coef.fitverse <- function(object, ...) {
  object$best_fit$params
}


# -- logLik.fitverse -----------------------------------------------------------

#' Extract the log-likelihood of the best-fit distribution
#'
#' @param object A \code{fitverse} object.
#' @param ...    Ignored.
#' @return Object of class \code{"logLik"}.
#' @export
#' @method logLik fitverse
logLik.fitverse <- function(object, ...) {
  fit <- object$best_fit
  ll  <- structure(fit$loglik,
                   df    = fit$n_params,
                   nobs  = fit$n,
                   class = "logLik")
  ll
}


# -- AIC.fitverse -------------------------------------------------------------

#' Extract the AIC of the best-fit distribution
#'
#' @param object A \code{fitverse} object.
#' @param ...    Ignored.
#' @param k      Penalty per parameter (default 2 for standard AIC).
#' @return Numeric AIC value.
#' @export
#' @method AIC fitverse
AIC.fitverse <- function(object, ..., k = 2) {
  fit <- object$best_fit
  k * fit$n_params - 2 * fit$loglik
}


# -- BIC.fitverse -------------------------------------------------------------

#' Extract the BIC of the best-fit distribution
#'
#' @param object A \code{fitverse} object.
#' @param ...    Ignored.
#' @return Numeric BIC value.
#' @export
#' @method BIC fitverse
BIC.fitverse <- function(object, ...) {
  object$best_fit$bic
}


# -- as.list.fitverse ----------------------------------------------------------

#' Convert a fitverse object to a plain list suitable for serialisation
#'
#' Strips non-serialisable elements (raw fit objects from fitdistrplus,
#' environments, etc.) and converts the result to a flat, JSON-friendly list.
#' This is the primary bridge to the Shiny web interface.
#'
#' @param x   A \code{fitverse} object.
#' @param ... Ignored.
#' @return A named list with the following top-level keys:
#'   \code{best_fit}, \code{ranking}, \code{data_summary}, \code{n},
#'   \code{xlab}, \code{criterion}, \code{method}.
#' @export
#' @method as.list fitverse
as.list.fitverse <- function(x, ...) {

  bf    <- x$best_fit
  entry <- get_dist_entry(bf$dist_key)

  # Serialise GoF table
  gof_list <- if (!is.null(bf$gof)) {
    lapply(seq_len(nrow(bf$gof)), function(i) as.list(bf$gof[i, ]))
  } else NULL

  # Serialise ranking table (drop dist_key internal column for output)
  rank_df <- x$ranking
  rank_out <- lapply(seq_len(nrow(rank_df)), function(i) as.list(rank_df[i, ]))

  list(
    best_fit = list(
      dist_key     = bf$dist_key,
      display_name = entry$display_name,
      family       = entry$family,
      method       = bf$method,
      params       = as.list(bf$params),
      loglik       = bf$loglik,
      aic          = bf$aic,
      bic          = bf$bic,
      n_params     = bf$n_params,
      converged    = bf$converged,
      gof          = gof_list
    ),
    ranking      = rank_out,
    data_summary = x$data_summary,
    n            = x$n,
    xlab         = x$xlab %||% "x",
    criterion    = x$criterion,
    method       = x$method
  )
}


# -- to_json -------------------------------------------------------------------

#' Serialise a fitverse object to a JSON string
#'
#' Converts the \code{fitverse} object to a JSON string via
#' \code{\link{as.list.fitverse}()}.  Requires the \pkg{jsonlite} package.
#'
#' @param x      A \code{fitverse} object.
#' @param pretty Logical: indent output for readability? Default \code{TRUE}.
#' @param ...    Ignored.
#' @return A JSON character string.
#' @export
to_json <- function(x, pretty = TRUE, ...) {
  UseMethod("to_json")
}

#' @export
to_json.fitverse <- function(x, pretty = TRUE, ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for to_json(). Install it with install.packages('jsonlite').",
         call. = FALSE)
  jsonlite::toJSON(as.list(x), pretty = pretty, auto_unbox = TRUE, na = "null")
}


# -- NULL coalescing helper ----------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b
