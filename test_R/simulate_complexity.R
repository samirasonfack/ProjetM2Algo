#' Simulation and Complexity Analysis for DTW
#'
#' Measures and compares the execution time of R and Rcpp DTW implementations
#' using a log-log regression to estimate the time complexity exponent.
#'
#' @param n_values Vector of sequence lengths (N) to test.
#' @param algo_r The R DTW function (e.g., dtw_sakoe_chiba).
#' @param algo_rcpp The Rcpp DTW function.
#' @param radius The fixed Sakoe-Chiba radius (W). Use a value >= 10 for better results.
#' @param plot Logical, whether to generate a log-log plot.
#'
#' @return A list containing the raw data, linear model fits, and complexity slopes.
#' @import microbenchmark ggplot2
#' @export
simulate_dtw <- function(n_values,
                         algo_r,
                         algo_rcpp,
                         radius = 10, # Increased radius for better signal-to-noise ratio
                         plot = TRUE) {

  if (!requireNamespace("microbenchmark", quietly = TRUE)) {
    stop("Package 'microbenchmark' needed for this function to work. Please install it.", call. = FALSE)
  }

  # Internal function to measure time using microbenchmark
  time_dtw <- function(n, algo, radius){
    # Generate random sequences for input
    set.seed(123)
    x <- runif(n, 0, 100)
    y <- runif(n, 0, 100)

    # Use microbenchmark for precise timing (10 runs, 5s max)
    t <- microbenchmark::microbenchmark(
      algo(x, y, radius),
      times = 10,
      unit = "s",
      control = list(maxit = 5)
    )

    # Use the median time in seconds for robustness
    t_median <- summary(t)$median

    # Avoid log(0)
    if(t_median == 0) t_median <- 1e-9
    t_median
  }

  # Measure times
  times_r <- sapply(n_values, time_dtw, algo = algo_r, radius = radius)
  times_rcpp <- sapply(n_values, time_dtw, algo = algo_rcpp, radius = radius)

  # Log-log regression: log(Time) ~ Slope * log(N) + Intercept
  # The slope of this line is the empirical complexity exponent.
  fit_r <- lm(log(times_r) ~ log(n_values))
  fit_rcpp <- lm(log(times_rcpp) ~ log(n_values))

  # Graphique optionnel
  p <- NULL
  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' needed for plotting. Skipping plot.", call. = FALSE)
    } else {
      df <- data.frame(
        n = rep(n_values, 2),
        time = c(times_r, times_rcpp),
        algo = factor(rep(c("R", "Rcpp"), each = length(n_values)))
      )

      p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = time, color = algo)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::geom_abline(intercept = coef(fit_r)[1], slope = coef(fit_r)[2], linetype = "dashed", alpha = 0.5, color = "red") +
        ggplot2::geom_abline(intercept = coef(fit_rcpp)[1], slope = coef(fit_rcpp)[2], linetype = "dashed", alpha = 0.5, color = "blue") +
        ggplot2::labs(x = "Taille n (log scale)",
                      y = "Temps (s) (log scale)",
                      color = "Algorithme",
                      title = paste("Comparaison DTW R vs Rcpp")) +
        ggplot2::theme_minimal()
    }
  }

  # Retour
  return(list(
    n_values = n_values,
    times_r = times_r,
    times_rcpp = times_rcpp,
    fit_r = fit_r,
    fit_rcpp = fit_rcpp,
    slope_r = coef(fit_r)[2],
    slope_rcpp = coef(fit_rcpp)[2],
    plot = p
  ))
}

