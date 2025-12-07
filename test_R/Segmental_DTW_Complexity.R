library(Rcpp)
source('R/Segmental_DTW.R')
library(ggplot2)

# Fonction pour mesurer la complexité empirique du Segmental DTW
segmental_dtw_complexity <- function(N = c(50, 100, 200, 400, 800, 1000, 2000, 3000, 4000),
                                     segment_length = 10,
                                     plot_graph = TRUE,
                                     seed = 42) {

  set.seed(seed)

  # Stocker les temps
  time_results <- numeric(length(N))

  for(i in seq_along(N)){
    n <- N[i]
    serie1 <- runif(n, 0, 100)
    serie2 <- runif(n, 0, 100)

    start_time <- Sys.time()
    result <- segmental_dtw_two_series(serie1, serie2, segment_length)
    end_time <- Sys.time()

    elapsed <- as.numeric(difftime(end_time, start_time, units="secs"))
    time_results[i] <- elapsed
    message("n = ", n, " -> Temps = ", round(elapsed, 4), " sec")
  }

  # log-log transformation
  log_n <- log(N)
  log_T <- log(time_results)

  # Estimation de la pente (exposant)
  fit <- lm(log_T ~ log_n)
  slope <- coef(fit)[2]

  # Création du graphique
  if(plot_graph){
    df <- data.frame(log_n = log_n, log_T = log_T)

    gg <- ggplot(df, aes(x = log_n, y = log_T)) +
      geom_point(color = "blue", size = 3) +
      geom_line(color = "blue") +
      geom_smooth(method = "lm", se = FALSE, color = "red", linetype="dashed") +
      labs(title = "Complexité empirique du Segmental DTW",
           x = "log(n)",
           y = "log(T) [sec]") +
      theme_minimal()

    print(gg)
  }

  return(list(
    N = N,
    time = time_results,
    slope = slope,
    lm_fit = fit
  ))
}



