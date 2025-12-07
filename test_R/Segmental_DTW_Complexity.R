library(Rcpp)
source('R/Segmental_DTW.R')
N <- c(50, 100, 200, 400, 800,1000,2000,3000,4000)
segment_length <- 10

# Storing time results
time_results <- numeric(length(N))


# Measure time for different sizes

set.seed(42)
for(i in seq_along(N)){
  n <- N[i]
  serie1 <- runif(n, 0, 100)
  serie2 <- runif(n, 0, 100)

  start_time <- Sys.time()
  result <- segmental_dtw_two_series(serie1, serie2, segment_length)
  end_time <- Sys.time()

  elapsed <- as.numeric(difftime(end_time, start_time, units="secs"))
  time_results[i] <- elapsed
  cat("n =", n, "-> Temps =", elapsed, "sec\n")
}


# log-log plot
log_n <- log(N)
log_T <- log(time_results)

install.packages("ggplot2")
library(ggplot2)

# data.frame for ggplot
df <- data.frame(
  log_n = log(N),
  log_T = log(time_results)
)

# Tracer la courbe log-log
gg <- ggplot(df, aes(x = log_n, y = log_T)) +
  geom_point(color = "blue", size = 3) +
  geom_line(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype="dashed") +
  labs(title = "Empeircal Complexity for Segmental DTW",
       x = "log(n)",
       y = "log(T)") +
  theme_minimal()

gg

# Estimate complexity slope
fit <- lm(log_T ~ log_n)
abline(fit, col="red", lwd=2)
slope <- coef(fit)[2]
print(paste("Slope (exponent) =", slope))
