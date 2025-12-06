#' DTW algorithm testing with a case example
x <- c(1, 3, 4, 9,2)
y <- c(1, 3, 7, 8, 9)
dtw_result <- dtw(x, y)

dtw_result

res <- simulation_dtw(debut = 0, fin = 100, pas = 1)

