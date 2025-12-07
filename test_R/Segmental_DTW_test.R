library(Rcpp)
source('R/Segmental_DTW.R')

serie1 <- c(1, 3, 4, 9, 8)
serie2 <- c(1, 2, 3, 4, 7, 8, 9, 10)

result <- segmental_dtw_two_series(serie1, serie2, segment_length = 2)

cat("Distance DTW finale segmentale :", result$dtw_distance, "\n")
cat("Chemin frame-level final :\n")
print(result$final_path)
