#' Dynamic Time Warping (DTW) with Traceback
#'
#' Computes the Dynamic Time Warping distance between two numeric sequences
#' and returns the DTW matrix, the traceback matrix, and the optimal alignment path.
#' In this file, we will applicate this to segmented time series in order to classify them using K-means clustering 
#' which uses DTW as a distance measure. 

#
#' Algorithm: 
#'  1. Segmentation of the time series using sliding windows.
#'  2. For each segment, compute the DTW distance to all other segments.
#'  3. Normalize the segments to ensure comparability.
#'  4. Calculate the DTW distance matrix for all segments.

# Segmentation 
sliding_windows <- function(serie, segment_size, step){
    n <- length(serie)
    segments <- list()
    idx <- 1
    for (i in seq(1, n - segment_size + 1, by = step)){
        segments[[idx]] <- serie[i:(i + segment_size - 1)]
        idx <- idx + 1
    }
    return ( segments )
} 

# Normalization function
normalize_segment <- function(segments) {
    normalize_segments <- list()
    for ( segment in segments){
        if (length(segment) == 0) {
            next
        }
        if (max(segment) == min(segment)) {
            segment[] <- 0  # If all values are the same, set to zero
        } else {
            segment[] <- (segment - min(segment)) / (max(segment) - min(segment))
        normalize_segments <- c(normalize_segments, list(segment))
        }
    }
    return (normalize_segments)
}

# DTW function

dtw_segment <- function(normalized_segments_1, normalized_segments_2){
    dtw_results <- list()
    
    #' #___dtw_results____ : 
    #' will store the dtw_result for each pair of segments from the two series.
    #' #___dtw_result___ contains :
    #' DTW matrix, traceback matrix, and optimal path for each segment pair.
    
    for (i in seq_along(normalized_segments_1)) {

        for (j in seq_along(normalized_segments_2)) {
            dtw_result <- dtw(normalized_segments_1[[i]], normalized_segments_2[[j]])
            dtw_results[[paste0("Segment_", i, "_", j)]] <- dtw_result
        }
        
    }

    return(dtw_results)
}
