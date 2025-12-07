#' Segmental Dynamic Time Warping (DTW) Implementation in R
#' It's a A Parallelizable Alternative to Dynamic Time Warping for Time Series
#' This implementation computes the Segmental DTW between two time series using a segment-level approach.
#' The algorithm involves segmenting the time series, computing DTW for each segment,
#' and then aggregating the results to obtain the final DTW distance and path.
#' pipelines the process of algorithm : segmentation -> DTW computation using dtw naive ->
#' segment-level cost list -> cumulative cost ->
#' backtrace segment-level -> reconstruct frame-level path




library(Rcpp)
# Import the dtw function from dtw.R
source("R/dtw.R")

#' @export
# ________Segmentation of the time series 1 into n=(segment_length) segments
#' Input : series, segment length
#' Output : list of segments (start, end) positions
segment_series <- function(series, segment_length){
  n <- length(series)
  segments <- list()
  start <- 1
  while(start <= n){
    end <- min(start + segment_length - 1, n)
    segments <- append(segments, list(c(start, end)))
    start <- end + 1
  }
  return(segments)
}

#' @export
# _______Calculate the DTW between (each segment of serie 2) and (serie 1) using dtw function
#' Input : serie1, serie2, list of segments (start, end) positions
#' Output : list of dtw results for each segment
compute_segment_dtw <- function(serie1, serie2, segments){
  dtw_list <- lapply(segments, function(seg){
    subserie <- serie2[seg[1]:seg[2]]
    dtw(serie1, subserie)
  })
  return(dtw_list)
}

#' @export
# _____ Segment-level cost list
# input : list of dtw results for each segment
# output : list of segment-level cost vectors

compute_segment_level_cost <- function(dtw_segments){
  n_seg <- length(dtw_segments)
  segment_cost_list <- vector("list", n_seg)

  for(i in 1:n_seg){
    dtw_i <- dtw_segments[[i]]
    last_row <- dtw_i$DTW_matrix[nrow(dtw_i$DTW_matrix), ] # Last row of DTW matrix between serie1 and segment i
    segment_cost_list[[i]] <- last_row
  }

  return(segment_cost_list)
}

#' @export
# Segmental DTW for the two series aand the frame-level path reconstruction
# input : serie1, serie2, segment length
# output : list with segments, dtw_segments, cum_cost, segment_path, final_path, dtw_distance
segmental_dtw_two_series <- function(serie1, serie2, segment_length){

  # Segmentation
  segments <- segment_series(serie2, segment_length)

  # DTW per segment
  dtw_segments <- compute_segment_dtw(serie1, serie2, segments)

  # Segment-level cost list
  seg_cost_list <- compute_segment_level_cost(dtw_segments)

  # Calcul cumulative cost of segment-level
  n_seg <- length(seg_cost_list)
  n_len <- length(serie1)
  cum_cost <- matrix(Inf, nrow=n_seg, ncol=n_len)
  cum_cost[1, 1:length(seg_cost_list[[1]])] <- seg_cost_list[[1]]

  for(i in 2:n_seg){
    seg_cost <- seg_cost_list[[i]]
    for(j in 1:length(seg_cost)){
      min_prev <- min(cum_cost[i-1, j:n_len])
      cum_cost[i,j] <- seg_cost[j] + min_prev
    }
  }

  # Final Distance DTW
  dtw_distance <- min(cum_cost[n_seg, ])

  # Backtrace segment-level
  path_seg <- list()
  i <- n_seg
  j <- which.min(cum_cost[n_seg, ])
  while(i > 0){
    path_seg <- append(list(c(i,j)), path_seg)
    if(i==1) break
    j <- which.min(cum_cost[i-1, j:n_len])
    i <- i-1
  }

  # fianl optimal path of frame-level
  final_path <- list()
  for(k in 1:length(path_seg)){
    seg_idx <- path_seg[[k]][1]
    serie1_pos <- path_seg[[k]][2]
    dtw_i <- dtw_segments[[seg_idx]]

    # Extraire le chemin jusqu'à la position serie1_pos
    path_vec <- dtw_i$path
    path_sub <- lapply(path_vec, function(p){
      if(p[1] <= serie1_pos){
        c(p[1], p[2] + segments[[seg_idx]][1] - 1)
      } else {
        NULL
      }
    })
    path_sub <- do.call(rbind, path_sub)
    if(!is.null(path_sub)){
      final_path <- append(final_path, list(path_sub))
    }
  }

  final_path_matrix <- do.call(rbind, final_path)

  return(list(
    segments = segments,
    dtw_segments = dtw_segments,
    cum_cost = cum_cost,
    segment_path = path_seg,
    final_path = final_path_matrix,
    dtw_distance = dtw_distance
  ))
}

