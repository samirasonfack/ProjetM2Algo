
#' k-NN Classification for Time Series Using DTW (Rcpp Implementation)
#'
#' This function performs time-series classification using a custom
#' Dynamic Time Warping (DTW) distance implemented via Rcpp.
#'
#' It computes the DTW distance between one test series and each
#' training series, selects the _k_ nearest neighbors, and returns the
#' predicted class based on majority vote.
#'
#' @param train_matrix A numeric matrix of dimension \eqn{N × T} where each row
#'   is a time series of length \eqn{T}.
#' @param train_labels A factor or character vector of length \eqn{N}
#'   containing the class labels associated with each training series.
#' @param test_series A numeric vector of length \eqn{T} representing the
#'   time series to classify.
#' @param k Integer. Number of neighbors to use for the k-NN vote
#'   (default: 5).
#' @param dtw_func A function taking two numeric vectors and
#'   \code{radius}, returning a list containing a \code{distance} field.
#'   Default is \code{dtw_sakoe_chiba_rcpp}.
#' @param radius Integer. Sakoe-Chiba window radius used by the DTW algorithm.
#'
#' @return A list containing:
#' \describe{
#'   \item{predicted_label}{The predicted class (character).}
#'   \item{neighbor_labels}{Labels of the \eqn{k} nearest neighbors.}
#'   \item{neighbor_distances}{DTW distances of the nearest neighbors.}
#'   \item{all_distances}{Vector of DTW distances to all training series.}
#'   \item{neighbor_index}{Indices of the nearest neighbors.}
#' }
#'
#' @details
#' This function is designed for high-performance classification with
#' DTW computed in C++ for speed.
#'
#' It is compatible with multivariate extensions if DTW is computed
#' accordingly (e.g., multidimensional DTW).
#'
#' @examples
#' \dontrun{
#'   result <- knn_dtw_rcpp(
#'     train_matrix = M,
#'     train_labels = y,
#'     test_series  = M[10, ],
#'     k = 3,
#'     radius = 30
#'   )
#'
#'   print(result$predicted_label)
#'   print(result$neighbor_distances)
#' }
#'
#' @references
#' 1. Berndt, D. J., & Clifford, J. (1994).
#'    *Using Dynamic Time Warping to Find Patterns in Time Series*.
#'    Proceedings of the AAAI Workshop on Knowledge Discovery.
#'
#' 2. Sakoe, H., & Chiba, S. (1978).
#'    *Dynamic Programming Algorithm Optimization for Spoken Word Recognition*.
#'    IEEE Transactions on Acoustics, Speech, and Signal Processing.
#'
#' 3. Tavenard, R. et al. (2020).
#'    *Tslearn: A Machine Learning Toolkit for Time Series Data*.
#'    Journal of Machine Learning Research.
#'
#' @export
knn_dtw_rcpp <- function(train_matrix,
                         train_labels,
                         test_series,
                         k = 5,
                         dtw_func = dtw_sakoe_chiba_rcpp,
                         radius = 30) {

  # 1. Checks
  if (!is.matrix(train_matrix)) stop("train_matrix must be a matrix")
  if (!is.vector(test_series)) stop("test_series must be a numeric vector")
  if (length(train_labels) != nrow(train_matrix))
    stop("train_labels length does not match number of rows in train_matrix")

  dtw_distances <- sapply(1:nrow(train_matrix), function(i) {
    train_series_vec <- as.numeric(train_matrix[i, ])
    test_series_vec  <- as.numeric(test_series)

    result <- dtw_func(test_series_vec, train_series_vec, radius = radius)


    # EXTRACTION SÉCURISÉE DU SCALAIRE
    return(as.numeric(unlist(result$dtw_dist))[1])
  })


  # 3. Select nearest neighbors
  neighbor_idx <- order(dtw_distances)[1:k]
  neighbor_labels <- train_labels[neighbor_idx]

  # 4. Majority vote
  predicted_label <- names(sort(table(neighbor_labels), decreasing = TRUE))[1]

  # 5. Output structure
  return(list(
    predicted_label = predicted_label,
    neighbor_labels = neighbor_labels,
    neighbor_distances = dtw_distances[neighbor_idx],
    all_distances = dtw_distances,
    neighbor_index = neighbor_idx
  ))
}
