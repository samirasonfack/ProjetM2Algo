#' Dynamic Time Warping (DTW) with Traceback
#'
#' Computes the Dynamic Time Warping distance between two numeric sequences
#' and returns the DTW matrix, the traceback matrix, and the optimal alignment path.
#'
#' @param x Numeric vector. The first sequence.
#' @param y Numeric vector. The second sequence.
#'
#' @return A list containing:
#' \describe{
#'   \item{DTW_matrix}{The cumulative DTW cost matrix.}
#'   \item{traceback_matrix}{A matrix storing the direction of the optimal path ("diag", "up", "left").}
#'   \item{path}{A list of coordinates representing the optimal alignment path.}
#' }
#' @examples
#' x <- c(1, 3, 4, 9)
#' y <- c(1, 3, 7, 8, 9)
#' dtw_result <- dtw(x, y)
#'
#' @export
dtw <- function(x, y) {
  n <- length(x)
  m <- length(y)

  # Cost matrix
  D <- matrix(0, nrow = n, ncol = m)
  for(i in 1:n){
    for(j in 1:m){
      D[i,j] <- abs(x[i] - y[j])  # using absolute distance
    }
  }

  # DTW matrix of cumulative cost for a given optimal path and traceback
  C <- matrix(Inf, nrow = n, ncol = m)
  traceback <- matrix("", nrow = n, ncol = m)

  # Initialize first cell
  C[1,1] <- D[1,1]

  # First row
  for(j in 2:m){
    C[1,j] <- C[1,j-1] + D[1,j]
    traceback[1,j] <- "left"
  }

  # First column
  for(i in 2:n){
    C[i,1] <- C[i-1,1] + D[i,1]
    traceback[i,1] <- "up"
  }

  # Fill the rest of the DTW matrix
  for(i in 2:n){
    for(j in 2:m){
      choices <- c(C[i-1,j-1], C[i-1,j], C[i,j-1])
      best <- min(choices)
      C[i,j] <- D[i,j] + best

      # Traceback
      if(best == C[i-1,j-1]){
        traceback[i,j] <- "diag"
      } else if(best == C[i-1,j]){
        traceback[i,j] <- "up"
      } else {
        traceback[i,j] <- "left"
      }
    }
  }

  # Traceback path
  i <- n
  j <- m
  path <- list()
  while(i > 0 & j > 0){
    path <- append(list(c(i,j)), path)
    if(traceback[i,j] == "diag"){
      i <- i-1
      j <- j-1
    } else if(traceback[i,j] == "up"){
      i <- i-1
    } else if(traceback[i,j] == "left"){
      j <- j-1
    } else {
      break
    }
  }

  # Return results as a list
  list(
    DTW_matrix = C,
    traceback_matrix = traceback,
    path = path
  )
}

