dtw_sakoe_chiba <- function(x, y, radius = 1) {
  n <- length(x)
  m <- length(y)

  # cumulative cost matrix (Inf outside the band)
  C <- matrix(Inf, nrow = n, ncol = m)
  traceback <- matrix("", nrow = n, ncol = m)

  # ----- INITIALISATION -----

  # first cell
  C[1,1] <- abs(x[1] - y[1])

  # first row (inside band)
  for (j in 2:m) {
    if (abs(1 - j) <= radius) {
      C[1,j] <- C[1,j-1] + abs(x[1] - y[j])
      traceback[1,j] <- "left"
    }
  }

  # first column (inside band)
  for (i in 2:n) {
    if (abs(i - 1) <= radius) {
      C[i,1] <- C[i-1,1] + abs(x[i] - y[1])
      traceback[i,1] <- "up"
    }
  }

  # ----- MAIN DP WITH BAND -----

  for (i in 2:n) {

    j_min <- max(2, i - radius)
    j_max <- min(m, i + radius)

    for (j in j_min:j_max) {

      # local cost computed on the fly (no full D matrix)
      d <- abs(x[i] - y[j])

      # only valid neighbors (they may be Inf outside band)
      choices <- c(
        C[i-1, j-1],  # diag
        C[i-1, j],    # up
        C[i, j-1]     # left
      )

      best <- min(choices)
      C[i,j] <- d + best

      # record direction
      if (best == C[i-1, j-1]) {
        traceback[i,j] <- "diag"
      } else if (best == C[i-1, j]) {
        traceback[i,j] <- "up"
      } else {
        traceback[i,j] <- "left"
      }
    }
  }

  # ----- TRACEBACK -----

  i <- n
  j <- m
  path <- list()

  while (i > 0 && j > 0) {

    path <- append(list(c(i,j)), path)

    direction <- traceback[i,j]

    if (direction == "diag") {
      i <- i - 1
      j <- j - 1
    } else if (direction == "up") {
      i <- i - 1
    } else if (direction == "left") {
      j <- j - 1
    } else {
      # reached the beginning (normally i=j=1)
      break
    }
  }

  return(list(
    DTW_matrix = C,
    traceback_matrix = traceback,
    path = path
  ))
}

