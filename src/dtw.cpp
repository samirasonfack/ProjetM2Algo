#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List dtw_rcpp(NumericVector x, NumericVector y) {
  int n = x.size();
  int m = y.size();

  // Cost and DTW matrices
  NumericMatrix D(n, m);
  NumericMatrix C(n, m);
  CharacterMatrix traceback(n, m);

  // Step 1: cost matrix
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      D(i,j) = std::abs(x[i] - y[j]);
    }
  }

  // Step 2: initialize first cell
  C(0,0) = D(0,0);

  // first row
  for(int j=1; j<m; j++){
    C(0,j) = C(0,j-1) + D(0,j);
    traceback(0,j) = "left";
  }

  // first column
  for(int i=1; i<n; i++){
    C(i,0) = C(i-1,0) + D(i,0);
    traceback(i,0) = "up";
  }

  // fill the rest
  for(int i=1; i<n; i++){
    for(int j=1; j<m; j++){
      double diag = C(i-1,j-1);
      double up = C(i-1,j);
      double left = C(i,j-1);
      double best = std::min({diag, up, left});
      C(i,j) = D(i,j) + best;

      if(best == diag) traceback(i,j) = "diag";
      else if(best == up) traceback(i,j) = "up";
      else traceback(i,j) = "left";
    }
  }

  // Step 3: traceback path
  int i = n-1, j = m-1;
  List path;
  while(i >= 0 && j >= 0){
    path.push_front(NumericVector::create(i+1,j+1)); // 1-based indexing
    std::string move = Rcpp::as<std::string>(traceback(i,j));
    if(move == "diag") { i--; j--; }
    else if(move == "up") { i--; }
    else if(move == "left") { j--; }
    else break;
  }

  // return results
  return List::create(
    Named("DTW_matrix") = C,
    Named("traceback_matrix") = traceback,
    Named("path") = path
  );
}
