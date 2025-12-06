#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List dtw_sakoe_chiba_rcpp(NumericVector x, NumericVector y, int radius = 1) {
  int n = x.size();
  int m = y.size();

  // cumulative cost matrix
  NumericMatrix C(n, m);
  CharacterMatrix traceback(n, m);

  // initialize all to +Inf
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      C(i,j) = R_PosInf;
      traceback(i,j) = "";
    }
  }

  // first cell
  C(0,0) = std::abs(x[0] - y[0]);

  // first row
  for(int j=1; j<m; j++){
    if(std::abs(0 - j) <= radius){
      C(0,j) = C(0,j-1) + std::abs(x[0] - y[j]);
      traceback(0,j) = "left";
    }
  }

  // first column
  for(int i=1; i<n; i++){
    if(std::abs(i - 0) <= radius){
      C(i,0) = C(i-1,0) + std::abs(x[i] - y[0]);
      traceback(i,0) = "up";
    }
  }

  // main DP
  for(int i=1; i<n; i++){
    int j_min = std::max(1, i - radius);
    int j_max = std::min(m-1, i + radius);

    for(int j=j_min; j<=j_max; j++){
      double d = std::abs(x[i] - y[j]);
      double best = C(i-1,j-1);
      std::string dir = "diag";

      if(C(i-1,j) < best){
        best = C(i-1,j);
        dir = "up";
      }
      if(C(i,j-1) < best){
        best = C(i,j-1);
        dir = "left";
      }

      C(i,j) = d + best;
      traceback(i,j) = dir;
    }
  }

  // traceback
  int i = n-1;
  int j = m-1;
  std::vector< std::vector<int> > path;

  while(i >=0 && j >=0){
    path.push_back({i+1, j+1}); // 1-based indices for R
    std::string dir = Rcpp::as<std::string>(traceback(i,j));

    if(dir == "diag"){
      i--; j--;
    } else if(dir == "up"){
      i--;
    } else if(dir == "left"){
      j--;
    } else {
      break;
    }
  }

  // reverse path
  std::reverse(path.begin(), path.end());

  // convert to R list
  List path_list(path.size());
  for(size_t k=0; k<path.size(); k++){
    path_list[k] = IntegerVector(path[k].begin(), path[k].end());
  }

  return List::create(
    Named("DTW_matrix") = C,
    Named("traceback_matrix") = traceback,
    Named("path") = path_list
  );
}
