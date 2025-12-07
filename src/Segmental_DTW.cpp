#include <Rcpp.h>
#include <algorithm> // pour std::reverse
using namespace Rcpp;

// ------------------ DTW ------------------
// [[Rcpp::export]]
List dtw_cpp(NumericVector x, NumericVector y) {
  int n = x.size();
  int m = y.size();

  NumericMatrix D(n, m);
  NumericMatrix C(n, m);
  CharacterMatrix traceback(n, m);

  // cost matrix
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      D(i,j) = std::abs(x[i] - y[j]);
      C(i,j) = R_PosInf;
      traceback(i,j) = "";
    }
  }

  // initialization
  C(0,0) = D(0,0);
  for(int j=1;j<m;j++){ C(0,j)=C(0,j-1)+D(0,j); traceback(0,j)="left"; }
  for(int i=1;i<n;i++){ C(i,0)=C(i-1,0)+D(i,0); traceback(i,0)="up"; }

  // filling
  for(int i=1;i<n;i++){
    for(int j=1;j<m;j++){
      double diag=C(i-1,j-1), up=C(i-1,j), left=C(i,j-1);
      double best = std::min({diag, up, left});
      C(i,j)=D(i,j)+best;
      if(best==diag) traceback(i,j)="diag";
      else if(best==up) traceback(i,j)="up";
      else traceback(i,j)="left";
    }
  }

  // path
  IntegerMatrix path(n+m, 2); // max size
  int pi = n+m-1;
  int i=n-1,j=m-1;
  while(i>=0 && j>=0){
    path(pi,0)=i+1;
    path(pi,1)=j+1;
    pi--;
    if(traceback(i,j)=="diag"){ i--; j--; }
    else if(traceback(i,j)=="up"){ i--; }
    else if(traceback(i,j)=="left"){ j--; }
    else break;
  }

  IntegerMatrix final_path(n+m-pi-1, 2);
  for(int k=0;k<n+m-pi-1;k++){
    final_path(k,0)=path(pi+k+1,0);
    final_path(k,1)=path(pi+k+1,1);
  }

  return List::create(
    _["DTW_matrix"]=C,
    _["traceback_matrix"]=traceback,
    _["path"]=final_path
  );
}

// ------------------ Segmental DTW ------------------
// [[Rcpp::export]]
List segmental_dtw_cpp(NumericVector serie1, NumericVector serie2, int segment_length) {
  int n = serie1.size();
  int m = serie2.size();
  int n_seg = (m + segment_length - 1)/segment_length;

  // cumulative cost matrix
  NumericMatrix cum_cost(n_seg, n);
  for(int i=0;i<n_seg;i++)
    for(int j=0;j<n;j++)
      cum_cost(i,j)=R_PosInf;

  std::vector<IntegerMatrix> dtw_paths; // vecteur vide
  std::vector<NumericMatrix> dtw_matrices;

  // calcul DTW par segment
  for(int s=0;s<n_seg;s++){
    int start = s*segment_length;
    int end = std::min(start+segment_length-1, m-1);
    NumericVector subserie = serie2[Range(start,end)];
    List res = dtw_cpp(serie1, subserie);

    dtw_paths.push_back(res["path"]); // ✅ push_back au lieu de dtw_paths[s]=
    dtw_matrices.push_back(res["DTW_matrix"]);

    // cumulative cost
    NumericMatrix C = dtw_matrices.back();
    for(int j=0;j<n;j++){
      double min_prev = 0;
      if(s>0){
        min_prev = R_PosInf;
        for(int k=j;k<n;k++)
          if(cum_cost(s-1,k)<min_prev) min_prev=cum_cost(s-1,k);
      }
      cum_cost(s,j) = C(n-1,j)+min_prev;
    }
  }

  // DTW distance finale
  double dtw_distance = R_PosInf;
  int last = n_seg-1;
  for(int j=0;j<n;j++)
    if(cum_cost(last,j)<dtw_distance) dtw_distance=cum_cost(last,j);

    // path segmental
    IntegerMatrix path_seg(n_seg,2);
    int jmin = 0;
    double min_val = cum_cost(last,0);
    for(int j=0;j<n;j++) if(cum_cost(last,j)<min_val){ min_val=cum_cost(last,j); jmin=j; }
    int j=jmin;
    for(int i=last;i>=0;i--){
      path_seg(i,0)=i+1;
      path_seg(i,1)=j+1;
      if(i>0){
        min_val = R_PosInf;
        int new_j = j;
        for(int k=j;k<n;k++){
          if(cum_cost(i-1,k)<min_val){ min_val=cum_cost(i-1,k); new_j=k; }
        }
        j=new_j;
      }
    }

    return List::create(
      _["cum_cost"]=cum_cost,
      _["dtw_distance"]=dtw_distance,
      _["path_seg"]=path_seg
    );
}

