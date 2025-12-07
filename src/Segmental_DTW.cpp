#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// ------------------ DTW simple pour deux séries ------------------
// [[Rcpp::export]]
List dtw_cpp(NumericVector x, NumericVector y) {
  int n = x.size();
  int m = y.size();

  NumericMatrix D(n, m);
  NumericMatrix C(n, m);
  CharacterMatrix traceback(n, m);

  // matrice de coût
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      D(i,j) = std::abs(x[i] - y[j]);
      C(i,j) = R_PosInf;
      traceback(i,j) = "";
    }
  }

  // initialisation
  C(0,0) = D(0,0);

  for(int j=1;j<m;j++){ C(0,j)=C(0,j-1)+D(0,j); traceback(0,j)="left"; }
  for(int i=1;i<n;i++){ C(i,0)=C(i-1,0)+D(i,0); traceback(i,0)="up"; }

  // remplissage
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

  // path frame-level
  int i=n-1,j=m-1;
  List path;
  while(i>=0 && j>=0){
    path.push_back(IntegerVector::create(i+1,j+1));
    if(traceback(i,j)=="diag"){ i--; j--; }
    else if(traceback(i,j)=="up"){ i--; }
    else if(traceback(i,j)=="left"){ j--; }
    else break;
  }
  std::reverse(path.begin(), path.end());

  return List::create(
    _["DTW_matrix"]=C,
    _["traceback_matrix"]=traceback,
    _["path"]=path
  );
}

// ------------------ Segmental DTW ------------------
// [[Rcpp::export]]
List segmental_dtw_cpp(NumericVector serie1, NumericVector serie2, int segment_length) {

  int n = serie1.size();
  int m = serie2.size();
  std::vector< std::pair<int,int> > segments;

  // segmentation
  for(int start=0; start<m; start+=segment_length){
    int end = std::min(start+segment_length-1, m-1);
    segments.push_back(std::make_pair(start,end));
  }

  // DTW par segment
  std::vector<List> dtw_segments;
  for(int s=0;s<segments.size();s++){
    int start = segments[s].first;
    int end = segments[s].second;
    NumericVector subserie = serie2[Range(start,end)];
    dtw_segments.push_back(dtw_cpp(serie1, subserie));
  }

  // Segment-level cost matrix
  int n_seg = dtw_segments.size();
  NumericMatrix cum_cost(n_seg, n);
  for(int i=0;i<n_seg;i++)
    for(int j=0;j<n;j++)
      cum_cost(i,j)=R_PosInf;

  // première ligne
  NumericVector first_seg = dtw_segments[0]["DTW_matrix"];
  NumericMatrix first_mat = dtw_segments[0]["DTW_matrix"];
  for(int j=0;j<n;j++)
    cum_cost(0,j) = first_mat(n-1,j);

  // cumulative cost
  for(int i=1;i<n_seg;i++){
    NumericMatrix seg_mat = dtw_segments[i]["DTW_matrix"];
    for(int j=0;j<n;j++){
      double min_prev = R_PosInf;
      for(int k=j;k<n;k++)
        if(cum_cost(i-1,k)<min_prev) min_prev=cum_cost(i-1,k);
      cum_cost(i,j) = seg_mat(n-1,j)+min_prev;
    }
  }

  // DTW distance finale
  double dtw_distance = R_PosInf;
  for(int j=0;j<n;j++)
    if(cum_cost(n_seg-1,j)<dtw_distance) dtw_distance=cum_cost(n_seg-1,j);

  // Backtrace segment-level
  std::vector< std::pair<int,int> > path_seg;
  int i = n_seg-1;
  int j = 0;
  double min_val = cum_cost(i,0);
  for(int k=0;k<n;k++) if(cum_cost(i,k)<min_val){ min_val=cum_cost(i,k); j=k; }

  while(i>=0){
    path_seg.push_back(std::make_pair(i,j));
    if(i==0) break;
    min_val = R_PosInf;
    int new_j = j;
    for(int k=j;k<n;k++){
      if(cum_cost(i-1,k)<min_val){ min_val=cum_cost(i-1,k); new_j=k; }
    }
    j=new_j;
    i--;
  }
  std::reverse(path_seg.begin(), path_seg.end());

  return List::create(
    _["cum_cost"]=cum_cost,
    _["dtw_distance"]=dtw_distance,
    _["path_seg"]=path_seg
  );
}
