#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List dtw_sakoe_chiba_rcpp(NumericVector x, NumericVector y, int radius = 1) {
  int n = x.size();
  int m = y.size();
  if(n == 0 || m == 0) stop("Input series cannot be empty.");

  int w = 2*radius + 1; // largeur de la bande

  // Matrice des coûts (n x w)
  NumericMatrix C(n, w);
  for(int i=0;i<n;i++)
    for(int j=0;j<w;j++)
      C(i,j) = R_PosInf;

  // Matrice de traceback (n x w)
  IntegerMatrix traceback(n, w);

  // Initialisation de la première cellule
  C(0, radius) = std::abs(x[0] - y[0]);
  traceback(0, radius) = 0;

  // Boucle principale
  for(int i=0;i<n;i++){
    int j_min = std::max(0, i - radius);
    int j_max = std::min(m-1, i + radius);

    for(int j=j_min;j<=j_max;j++){
      int jb = j - (i - radius); // index dans la bande
      if(i==0 && j==0) continue;

      double cost = std::abs(x[i]-y[j]);
      double c_diag = R_PosInf, c_up=R_PosInf, c_left=R_PosInf;

      // Diagonal
      if(i>0 && j>0){
        int jb_diag = j -1 - ((i-1) - radius);
        if(jb_diag>=0 && jb_diag<w) c_diag = C(i-1, jb_diag);
      }
      // Up
      if(i>0){
        int jb_up = j - (i-1 - radius);
        if(jb_up>=0 && jb_up<w) c_up = C(i-1, jb_up);
      }
      // Left
      if(j>0){
        int jb_left = j-1 - (i - radius);
        if(jb_left>=0 && jb_left<w) c_left = C(i, jb_left);
      }

      double best = c_diag;
      int dir = 0; // diag
      if(c_up < best){ best = c_up; dir = 1; }
      if(c_left < best){ best = c_left; dir = 2; }

      C(i,jb) = cost + best;
      traceback(i,jb) = dir;
    }
  }

  // Traceback
  std::vector< std::vector<int> > path;
  int i = n-1;
  int j = m-1;
  while(i>=0 && j>=0){
    int jb = j - (i - radius);
    if(jb<0 || jb>=w) break;
    path.push_back({i+1,j+1});
    int dir = traceback(i,jb);
    if(dir==0){ i--; j--; }
    else if(dir==1){ i--; }
    else if(dir==2){ j--; }
    else break;
  }
  std::reverse(path.begin(), path.end());

  // Conversion vers R
  List path_list(path.size());
  for(size_t k=0;k<path.size();k++)
    path_list[k] = IntegerVector(path[k].begin(), path[k].end());

  // Distance finale
  int jb_end = m-1 - (n-1 - radius);
  double dtw_dist = (jb_end>=0 && jb_end<w) ? C(n-1,jb_end) : R_PosInf;

  return List::create(
    Named("DTW_distance") = dtw_dist,
    Named("DTW_matrix_band") = C,
    Named("traceback_matrix_band") = traceback,
    Named("path") = path_list
  );
}

