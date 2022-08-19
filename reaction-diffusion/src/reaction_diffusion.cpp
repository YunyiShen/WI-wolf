#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double subm_index(const arma::mat& X,int M, int x,int N, int y,int x0,int y0)
{
  if (x<0 || y<0 || x>=M || y>= N)
      return ((double)((x<0)+(y<0)+(x>=M)+(y>=N))) * X(x0,y0); // anywhere at boundary
  return X(x,y);
}

// This function performs a single iteration of Gray-Scott algorithm
arma::mat update_density(const arma::mat & X, 
                         const arma::mat&L,
                          const arma::mat & Landscape, 
                          double r
                                 ){
  
  double laplace;
  int m = X.n_rows;
  int n = X.n_cols;
  
  arma::mat X_new(m, n);
  
  for(int x = 0; x < m; x++) {
    for(int y = 0; y < n; y++){
      laplace = 0.0;
      for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
          laplace +=  L(i+1, j+1) * subm_index(X,m,x-i,n,y-j,x,y);
        }
      }
      X_new(x,y) = X(x,y)+r*laplace+X(x,y)*(1-X(x,y)/Landscape(x,y));
      X_new(x,y) = X_new(x,y) >=0 ? X_new(x,y) : 0;
    }
  }
  return X_new;
}

// This function performs Gray-Scott a number of times n
// [[Rcpp::export]]
arma::cube simulate_reaction_diffusion(const arma::mat& X_init, 
                              const arma::mat& Landscape, 
                              double r, 
                              int n){
  int M = X_init.n_rows;
  int N = X_init.n_cols;
  arma::cube res(M,N,n+1);
  arma::mat L = {{0,1,0},{1,-4,1},{0,1,0}}; // discrete Laplace operator
  res.slice(0) = X_init;
  for(int i = 1; i <= n; i++){
    res.slice(i) = update_density(res.slice(i-1),L,Landscape, r);
  }
  return res;
}
