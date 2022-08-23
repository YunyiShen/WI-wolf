#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

double subm_index(const arma::mat& X,int M, int x,int N, int y,int x0,int y0)
{
  if (x<0 || y<0 || x>=M || y>= N)
      return 0; //((double)((x<0)+(y<0)+(x>=M)+(y>=N))) * X(x0,y0); // anywhere at boundary
  return X(x,y);
}

// This function performs a single iteration of Gray-Scott algorithm
arma::cube update_density(const arma::mat & X, 
                         const arma::mat&L,
                          const arma::mat & Landscape, 
                          double r
                                 ){
  
  double laplace;
  int m = X.n_rows;
  int n = X.n_cols;
  
  arma::cube X_new(m, n, 3);
  
  for(int x = 0; x < m; x++) {
    for(int y = 0; y < n; y++){
      laplace = 0.0;
      for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
          laplace +=  L(i+1, j+1) * subm_index(X,m,x+i,n,y+j,x,y);
        }
      }
      X_new(x,y,0) = X(x,y)+r*laplace+X(x,y)*(1-X(x,y)/Landscape(x,y));
      X_new(x,y,0) = X_new(x,y,0) >=0 ? X_new(x,y,0) : 0;
      X_new(x,y,1) = r*laplace; //(abs(r*laplace) >  0.95*abs(X(x,y)*(1-X(x,y)/Landscape(x,y)))) && (abs(r*laplace) <  1.05*abs(X(x,y)*(1-X(x,y)/Landscape(x,y)))) ? 1:0;
      X_new(x,y,2) = X(x,y)*(1-X(x,y)/Landscape(x,y));
    }
  }
  return X_new;
}

// This function performs Gray-Scott a number of times n
// [[Rcpp::export]]
Rcpp::List simulate_reaction_diffusion(const arma::mat& X_init, 
                              const arma::mat& Landscape, 
                              double r, 
                              int n){
  int M = X_init.n_rows;
  int N = X_init.n_cols;
  arma::cube res(M,N,n+1);
  arma::cube frontier(M,N,n+1, arma::fill::zeros);
  arma::cube growth(M,N,n+1, arma::fill::zeros);
  arma::cube tmp;
  //arma::mat L = {{0.25,.5,0.25},{.5,-3,.5},{0.25,.5,0.25}}; // discrete Laplace operator
  arma::mat L = {{0,1,0.},{1,-4,1},{0,1,0}}; // discrete Laplace operator
  
  res.slice(0) = X_init;
  for(int i = 1; i <= n; i++){
    tmp= update_density(res.slice(i-1),L,Landscape, r);
    res.slice(i) = tmp.slice(0);
    frontier.slice(i) = tmp.slice(1);
    growth.slice(i) = tmp.slice(2);
  }
  Rcpp::List result;
  result["density"]=res;
  result["frontier"] = frontier;
  result["growth"] = growth;
  return result;
}
