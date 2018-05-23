#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat imp0clC (const arma::mat X, arma::mat cls) {
  const int n_gene = X.n_rows;
  const int n_cell = X.n_cols;
  const int n_cls = cls.n_rows;
  //  const int npow = powV.n_rows;
  // containers
  //  printf("\n-----------------------------------------------------------\n");
  //  printf("Number of genes : %d \nNumber of cells : %d \nNumber of clusterings to average : %d\n", n_gene, n_cell, n_cls);
  arma::mat tpX(n_gene, n_cell);
  tpX.fill(0);

  arma::vec indx(n_cell);
  for(int ii =0 ; ii < n_cell; ii++) {
    indx(ii) = ii;
  }

  for(int k=0 ; k < n_cls ; k ++ ) {	// each clustering setting
    //    printf("processing cluster set : %d \n", k+1 );
    arma::mat Xn = X;    
    int km = cls.row(k).max() ;
    for(int i = 0 ; i < n_gene ; i ++ ) { // for each gene
      for(int kk = 0 ; kk < km ; kk ++ ) {	// for each cluster under current clustering results

        arma::vec ind_grp = indx(find(cls.row(k) == kk+1)) ;	// cell clusters
        arma::vec ind_grp_zeros = indx(find(cls.row(k) == kk+1 && Xn.row(i) == 0)) ;
        arma::vec ind_grp_nonzeros = indx(find(cls.row(k) == kk+1 && Xn.row(i) != 0)) ;

        int tot_grp = ind_grp.size();
        int tot_grp_zeros = ind_grp_zeros.size();
        int tot_grp_nonzeros = ind_grp_nonzeros.size();

        if( tot_grp_zeros > 0 && tot_grp_zeros < tot_grp && tot_grp > 1 ) {
          
          float imp_val = 0;
          for(int j = 0 ; j < tot_grp_nonzeros ; j++ ) {
            imp_val += Xn(i, ind_grp_nonzeros(j) ); 
          }
            
          imp_val = imp_val / tot_grp;
          
          for(int j = 0 ; j < tot_grp_zeros ; j++ ) {
            Xn(i, ind_grp_zeros(j)) = imp_val ;
          }
          
        }
      }
    }
    tpX = tpX + Xn;        
  }
            
  return(tpX / n_cls);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::mat imp0clC2 (const arma::mat X, arma::mat cls) {
  const int n_gene = X.n_rows;
  const int n_cell = X.n_cols;
  const int n_cls = cls.n_rows;
  //  const int npow = powV.n_rows;
  // containers
  //  printf("\n-----------------------------------------------------------\n");
  //  printf("Number of genes : %d \nNumber of cells : %d \nNumber of clusterings to average : %d\n", n_gene, n_cell, n_cls);
  arma::mat tpX(n_gene, n_cell);
  tpX.fill(0);

  arma::vec indx(n_cell);
  for(int ii =0 ; ii < n_cell; ii++) {
    indx(ii) = ii;
  }

  for(int k=0 ; k < n_cls ; k ++ ) {
    //    printf("processing cluster set : %d \n", k+1 );
    arma::mat Xn = X;    
    int km = cls.row(k).max() ;
    for(int i = 0 ; i < n_gene ; i ++ ) {
      for(int kk = 0 ; kk < km ; kk ++ ) {

        arma::vec ind_grp = indx(find(cls.row(k) == kk+1)) ;
        arma::vec ind_grp_zeros = indx(find(cls.row(k) == kk+1 && Xn.row(i) == 0)) ;
        arma::vec ind_grp_nonzeros = indx(find(cls.row(k) == kk+1 && Xn.row(i) != 0)) ;

        int tot_grp = ind_grp.size();
        int tot_grp_zeros = ind_grp_zeros.size();
        //        int tot_grp_nonzeros = ind_grp_nonzeros.size();

        if( tot_grp_zeros > 0 && tot_grp_zeros < tot_grp && tot_grp > 4 ) {

          arma::vec imp_vec(ind_grp);
          for(int j = 0; j < tot_grp ; j++) {
            imp_vec(j) = Xn(i, ind_grp(j) );
          }
          float imp_val = median(imp_vec);
          
          //          float imp_val = 0;
          //for(int j = 0 ; j < tot_grp_nonzeros ; j++ ) {
          //  imp_val += Xn(i, ind_grp_nonzeros(j) ); 
          //}
            
          //imp_val = imp_val / tot_grp;
          
          for(int j = 0 ; j < tot_grp_zeros ; j++ ) {
            Xn(i, ind_grp_zeros(j)) = imp_val ;
          }
          
        }
      }
    }
    tpX = tpX + Xn;        
  }
            
  return(tpX / n_cls);
}

