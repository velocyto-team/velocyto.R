// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// given a distance matrix, form balanced kNN assignment using greedy approach
// [[Rcpp::export]]
arma::sp_mat  balanced_knn(const arma::mat& d,int k,int maxl,bool returnDistanceValues=false, int nthreads=1) {
  arma::uvec l(d.n_cols,arma::fill::zeros); // reciprocal neighbor count
  arma::umat dsi(d.n_rows,d.n_cols); // will store column sort indices of the d matrix
#pragma omp parallel for shared(l) num_threads(nthreads)
  // run regular knn, record number of reciprocal neighbors
  for(int i=0;i<d.n_cols;i++) {
    arma::uvec si=sort_index(d.col(i));
    dsi.col(i)=si;
    l.elem(si.subvec(0,k-1))+=1;
  }

  arma::uvec lsi=sort_index(l,"descend"); // greedy order of column considerations
  l.zeros(); // reset so that l can be used to keep track of reciprocal counts in the kNN being constucted
  // greedy knn (construct row index vector)
  arma::uvec rowind(k*d.n_cols);
  arma::vec vals(k*d.n_cols,arma::fill::zeros); // values
  vals+=1.0; // default to ones

  // run regular knn, record number of reciprocal neighbors
  for(int i=0;i<d.n_cols;i++) {
    int el=lsi[i]; // element for which the neighbors are being found
    arma::uvec si=dsi.col(el);
    int p=0;
    int j;
    for(j=0;j<dsi.n_rows && p<k;j++) {
      // consider si[j] for p-th neighbor position
      int m=si[j];
      if(el==m) { continue; } // dont record or count self-relationships
      if(l[m]>= maxl) { continue; } // element has already maxed out its neighbors
      rowind[el*k+p]=m; l[m]++; p++; // record neighbor
      if(returnDistanceValues) { vals[el*k+p]=d(m,el); }
    }
    if(j==dsi.n_rows && p<k) {
      rowind[el*k+p]=el; p++; // fill in last element(s) with self-identity if there were not enough spare elements
      while(p<k) {
	std::stringstream es;
	es<<"unable to find unfilled neighbors for i="<<i<<" el="<<el<<" p="<<p;
	Rf_warning(es.str().c_str());
	rowind[el*k+p]=el; p++;
      }
    }
  }


  // construct a column pointer index - every coulmn has only k entries, so its easy

  arma::uvec colptr(d.n_cols+1);
  int cc=0;
  for(int i=0;i<colptr.n_elem;i++) { colptr[i]=cc; cc+=k; }
  // make return matrix
  arma::sp_mat knn(rowind,colptr,vals,d.n_rows,d.n_cols); // kNN matrix
  return(knn);
}

// quick matrix correlation function
// [[Rcpp::export]]
arma::mat arma_mat_cor(const arma::mat& m) {
  return(cor(m));
}


// calculates correlation matrix between e-e_i and d_i
// [[Rcpp::export]]
arma::mat  colDeltaCor(const arma::mat& e, const arma::mat& d, int nthreads=1) {
  arma::mat rm(e.n_cols,e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for(int i=0;i<e.n_cols;i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    rm.col(i)=cor(t,d.col(i));
  }

  return(rm);
}

// [[Rcpp::export]]
arma::mat  colDeltaCorSqrt(const arma::mat& e, const arma::mat& d, int nthreads=1) {
  arma::mat rm(e.n_cols,e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for(int i=0;i<e.n_cols;i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    t=sqrt(abs(t)) % sign(t);
    rm.col(i)=cor(t,d.col(i));
  }

  return(rm);
}

// [[Rcpp::export]]
arma::mat  colDeltaCorLog10(const arma::mat& e, const arma::mat& d, double pseudocount=1.0, int nthreads=1) {
  arma::mat rm(e.n_cols,e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for(int i=0;i<e.n_cols;i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    t=log10(abs(t)+pseudocount) % sign(t);
    rm.col(i)=cor(t,d.col(i));
  }

  return(rm);
}


// eucledian distance of m columns to p columns
// [[Rcpp::export]]
arma::mat  colEuclid(const arma::mat& e, const arma::mat& p, int nthreads=1) {
  arma::mat rm(e.n_cols,e.n_cols);

#pragma omp parallel for shared(rm) num_threads(nthreads)
  for(int i=0;i<e.n_cols;i++) {
    arma::mat t(e); t.each_col() -= p.col(i);
    t%=t;
    arma::rowvec v=sqrt(sum(t,0));
    rm.col(i)=trans(v);
  }
  return(rm);
}


// calculates arrows on a provided embedding (emb) given cell-cell transition probabilities (tp)
// returns arrow deltas
// [[Rcpp::export]]
arma::mat  embArrows(const arma::mat& emb, const arma::mat& tp, double arrowScale=1.0, int nthreads=1) {

  arma::mat dm(emb.n_cols,emb.n_rows);
  arma::mat tpb(ceil(tp)); // binarized version to give equal weights to each cell in the neighborhood
  tpb.each_row() /= sum(tpb,0); // normalize
  arma::colvec zv(emb.n_cols,arma::fill::zeros);
  arma::mat temb=trans(emb);
#pragma omp parallel for shared(dm) num_threads(nthreads)
  for(int i=0;i<emb.n_rows;i++) {
    arma::mat di(temb);
    di.each_col()-=di.col(i); // coordinate difference to every cell
    di=arma::normalise(di,2,0) * arrowScale; // normalized, scaled
    di.col(i)=zv; // no distance to itself
    arma::vec ds=di * tp.col(i) - di * tpb.col(i);
    dm.col(i)=ds;
  }
  return(dm);
}
