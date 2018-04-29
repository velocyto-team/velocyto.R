// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

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
arma::mat  embArrows(const arma::mat& emb, const arma::sp_mat& tp, double arrowScale=1.0, int nthreads=1) {
  arma::mat dm(emb.n_cols,emb.n_rows);
  arma::sp_mat tpb(tp); // make binarized version to give equal weights to each cell in the neighborhood
  arma::vec tprs(tp.n_cols,arma::fill::zeros);
  arma::sp_mat::iterator ei=tpb.end();
  for(arma::sp_mat::iterator ci=tpb.begin(); ci!=ei; ++ci) {
    tprs[ci.col()]++; 
  }
  for(arma::sp_mat::iterator ci=tpb.begin(); ci!=ei; ++ci) {
    (*ci)=1.0/tprs[ci.col()];
  }

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

// based on the transition probabilities to each cell, and the diference of expression between them calculates the expected 
// scale 1- linear; 2-sqrt; 3 log
// [[Rcpp::export]]
arma::mat  expectedExpressionShift(const arma::mat& e, const arma::sp_mat& tp, int scale=1, double pseudocount=1.0, int nthreads=1) {
  arma::mat rm(e.n_rows,e.n_cols);
#pragma omp parallel for shared(rm) num_threads(nthreads)
  for(int i=0;i<e.n_cols;i++) {
    arma::mat t(e); t.each_col() -= e.col(i);
    if(scale==2) { 
      t=sqrt(abs(t)) % sign(t);
    } else if(scale==3) {
      t=log10(abs(t)+pseudocount) % sign(t);
    }
    rm.col(i)= t*tp.col(i);
  }
      
  return(rm);
}

// calculates total length of the segments (m has two columns - start and end positions of the segments)
// discounting the overlap. m must be sorted with respect to the first (start) column position.
// [[Rcpp::export]]
double flatLength(NumericMatrix m) {
  double totsize=0;
  double lastEnd=0;
  double lastBegin=1;
  for(int i=0;i<m.nrow();i++) {
    if(m(i,0) > lastEnd) { // new block
      totsize+=lastEnd-lastBegin+1;
      lastEnd=m(i,1);
      lastBegin=m(i,0);
    } else { // continuing old block
      if(m(i,1)>lastEnd) { lastEnd=m(i,1); }
    }
  }
  totsize+=lastEnd-lastBegin+1; // push in final block
  return(totsize);
}


// group overlapping motifs, keeping track of the scores
// note: df must be sorted
// [[Rcpp::export]]
DataFrame groupMotifs(DataFrame df,int msize=1) {
  CharacterVector chr = df["chr"];
  IntegerVector s = df["start"];
  IntegerVector e = df["end"];
  CharacterVector name = df["name"];
  NumericVector score = df["score"];
  CharacterVector strand = df["strand"];
  std::vector<int> cs,ce,ci;
  std::vector<double> cscore;
  std::string rchr(chr[0]);
  int rs=0;
  std::string ichr;
  double maxScore=-1;
  for(int i=1;i<s.size();i++) {
    ichr=chr[i];
    if(rchr!=ichr || s[i]-e[i-1] >=msize) { // end block
      cs.push_back(s[rs]); ce.push_back(e[i-1]);
      ci.push_back(rs);
      cscore.push_back(maxScore);
      rs=i; rchr=ichr; maxScore=-1;
    }
    if(score[i]>maxScore) { maxScore=score[i]; }
  }
  // final check
  if(rs!=s.size()-1) {
    cs.push_back(s[rs]); ce.push_back(e[s.size()-1]); ci.push_back(rs); cscore.push_back(maxScore);
  }
  CharacterVector cchr(ci.size());
  CharacterVector cstrand(ci.size());
  CharacterVector cname(ci.size());
  for(int i=0;i<ci.size();i++) {
    cchr[i]=chr[ci[i]];
    cstrand[i]=strand[ci[i]];
    cname[i]=name[ci[i]];
  }
  //DataFrame rdf = DataFrame::create( Named("chr")=cchr, Named("start") = cs , Named("end") = ce, Named("name")=cname, Named("score")=cscore,Named("strand")=cstrand);
  DataFrame rdf = DataFrame::create(Named("chr")=cchr, Named("start") = cs , Named("end") = ce, Named("name")=cname, Named("score")=cscore, Named("strand")=cstrand);
  return(rdf);
}
