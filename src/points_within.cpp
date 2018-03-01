#include <R.h>
#include <Rinternals.h>
#include <set>
using namespace std;

/* -- points_witin from the spp package
# determines a set of points within a set of fragments
# note: all vectors sorted in ascending order
# note: all vectors are integers
# x_R - vector of point positions
# se_R - vector of start and end positions
# fi_R - vector of signed fragment indecies
# return_list_R - whether a list of fragments should be returned for each point
# return_unique_R - whether points in multiple fragments should be omitted
*/
// [[Rcpp::export]]
SEXP points_within2(SEXP x_R,SEXP se_R,SEXP fi_R,SEXP return_list_R,SEXP return_unique_R,SEXP return_point_counts_R) {

  int* x=INTEGER(x_R);
  int nx=LENGTH(x_R);
  int* se=INTEGER(se_R);
  int* fi=INTEGER(fi_R);
  int nf=LENGTH(se_R);

  int return_list=*(INTEGER(return_list_R));
  int return_unique=*(INTEGER(return_unique_R));
  int return_point_counts=*(INTEGER(return_point_counts_R));

  set<int> fset;


  SEXP nv; int *i_nv;
  int np=0;
  if(return_point_counts) {
    PROTECT(nv = allocVector(INTSXP, nf/2)); np++;
    i_nv=INTEGER(nv);
    for(int i=0;i<nf/2;i++) { i_nv[i]=0; }
  } else if(return_list) {
    PROTECT(nv = allocVector(VECSXP, nx)); np++;
  } else {
    PROTECT(nv=allocVector(INTSXP,nx));  np++;
    i_nv=INTEGER(nv);
  }

  int j=0;
  for(int i=0;i<nx;i++) {
    // advance j
    while(j<nf && se[j]<x[i]) {
      int frag=fi[j];
      if(frag>0) { // insert
	fset.insert(frag);

      } else { // remove
	fset.erase(-frag);

      }
      j++;
    }

    if(return_list) {
      if(fset.empty() || (return_unique && fset.size()>1)) {
	// assign null list?
      } else {
	SEXP fil_R;
	PROTECT(fil_R=allocVector(INTSXP,fset.size()));  np++;
	int* fil=INTEGER(fil_R);
	int k=0;
	for(set<int>::const_iterator ki=fset.begin();ki!=fset.end();++ki) {
	  fil[k]=*ki; k++;
	}
	SET_VECTOR_ELT(nv, i, fil_R);
	UNPROTECT(1); np--;
      }
    } else {
      if(return_point_counts) {
	for(set<int>::const_iterator ki=fset.begin();ki!=fset.end();++ki) {
	  i_nv[*ki-1]++;
	}
      } else {
	if(fset.empty() || (return_unique && fset.size()>1)) {
	  i_nv[i]=-1;
	} else {
	  i_nv[i]=*fset.begin();
	}
      }
    }
  }

  UNPROTECT(np);
  return nv;
}
