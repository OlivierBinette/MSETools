#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec colSumsSub(const arma::mat& mat,
                        const arma::uvec& rowIDs) {

  return arma::sum(mat.rows(rowIDs-1), 0);
}

// [[Rcpp::export]]
double dec(const Rcpp::IntegerVector I, int p) {
  double s = 0;
  for (int i = 0; i < I.size(); i++) {
    s += pow(2, p - I(i));
  }
  return s;
}

// [[Rcpp::export]]
arma::uvec decClique(const Rcpp::List cliques, int p) {
  int n = cliques.size();
  arma::uvec res(n);

  for(int i = 0; i < n; i++) {
    res(i) = dec(as<IntegerVector>(cliques[i]), p);
  }

  return res;
}

// [[Rcpp::export]]
void computeML(arma::mat& inPlace,
                        int j,
                        const arma::mat& compMat,
                        const Rcpp::List& C,
                        const Rcpp::List& S,
                        const arma::rowvec& denominator,
                        int p) {
  arma::uvec decC = decClique(C, p);
  arma::uvec decS = decClique(S, p);

  double nsubgraphs = decC.n_elem - decS.n_elem;

  inPlace.row(j-1) = colSumsSub(compMat, decC)
    - colSumsSub(compMat, decS)
    + nsubgraphs*denominator;
}

// [[Rcpp::export]]
arma::mat computeLogPostProbs(const arma::mat& compMat,
                              const Rcpp::List& graphs,
                              const arma::rowvec& denominator,
                              int p) {
  int N = graphs.size();
  int Nmissing = compMat.n_cols;

  arma::mat weights(N, Nmissing);
  Rcpp::List graph;
  Rcpp::List emptyList = Rcpp::List::create();
  for (int j = 1; j <= N; j++) {
    graph = as<List>(graphs[j-1]);
    if (graph.size() == 2) { // Both cliques and separators
      computeML(weights, j,
                compMat,
                as<List>(graph)[0],
                as<List>(graph)[1],
                denominator,
                p);
    } else { // Only cliques
      computeML(weights, j,
                compMat,
                as<List>(graph)[0],
                emptyList,
                denominator,
                p);
    }
  }

  return weights;
}

// [[Rcpp::export]]
void rowAdd(arma::mat& mat,
            const arma::rowvec& v) {

  mat.each_row() += v;
}

// [[Rcpp::export]]
void colAdd(arma::mat& mat,
            const arma::colvec& v) {

  mat.each_col() += v;
}

// [[Rcpp::export]]
void expNormalize(arma::mat& mat) {
  mat = exp(mat);
  mat = mat/sum(sum(mat));
}


