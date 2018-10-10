#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double ARImp(const arma::ivec&  Q, const arma::Mat<double>& M) {
    
    arma::ivec clusters = unique(Q);
    int K = clusters.size();
    int N = Q.size();
    arma::vec t0 = arma::zeros<arma::vec>(2);
    t0[1] = arma::accu(M);

    arma::uvec kidx;
    for (int kk=0; kk<K; kk++)
    {
        kidx = arma::find(Q == clusters[kk]);

        t0[0] += arma::accu(M.submat(kidx, kidx));
    }
    
    t0 -= arma::accu(M.diag());
    t0 /= 2.0;
    
    // t2
    double t2 = sum(choose(as<NumericVector>(table(IntegerVector(Q.begin(),Q.end()))),2));

    // t3    
    double t3 = 2*t0[1]*t2 / (N*(N-1));

    // ARImp
    return (t0[0]-t3)/((t0[1] + t2)/2 - t3);
}

// [[Rcpp::export]]
arma::Mat<double> update_coclus(arma::uvec indii, arma::uvec indjj, arma::Mat<double> M, arma::Mat<double> update) {
    M(indii - 1,indjj - 1) = update;
    M(indjj - 1,indii - 1) = update.t();
    return M;
}

// [[Rcpp::export]]
double update_arimp(const arma::ivec& Q, arma::uvec indii, arma::uvec indjj, arma::Mat<double> M, arma::Mat<double> update) {
    M(indii - 1,indjj - 1) = update;
    M(indjj - 1,indii - 1) = update.t();
    return ARImp(Q, M);
}

// [[Rcpp::export]]
arma::dvec test_pairs(const arma::ivec Q, const arma::Mat<int> clusPairs, arma::Mat<int> num, arma::Mat<double> denom) {
    
    arma::dvec arimps(clusPairs.n_cols);
    int k1p = Q.max() + 1;
    
    arma::uvec indii;
    arma::uvec indjj;
    for (int pp=0; pp<clusPairs.n_cols; pp++) {
        Rcout << clusPairs(0,pp) << arma::endl;
        arma::ivec Qm = Q;
        indii = arma::find(Q == clusPairs(0,pp));
        indjj = arma::find(Q == clusPairs(1,pp));

        Qm.replace(clusPairs(0,pp),k1p);
        Qm.replace(clusPairs(1,pp),k1p);

        num(indii, indjj) += 1;
        num(indjj, indii) += 1;

        arimps[pp] = ARImp(Qm, num/denom);
        
        num(indii, indjj) -= 1;
        num(indjj, indii) -= 1;
    }
    
    return arimps;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
#test_pairs(clus, clusPairs, num, denom)
*/
