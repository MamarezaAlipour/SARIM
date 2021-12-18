// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


// [[Rcpp::export]]
Rcpp::List icholCpp(const Eigen::MappedSparseMatrix<double>& Q)
{
	Eigen::IncompleteCholesky<double> icholOfQ(Q);
	return Rcpp::List::create(Rcpp::Named("ichol") = icholOfQ.matrixL());
};

