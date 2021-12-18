// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "response_function.hpp"


// calculation of response function for log or logit
Eigen::VectorXd response_function(const Eigen::VectorXd& eta,
	const std::string& link) {
	Eigen::VectorXd out;

	// log-function for poisson model
	if (link == "log") {
		out = eta.array().exp();
	}

	// logit-function for binomial model
	if (link == "logit") {
		out = 1 / (1 + (-eta).array().exp());
	}

	return out;
};