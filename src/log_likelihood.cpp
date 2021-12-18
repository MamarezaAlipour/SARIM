// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "log_likelihood.hpp"
#include "iwls.hpp"
#include "response_function.hpp"


// Compute the loglikelihood depending on the family and link function 
Eigen::VectorXd loglike(const Eigen::VectorXd& y,
	const Eigen::VectorXd& eta,
	const std::string& family,
	const std::string& link,
	const double& Ntrials) {

	// Intialise output vectors
	Eigen::VectorXd out;

	Eigen::VectorXd h = response_function(eta, link);

	if (family == "poisson") {
		out = (y.array() * eta.array() - h.array()).colwise().sum();
	}

	if (family == "binomial") {

		/*  Try to make log-likelihood independent of starting values in gamma
		 *  for (int i = 0; i < nrows; ++i) {
		 *      if (h(i) == 1) h(i) = 1 - 1e-10;
		 *      if (h(i) == 0) h(i) = 1e-10;
		 *  }
		 */
		out = (y.array() * h.array().log() + (Ntrials - y.array()).array() * (1 - h.array()).array().log()).colwise().sum();
	}

	return out;
}




