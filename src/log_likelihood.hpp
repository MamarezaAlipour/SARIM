#ifndef LOG_LIKELIHOOD_H_
#define LOG_LIKELIHOOD_H_

// Headerfunction for log-likelihood

#include <RcppEigen.h>

Eigen::VectorXd loglike(const Eigen::VectorXd& y,
	const Eigen::VectorXd& eta,
	const std::string& family,
	const std::string& link,
	const double& Ntrials);



#endif // LOG_LIKELIHOOD_H_
