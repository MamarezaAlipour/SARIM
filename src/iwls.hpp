#ifndef IWLS_H_
#define IWLS_H_

#include <RcppEigen.h>

struct IWLS {
	Eigen::SparseMatrix<double> W;
	Eigen::VectorXd y_tilde;
};


IWLS compute(const Eigen::VectorXd& y,
	const Eigen::VectorXd& eta,
	const std::string& fam,
	const std::string& link,
	const double& Ntrials);

#endif // IWLS_H_