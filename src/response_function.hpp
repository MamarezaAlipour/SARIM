#ifndef RESPONSEFUNC_H_
#define RESPONSEFUNC_H_

#include <RcppEigen.h>

Eigen::VectorXd response_function(const Eigen::VectorXd& eta,
	const std::string& link);


#endif // IWLS_H_