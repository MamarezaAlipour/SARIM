#ifndef MISC_H_
#define MISC_H_

#include <RcppEigen.h>
#include <random>

Eigen::SparseMatrix<double> ichol(const Eigen::SparseMatrix<double>& Q);

Eigen::VectorXd random_gamma(const int& n, const double& shape, const double& scale);

Eigen::VectorXd random_gauss(const int& n);

Eigen::VectorXd random_uniform(const int& n);

#endif // MISC_H_
