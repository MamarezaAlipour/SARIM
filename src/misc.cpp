#ifndef RUE_H_
#define RUE_H_

#include <RcppEigen.h>
#include <random>
#include "misc.hpp"


// function to calculate the incomplete cholesky decomposition, return a sparse matrix
// need for lanczos algorithm
Eigen::SparseMatrix<double> ichol(const Eigen::SparseMatrix<double>& Q) {
	Eigen::IncompleteCholesky<double> icholOfQ(Q);
	return icholOfQ.matrixL();
};


// faster C++-sampling from a gamma distribution
Eigen::VectorXd random_gamma(const int& n, const double& shape, const double& scale) {
	Eigen::VectorXd e(n);

	std::random_device device_random;
	std::default_random_engine generator(device_random());
	std::gamma_distribution<> distribution(shape, scale);

	for (int i = 0; i < n; ++i)
	{
		e(i) = distribution(generator);
	}
	return e;
};


// faster C++-sampling from a gaussian distribution
Eigen::VectorXd random_gauss(const int& n)
{
	Eigen::VectorXd e(n);

	std::random_device device_random;
	std::default_random_engine generator(device_random());
	std::normal_distribution<> distribution(0, 1);

	for (int i = 0; i < n; ++i)
	{
		e(i) = distribution(generator);
	}
	return e;
};


// faster C++-sampling from a uniform distribution
Eigen::VectorXd random_uniform(const int& n) {
	Eigen::VectorXd e(n);

	std::random_device device_random;
	std::default_random_engine generator(device_random());
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	for (int i = 0; i < n; ++i)
	{
		e(i) = distribution(generator);
	}
	return e;
};

#endif // RUE_H_