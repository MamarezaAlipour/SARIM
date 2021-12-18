#ifndef LANCZOS_H_
#define LANCZOS_H_

#include <RcppEigen.h>


Eigen::MatrixXd tri_mat_creation(const int& n,
	const Eigen::VectorXd& a,
	const Eigen::VectorXd& b);


struct Lanczos
{
	Eigen::VectorXd x;
	Eigen::VectorXd error;
	int Iteration;
};


Lanczos algorithm(const Eigen::SparseMatrix<double>& Q,
	const int& m,
	const Eigen::SparseMatrix<double>& F1,
	const Eigen::SparseMatrix<double>& F2,
	const double& thr);

#endif // LANCZOS_H_