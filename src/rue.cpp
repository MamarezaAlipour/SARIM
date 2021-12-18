// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "rue.hpp"
#include "misc.hpp"


// Rue-algorithm
//
// first decompose Q by Cholesky, i.e. Q = LL'
// a sample x ~ N(0, Q^{-1}) can then be obtained by sampling z from N(0, I)
// and solving L'x = z
// A sample from N(mu, Q^{-1}) is solving Q mu = b; first by solving Lw = b 
// and subsequently solving L'mu = v
// Apply a linear constraint with  Q^{-1}*A'*(A*Q^{-1}*A')^{-1} * (Ax - e)
// where A is a row vector of ones and e = 0

RueSolv algorithm(const Eigen::SparseMatrix<double>& Q,
	const Eigen::VectorXd& b,
	const std::string& lin_con)
{
	RueSolv rue_solver;
	// number of rows in Q, necessary for sampling from z
	int n = Q.rows();

	// initialise vectors for solving system
	Eigen::VectorXd x(n), w(n), mu(n), gamma(n);

	// sample z from N(0, I)
	Eigen::VectorXd z = random_gauss(n);

	// compute Cholesky decomposition of Q and save as "Lt" for upper matrix
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > lltOfQ(Q);
	Eigen::SparseMatrix<double> Lt = lltOfQ.matrixU();

	// solving "L'x = z" to get a sample from N(0, Q^{-1})
	Eigen::SparseLU<Eigen::SparseMatrix<double> > LUtsolve(Lt);
	x = LUtsolve.solve(z);

	// solving "Lw = b"
	Eigen::SparseLU<Eigen::SparseMatrix<double> > LUsolve(Lt.transpose());
	w = LUsolve.solve(b);

	// solving "L'mu = w"
	mu = LUtsolve.solve(w);

	// gamma sample from N(mu, Q^{-1}), i.e. ga = x + mu
	gamma = x + mu;

	// apply a linear constraint if needed
	if (lin_con == "TRUE") {
		Eigen::VectorXd v(n), u(n);
		Eigen::VectorXd At = Eigen::VectorXd::Ones(n);
		Eigen::VectorXd e = Eigen::VectorXd::Zero(1);
		u = LUsolve.solve(At);
		v = LUtsolve.solve(u);

		gamma = gamma - v * ((At.transpose() * v).cwiseInverse()) * (At.transpose() * gamma - e);
		mu = mu - v * ((At.transpose() * v).cwiseInverse()) * (At.transpose() * mu - e);
	}

	// return of gamma and mu
	rue_solver.ga = gamma;
	rue_solver.mu = mu;
	return rue_solver;
};