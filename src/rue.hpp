#ifndef RUE_H_
#define RUE_H_

#include <RcppEigen.h>


// structure for multiple returns
struct RueSolv
{
	Eigen::VectorXd ga;
	Eigen::VectorXd mu;
};


// Rue-algorithm
//
// first decompose Q by Cholesky, i.e. Q = LL'
// a sample x ~ N(0, Q^{-1}) can then be obtained by sampling z from N(0, I)
// and solving L'x = z
// A sample from N(mu, Q^{-1}) is solving Q mu = b; first by solving Lw = b 
// and subsequently solving L'mu = v
// Apply a linear constraint with  Q^{-1}*A'*(A*Q^{-1}*A')^{-1} * (Ax - e)
// where A a row vector of ones and e = 0
RueSolv algorithm(const Eigen::SparseMatrix<double>& Q,
	const Eigen::VectorXd& b,
	const std::string& lin_con);


#endif // RUE_H_