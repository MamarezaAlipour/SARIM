// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

inline Eigen::MatrixXd tri_mat_creation(const int& n,
	const Eigen::VectorXd& a,
	const Eigen::VectorXd& b)
{
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
	Eigen::MatrixXd BU = A;
	Eigen::MatrixXd BL = A;

	A.diagonal() = a.segment(0, n);
	BU.diagonal<1>() = b.segment(1, n);
	BL.diagonal<-1>() = b.segment(1, n);

	return A + BU + BL;
};


// [[Rcpp::export]]
Rcpp::List lanczosCpp(const Eigen::MappedSparseMatrix<double>& Q,
	const int& m,
	const Eigen::MappedSparseMatrix<double>& F1,
	const Eigen::MappedSparseMatrix<double>& F2,
	const double& thr) {
	int iter;
	int n = Q.cols();

	Eigen::VectorXd z = Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(n));

	Eigen::MatrixXd V(n, m + 1);

	V.col(0) = z / z.norm();

	Eigen::VectorXd beta = Eigen::VectorXd::Zero(m + 1);
	Eigen::VectorXd alpha = Eigen::VectorXd::Zero(m);

	Eigen::VectorXd w, x = Eigen::VectorXd::Zero(n);
	Eigen::VectorXd x_old = x;
	Eigen::MatrixXd T;

	Eigen::VectorXd err = alpha;

	Eigen::VectorXd a = x;
	Eigen::VectorXd b = a;


	Eigen::SparseLU<Eigen::SparseMatrix<double> > LUF1(F1);
	Eigen::SparseLU<Eigen::SparseMatrix<double> > LUF2(F2);



	for (int j = 0; j < m; j++) {
		iter = j + 1;
		x_old = x;

		a = LUF2.solve(V.col(j));
		b = Q * a;
		w = LUF1.solve(b);

		if (j > 0) {
			w = w - beta(j) * V.col(j - 1);
		}

		alpha.row(j) = w.transpose() * V.col(j);

		w = w - alpha(j) * V.col(j);

		beta(j + 1) = w.norm();

		V.col(j + 1) = w / w.norm();

		if (j > 1) {
			T = tri_mat_creation(j + 1, alpha, beta);
			Eigen::VectorXd e1 = Eigen::VectorXd::Zero(j + 1);
			e1(0) = 1;
			x = z.norm() * V.leftCols(j) * T.sqrt().lu().solve(e1);
			err(j) = (x - x_old).norm() / x.norm();
			if (err(j) < thr) {
				break;
			}
		}

	}
	x = LUF2.solve(x);

	return Rcpp::List::create(Rcpp::Named("x") = x,
		Rcpp::Named("error") = err,
		Rcpp::Named("No. of Iteration") = iter);

}

