#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// icholCpp
Rcpp::List icholCpp(const Eigen::MappedSparseMatrix<double>& Q);
RcppExport SEXP _Sarim_icholCpp(SEXP QSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type Q(QSEXP);
	rcpp_result_gen = Rcpp::wrap(icholCpp(Q));
	return rcpp_result_gen;
	END_RCPP
}
// lanczosCpp
Rcpp::List lanczosCpp(const Eigen::MappedSparseMatrix<double>& Q, const int& m, const Eigen::MappedSparseMatrix<double>& F1, const Eigen::MappedSparseMatrix<double>& F2, const double& thr);
RcppExport SEXP _Sarim_lanczosCpp(SEXP QSEXP, SEXP mSEXP, SEXP F1SEXP, SEXP F2SEXP, SEXP thrSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type Q(QSEXP);
	Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
	Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type F1(F1SEXP);
	Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type F2(F2SEXP);
	Rcpp::traits::input_parameter< const double& >::type thr(thrSEXP);
	rcpp_result_gen = Rcpp::wrap(lanczosCpp(Q, m, F1, F2, thr));
	return rcpp_result_gen;
	END_RCPP
}
// sarim_gibbs
Rcpp::List sarim_gibbs(const Eigen::Map<Eigen::VectorXd>& y, const Rcpp::List& Z, const Rcpp::List& K, const Rcpp::List& K_rank, const Rcpp::List& gamma, const Rcpp::List& ka_start, const Rcpp::List& ka_values, const Rcpp::List& solver, const Rcpp::List& lin_constraint, const double& sigma, const Eigen::Map<Eigen::VectorXd>& sigma_values, const int& nIter, const int& m, const double& thr, const bool& display_progress);
RcppExport SEXP _Sarim_sarim_gibbs(SEXP ySEXP, SEXP ZSEXP, SEXP KSEXP, SEXP K_rankSEXP, SEXP gammaSEXP, SEXP ka_startSEXP, SEXP ka_valuesSEXP, SEXP solverSEXP, SEXP lin_constraintSEXP, SEXP sigmaSEXP, SEXP sigma_valuesSEXP, SEXP nIterSEXP, SEXP mSEXP, SEXP thrSEXP, SEXP display_progressSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type y(ySEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type Z(ZSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type K(KSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type K_rank(K_rankSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type gamma(gammaSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type ka_start(ka_startSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type ka_values(ka_valuesSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type solver(solverSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type lin_constraint(lin_constraintSEXP);
	Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
	Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type sigma_values(sigma_valuesSEXP);
	Rcpp::traits::input_parameter< const int& >::type nIter(nIterSEXP);
	Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
	Rcpp::traits::input_parameter< const double& >::type thr(thrSEXP);
	Rcpp::traits::input_parameter< const bool& >::type display_progress(display_progressSEXP);
	rcpp_result_gen = Rcpp::wrap(sarim_gibbs(y, Z, K, K_rank, gamma, ka_start, ka_values, solver, lin_constraint, sigma, sigma_values, nIter, m, thr, display_progress));
	return rcpp_result_gen;
	END_RCPP
}
// sarim_mcmc
Rcpp::List sarim_mcmc(const Eigen::Map<Eigen::VectorXd>& y, const Rcpp::List& Z, const Rcpp::List& K, const Rcpp::List& K_rank, const Rcpp::List& gamma, const Rcpp::List& ka_start, const Rcpp::List& ka_values, const Rcpp::List& solver, const Rcpp::List& lin_constraint, const Rcpp::String& family, const Rcpp::String& link, const double& Ntrials, const int& nIter, const int& m, const double& thr, const bool& display_progress);
RcppExport SEXP _Sarim_sarim_mcmc(SEXP ySEXP, SEXP ZSEXP, SEXP KSEXP, SEXP K_rankSEXP, SEXP gammaSEXP, SEXP ka_startSEXP, SEXP ka_valuesSEXP, SEXP solverSEXP, SEXP lin_constraintSEXP, SEXP familySEXP, SEXP linkSEXP, SEXP NtrialsSEXP, SEXP nIterSEXP, SEXP mSEXP, SEXP thrSEXP, SEXP display_progressSEXP) {
	BEGIN_RCPP
		Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type y(ySEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type Z(ZSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type K(KSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type K_rank(K_rankSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type gamma(gammaSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type ka_start(ka_startSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type ka_values(ka_valuesSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type solver(solverSEXP);
	Rcpp::traits::input_parameter< const Rcpp::List& >::type lin_constraint(lin_constraintSEXP);
	Rcpp::traits::input_parameter< const Rcpp::String& >::type family(familySEXP);
	Rcpp::traits::input_parameter< const Rcpp::String& >::type link(linkSEXP);
	Rcpp::traits::input_parameter< const double& >::type Ntrials(NtrialsSEXP);
	Rcpp::traits::input_parameter< const int& >::type nIter(nIterSEXP);
	Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
	Rcpp::traits::input_parameter< const double& >::type thr(thrSEXP);
	Rcpp::traits::input_parameter< const bool& >::type display_progress(display_progressSEXP);
	rcpp_result_gen = Rcpp::wrap(sarim_mcmc(y, Z, K, K_rank, gamma, ka_start, ka_values, solver, lin_constraint, family, link, Ntrials, nIter, m, thr, display_progress));
	return rcpp_result_gen;
	END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
	{"_Sarim_icholCpp", (DL_FUNC)&_Sarim_icholCpp, 1},
	{"_Sarim_lanczosCpp", (DL_FUNC)&_Sarim_lanczosCpp, 5},
	{"_Sarim_sarim_gibbs", (DL_FUNC)&_Sarim_sarim_gibbs, 15},
	{"_Sarim_sarim_mcmc", (DL_FUNC)&_Sarim_sarim_mcmc, 16},
	{NULL, NULL, 0}
};

RcppExport void R_init_Sarim(DllInfo* dll) {
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
