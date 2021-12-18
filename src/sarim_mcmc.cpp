// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

#include <RcppEigen.h>
#include "rue.hpp"
#include "lanczos.hpp"
#include "misc.hpp"
#include "iwls.hpp"
#include "log_likelihood.hpp"

//' MCMC sampler with Metropolis-Hasting step in Gibbs-sampler for use in sarim()-function
//' 
//' This generates samples for the coefficients using a Metropolis-Hasting step in
//' the Gibbs-sampling for the precision parameters and further assumed that the 
//' response variable is either binomial or poisson. So family = "binomial"/"poisson" with
//' link = "logit"/"log" is choosen respectively.
//' 
//' @param y Response variable, given as a vector, use as.numeric() if error occur.
//' @param Z Design matrices of the covariates, given as list with sparse matrix, 
//'          use e.g. as(matrix, "dgCMatrix") from library(Matrix). Use for example the
//'          useful sx()-function for smoothing.
//' @param K Structure/penalty matrices for the coefficients, given as list, 
//'          also with sparse matrix. Can be choosen in sx()-function.
//' @param K_rank List of ranks of the structure/penalty matrix.
//' @param gamma List of coefficient, given as vector. Row-length need to be the 
//'          same as the columns of Z. Per starting default from uniform distribution 
//'          is sampled, but a specific starting value can be given, using the sx()-function
//'          in the formular, e.g. y ~ sx(x1, gamma = c(rep(1, 5))).
//' @param ka_start Start value for kappa, given as a list and double/float value.
//' @param ka_values Coefficients for kappa, given as list within as vector c(kappa_a, kappa_b).
//' @param solver List of the solvers ("rue" or "lanczos") for sampling from a gaussian distribution, i.e. gamma ~ N(eta, Q)
//'          with Q as precision matrix. Can be choosen in sx()-function.
//' @param lin_constraint Specify if a linear constraint for the coefficient is needed. 
//'          Given as list with "TRUE" or "FALSE" values. Can be choosen in sx()-function.
//' @param family Can currently only be "binomial" or "poisson".
//' @param link Link need to be chooses correspondingly to the family function. If
//'          family = "binomial", please choose link = "logit" or if family = "poisson", 
//'          please choose link = "log". Currently no other link function is present.
//' @param nIter Number of iterations for MCMC-algorithm
//' @param Ntrials Number of trails, only interesting for binomial distribution
//' @param m Number of maximal Lanczos-iterations
//' @param thr threshold when the Lanczos-algorithm or conjugate gradient-algorithm should stop
//' 
//' @return Return a list of values:
//' "coef_results" = result of the estimated coefficient, output given as matrix;
//' "kappa_result" = result of the estimated kappa (precision) parameters, output given as vector;
//' "accept_rate" = acceptance rate for each coefficient;
//' "lanzcos_iterations" = number of the lanczos-iteration in each step.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List sarim_mcmc(const Eigen::Map<Eigen::VectorXd> & y,
                      const Rcpp::List & Z,
                      const Rcpp::List & K,
                      const Rcpp::List & K_rank,
                      const Rcpp::List & gamma,
                      const Rcpp::List & ka_start,
                      const Rcpp::List & ka_values,
                      const Rcpp::List & solver,
                      const Rcpp::List & lin_constraint,
                      const Rcpp::String & family,
                      const Rcpp::String & link,
                      const double & Ntrials,
                      const int & nIter,
                      const int & m,
                      const double & thr,
                      const bool & display_progress = true) {
    
    
    // display progress bar
    Progress pro(nIter, display_progress);
    
    // number of covariates (p), get from R:List Z
    int p = Z.size();
    
    // number of observations
    int n = y.rows();
    
    
    ////////////////////////////////////////////////////////////////////////////
    // generate list for further calculations and better output for R
    Rcpp::List coef_results(p);         // list for coefficients
    Rcpp::List mu_results(p);           // list for temporary mu for faster iterations steps
    Rcpp::List v_history(p);            // list for temporary v for faster iterations steps
    Rcpp::List iterative_sampling(p);   // list for iterative samples in lanczos-algo
    Rcpp::List kappa_results(p);        // list for kappa values
    Rcpp::List m_iter(p);               // list for max-lanczos-iterations
    
    for (int i = 0; i < p; ++i) {
        // gamma matrix 
        Eigen::VectorXd gamma_tmp = gamma(i);
        int k_size = gamma_tmp.size();
        Eigen::MatrixXd gamma_results = Eigen::MatrixXd::Zero(k_size, nIter + 1);
        gamma_results.col(0) = gamma_tmp;
        coef_results[i] = gamma_results;
        
        // kappa-value vector
        double kappa_tmp = ka_start(i);
        Eigen::VectorXd kappa_results_tmp(nIter + 1);
        kappa_results_tmp(0) = kappa_tmp;
        kappa_results[i] = kappa_results_tmp;
        
        // mu, for eventually faster calculation of mean form gamma ~ N(mu, Q)
        mu_results[i] = 0 * gamma_tmp;
        
        // history for eventually linear constraint, for faster calculation
        v_history[i] = 0 * gamma_tmp;
        
        // number of iterations in lanczos-algo
        iterative_sampling[i] = Eigen::VectorXd::Zero(nIter);
        
        // starting value for lanczos-algo
        m_iter[i] = m;
    };
    

    // list for acceptance rate
    int ac = 0;
    Rcpp::List ac_list(p);
    for (int i = 0; i < p; ++i) {
        ac_list[i] = ac;
    }
    
    // initialise eta the first time, eta = Z_1 * gamma_1 + ... + Z_p * gamma_p
    Eigen::VectorXd eta = Eigen::VectorXd::Zero(n);
    // eta_tmp for calculating eta_{-k}
    Eigen::VectorXd eta_tmp = eta;
    /*
    for (int i = 0; i < p; ++i) {
        Eigen::SparseMatrix<double> Z_tmp = Z(i);
        Eigen::VectorXd gamma_tmp = gamma(i);
        eta += Z_tmp * gamma_tmp;
    };
    */
    
    // initialise values, required for computation
    Eigen::SparseMatrix<double> Z_k, K_k, Q, W, ZtW, M, Mt;
    Eigen::MatrixXd gamma_matrix, gamma_current;
    
    Eigen::VectorXd K_rk, // Rank for penalty matrix
        b,              // response vector, i.e. Z' * W (y_tilde - eta_{-k})
        ka_vector,      // vector for 
        y_tilde,        // working observations
        gamma_proposal, // proposal for gamma, could be or not rejected
        mu,             // current mean for gamma ~ N(mu, Q)
        mu_tmp,         // proposal mean
        ka_tmp,         // vector for kappa-values, as temporary value 
        x,              // proposal for sampling from N(0, Q^{-1})
        v_hist,         // history of v, needed linear constraint and for eventually faster solving 
        it_sampling,    // vector for numbers of lanczos-iterations
        ll,             // log-likelihood with current gamma
        ll_proposal,    // log-likelihood with proposal gamma
        proposal_c_given_p,     // conditional likelihood given 
        proposal_p_given_c, 
        prior_c, 
        prior_p, 
        alpha,          // acceptance probability
        u;              // u ~ U[0, 1]
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Initialise working observations and weights
    for (int k = 0; k < p; ++k) {
        Z_k = Z(k);                     // design matrix
        K_k = K(k);                     // penalty/structure matrix
        gamma_matrix = coef_results(k); // current 
        Eigen::VectorXd gamma_old = gamma_matrix.col(0);
        ka_vector = kappa_results(k);
        std::string solv = solver(k);
        std::string lin_con = lin_constraint(k);
        
        
        // compute weights and working response from file "iwls.cpp"
        IWLS iwls = compute(y, eta, family, link, Ntrials);
        W = iwls.W;
        y_tilde = iwls.y_tilde;
        
        // compute Z' * W once for efficiency
        ZtW = Z_k.transpose() * W;
        
        // calculate Q = Z' * W * Z + kappa * K  and  b = Z' * W * (y - eta_{-k})
        b = ZtW * ( y_tilde - (eta - Z_k * gamma_matrix.col(0)) );
        Q = ZtW * Z_k + ka_vector(0) * K_k;
            
        // initalise temporary gamma
        Eigen::VectorXd ga_tmp; 
        
        // set start values for gamma
        if (solv == "rue") {
            // solve linear system Q * ga = b by Rue's method
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > lltOfQ(Q); 
            Eigen::SparseMatrix<double> Lt = lltOfQ.matrixU();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > LUtsolve(Lt);
            Eigen::SparseLU<Eigen::SparseMatrix<double> > LUsolve(Lt.transpose());
            ga_tmp = LUtsolve.solve(LUsolve.solve(b));
            
            if (lin_con == "TRUE") {
                unsigned int n = ga_tmp.rows();
                Eigen::VectorXd At = Eigen::VectorXd::Ones(n);
                Eigen::VectorXd e = Eigen::VectorXd::Zero(1);
                Eigen::VectorXd v = LUtsolve.solve(LUsolve.solve(At));
                
                ga_tmp = ga_tmp - v * ( (At.transpose() * v).cwiseInverse() ) * (At.transpose() * ga_tmp - e);
            }
            
        } else {
            // solve linear system Q * ga = b by conjugate gradient method
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 
                                     Eigen::Lower|Eigen::Upper, 
                                     Eigen::IncompleteCholesky<double> > iccg(Q);
            iccg.setTolerance(thr); // set tolerance for convergence
            ga_tmp = iccg.solve(b); // solve Q * ga = b
            
            // apply linear constraint if needed
            if (lin_con == "TRUE") {
                // get number of parameters
                unsigned int n = ga_tmp.rows();
                Eigen::VectorXd At = Eigen::VectorXd::Ones(n); // column vector of ones
                Eigen::VectorXd e = Eigen::VectorXd::Zero(1);  // single value of zero   
                Eigen::VectorXd v = iccg.solve(At);
                
                ga_tmp = ga_tmp - v * ( (At.transpose() * v).cwiseInverse() ) * (At.transpose() * ga_tmp - e);
            }
            
        };
        
        eta += Z_k * (ga_tmp - gamma_old);
        gamma_matrix.col(0) = ga_tmp;
        coef_results[k] = gamma_matrix;
        mu_results[k] = ga_tmp;
    }
    // end initalisation
    
    
    ////////////////////////////////////////////////////////////////////////////
    // ITERATION
    // loop for number of iterations "nIter" 
    // C++ begin in 0 !!!   important for loops
    for (int n_mcmc = 1; n_mcmc <= nIter; ++n_mcmc) {
        // if (Progress::check_abort() )
        //     return -1.0;
        pro.increment();
        
        // Update regression coefficients, one at a time
        // loop for covariates p 
        for (int k = 0; k < p; ++k) {
            
            // initialise required matrices Z and K and coefficient gamma and kappa
            // set ._k for required calculations
            Z_k = Z(k);                     // design matrix
            K_k = K(k);                     // penalty/structure matrix
            K_rk = K_rank(k);               // rank of penalty/structure matrix
            std::string solv = solver(k);   // which solver should be used
            std::string lin_con = lin_constraint(k); // linear constraint necessary?
            gamma_matrix = coef_results(k); // // coefficient matrix
            gamma_current = gamma_matrix.col(n_mcmc - 1); // current gamma for simpler calculation
            ka_vector = kappa_results(k);   // vector for precision parameter
            mu = mu_results(k);             // previous value for expectation
            v_hist = v_history(k);          // history of v, linear constraint
            it_sampling = iterative_sampling(k);    // number of iterations in lanczos
            int max_iter_of_lanczos = m_iter(k);    // max-lanczos-iterations
            ////////////////////////////////////////////////////////////////////
            
            
            // First calculate the log-likelihood and 
            // eta_tmp = eta + Z_k * (mu - gamma), i.e. no requirement of log|Q|
            ll = loglike(y, eta, family, link, Ntrials);
            eta_tmp = eta + Z_k * (mu - gamma_current);
            
            // calculate weight matrix W and working observations
            IWLS iwls = compute(y, eta_tmp, family, link, Ntrials);
            W = iwls.W;
            y_tilde = iwls.y_tilde;
            
            // calculate one time for efficiency Z'W
            ZtW = Z_k.transpose() * W;
            
            // calculate Q (precision) matrix 
            Q = ZtW * Z_k + ka_vector(n_mcmc - 1) * K_k;
            
            // calculate b;  b =  Z'_{k} W * (y_tilde - eta_{-k})
            b = ZtW * (y_tilde - (eta_tmp - Z_k * mu) ); 
            
            
            // sample from gaussian distribution
            // depending on solver 'rue' or 'lanczos'
            if (solv == "rue") {
                // simple Rue-block-algorithm; see file "rue.cpp"/"rue.hpp"
                // include linear-constraint if needed
                RueSolv rue_solver = algorithm(Q, b, lin_con);
                gamma_proposal = rue_solver.ga;
                mu_tmp = rue_solver.mu;
                
            } else {
                // Lanczos approximation of x ~ N(0, Q^{-1})
                // first calculate incomplete cholesky; see file "misc.cpp"/"misc.hpp"
                M = ichol(Q);
                Mt = M.transpose();
                
                // Lanczos algorithm; see file "lanczos.cpp"/"lanczos.hpp"
                Lanczos lanczos_solver = algorithm(Q, max_iter_of_lanczos, M, Mt, thr);
                
                // save number of lanczos-iteration
                it_sampling(n_mcmc - 1) = lanczos_solver.Iteration;
                iterative_sampling[k] = it_sampling;
                m_iter[k] = it_sampling(n_mcmc - 1) + 10;
                if (lanczos_solver.Iteration >= m) {
                    Rcpp::Rcout << "Number of lanczos-iteration <= " << m <<" in iteration " 
                                << n_mcmc << ". Set m higher (>" << m << ")" << std::endl;
                }
                
                // sample from N(0, Q^{-1})
                x = lanczos_solver.x;
                
                // solve linear system Q*mu = b
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 
                                         Eigen::Lower|Eigen::Upper, 
                                         Eigen::IncompleteCholesky<double> > iccg(Q);
                iccg.setTolerance(thr);     // set tolerance threshold for convergence
                iccg.solveWithGuess(b, mu); // guess initial solution for eventually faster solving
                mu_tmp = iccg.solve(b);
                
                // proposal for gamma ~ N(mu, Q^{-1}), i.e. gamma = x + mu
                gamma_proposal = x + mu_tmp;
                
                // apply linear constraint if needed
                if (lin_con == "TRUE") {
                    // get number of parameters
                    unsigned int n = gamma_proposal.rows();
                    Eigen::VectorXd At = Eigen::VectorXd::Ones(n);
                    Eigen::VectorXd e = Eigen::VectorXd::Zero(1);
                    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 
                                             Eigen::Lower|Eigen::Upper, 
                                             Eigen::IncompleteCholesky<double> > iccgOfA(Q);
                    iccgOfA.setTolerance(thr);          // set tolerance threshold for convergence
                    iccgOfA.solveWithGuess(b, v_hist);  // guess initial solution for eventually faster solving
                    Eigen::VectorXd v = iccgOfA.solve(At);
                    
                    v_history[k] = v;
                    
                    gamma_proposal = gamma_proposal - v * ( (At.transpose() * v).cwiseInverse() ) * (At.transpose() * gamma_proposal - e);
                    mu_tmp = mu_tmp - v * ( (At.transpose() * v).cwiseInverse() ) * (At.transpose() * mu_tmp - e);
                }
                
            };
            
            // compute eta and loglikelihood
            Eigen::VectorXd eta_tmp_proposal = eta_tmp + Z_k * (gamma_proposal - mu);
            ll_proposal = loglike(y, eta_tmp_proposal, family, link, Ntrials);
            
            // compute p(ga_c | \mu^p, Q^p)
            proposal_c_given_p = -0.5 * ( (gamma_current - mu).transpose() * Q * (gamma_current - mu) );
            
            // compute p(ga_p | \mu^c, Q^c)  
            proposal_p_given_c = -0.5 * ( (gamma_proposal - mu).transpose() * Q * (gamma_proposal - mu) );
            
            // compute p(ga_c | \ka^c)
            prior_c = - (ka_vector.row(n_mcmc - 1) * 0.5) * gamma_current.transpose() * K_k * gamma_current;
            
            // compute p(ga_p | \ka^c)
            prior_p = - (ka_vector.row(n_mcmc - 1) * 0.5) * gamma_proposal.transpose() * K_k * gamma_proposal;
            
            // acceptance probability
            alpha = (ll_proposal + prior_p + proposal_c_given_p) - (ll + prior_c + proposal_p_given_c); 
            
            // accept or reject the proposal?
            u = (random_uniform(1)).array().log(); // sample from file "misc.cpp"/"misc.hpp"
            if (alpha(0) > u(0)) {
                gamma_matrix.col(n_mcmc) = gamma_proposal;
                coef_results[k] = gamma_matrix;
                mu_results[k] = mu_tmp;
                eta = eta_tmp_proposal;   
                // increase accept
                int accept = ac_list(k);
                accept += 1;
                ac_list[k] = accept;
                
            } else {
                gamma_matrix.col(n_mcmc) = gamma_current;
                coef_results[k] = gamma_matrix;
                mu_results[k] = mu;
            }
            
            
            // update kappa_k by sampling from  
            // ~ Ga(kappa_alpha + rk(K_i)/2 , kappa_beta + gamma_k ' * K_k * gamma_k / 2)
            ka_tmp = ka_values(k);
            double ka_alpha = ka_tmp.coeff(0, 0) + 0.5 * K_rk(0);
            double ka_beta;
            ka_beta = ka_tmp.coeff(1, 0) + 
                0.5 * (gamma_matrix.col(n_mcmc)).transpose() * K_k * gamma_matrix.col(n_mcmc);
            ka_vector.row(n_mcmc) = random_gamma(1, ka_alpha, 1/ka_beta);
            kappa_results[k] = ka_vector;
        };
        
    };
    
    
    // calculate the acceptance rate
    Rcpp::List ac_rate(p);
    for (int i = 0; i < p; ++i) {
        double accept = ac_list(i);
        double acceptrate = accept / (1.0 * nIter);
        ac_rate[i] = acceptrate;
    }
    
    // return lists for gamma, kappa and accept_rate 
    return Rcpp::List::create(Rcpp::Named("coef_results") = coef_results, 
                              Rcpp::Named("kappa_results") = kappa_results,
                              Rcpp::Named("accept_rate") = ac_rate,
                              Rcpp::Named("lanzcos_iterations") = iterative_sampling);

}
