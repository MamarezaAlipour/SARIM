#' Defining smooth in sarim formulae
#' 
#' @description Function used in definition of smooth terms within sarim model formulae.
#' 
#' @param x covariates to be smooth, similar to bs form spline package.
#' @param knots number of knots, only equidistant is acutual available.
#' @param degree the order of the Basis functions, default = 3 (cubic splines).
#' @param penalty "identity", "rw1", "rw2" for structure/penalty matrix 
#'                or "gmrf" if your own design matrix Z is given.
#' @param solver "rue" or "lanczos" for sampling from a normal distribution. "rue" is
#'               default and a good choose for normal smoothing. "lanczos" is an iterative methods
#'               especially suitable for BIG STRUCTURE MATRICES, i.e. dim(Z) > 5000.
#' @param ka_start starting value for kappa .
#' @param ka_a a value for sampling form priori kappa~Ga(ka_a, ka_b), see kappa.
#' @param ka_b a value for sampling form priori kappa~Ga(ka_a, ka_b), see kappa.
#' @param linear_constraint apply a linear constraint if it needed.
#' @param gamma optional starting values for coefficent of gamma, the length(gamma)  
#' should be the same as the columns of x, i.e. the number knots, which are  
#' choosen in this sx-function.
#' @param Z optional design matrix, given for spatial or temporal effects, in this case
#' the penalty need to be set to "gmrf". If Z is given to the function, an additional
#' structure/penalty matrix K need to be set up.
#' @param K structure/penalty matrix, given as the penalty for the prior parameters,
#' could be choose for example as intrinsic Gaussian Markov random field (IGMRF).
#' 
#' @export
sx <- function(x = x, knots = 5, penalty = "rw2", degree = 3, solver = "rue", 
               ka_start = 1, ka_a = 1, ka_b = 0.0001, linear_constraint = "FALSE",  
               gamma = NULL, Z = NULL, K = NULL) {
    require(Matrix)
    penlist <- list("identity", "rw1", "rw2", "gmrf")
    if (!is.element(penalty, penlist)) 
        stop('The penalty is unknown, use "identity", "rw1", "rw2" or "gmrf" ')
    
    solvlist <- list("rue", "lanczos")
    if (!is.element(solver, solvlist)) 
        stop('The solver is unknown, use "rue" (as default) or "lanczos", iff Z is a BIG matrix ')
    
    constraintlist <- list("TRUE", "FALSE")
    if (!is.element(linear_constraint, constraintlist)) 
        stop('Linear constraint can only be "TRUE" or "FALSE" ')
    
    if ( (knots == 1) && (is.null(Z)) ) {
        mat <- as.matrix(x)
        if (penalty != "identity")
            warning('The penalty "identity" is choosen per default ')
        attr(mat, "penalty") <- "identity"
    }
    
    if ( (knots > 1) && (is.null(Z)) ) {
        require(splines)
        if (knots < 3) {
            stop('Choose more knots, default is "knots = 5", recommended is also more 
                             for example "knots = 30", or choose "knots = 1" for one column covariate ')
        }
        mat <- splines::bs(x = x, df = knots, degree = degree)
        attr(mat, "penalty") <- penalty
    }
    
    
    if (!is.null(Z)) {
        mat <- as.matrix(rep(0, nrow(Z)))
        if (penalty != "gmrf") 
            stop('If Z is given, choose penalty as gmrf, i.e. penalty = "gmrf" ')
        
        attr(mat, "penalty") <- penalty
        attr(mat, "K") <- K
        attr(mat, "Z") <- Z
    }
    
    attr(mat, "ka_start") <- ka_start
    attr(mat, "ka_a") <- ka_a
    attr(mat, "ka_b") <- ka_b
    attr(mat, "solver") <- solver
    attr(mat, "linear_const") <- linear_constraint
    
    if (!is.null(gamma)) {
        if (is.null(attr(mat, "Z"))) {
            if (length(gamma) != ncol(mat)) {
                stop('The gamma length does not match with columns of x ')
            } else {
                attr(mat, "gamma") <- gamma
            }
        } else {
            if (length(gamma) != ncol(attr(mat, "Z"))) {
                stop('The gamma length does not match with columns of x ')
            } else {
                attr(mat, "gamma") <- gamma
            }
        }
    }
    
    return(mat)
}