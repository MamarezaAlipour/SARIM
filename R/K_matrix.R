#' Generate the structure/penalty matrix for sarim formulae
#' 
#' @description Function used in definition of smooth terms within sarim model formulae.
#' 
#' @param dimension Generate a squared matrix by (dimension x dimension).
#' @param penalty "identity", "rw1" or "rw2" for structure/penalty matrix.
#' @param data only interesting if penalty "gmrf" is choosen, in this case the data is returned.
#' 
#' @export
K_matrix <- function(dimension, penalty, data) {
    require(Matrix)
    # Identity
    if (penalty == "identity") {
        S = matrix(c(1))
    }
    
    # Random Walk 1st order
    if (penalty == "rw1") {
        P <- diff(.symDiagonal(dimension), differences = 1)
        S <- Matrix::crossprod(P, P)
    }
    
    # Random Walk 2nd order
    if (penalty == "rw2") {
        P <- diff(.symDiagonal(dimension), differences = 2)
        S <- Matrix::crossprod(P, P)
    } 
    
    # Gaussian Markov random field
    if (penalty == "gmrf") {
        S <- data
    }
    
    return(S)
}
