#######################################
## CONCORD - ISTA
#######################################


#' @useDynLib gconcordopt
#' @export
concordista_e <- function(D, lam, pMat, penalize.diagonal = 0L, BBstep = 0L,
                       tol = 1e-5, maxit = 100L, info = FALSE, 
                       trkeigs = FALSE, DisS = FALSE, en = 1.0, ...) {
    
    soln <- .Call('gconcordopt_cc_elastic_ista', PACKAGE = 'gconcordopt',
                  D, lam, pMat, as.integer(penalize.diagonal), as.integer(BBstep),
                  tol, maxit, as.integer(info), as.integer(trkeigs), as.integer(DisS),
                  as.double(en))

    soln$omega <- as(soln$omega, 'symmetricMatrix')
    output <- structure(soln$omega,
                        lambda = lam,
                        iter = soln$itr,
                        misc = soln$info)

    output
    
}

#' @useDynLib gconcordopt
#' @export
concordista <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                        penalize.diagonal = 0L, BBstep = 0L,
                        tol = 1e-5, maxit = 100L, info = FALSE, 
                        trkeigs = FALSE, DisS = FALSE, ...) {
    p <- ncol(D)
    soln <- .Call('gconcordopt_cc_ista', PACKAGE = 'gconcordopt',
                  D, lam, pMat, as(X0, 'dgCMatrix'),
                  as.integer(penalize.diagonal), as.integer(BBstep),
                  tol, maxit, as.integer(info), as.integer(trkeigs), as.integer(DisS))

    soln$omega <- as(soln$omega, 'symmetricMatrix')
    output <- structure(soln$omega,
                        lambda = lam,
                        iter = soln$itr,
                        misc = soln$info)

    output
    
}

#' @useDynLib gconcordopt
#' @export
concordista_0_0 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                            penalize.diagonal = 0L, BBstep = 0L, ...) {

    concordista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(BBstep), ...)
    
}

#' @useDynLib gconcordopt
#' @export
concordista_1_0 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                            penalize.diagonal = 1L, BBstep = 0L, ...) {

    concordista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(BBstep), ...)
    
}

#' @useDynLib gconcordopt
#' @export
concordista_0_1 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                            penalize.diagonal = 0L, BBstep = 1L, ...) {
    
    concordista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(BBstep), ...)
    
}

#' @useDynLib gconcordopt
#' @export
concordista_1_1 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                            penalize.diagonal = 1L, BBstep = 1L, ...) {
    
    concordista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(BBstep), ...)
    
}

#######################################
## CONCORD - FISTA
#######################################

#' @useDynLib gconcordopt
#' @export
concordfista <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                         penalize.diagonal = 0L, steptype = 0L,
                         tol = 1e-5, maxit = 100L, info = FALSE, 
                         trkeigs = FALSE, DisS = FALSE, ...) {
    
    soln <- .Call('gconcordopt_cc_fista', PACKAGE = 'gconcordopt',
                  D, lam, pMat, as(X0, 'dgCMatrix'), as.integer(penalize.diagonal), as.integer(steptype),
                  tol, maxit, as.integer(info), as.integer(trkeigs),
                  as.integer(DisS))

    soln$omega <- as(soln$omega, 'symmetricMatrix')
    output <- structure(soln$omega,
                        lambda = lam,
                        convcond = soln$convcond,
                        iter = soln$itr,
                        misc = soln$info)

    output
    
}

#' @useDynLib gconcordopt
#' @export
concordfista_0_0 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                             penalize.diagonal = 0L, steptype = 0L, ...) {    

    concordfista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(steptype), ...)
    
}

#' @useDynLib gconcordopt
#' @export
concordfista_0_1 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                             penalize.diagonal = 0L, steptype = 1L, ...) {
    
    concordfista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(steptype), ...)
    
}

#' @useDynLib gconcordopt
#' @export
concordfista_0_2 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                             penalize.diagonal = 0L, steptype = 2L, ...) {
    
    concordfista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(steptype), ...)
    
}


#' @useDynLib gconcordopt
#' @export
concordfista_1_0 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                             penalize.diagonal = 1L, steptype = 0L, ...) {
    
    concordfista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(steptype), ...)
    
}

#' @useDynLib gconcordopt
#' @export
concordfista_1_1 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                             penalize.diagonal = 1L, steptype = 1L, ...) {
    
    concordfista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(steptype), ...)

}

#' @useDynLib gconcordopt
#' @export
concordfista_1_2 <- function(D, lam, pMat, X0=as(diag(ncol(D)), 'dgCMatrix'),
                             penalize.diagonal = 1L, steptype = 2L, ...) {
    
    concordfista(D, lam, pMat, X0, as.integer(penalize.diagonal), as.integer(steptype), ...)
    
}
