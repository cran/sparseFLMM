###################################################################################
# author: Jona Cederbaum (with thanks to Fabian Scheipl)
# NOTE: this constructor builds on the wrapper function smoothCon provided
# by Simon Wood in package mgcv.
###################################################################################
# description: smooth construct class for smoothing with our symmetry constraint
# NOTE: this class is so far applicable to auto-covariances only.
# It is implemented for tensor product P-splines and
# allows for two different penalty types.
# So far, it assumes the same number and type of basis functions in each direction.
###################################################################################

######################
# underlying procedure
######################
# 1.) For each auto-covariance, we first build the marginal spline design matrices and the corresponding
# marginal difference penalties.
# 2.) The tensor product of the marginal design matrices is built and the bivariate penalty matrix is set up.
# 3.) The constraint matrix is applied to the tensor product design matrix and to the penalty matrix.

##############
# what is what
##############
# k: number of basis functions
# bsmargin: type of penalty for both directions
# m: splines and difference order
# kroneckersum: which penalty matrix should be used
## TRUE to specify a Kronecker sum penalty of the form: S1 \otimes I + I \otimes S2
## FALSE to specify a Kronecker product penalty of the form: S1 \otimes S2

###################
# constraint matrix
###################

#' Construct symmetry constraint matrix for bivariate symmetric smoothing.
#'
#' This function can be used to construct a symmetry constraint matrix that imposes
#' a symmetry constraint on spline coefficients in symmetric bivariate smoothing problems and is especially
#' designed for constructing objects of the class "symm.smooth", see \code{\link[sparseFLMM]{smooth.construct.symm.smooth.spec}}.
#'
#' @details Imposing a symmetry constraint to the spline coefficients in order to obtain a reduced coefficient vector is
#' equivalent to right multiplication of the bivariate design matrix
#' with the symmetry constraint matrix obtained with function \code{make_summation_matrix}.
#' The penalty matrix of the bivariate smooth needs to be adjusted to the reduced coefficient vector
#' by left and right multiplication with the symmetry constraint matrix.
#' This function is used in the constructor function \code{\link[sparseFLMM]{smooth.construct.symm.smooth.spec}}.
#'
#'
#' @param F number of marginal basis functions.
#' @seealso \code{\link[mgcv]{smooth.construct}} and \code{\link[mgcv]{smoothCon}} for details on constructors
#' @export
#' @return A symmetry constraint matrix of dimension \eqn{F^2 x F(F+1)/2}.
#' @references Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
#' Submitted on arXiv.
make_summation_matrix <- function(F){
  ind_mat <- matrix(1:F^2, nrow = F, ncol = F) # index square
  pairs <- cbind(c(ind_mat), c(t(ind_mat))) # all pairs using transposed = mirror
  cons <- pairs[pairs[, 1]<pairs[, 2], , drop = FALSE] # pairs to use
  C <- diag(F^2) # initialize matrix
  C[, apply(cons, 1, min)] <-C[, cons[, 1]] + C[, cons[, 2]] # add up paired columns
  C <- C[, -apply(cons, 1, max)] # remove unnecessary columns
  C
}

######################
# constructor function
######################
#' Symmetric bivariate smooths constructor
#'
#' The \code{symm} class is a new smooth class that is appropriate for symmetric bivariate smooths, e.g. of covariance functions,
#' using tensor-product smooths in a \code{gam} formula. A symmetry constraint matrix is constructed
#' (see \code{\link[sparseFLMM]{make_summation_matrix}}) to impose
#' a symmetry constraint on the spline coefficients, which considerably reduces the number of coefficients that have to be estimated.
#'
#' @details The underlying procedure is the following: First, the marginal spline design matrices and the corresponding
#' marginal difference penalties are built. Second, the tensor product of the marginal design matrices is built
#' and the bivariate penalty matrix is set up. Third, the constraint matrix is applied
#' to the tensor product design matrix and to the penalty matrix.
#'
#' @param object is a smooth specification object or a smooth object.
#' @param data a data frame, model frame or list containing the values
#'  of the (named) covariates at which the smooth term is to be evaluated.
#' @param knots an optional data frame supplying any knot locations
#'  to be supplied for basis construction.
#' @seealso \code{\link[mgcv]{smooth.construct}} and \code{\link[mgcv]{smoothCon}} for details on constructors
#' @export
#' @return An object of class "symm.smooth". See \code{\link[mgcv]{smooth.construct}} for the elements it will contain.
#' @references Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
#' Submitted on arXiv.
smooth.construct.symm.smooth.spec <- function(object, data, knots){
  ##############
  # check inputs
  ##############
  if(length(object$term) != 2) stop("basis only handels 2D smooths") # check if two marginal smooths

  x <- data[[object$term[1]]]
  y <- data[[object$term[2]]]

  if(length(unique(x)) < object$bs.dim) warning("basis dimension is larger than number of unique covariates")

  #############################
  # set defaults if no optional
  # arguments are given
  #############################
  if(is.null(object$xt))
    object$xt <- list(bsmargin = "ps", kroneckersum = TRUE) # set defaults

  if(is.null(object$xt$kroneckersum)) # if only kroneckersum is missing in xt
    object$xt$kroneckersum <- TRUE

  if(is.null(object$xt$bsmargin))  # if only bsmargin is missing in xt
    object$xt$bsmargin <- "ps"

  if(object$xt$bsmargin != "ps") stop("marginal smooth class need to be 'ps'") # only allow marginal b-splines

  #########################
  # check input for margins
  #########################
  if (length(object$p.order) == 1){ # if e.g. m = c(1) -> m = c(1, 1)
    m <- rep(object$p.order, 2)
  }else{
    m <- object$p.order  # m[1] - basis order, m[2] - penalty order
  }
  m[is.na(m)] <- 2 # default if object$p.order is missing -> m = c(2, 2)
  object$p.order <- m

  if (object$bs.dim<0) object$bs.dim <- max(10, m[1]) # default
  nk <- object$bs.dim - m[1] # number of interior knots
  if (nk <= 0) stop("basis dimension too small for b-spline order")

  #############
  # check knots
  #############
  k1 <- knots[[object$term[1]]]
  k2 <- knots[[object$term[2]]]
  if(!is.null(k1) & !is.null(k2)){
    if((k1 != k2)) stop("number of specified knots is not equal for both margins")
  }

  Sm <- list()

  ##############################
  # build marginal design matrix
  # and marginal penalties
  ##############################
  smooth1 <- smooth.construct(eval(as.call(list(as.symbol("s"), as.symbol(object$term[1]), bs = object$xt$bsmargin,
                                              k = object$bs.dim, m = object$p.order))), data = data, knots = knots)

  smooth2 <- smooth.construct(eval(as.call(list(as.symbol("s"), as.symbol(object$term[2]), bs = object$xt$bsmargin,
                                              k = object$bs.dim, m = object$p.order))), data = data, knots = knots)
  ############################
  # build tensor product model
  # matrix and penalty matrix
  ############################
  X <- tensor.prod.model.matrix(X = list(smooth1$X, smooth2$X))

  Sm[[1]] <- smooth1$S[[1]]
  Sm[[2]] <- smooth2$S[[1]]

  if(object$xt$kroneckersum){
    S <- tensor.prod.penalties(list(Sm[[1]], Sm[[2]]))
    S <- S[[1]] + S[[2]]
  } else{
    S <- Sm[[1]]%x%Sm[[2]]
  }

  ################################################
  # constraint equal coefficients by summation
  # of columns of X and adaption of penalty matrix
  ################################################

  Z <- make_summation_matrix(F = object$bs.dim)

  X_tri <- X %*% Z
  S_tri <- t(Z) %*% S %*% Z

  # rank and null space dimension of penalty matrix
  r <- qr(S_tri)$rank
  nsd <- nrow(S_tri)-r

  #########################
  # make symm.smooth object
  #########################
  object$S <- list(S_tri) # penalty
  object$X <- X_tri  # design matrix
  object$rank <- r
  object$null.space.dim <- nsd # dimension of unpenalized space
  object$m <- m # store p-splines specific info
  object$knots <- k1

  object$margin<list()
  object$margin[[1]] <- smooth1
  object$margin[[2]] <- smooth2

  class(object) <- "symm.smooth" # gives object a class
  object
}


##########################
# predict method function
##########################
# needed for functions plot.gam(), model.matrix()
# also needed when bam() is used instead of gam()
# NOTE: the object here is: gam$smooth[[i]] of class symm.smooth which
# can also be generated using smooth.construct()

#' Predict matrix method for symmetric bivariate smooths.
#'
#' @param object is a \code{symm.smooth} object created by \code{\link{smooth.construct.symm.smooth.spec}},
#' see \code{\link[mgcv]{smooth.construct}}.
#' @param data see \code{\link[mgcv]{smooth.construct}}.
#' @seealso \code{\link[mgcv]{Predict.matrix}} and \code{\link[mgcv]{smoothCon}} for details on constructors.
#' @export
Predict.matrix.symm.smooth <- function(object, data){
  m <- length(object$margin)
  X <- list()
  for(i in 1:m){
    term <- object$margin[[i]]$term
    dat <- list()
    for(j in 1:length(term)){
      dat[[term[j]]] <- data[[term[j]]]
    }
    X[[i]] <- PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
  }
  X <- tensor.prod.model.matrix(X)

  ############################################
  # constraint equal coefficients by summation
  # of columns of X and adaption of penalty
  ############################################
  Z <- make_summation_matrix(F = object$bs.dim)
  X_tri <- X %*% Z
  X_tri
}

###########################################################################
