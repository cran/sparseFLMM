###################################################################################
# author: Jona Cederbaum (creator) and Almond Stoecker (with thanks to Fabian Scheipl)
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
#' a (skew-)symmetry constraint on (cyclic) spline coefficients in symmetric bivariate smoothing problems and is especially
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
#' @param skew logical, should the basis be constraint to skew-symmetry instead
#' of symmetry.
#' @param cyclic.degree integer, specifying the number of basis functions identified
#' with each other at the boundaries in order to implement periodicity. Should
#' be specified to match the degree of the utilized B-spline basis.
#' @seealso \code{\link[mgcv]{smooth.construct}} and \code{\link[mgcv]{smoothCon}} for details on constructors
#' @export
#' @author Jona Cederbaum, Almond Stoecker
#' @return A basis transformation matrix of dimension \eqn{F^2 \times G} with
#' \eqn{G<F^2} depending on the specified constraint.
#' @references Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
#' Submitted on arXiv.
make_summation_matrix <- function(F, skew = FALSE, cyclic.degree = 0){
  ind_mat <- matrix(1:F^2, nrow = F, ncol = F) # index square
  pairs <- cbind(c(ind_mat), c(t(ind_mat))) # all pairs using transposed = mirror
  cons <- pairs[pairs[, 1]<pairs[, 2], , drop = FALSE] # pairs to use
  C <- diag(F^2) # initialize matrix

  # generate summation matrix for skew-symmetric / symmetric case
  if(skew) {
    C[, cons[, 1]] <- C[, cons[, 1]] - C[, cons[, 2]]
  } else {
    C[, cons[, 1]] <- C[, cons[, 1]] + C[, cons[, 2]]
  }
  if(skew) ind_vec <- cons[,1] else
    ind_vec <- pairs[pairs[, 1] <= pairs[, 2], 1, drop = FALSE]
  C <- C[, ind_vec]
  # => done in the non-cyclic case

  if(cyclic.degree>0) {
    if(F < 2*cyclic.degree) stop("For F<2*cyclic.degree not implemented, yet.")
    # helper function for mapping values of ind_mat to column idx of C
    match_ind <- function(ind)
      sapply(ind, function(x) which.max(x == ind_vec))

    # match edges
    if(F > 2*cyclic.degree) {
      edge_pairs <- cbind(
        c(ind_mat[(cyclic.degree+1):(F-cyclic.degree), 1:cyclic.degree]),
        c(t(ind_mat[F-cyclic.degree + 1:cyclic.degree, (cyclic.degree+1):(F-cyclic.degree)]))
      )
      edge_pairs <- apply(edge_pairs, 2, match_ind)
      if(skew) {
        C[, edge_pairs[, 1]] <- C[, edge_pairs[, 1]] - C[, edge_pairs[, 2]]
      } else {
        C[, edge_pairs[, 1]] <- C[, edge_pairs[, 1]] + C[, edge_pairs[, 2]]
      }
      C <- C[, -edge_pairs[, 2]]
      ind_vec <- ind_vec[-edge_pairs[, 2]]
    }
    # match corners
    corner_loc <- subset(expand.grid(row = 1:cyclic.degree,
                                     col = 1:cyclic.degree),
                         if(skew) row > col else row >= col)
    if(nrow(corner_loc)>0) {
      corner_pairs <- cbind(
        mapply(function(x,y) ind_mat[x,y], corner_loc$row + F-cyclic.degree, corner_loc$col),
        mapply(function(x,y) ind_mat[x,y], corner_loc$row, corner_loc$col),
        mapply(function(x,y) ind_mat[x,y], corner_loc$row + F-cyclic.degree, corner_loc$col + F-cyclic.degree)
      )

      corner_pairs <- matrix(apply(corner_pairs, 2, match_ind), ncol = ncol(corner_pairs))

      C[, corner_pairs[, 1]] <- C[, corner_pairs[, 1]] +
        C[, corner_pairs[, 2]] + C[, corner_pairs[, 3]]
      C <- C[, -c(corner_pairs[, 2:3]), drop = FALSE]
      ind_vec <- ind_vec[-c(corner_pairs[, 2:3])]
    }

    # symmetrize lower corner square
    lower_square <- ind_mat[F-cyclic.degree + 1:cyclic.degree,
                            1:cyclic.degree ]
    if(length(lower_square)>1) {
      lower_pairs <- cbind(c(lower_square), c(t(lower_square)))
      lower_pairs <- lower_pairs[lower_pairs[, 1] < lower_pairs[, 2], , drop = FALSE]
      lower_pairs <- matrix(apply(lower_pairs, 2, match_ind), ncol = ncol(lower_pairs))
      if(skew) {
        C[, lower_pairs[, 1]] <- C[, lower_pairs[, 1]] - C[, lower_pairs[, 2]]
      } else {
        C[, lower_pairs[, 1]] <- C[, lower_pairs[, 1]] + C[, lower_pairs[, 2]]
      }
      C <- C[,  -lower_pairs[,2], drop = FALSE]
      ind_vec <- ind_vec[-lower_pairs[,2]]
    }
    if(skew) C <- C[, -match_ind(diag(as.matrix(lower_square)))]
  }

  C
}

######################
# constructor function
######################
#' Symmetric bivariate smooths constructor
#'
#' The \code{symm} class is a smooth class that is appropriate for symmetric bivariate smooths, e.g. of covariance functions,
#' using tensor-product smooths in a \code{gam} formula. A constraint matrix is constructed
#' (see \code{\link[sparseFLMM]{make_summation_matrix}}) to impose
#' a (skew-)symmetry constraint on the (cyclic) spline coefficients,
#' which considerably reduces the number of coefficients that have to be estimated.
#'
#' @details By default a symmetric bivariate B-spline smooth \eqn{g} is specified,
#' in the sense that \eqn{g(s, t) = g(t, s)}. By setting
#' \code{s(..., bs = "symm", xt = list(skew = TRUE))}, a skew-symmetric (or anti-smmetric)
#' smooth with \eqn{g(s, t) = -g(t, s)} can be specified instead.
#' In both cases, the smooth can also be constraint to be cyclic
#' with the property \eqn{g(s, t) = g(s + c, t) = g(s, t + c)}
#' for some fixed constant \eqn{c} via specifying \code{xt = list(cyclic = TRUE)}.
#' Note that this does not correspond to specifying a tensor-product smooth from
#' cyclic marginal B-splines as given by the \code{cp}-smooth.
#' In the cyclic case, it is recommended to explicitly specify the range of the domain
#' of the smooth via the \code{knots} argument, as this determines the period and
#' often deviates from the observed range.
#'
#' The underlying procedure is the following: First, the marginal spline design matrices and the corresponding
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
#' @author Jona Cederbaum, Almond Stoecker
#' @return An object of class "symm.smooth". See \code{\link[mgcv]{smooth.construct}} for the elements it will contain.
#' @references Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
#' Submitted on arXiv.
#' @example tests/smooth.construct.symm.smooth.spec_example.R
smooth.construct.symm.smooth.spec <- function(object, data, knots){
  ##############
  # check inputs
  ##############
  if(length(object$term) > 2)
    stop("basis only handels 1D and 2D smooths")

  #############################
  # set defaults if no optional
  # arguments are given
  #############################
  if (is.null(object$xt))
    object$xt <- list(skew = FALSE, cyclic = FALSE)
  if(is.null(object$xt$skew))
    object$xt$skew <- FALSE
  if(is.null(object$xt$cyclic))
    object$xt$cyclic <- FALSE
  if(is.null(object$xt$bsmargin))
    object$xt$bsmargin <- "ps"

  # __ 1D case _____________________________________________________
  # determine design mat X, penalty mat S and Z transformation matrix

  if(length(object$term) == 1) {

    if(object$xt$cyclic)
      warning("Only 2D splines can be cross-cyclic.
              Hence, cyclic = TRUE is ignored.
              You might want to specify bsmargin = 'cp' instead
              to get cyclic B-splines.")

    # borrow form pspline smooth
    object <- smooth.construct(eval(as.call(list(as.symbol("s"),
                                                 as.symbol(object$term[1]),
                                                 bs = object$xt$bsmargin,
                                                 pc = object$point.con, xt = object$xt,
                                                 k = object$bs.dim, m = object$p.order))),
                               data = data,
                               knots = knots)
    # make (skew)-symmetric coefficient basis
    if(object$xt$skew) {
      bs.dim <- floor(object$bs.dim/2)
      Z <- rbind( diag(nrow = bs.dim),
                  if(object$bs.dim %% 2) 0,
                  - diag(nrow = bs.dim)[, bs.dim:1] )
    } else {
      bs.dim <- ceiling(object$bs.dim/2)
      Z <- rbind( diag(nrow = bs.dim),
                  diag(nrow = bs.dim)[,bs.dim:1] )
      if(object$bs.dim %% 2) Z <- Z[-bs.dim, ]
    }
    S <- object$S[[1]]
  }

  # __ 2D case _____________________________________________________
  # determine designmat X, penaltymat S and Z transformation matrix

  if(length(object$term) == 2) {

    x <- data[[object$term[1]]]
    y <- data[[object$term[2]]]

    if(length(unique(x)) < object$bs.dim)
      warning("basis dimension is larger than number of unique covariates")

    #############################
    # set defaults if no optional
    # arguments are given
    #############################
    if(is.null(object$xt))
      object$xt <- list(bsmargin = "ps", kroneckersum = TRUE) # set defaults

    if(is.null(object$xt$kroneckersum)) # if only kroneckersum is missing in xt
      object$xt$kroneckersum <- TRUE

    if (!all(sapply(object$xt$bsmargin, '%in%', c("ps", "cp"))))
      stop("marginal smooth classes need to be 'ps' or 'cp'.")
    # only allow marginal (cyclic) b-splines

    #########################
    # check input for margins
    #########################
    if(length(object$xt$bsmargin) == 1)
      object$xt$bsmargin <- rep(object$xt$bsmargin, 2)
    if (length(object$p.order) == 1) {
      m <- rep(object$p.order, 2)
    }
    else {
      m <- object$p.order
    }
    m[is.na(m)] <- 2
    object$p.order <- m
    if (object$bs.dim < 0)
      object$bs.dim <- max(10, m[1])
    nk <- object$bs.dim - m[1]
    if (nk <= 0)
      stop("basis dimension too small for b-spline order")

    #############
    # check knots
    #############
    k1 <- if(is.null(knots[[object$term[1]]]))
      knots[[object$term[2]]] else knots[[object$term[1]]]
    k2 <- knots[[object$term[2]]]
    if(!is.null(k2)) {
      if(!identical(k1, k2))
        stop("number of specified knots is not equal for both margins")
    }
    if(is.null(k1)) k1 <- range(data[object$term])
    object$knots <- list(k1, k1)
    names(object$knots) <- object$term

    Sm <- list()

    ##############################
    # build marginal design matrix
    # and marginal penalties
    ##############################
    smooth1 <- smooth.construct(eval(as.call(list(as.symbol("s"),
                                                  as.symbol(object$term[1]), bs = object$xt$bsmargin[1],
                                                  k = object$bs.dim, m = object$p.order))), data = data,
                                knots = object$knots[object$term[1]])
    smooth2 <- smooth.construct(eval(as.call(list(as.symbol("s"),
                                                  as.symbol(object$term[2]), bs = object$xt$bsmargin[2],
                                                  k = object$bs.dim, m = object$p.order))), data = data,
                                knots = object$knots[object$term[2]])
    ############################
    # build tensor product model
    # matrix and penalty matrix
    ############################
    object$X <- tensor.prod.model.matrix(X = list(smooth1$X, smooth2$X))

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

    Z <- make_summation_matrix(F = object$bs.dim,
                               skew = object$xt$skew,
                               cyclic.degree = object$xt$cyclic * (m[1]+1))
    object$margin < list()
    object$margin[[1]] <- smooth1
    object$margin[[2]] <- smooth2
    object$knots <- k1
    object$m <- m
    bs.dim <- ncol(Z)
  }

  # __ general _____________________________________________________
  # apply Z trafo and prepare and return object

  #########################
  # make symm.smooth object
  #########################

  object$X <- object$X %*% Z
  object$S <- list(crossprod(Z, S) %*% Z)
  object$Z <- Z
  # object$bs.dim <- bs.dim
  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- bs.dim - object$rank
  # no sum-to-zero constraint for skew-symm bases:
  if(object$xt$skew) object$C <- matrix(0, 0, bs.dim)

  class(object) <- "symm.smooth"
  object
}


##########################
# predict method function
##########################
# needed for functions plot.gam(), model.matrix()
# also needed when bam() is used instead of gam()
# NOTE: the object here is: gam$smooth[[i]] of class symm.smooth which
# can also be generated using smooth.construct()

#' Predict matrix method for (skew-)symmetric bivariate smooths.
#'
#' @param object is a \code{symm.smooth} object created by \code{\link{smooth.construct.symm.smooth.spec}},
#' see \code{\link[mgcv]{smooth.construct}}.
#' @param data see \code{\link[mgcv]{smooth.construct}}.
#' @seealso \code{\link[mgcv]{Predict.matrix}} and \code{\link[mgcv]{smoothCon}} for details on constructors.
#' @export
#' @author Jona Cederbaum, Almond Stoecker
Predict.matrix.symm.smooth <- function (object, data) {

  # __ 1D case _______________________________________________
  # determine design mat X and apply Z transformation matrix

  # almost identical to Predict.matrix.pspline.smooth
  # only with (skew)-symmetric basis in the end

  if(length(object$term) == 1) {
    m <- object$m[1] + 1
    ll <- object$knots[m + 1]
    ul <- object$knots[length(object$knots) - m]
    m <- m + 1
    x <- data[[object$term]]
    n <- length(x)
    ind <- x <= ul & x >= ll
    if (is.null(object$deriv))
      object$deriv <- 0
    if (sum(ind) == n) {
      X <- splines::spline.des(object$knots, x, m, rep(object$deriv,
                                                       n))$design
    }
    else {
      D <- splines::spline.des(object$knots, c(ll, ll, ul,
                                               ul), m, c(0, 1, 0, 1))$design
      X <- matrix(0, n, ncol(D))
      nin <- sum(ind)
      if (nin > 0)
        X[ind, ] <- splines::spline.des(object$knots, x[ind],
                                        m, rep(object$deriv, nin))$design
      if (object$deriv < 2) {
        ind <- x < ll
        if (sum(ind) > 0)
          X[ind, ] <- if (object$deriv == 0)
            cbind(1, x[ind] - ll) %*% D[1:2, ]
        else matrix(D[2, ], sum(ind), ncol(D), byrow = TRUE)
        ind <- x > ul
        if (sum(ind) > 0)
          X[ind, ] <- if (object$deriv == 0)
            cbind(1, x[ind] - ul) %*% D[3:4, ]
        else matrix(D[4, ], sum(ind), ncol(D), byrow = TRUE)
      }
    }
    # apply (skew-)symmetry constraint
    if(object$xt$skew) {
      bs.dim <- floor(object$bs.dim/2)
      X <- X[, 1:bs.dim] -
        X[, ncol(X)+1 - (1:bs.dim)]
    } else {
      bs.dim <- ceiling(object$bs.dim/2)
      X <- X[, 1:bs.dim] +
        X[, ncol(X)+1 - (ifelse(ncol(X)%%2, 2, 1):bs.dim)]
    }

    if (object$mono == 0)
      return(X)
    else return(X %*% object$Bs)
  }

  # __ 2D case _______________________________________________
  # determine designmat X and apply Z transformation matrix

  # almost identical to earlier version of Predict.matrix.symm.smooth
  # only also allowing for the skew-symmetric option
  # in make_summation_matrix

  if(length(object$term) == 2) {
    m <- length(object$margin)
    X <- list()
    for (i in 1:m) {
      term <- object$margin[[i]]$term
      dat <- list()
      for (j in 1:length(term)) {
        dat[[term[j]]] <- data[[term[j]]]
      }
      X[[i]] <- PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
    }
    X <- tensor.prod.model.matrix(X)
    if(is.null(object$Z)) {
      Z <- make_summation_matrix(F = object$bs.dim, skew = object$xt$skew,
                                 cyclic.degree = object$xt$cyclic *
                                   (object$m[1]+1) )
    } else {
      Z <- object$Z
    }

    X %*% Z
  }

}

###########################################################################
