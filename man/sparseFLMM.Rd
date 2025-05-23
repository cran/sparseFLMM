% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/call_all_functions.R
\name{sparseFLMM}
\alias{sparseFLMM}
\title{Functional Linear Mixed Models for Irregularly or Sparsely Sampled Data}
\usage{
sparseFLMM(
  curve_info,
  use_RI = FALSE,
  use_simple = FALSE,
  method = "fREML",
  use_bam = TRUE,
  bs = "ps",
  d_grid = 100,
  min_grid = 0,
  max_grid = 1,
  my_grid = NULL,
  bf_mean = 8,
  bf_covariates = 8,
  m_mean = c(2, 3),
  covariate = FALSE,
  num_covariates,
  covariate_form,
  interaction,
  which_interaction = matrix(NA),
  save_model_mean = FALSE,
  para_estim_mean = FALSE,
  para_estim_mean_nc = 0,
  bf_covs,
  m_covs,
  use_whole = FALSE,
  use_tri = FALSE,
  use_tri_constr = TRUE,
  use_tri_constr_weights = FALSE,
  np = TRUE,
  mp = TRUE,
  use_discrete_cov = FALSE,
  para_estim_cov = FALSE,
  para_estim_cov_nc = 0,
  var_level = 0.95,
  N_B = NA,
  N_C = NA,
  N_E = NA,
  use_famm = FALSE,
  use_bam_famm = TRUE,
  bs_int_famm = list(bs = "ps", k = 8, m = c(2, 3)),
  bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3)),
  save_model_famm = FALSE,
  use_discrete_famm = FALSE,
  para_estim_famm = FALSE,
  para_estim_famm_nc = 0,
  nested = FALSE
)
}
\arguments{
\item{curve_info}{data table in which each row represents a single observation
point. \code{curve_info} needs to contain the following columns:
\itemize{
\item \code{y_vec} (numeric): the response values for each observation point
\item \code{t}  (numeric): the observations point locations, i.e. \eqn{t_{ij}}
\item \code{n_long} (integer): unique identification number for each curve
\item \code{subject_long} (integer): unique identification number for each
level of the first grouping variable
(e.g. speakers for the phonetics data in the example below).
In the case of independent functions, \code{subject_long} should be set equal to \code{n_long}.
}
For models with two crossed functional random intercepts, the data table additionally needs to have columns:
\itemize{
\item \code{word_long} (integer): unique identification number for each level of
the second grouping variable (e.g. words for the phonetics data in the example below)
\item \code{combi_long} (integer): number of the repetition of the combination of the
corresponding level of the first and of the second grouping variable.
}
For models with two nested functional random intercepts, the data table
additionally needs to have columns:
#' \itemize{
\item \code{word_long} (integer): unique identification number for each level
of the second grouping variable (e.g. phases of a randomized controled
trial). Note that the nested model is only implemented for two levels in the
second grouping variable.
\item \code{combi_long} (integer): number of the repetition of the combination of the
corresponding level of the first and of the second grouping variable.
}
For models with covariates as part of the mean function \eqn{\mu(t_{ij},x_i)},
the covariate values (numeric) need to be in
separate columns with names: \code{covariate.1}, \code{covariate.2}, etc.}

\item{use_RI}{TRUE to specify a model with one functional random intercept
for the first grouping variable (\code{subject_long}) and a smooth random error curve.
Defaults to \code{FALSE}, which specifies a model with crossed functional random
intercepts for the
first and second grouping variable and a smooth error curve.}

\item{use_simple}{\code{TRUE} to specify a model with only a smooth random error function,
\code{use_RI} should then also be set to \code{TRUE}. Defaults to \code{FALSE}.}

\item{method}{estimation method for \code{gam} or \code{bam}, see \code{\link[mgcv]{mgcv}} for more details.
Defaults to \code{"fREML"}.}

\item{use_bam}{\code{TRUE} to use function bam instead of function \code{gam} (syntax is the same, bam is faster for large data sets).
\code{bam} is recommended and set as default.}

\item{bs}{spline basis function type for the estimation of the mean function and
the auto-covariance, see \code{\link[mgcv]{s}} and \code{\link[mgcv]{te}} for more details.
Defaults to penalized B-splines, i.e. \code{bs = "ps"}. This choice is recommended as others have not been
tested yet.}

\item{d_grid}{pre-specified grid length for equidistant grid on which the mean, the auto-covariance surfaces, the eigenfunctions
and the functional random effects are evaluated. NOTE: the length of the grid can be important for computation time (approx. quadratic influence).
Defaults to \code{d_grid = 100}.}

\item{min_grid}{minimum value of equidistant grid (should approx. correspond to minimum value of time interval). Defaults to \code{min_grid = 0}.}

\item{max_grid}{maximum value of equidistant grid (should approx. correspond to maximum value of time interval). Defaults to \code{max_grid = 1}.}

\item{my_grid}{optional evaluation grid, which can be specified and used instead of \code{d_grid}, \code{min_grid}, \code{max_grid}.
NOTE: the grid should be equidistant.}

\item{bf_mean}{basis dimension (number of basis functions) used for the functional intercept
\eqn{f_0(t_{ij})} in the mean estimation via \code{bam/gam}. Defaults to \code{bf_mean = 8}.}

\item{bf_covariates}{basis dimension (number of basis functions) used for the functional effects
of covariates in the mean estimation via \code{bam/gam}. Defaults to \code{bf_covariates = 8}.
NOTE: in the current implementation, the same basis dimension for all covariates is used.}

\item{m_mean}{order of the penalty for this term in \code{bam/gam} of mean estimation, for \code{bs = "ps"} spline and penalty order,
defaults to \code{m_mean = c(2, 3)}, i.e., cubic B-splines with third order difference penalty, see \code{\link[mgcv]{s}} for details.}

\item{covariate}{\code{TRUE} to estimate covariate effects (as part of the mean function).}

\item{num_covariates}{number of covariates that are included in the model.
NOTE: not number of effects in case interactions of covariates are specified.}

\item{covariate_form}{vector with entries for each covariate that specify the form in which the
respective covariate enters the mean function. Possible forms are \code{"by"} for varying-coefficient \eqn{(f(t_{ij})*covariate)}, which is
possible for dummy coded covariates and metric covariates and
\code{"smooth"} for smooth effect in t and in covariate \eqn{(f(t_{ij}, covariate))}, which is only
possible for metric covariates! NOTE: metric covariates should be centered such that the global functional intercept \eqn{f_0(t_{ij})} can be interpreted as
global mean function and the effect can be interpreted as difference from the global mean.}

\item{interaction}{\code{TRUE} to estimate interaction effects of covariates, which interactions, see \code{which_interaction} (below).
Interactions are possible for dummy-coded covariates that act as varying coefficients.}

\item{which_interaction}{symmetric matrix that specifies which interactions should be considered in case \code{covariate = TRUE} and \code{interaction = TRUE}.
Entry \code{which_interaction}[k, l] specifies that the interaction between \code{covariate.k} and \code{covariate.l} is modeled (example below).
NOTE: entries are redundant, \code{which_interaction}[l, k] should be set to the
same as \code{which_interaction}[k, l] (symmetric). Defaults to \code{which_interaction = matrix(NA)} which should be specified when \code{interaction = FALSE}.}

\item{save_model_mean}{\code{TRUE} to give out \code{gam/bam} object (attention: can be large!), defaults to \code{FALSE}.}

\item{para_estim_mean}{\code{TRUE} to parallelize mean estimation (only possible using \code{bam}), defaults to \code{FALSE}.}

\item{para_estim_mean_nc}{number of cores for parallelization of mean estimation (only possible using \code{bam},
only active if \code{para_estim_mean = TRUE}). Defaults to 0.}

\item{bf_covs}{vector of marginal basis dimensions (number of basis functions) used for covariance estimation via \code{bam/gam}
for each functional random effect (including the smooth error curve).
In the case of multiple grouping variables, the first entry corresponds to the first grouping variable,
the second vector entry corresponds to the second grouping variable, and the third to the smooth error curve.}

\item{m_covs}{list of marginal orders of the penalty for \code{bam/gam} for covariance estimation, for \code{bs = "ps"} marginal spline and
penalty order. As only symmetric surfaces are considered: same for both directions. \cr
For crossed fRIs: list of three vectors, e.g. \code{m_covs = list(c(2, 3), c(2, 3), c(2, 3))}, where first and second entry correspond to first
and second grouping variable, respectively and third entry corresponds to smooth error.
For one fRI: list of two vectors, e.g. \code{m_covs = list(c(2, 3), c(2, 3))}, where first entry corresponds to (first) grouping variable
and second entry corresponds to smooth error.
For independent curves: list of one vector, e.g. \code{m_covs = list(c(2,3))} corresponding to smooth error.}

\item{use_whole}{\code{TRUE} to estimate the whole auto-covariance surfaces without symmetry constraint.
Defaults to \code{FALSE} as is much slower than \code{use_tri_constr} and \code{use_tri_constr_weights}. For more details, see references below.}

\item{use_tri}{\code{TRUE} to estimate only the upper triangle of the auto-covariance surfaces without
symmetry constraint. Defaults to \code{FALSE} and not recommended. For more details, see references below.}

\item{use_tri_constr}{\code{TRUE} to estimate only the upper triangle of the auto-covariance surfaces with symmetry constraint using the
smooth class \code{'symm'}. Defaults to \code{TRUE}. For more details, see references below.}

\item{use_tri_constr_weights}{\code{TRUE} to estimate only the upper triangle of the auto-covariances with symmetry constraint, using the
smooth class \code{'symm'} and weights of 0.5 on the diagonal to use the same weights as for estimating the whole auto-covariance surfaces.
Defaults to \code{FALSE}. For more details, see references below.}

\item{np}{\code{TRUE} to use 'normal parameterization' for a tensor product smooth, see \code{\link[mgcv]{te}} for more details.
Defaults to \code{TRUE}.}

\item{mp}{\code{FALSE} to use Kronecker product penalty instead of Kronecker sum penalty
with only one smoothing parameter (\code{use_whole = TRUE} and \code{use_tri = TRUE}), for details see \code{\link[mgcv]{te}}.
For \code{use_tri_constr = TRUE} and \code{use_tri_constr_weights = TRUE}, only one smoothing parameter is estimated anyway.
Defaults to \code{TRUE}.}

\item{use_discrete_cov}{\code{TRUE} to further speed up the auto-covariance computation by discretization of
covariates for storage and efficiency reasons, includes parallelization controlled by \code{para_estim_cov_nc} (below),
see \code{\link[mgcv]{bam}} for more details. Defaults to \code{FALSE}.}

\item{para_estim_cov}{\code{TRUE} to parallelize auto-covariance estimation (only possible using \code{bam}), defaults to \code{FALSE}.}

\item{para_estim_cov_nc}{number of cores (if \code{use_discrete_cov = FALSE}) or number of threads (if \code{use_discrete_cov = TRUE})
for parallelization of auto-covariance estimation (only possible using \code{bam}, only active if \code{para_estim_cov = TRUE}).
Defaults to 0.}

\item{var_level}{pre-specified level of explained variance used for the choice of the number of the functional principal
components (FPCs). Alternatively, a specific number of FPCs can be specified (see below). Defaults to \code{var_level = 0.95}.}

\item{N_B}{number of components for B (fRI for first grouping variable) to keep, overrides \code{var_level} if not \code{NA}.}

\item{N_C}{number of components for C (fRI for second grouping variable) to keep, overrides \code{var_level} if not \code{NA}.}

\item{N_E}{number of components for E (smooth error) to keep, overrides \code{var_level} if not \code{NA}.}

\item{use_famm}{\code{TRUE} to embed the model into the framework of functional additive mixed
models (FAMMs) using re-estimation of the mean function together with the prediction of the FPC weights (scores).
This allows for point-wise confidence bands for the covariate effects. Defaults to \code{FALSE}.}

\item{use_bam_famm}{\code{TRUE} to use function \code{bam} instead of function \code{gam} in FAMM estimation
(reduces computation time for large data sets),
highly recommended. Defaults to \code{TRUE}.}

\item{bs_int_famm}{specification of the estimation of the functional intercept \eqn{f_0(t_{ij})}
(as part of the mean function), see \code{\link[refund]{pffr}} for details.
Defaults to \code{bs_int = list(bs = "ps", k = 8, m = c(2, 3))}, where
\code{bs}: type of basis functions, \code{k}: number of basis functions, \code{m}: order of the spline and order of the penalty.}

\item{bs_y_famm}{specification of the estimation of the covariates effects (as part of the mean function), see \code{\link[refund]{pffr}} for details.
Defaults to \code{bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3))}, where
\code{bs}: type of basis functions, \code{k}: number of basis functions, \code{m}: order of the spline and order of the penalty.}

\item{save_model_famm}{\code{TRUE} to give out the FAMM model object (attention: can be very large!).
Defaults to FALSE.}

\item{use_discrete_famm}{\code{TRUE} to further speed up the fpc-famm computation by discretization of
# covariates for storage and efficiency reasons, includes parallelization controlled by \code{para_estim_famm_nc} (below),
see \code{\link[mgcv]{bam}} for more details. Defaults to \code{FALSE}.}

\item{para_estim_famm}{\code{TRUE} to parallelize FAMM estimation. Defaults to \code{FALSE}.}

\item{para_estim_famm_nc}{number of cores (if \code{use_discrete_famm = FALSE}) or number of
threads (if \code{use_discrete_famm = TRUE})
for parallelization of FAMM estimation (only possible using \code{bam}, only active if \code{para_estim_famm = TRUE}).
Defaults to 0.}

\item{nested}{\code{TRUE} to specify a model with nested functional random
intercepts for the first and second grouping variable and a smooth error
curve. Defaults to \code{FALSE}.}
}
\value{
The function returns a list of two elements: \code{time_all} and \code{results}. \cr
\code{time_all} contains the total system.time() for calling function \code{sparseFLMM()}.\cr
\code{results} is a list of all estimates, including:
\itemize{
\item \code{mean_hat}: includes the components of the estimated mean function.
\itemize{
\item \code{mean_pred} contains effects of dummy covariates or metric covariates with a linear effect (varying coefficients).
\item \code{mean_pred_smooth} contains effects of metric covariates with a smooth effect.
\item \code{intercept} is the estimated intercept, which is part of \eqn{f_0(t_{ij})}.
}}
For each auto-covariance smoothing alternative \code{X} (\code{use_whole}, \code{use_tri},
\code{use_tri_constr}, \code{use_tri_constr_weights}):
\itemize{
\item \code{cov_hat_X}: includes
\itemize{
\item \code{sigmasq}: the estimated error variance
\item \code{sigmasq_int}: the integral of the estimated error variance over the domain
\item \code{grid_mat_B/C/E}: the estimated auto-covariance(s) evaluated on the pre-specified grid
\item \code{sp}: the smoothing parameter(s) for smoothing the auto-covariance(s)
\item \code{time_cov_estim}: the time for the smoothing the auto-covariance(s) only
\item \code{time_cov_pred_grid}: the time for evaluating the estimated auto-covariance(s) on the pre-specified grid.
}
\item \code{time_cov_X}: the total system.time() for the auto-covariance estimation
\item \code{fpc_hat_X}: including
\itemize{
\item \code{phi_B/C/E_hat_grid}: the estimated rescaled eigenfunctions evaluated on the pre-specified grid
\item \code{nu_B/C/E_hat}: the estimated rescaled eigenvalues
\item \code{N_B/C/E}: the estimated truncation numbers, i.e., number of FPCs which are chosen
\item \code{total_var}: the estimated total variance
\item \code{var_explained}: the estimated explained variance
\item \code{xi_B/C/E_hat}: the predicted FPC weights (scores).
}
\item \code{time_fpc_X}: the total system.time() for the eigen decompositions
 and prediction on the FPC weights (scores)
If \code{use_famm = TRUE}, the list \code{results} additionally contains:
\itemize{
\item \code{fpc_famm_hat_X}: including
\itemize{
\item \code{intercept}: the estimated intercept, which is part of \eqn{f_0(t_{ij})}
\item \code{residuals}: the residuals of the FAMM estimation
\item \code{xi_B/C/E_hat_famm}: the predicted basis weights
\item \code{famm_predict_B/C/E}: the predicted functional processes evaluated on the pre-specified grid
\item \code{famm_cb_mean}: the re-estimated functional intercept \eqn{f_0(t_{ij})}
\item \code{famm_cb_covariate.1}, \code{famm_cb_covariate.1}, etc: possible re-estimated covariate effects
\item \code{famm_cb_inter_1_2}, \code{famm_cb_inter_1_3}, etc: possible interaction effects
\item \code{time_fpc_famm_X}: the total system.time() for the FAMM estimation.
}
}
The unique identification numbers for the levels of the grouping variables and curves are
renumbered for convenience during estimation from 1 in ascending order.
The original identification numbers are returned in the list \code{results}:
\itemize{
\item \code{n_orig}: curve levels as they entered the estimation
\item \code{subject_orig}: levels of the first grouping variable as they entered the estimation
\item \code{word_orig}: levels of the second grouping variable (if existent) as they entered the estimation
\item \code{my_grid}: pre-specified grid.
}
}
}
\description{
Estimation of functional linear mixed models (FLMMs) for irregularly or sparsely
sampled data based on functional principal component analysis (FPCA).
The implemented models are special cases of the general FLMM
\deqn{Y_i(t_{ij}) = \mu(t_{ij},x_i) + z_i^T U(t_{ij}) + \epsilon_i(t_{ij}), i = 1,...,n, j = 1,...,D_i,}
with \eqn{Y_i(t_{ij})} the value of the response of curve \eqn{i} at observation point
\eqn{t_{ij}}, \eqn{\mu(t_{ij},x_i)} is a mean function, which may depend on covariates
\eqn{x_i = (x_{i1},\ldots,x_{ip})^T}. \eqn{z_i} is a covariate vector,
 which is multiplied with the vector of functional random effects \eqn{U(t_{ij})}.
\eqn{\epsilon_i(t_{ij})} is independent and identically distributed white noise
measurement error with homoscedastic, constant variance. For more details, see references below.\cr \cr
The current implementation can be used to fit four special cases
of the above general FLMM:
\itemize{
\item a model for independent functional data (e.g. longitudinal data),
for which \eqn{z_i^T U(t_{ij})} only consists of
a smooth curve-specific deviation (smooth error curve)
\item a model for correlated functional data with one
functional random intercept (fRI) for one grouping variable in addition
to a smooth curve-specific error
\item a model for correlated functional data with two crossed
fRIs for two grouping variables in addition to a smooth curve-specific error
\item a model for correlated functional data with two nested fRIs for
two grouping variables in addition to a smooth curve-specific error.}
}
\details{
The code can handle irregularly and possibly sparsely sampled
data. Of course, it can also be used to analyze regular grid data,
but as it is especially designed for the irregular case and there may
be a more efficient way to analyze regular grid data. \cr \cr
The mean function is of the form \deqn{\mu(t_{ij},x_i) = f_0(t_{ij}) +
\sum_{k=1}^r f_k(t_{ij},x_{ik}),} where \eqn{f_0(t_{ij})} is a functional
intercept.
Currently implemented are effects of dummy-coded and metric covariates which act as
varying-coefficients of the
form \eqn{f_k(t_{ij})*x_{ik}} and smooth effects of metric covariates (smooth in t and in the covariate)
of the form \eqn{f(t_{ij}, x_{ik})}. NOTE: metric covariates should be centered such that the global functional intercept can be interpreted as global mean function and
the effect can be interpreted as difference from the global mean. Interaction effects of dummy-coded
covariates acting as varying coefficients are possible.
\cr
\cr
The estimation consists of four main steps:
\enumerate{
\item estimation of the smooth mean function (including covariate effects)
under independence assumption using splines.
\item estimation of the smooth auto-covariances of the functional random effects.
A fast bivariate symmetric smoother implemented in the smooth class 'symm' can be used to speed up estimation (see below).
\item eigen decomposition of the estimated auto-covariances, which are evaluated on a
pre-specified equidistant grid. This yields estimated eigenvalues and eigenfunctions, which are
rescaled to ensure orthonormality with respect to the L2-scalar product.
\item prediction of the functional principal component weights (scores) yielding predictions for
the functional random effects.
}
The estimation of the mean function and auto-covariance functions is based on package \pkg{mgcv}.
\cr
The functional principal component weights (scores) are predicted as best (linear)
unbiased predictors. In addition, this implementation allows to embed the model
in the general framework of functional additive mixed models (FAMM) based on package \pkg{refund}, which allows for the construction of
point-wise confidence bands for covariate effects (in the mean function) conditional on the FPCA.
Note that the estimation as FAMM may be computationally expensive as the model
is re-estimated in a mixed model framework.

The four special cases of the general FLMM (two nested fRIs, two
crossed fRIs, one fRI, independent curves) are implemented as follows:
\itemize{
\item In the special case with two nested fRIs, three random processes B, C,
and E are considered, where B is the fRI for the first grouping variable (e.
g. patient in a random controlled trial), C denotes the fRI for the second
grouping variable (e.g. individual specific effect in the follow-up) and E
denotes the smooth error. For this special
case, arguments \code{use_RI} and \code{use_simple} are both set to
 \code{FALSE} and argument \code{nested} is set to \code{TRUE}. Note that
this implementation only allows for a simple before/after study design.
\item In the special case with two crossed fRIs,  three
random processes B, C, and E are considered, where B is the
fRI for the first grouping variable (e.g. speakers in the phonetics example below),
C denotes the fRI for the second grouping variable
(e.g. target words in the phonetics example below) and
E denotes the smooth error. For this special case, arguments
\code{use_RI} and \code{use_simple} are both set to \code{FALSE}.
\item In the special case with only one fRI, only B and E are considered
and the number of levels for the second grouping variable is to zero.
For this special case, argument \code{use_RI} is set to \code{TRUE} and argument
\code{use_simple} is set to \code{FALSE}.
\item The special case with independent curves is internally seen as a special case
of the model with one fRI for the first grouping variable, with the number of
levels for this grouping variable corresponding to the number of curves.
Thus, for each level of the first grouping variable there is one curve.
Therefore, for the special case of independent curves, the estimation returns an
estimate for the auto-covariance of B (instead of E) and all corresponding results are indicated with \code{'_B'}, although they
correspond to the smooth error. For this special case, arguments
\code{use_RI} and \code{use_simple} are both set to \code{TRUE}.}
}
\examples{
\dontrun{
# subset of acoustic data (very small subset, no meaningful results can be expected and
# FAMM estimation does not work for this subset example. For FAMM estimation, see below.)
data("acoustic_subset")

acoustic_results <- sparseFLMM(curve_info = acoustic_subset, use_RI = FALSE, use_simple = FALSE,
              method = "fREML", use_bam = TRUE, bs = "ps", d_grid = 100, min_grid = 0,
              max_grid = 1, my_grid = NULL, bf_mean = 8, bf_covariates = 8, m_mean = c(2,3),
              covariate = TRUE, num_covariates = 4, covariate_form = rep("by", 4),
              interaction = TRUE,
              which_interaction = matrix(c(FALSE, TRUE, TRUE, TRUE, TRUE,
              FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
              FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
              byrow = TRUE, nrow = 4, ncol = 4),
              save_model_mean = FALSE, para_estim_mean = FALSE, para_estim_mean_nc = 0,
              bf_covs = c(5, 5, 5), m_covs = list(c(2, 3), c(2, 3), c(2, 3)),
              use_whole = FALSE, use_tri = FALSE, use_tri_constr = TRUE,
              use_tri_constr_weights = FALSE, np = TRUE, mp = TRUE,
              use_discrete_cov = FALSE,
              para_estim_cov = FALSE, para_estim_cov_nc = 5,
              var_level = 0.95, N_B = NA, N_C = NA, N_E = NA,
              use_famm = FALSE, use_bam_famm = TRUE,
              bs_int_famm = list(bs = "ps", k = 8, m = c(2, 3)),
              bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3)),
              save_model_famm = FALSE, use_discrete_famm = FALSE,
              para_estim_famm = FALSE, para_estim_famm_nc = 0)}

\dontrun{
# whole data set with estimation in the FAMM framework

data("acoustic")
acoustic_results <- sparseFLMM(curve_info = acoustic, use_RI = FALSE, use_simple = FALSE,
              method = "fREML", use_bam = TRUE, bs = "ps", d_grid = 100, min_grid = 0,
              max_grid = 1, my_grid = NULL, bf_mean = 8, bf_covariates = 8, m_mean = c(2,3),
              covariate = TRUE, num_covariates = 4, covariate_form = rep("by", 4),
              interaction = TRUE,
              which_interaction = matrix(c(FALSE, TRUE, TRUE, TRUE, TRUE,
              FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
              FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
              byrow = TRUE, nrow = 4, ncol = 4),
              save_model_mean = FALSE, para_estim_mean = FALSE, para_estim_mean_nc = 0,
              bf_covs = c(5, 5, 5), m_covs = list(c(2, 3), c(2, 3), c(2, 3)),
              use_whole = FALSE, use_tri = FALSE, use_tri_constr = TRUE,
              use_tri_constr_weights = FALSE, np = TRUE, mp = TRUE,
              use_discrete_cov = FALSE,
              para_estim_cov = TRUE, para_estim_cov_nc = 5,
              var_level = 0.95, N_B = NA, N_C = NA, N_E = NA,
              use_famm = TRUE, use_bam_famm = TRUE,
              bs_int_famm = list(bs = "ps", k = 8, m = c(2, 3)),
              bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3)),
              save_model_famm = FALSE, use_discrete_famm = FALSE,
              para_estim_famm = TRUE, para_estim_famm_nc = 5)}

}
\references{
Cederbaum, Pouplier, Hoole, Greven (2016): Functional Linear Mixed Models
for Irregularly or Sparsely Sampled Data. Statistical Modelling, 16(1), 67-88.

Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
Submitted on arXiv.

Scheipl, F., Staicu, A.-M. and Greven, S. (2015):
Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501.
}
\seealso{
Note that \code{\link[sparseFLMM]{sparseFLMM}} calls \code{\link[mgcv]{bam}} or \code{\link[mgcv]{gam}} directly.

For functional linear mixed models with complex correlation structures
for data sampled on equal grids based on functional principal component analysis,
see function \code{denseFLMM} in package \code{denseFLMM}.
}
\author{
Jona Cederbaum
}
\keyword{FPCA}
\keyword{models}
