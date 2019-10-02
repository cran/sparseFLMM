#' Functional Linear Mixed Models for Irregularly or Sparsely Sampled Data
#'
#' Estimation of functional linear mixed models (FLMMs) for irregularly or sparsely
#' sampled data based on functional principal component analysis (FPCA).
#' The implemented models are special cases of the general FLMM
#' \deqn{Y_i(t_{ij}) = \mu(t_{ij},x_i) + z_i^T U(t_{ij}) + \epsilon_i(t_{ij}), i = 1,...,n, j = 1,...,D_i,}
#' with \eqn{Y_i(t_{ij})} the value of the response of curve \eqn{i} at observation point
#' \eqn{t_{ij}}, \eqn{\mu(t_{ij},x_i)} is a mean function, which may depend on covariates
#' \eqn{x_i = (x_{i1},\ldots,x_{ip})^T}. \eqn{z_i} is a covariate vector,
#'  which is multiplied with the vector of functional random effects \eqn{U(t_{ij})}.
#' \eqn{\epsilon_i(t_{ij})} is independent and identically distributed white noise
#' measurement error with homoscedastic, constant variance. For more details, see references below.\cr \cr
#' The current implementation can be used to fit three special cases
#' of the above general FLMM:
#' \itemize{
#' \item a model for independent functional data (e.g. longitudinal data),
#' for which \eqn{z_i^T U(t_{ij})} only consists of
#' a smooth curve-specific deviation (smooth error curve)
#' \item a model for correlated functional data with one
#' functional random intercept (fRI) for one grouping variable in addition
#' to a smooth curve-specific error
#' \item a model for correlated functional data with two crossed
#' fRIs for two grouping variables in addition to a smooth curve-specific error.}
#'
#' The code can handle irregularly and possibly sparsely sampled
#' data. Of course, it can also be used to analyze regular grid data,
#' but as it is especially designed for the irregular case and there may
#' be a more efficient way to analyze regular grid data. \cr \cr
#' The mean function is of the form \deqn{\mu(t_{ij},x_i) = f_0(t_{ij}) +
#' \sum_{k=1}^r f_k(t_{ij},x_{ik}),} where \eqn{f_0(t_{ij})} is a functional
#' intercept.
#' Currently implemented are effects of dummy-coded and metric covariates which act as
#' varying-coefficients of the
#' form \eqn{f_k(t_{ij})*x_{ik}} and smooth effects of metric covariates (smooth in t and in the covariate)
#' of the form \eqn{f(t_{ij}, x_{ik})}. NOTE: metric covariates should be centered such that the global functional intercept can be interpreted as global mean function and
#' the effect can be interpreted as difference from the global mean. Interaction effects of dummy-coded
#' covariates acting as varying coefficients are possible.
#' \cr
#' \cr
#' The estimation consists of four main steps:
#' \enumerate{
#' \item estimation of the smooth mean function (including covariate effects)
#' under independence assumption using splines.
#' \item estimation of the smooth auto-covariances of the functional random effects.
#' A fast bivariate symmetric smoother implemented in the smooth class 'symm' can be used to speed up estimation (see below).
#' \item eigen decomposition of the estimated auto-covariances, which are evaluated on a
#' pre-specified equidistant grid. This yields estimated eigenvalues and eigenfunctions, which are
#' rescaled to ensure orthonormality with respect to the L2-scalar product.
#' \item prediction of the functional principal component weights (scores) yielding predictions for
#' the functional random effects.
#' }
#' The estimation of the mean function and auto-covariance functions is based on package \pkg{mgcv}.
#' \cr
#' The functional principal component weights (scores) are predicted as best (linear)
#' unbiased predictors. In addition, this implementation allows to embed the model
#' in the general framework of functional additive mixed models (FAMM) based on package \pkg{refund}, which allows for the construction of
#' point-wise confidence bands for covariate effects (in the mean function) conditional on the FPCA.
#' Note that the estimation as FAMM may be computationally expensive as the model
#' is re-estimated in a mixed model framework.
#'
#'
#' @param curve_info data table in which each row represents a single observation
#' point. \code{curve_info} needs to contain the following columns:
#' \itemize{
#' \item \code{y_vec} (numeric): the response values for each observation point
#' \item \code{t}  (numeric): the observations point locations, i.e. \eqn{t_{ij}}
#' \item \code{n_long} (integer): unique identification number for each curve
#' \item \code{subject_long} (integer): unique identification number for each
#' level of the first grouping variable
#' (e.g. speakers for the phonetics data in the example below).
#' In the case of independent functions, \code{subject_long} should be set equal to \code{n_long}.
#' }
#' For models with two crossed functional random intercepts, the data table additionally needs to have columns:
#' \itemize{
#' \item \code{word_long} (integer): unique identification number for each level of
#' the second grouping variable (e.g. words for the phonetics data in the example below)
#' \item \code{combi_long} (integer): number of the repetition of the combination of the
#' corresponding level of the first and of the second grouping variable.
#' }
#' For models with covariates as part of the mean function \eqn{\mu(t_{ij},x_i)},
#' the covariate values (numeric) need to be in
#' separate columns with names: \code{covariate.1}, \code{covariate.2}, etc.
#' @param use_RI TRUE to specify a model with one functional random intercept
#' for the first grouping variable (\code{subject_long}) and a smooth random error curve.
#' Defaults to \code{FALSE}, which specifies a model with crossed functional random
#' intercepts for the
#' first and second grouping variable and a smooth error curve.
#' @param use_simple \code{TRUE} to specify a model with only a smooth random error function,
#' \code{use_RI} should then also be set to \code{TRUE}. Defaults to \code{FALSE}.
#' @param method estimation method for \code{gam} or \code{bam}, see \code{\link{mgcv}} for more details.
#' Defaults to \code{"fREML"}.
#' @param use_bam \code{TRUE} to use function bam instead of function \code{gam} (syntax is the same, bam is faster for large data sets).
#' \code{bam} is recommended and set as default.
#' @param bs spline basis function type for the estimation of the mean function and
#' the auto-covariance, see \code{\link{s}} and \code{\link{te}} for more details.
#' Defaults to penalized B-splines, i.e. \code{bs = "ps"}. This choice is recommended as others have not been
#' tested yet.
#' @param d_grid pre-specified grid length for equidistant grid on which the mean, the auto-covariance surfaces, the eigenfunctions
#' and the functional random effects are evaluated. NOTE: the length of the grid can be important for computation time (approx. quadratic influence).
#' Defaults to \code{d_grid = 100}.
#' @param min_grid         minimum value of equidistant grid (should approx. correspond to minimum value of time interval). Defaults to \code{min_grid = 0}.
#' @param max_grid         maximum value of equidistant grid (should approx. correspond to maximum value of time interval). Defaults to \code{max_grid = 1}.
#' @param my_grid          optional evaluation grid, which can be specified and used instead of \code{d_grid}, \code{min_grid}, \code{max_grid}.
#' NOTE: the grid should be equidistant.
#' @param bf_mean          basis dimension (number of basis functions) used for the functional intercept
#' \eqn{f_0(t_{ij})} in the mean estimation via \code{bam/gam}. Defaults to \code{bf_mean = 8}.
#' @param bf_covariates   basis dimension (number of basis functions) used for the functional effects
#' of covariates in the mean estimation via \code{bam/gam}. Defaults to \code{bf_covariates = 8}.
#' NOTE: in the current implementation, the same basis dimension for all covariates is used.
#' @param m_mean           order of the penalty for this term in \code{bam/gam} of mean estimation, for \code{bs = "ps"} spline and penalty order,
#' defaults to \code{m_mean = c(2, 3)}, i.e., cubic B-splines with third order difference penalty, see \code{\link{s}} for details.
#' @param covariate        \code{TRUE} to estimate covariate effects (as part of the mean function).
#' @param num_covariates   number of covariates that are included in the model.
#' NOTE: not number of effects in case interactions of covariates are specified.
#' @param covariate_form   vector with entries for each covariate that specify the form in which the
#' respective covariate enters the mean function. Possible forms are \code{"by"} for varying-coefficient \eqn{(f(t_{ij})*covariate)}, which is
#' possible for dummy coded covariates and metric covariates and
#' \code{"smooth"} for smooth effect in t and in covariate \eqn{(f(t_{ij}, covariate))}, which is only
#' possible for metric covariates! NOTE: metric covariates should be centered such that the global functional intercept \eqn{f_0(t_{ij})} can be interpreted as
#' global mean function and the effect can be interpreted as difference from the global mean.
#' @param interaction      \code{TRUE} to estimate interaction effects of covariates, which interactions, see \code{which_interaction} (below).
#' Interactions are possible for dummy-coded covariates that act as varying coefficients.
#' @param which_interaction  symmetric matrix that specifies which interactions should be considered in case \code{covariate = TRUE} and \code{interaction = TRUE}.
#' Entry \code{which_interaction}[k, l] specifies that the interaction between \code{covariate.k} and \code{covariate.l} is modeled (example below).
#' NOTE: entries are redundant, \code{which_interaction}[l, k] should be set to the
#' same as \code{which_interaction}[k, l] (symmetric). Defaults to \code{which_interaction = matrix(NA)} which should be specified when \code{interaction = FALSE}.
#' @param save_model_mean \code{TRUE} to give out \code{gam/bam} object (attention: can be large!), defaults to \code{FALSE}.
#' @param para_estim_mean \code{TRUE} to parallelize mean estimation (only possible using \code{bam}), defaults to \code{FALSE}.
#' @param para_estim_mean_nc  number of cores for parallelization of mean estimation (only possible using \code{bam},
#' only active if \code{para_estim_mean = TRUE}). Defaults to 0.
#' @param bf_covs  vector of marginal basis dimensions (number of basis functions) used for covariance estimation via \code{bam/gam}
#' for each functional random effect (including the smooth error curve).
#' In the case of multiple grouping variables, the first entry corresponds to the first grouping variable,
#' the second vector entry corresponds to the second grouping variable, and the third to the smooth error curve.
#' @param m_covs list of marginal orders of the penalty for \code{bam/gam} for covariance estimation, for \code{bs = "ps"} marginal spline and
#' penalty order. As only symmetric surfaces are considered: same for both directions. \cr
#' For crossed fRIs: list of three vectors, e.g. \code{m_covs = list(c(2, 3), c(2, 3), c(2, 3))}, where first and second entry correspond to first
#' and second grouping variable, respectively and third entry corresponds to smooth error.
#' For one fRI: list of two vectors, e.g. \code{m_covs = list(c(2, 3), c(2, 3))}, where first entry corresponds to (first) grouping variable
#' and second entry corresponds to smooth error.
#' For independent curves: list of one vector, e.g. \code{m_covs = list(c(2,3))} corresponding to smooth error.
#' @param use_whole \code{TRUE} to estimate the whole auto-covariance surfaces without symmetry constraint.
#' Defaults to \code{FALSE} as is much slower than \code{use_tri_constr} and \code{use_tri_constr_weights}. For more details, see references below.
#' @param use_tri \code{TRUE} to estimate only the upper triangle of the auto-covariance surfaces without
#' symmetry constraint. Defaults to \code{FALSE} and not recommended. For more details, see references below.
#' @param use_tri_constr \code{TRUE} to estimate only the upper triangle of the auto-covariance surfaces with symmetry constraint using the
#' smooth class \code{'symm'}. Defaults to \code{TRUE}. For more details, see references below.
#' @param use_tri_constr_weights \code{TRUE} to estimate only the upper triangle of the auto-covariances with symmetry constraint, using the
#' smooth class \code{'symm'} and weights of 0.5 on the diagonal to use the same weights as for estimating the whole auto-covariance surfaces.
#' Defaults to \code{FALSE}. For more details, see references below.
#' @param np \code{TRUE} to use 'normal parameterization' for a tensor product smooth, see \code{\link{te}} for more details.
#' Defaults to \code{TRUE}.
#' @param mp \code{FALSE} to use Kronecker product penalty instead of Kronecker sum penalty
#' with only one smoothing parameter (\code{use_whole = TRUE} and \code{use_tri = TRUE}), for details see \code{\link{te}}.
#' For \code{use_tri_constr = TRUE} and \code{use_tri_constr_weights = TRUE}, only one smoothing parameter is estimated anyway.
#' Defaults to \code{TRUE}.
#' @param use_discrete_cov \code{TRUE} to further speed up the auto-covariance computation by discretization of
#' covariates for storage and efficiency reasons, includes parallelization controlled by \code{para_estim_cov_nc} (below),
#' see \code{\link{bam}} for more details. Defaults to \code{FALSE}.
#' @param para_estim_cov \code{TRUE} to parallelize auto-covariance estimation (only possible using \code{bam}), defaults to \code{FALSE}.
#' @param para_estim_cov_nc number of cores (if \code{use_discrete_cov = FALSE}) or number of threads (if \code{use_discrete_cov = TRUE})
#' for parallelization of auto-covariance estimation (only possible using \code{bam}, only active if \code{para_estim_cov = TRUE}).
#' Defaults to 0.
#' @param var_level  pre-specified level of explained variance used for the choice of the number of the functional principal
#' components (FPCs). Alternatively, a specific number of FPCs can be specified (see below). Defaults to \code{var_level = 0.95}.
#' @param N_B  number of components for B (fRI for first grouping variable) to keep, overrides \code{var_level} if not \code{NA}.
#' @param N_C  number of components for C (fRI for second grouping variable) to keep, overrides \code{var_level} if not \code{NA}.
#' @param N_E  number of components for E (smooth error) to keep, overrides \code{var_level} if not \code{NA}.
#' @param use_famm \code{TRUE} to embed the model into the framework of functional additive mixed
#' models (FAMMs) using re-estimation of the mean function together with the prediction of the FPC weights (scores).
#' This allows for point-wise confidence bands for the covariate effects. Defaults to \code{FALSE}.
#' @param use_bam_famm \code{TRUE} to use function \code{bam} instead of function \code{gam} in FAMM estimation
#' (reduces computation time for large data sets),
#' highly recommended. Defaults to \code{TRUE}.
#' @param bs_int_famm  specification of the estimation of the functional intercept \eqn{f_0(t_{ij})}
#' (as part of the mean function), see \code{\link{pffr}} for details.
#' Defaults to \code{bs_int = list(bs = "ps", k = 8, m = c(2, 3))}, where
#' \code{bs}: type of basis functions, \code{k}: number of basis functions, \code{m}: order of the spline and order of the penalty.
#' @param bs_y_famm  specification of the estimation of the covariates effects (as part of the mean function), see \code{\link{pffr}} for details.
#' Defaults to \code{bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3))}, where
#' \code{bs}: type of basis functions, \code{k}: number of basis functions, \code{m}: order of the spline and order of the penalty.
#' @param save_model_famm  \code{TRUE} to give out the FAMM model object (attention: can be very large!).
#' Defaults to FALSE.
#' @param use_discrete_famm  \code{TRUE} to further speed up the fpc-famm computation by discretization of
#' # covariates for storage and efficiency reasons, includes parallelization controlled by \code{para_estim_famm_nc} (below),
#' see \code{\link{bam}} for more details. Defaults to \code{FALSE}.
#' @param para_estim_famm \code{TRUE} to parallelize FAMM estimation. Defaults to \code{FALSE}.
#' @param para_estim_famm_nc  number of cores (if \code{use_discrete_famm = FALSE}) or number of
#' threads (if \code{use_discrete_famm = TRUE})
#' for parallelization of FAMM estimation (only possible using \code{bam}, only active if \code{para_estim_famm = TRUE}).
#' Defaults to 0.
#'
#' @details The three special cases of the general FLMM (two
#' crossed fRIs, one fRI, independent curves) are implemented as follows:
#' \itemize{
#' \item In the special case with two crossed fRIs,  three
#' random processes B, C, and E are considered, where B is the
#' fRI for the first grouping variable (e.g. speakers in the phonetics example below),
#' C denotes the fRI for the second grouping variable
#' (e.g. target words in the phonetics example below) and
#' E denotes the smooth error. For this special case, arguments
#' \code{use_RI} and \code{use_simple} are both set to \code{FALSE}.
#' \item In the special case with only one fRI, only B and E are considered
#' and the number of levels for the second grouping variable is to zero.
#' For this special case, argument \code{use_RI} is set to \code{TRUE} and argument
#' \code{use_simple} is set to \code{FALSE}.
#' \item The special case with independent curves is internally seen as a special case
#' of the model with one fRI for the first grouping variable, with the number of
#' levels for this grouping variable corresponding to the number of curves.
#' Thus, for each level of the first grouping variable there is one curve.
#' Therefore, for the special case of independent curves, the estimation returns an
#' estimate for the auto-covariance of B (instead of E) and all corresponding results are indicated with \code{'_B'}, although they
#' correspond to the smooth error. For this special case, arguments
#' \code{use_RI} and \code{use_simple} are both set to \code{TRUE}.}
#'
#' @return The function returns a list of two elements: \code{time_all} and \code{results}. \cr
#' \code{time_all} contains the total system.time() for calling function \code{sparseFLMM()}.\cr
#' \code{results} is a list of all estimates, including:
#' \itemize{
#' \item \code{mean_hat}: includes the components of the estimated mean function.
#' \itemize{
#' \item \code{mean_pred} contains effects of dummy covariates or metric covariates with a linear effect (varying coefficients).
#' \item \code{mean_pred_smooth} contains effects of metric covariates with a smooth effect.
#' \item \code{intercept} is the estimated intercept, which is part of \eqn{f_0(t_{ij})}.
#' }}
#' For each auto-covariance smoothing alternative \code{X} (\code{use_whole}, \code{use_tri},
#' \code{use_tri_constr}, \code{use_tri_constr_weights}):
#' \itemize{
#' \item \code{cov_hat_X}: includes
#' \itemize{
#' \item \code{sigmasq}: the estimated error variance
#' \item \code{sigmasq_int}: the integral of the estimated error variance over the domain
#' \item \code{grid_mat_B/C/E}: the estimated auto-covariance(s) evaluated on the pre-specified grid
#' \item \code{sp}: the smoothing parameter(s) for smoothing the auto-covariance(s)
#' \item \code{time_cov_estim}: the time for the smoothing the auto-covariance(s) only
#' \item \code{time_cov_pred_grid}: the time for evaluating the estimated auto-covariance(s) on the pre-specified grid.
#' }
#' \item \code{time_cov_X}: the total system.time() for the auto-covariance estimation
#' \item \code{fpc_hat_X}: including
#' \itemize{
#' \item \code{phi_B/C/E_hat_grid}: the estimated rescaled eigenfunctions evaluated on the pre-specified grid
#' \item \code{nu_B/C/E_hat}: the estimated rescaled eigenvalues
#' \item \code{N_B/C/E}: the estimated truncation numbers, i.e., number of FPCs which are chosen
#' \item \code{total_var}: the estimated total variance
#' \item \code{var_explained}: the estimated explained variance
#' \item \code{xi_B/C/E_hat}: the predicted FPC weights (scores).
#' }
#' \item \code{time_fpc_X}: the total system.time() for the eigen decompositions
#'  and prediction on the FPC weights (scores)
#' If \code{use_famm = TRUE}, the list \code{results} additionally contains:
#' \itemize{
#' \item \code{fpc_famm_hat_X}: including
#' \itemize{
#' \item \code{intercept}: the estimated intercept, which is part of \eqn{f_0(t_{ij})}
#' \item \code{residuals}: the residuals of the FAMM estimation
#' \item \code{xi_B/C/E_hat_famm}: the predicted basis weights
#' \item \code{famm_predict_B/C/E}: the predicted functional processes evaluated on the pre-specified grid
#' \item \code{famm_cb_mean}: the re-estimated functional intercept \eqn{f_0(t_{ij})}
#' \item \code{famm_cb_covariate.1}, \code{famm_cb_covariate.1}, etc: possible re-estimated covariate effects
#' \item \code{famm_cb_inter_1_2}, \code{famm_cb_inter_1_3}, etc: possible interaction effects
#' \item \code{time_fpc_famm_X}: the total system.time() for the FAMM estimation.
#' }
#' }
#' The unique identification numbers for the levels of the grouping variables and curves are
#' renumbered for convenience during estimation from 1 in ascending order.
#' The original identification numbers are returned in the list \code{results}:
#' \itemize{
#' \item \code{n_orig}: curve levels as they entered the estimation
#' \item \code{subject_orig}: levels of the first grouping variable as they entered the estimation
#' \item \code{word_orig}: levels of the second grouping variable (if existent) as they entered the estimation
#' \item \code{my_grid}: pre-specified grid.
#' }
#'}
#'
#' @author Jona Cederbaum
#'
#' @seealso Note that \code{\link{sparseFLMM}} calls \code{\link[mgcv]{bam}} or \code{\link[mgcv]{gam}} directly.
#' @seealso For functional linear mixed models with complex correlation structures
#' for data sampled on equal grids based on functional principal component analysis,
#' see function \code{denseFLMM} in package \code{denseFLMM}.
#'
#' @keywords models FPCA
#'
#' @references
#' Cederbaum, Pouplier, Hoole, Greven (2016): Functional Linear Mixed Models
#' for Irregularly or Sparsely Sampled Data. Statistical Modelling, 16(1), 67-88.
#'
#' Cederbaum, Scheipl, Greven (2016): Fast symmetric additive covariance smoothing.
#' Submitted on arXiv.
#'
#' Scheipl, F., Staicu, A.-M. and Greven, S. (2015):
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501.
#'
#' @examples
#' \dontrun{
#' # subset of acoustic data (very small subset, no meaningful results can be expected and
#' # FAMM estimation does not work for this subset example. For FAMM estimation, see below.)
#' data("acoustic_subset")
#'
#' acoustic_results <- sparseFLMM(curve_info = acoustic_subset, use_RI = FALSE, use_simple = FALSE,
#'               method = "fREML", use_bam = TRUE, bs = "ps", d_grid = 100, min_grid = 0,
#'               max_grid = 1, my_grid = NULL, bf_mean = 8, bf_covariates = 8, m_mean = c(2,3),
#'               covariate = TRUE, num_covariates = 4, covariate_form = rep("by", 4),
#'               interaction = TRUE,
#'               which_interaction = matrix(c(FALSE, TRUE, TRUE, TRUE, TRUE,
#'               FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
#'               FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#'               byrow = TRUE, nrow = 4, ncol = 4),
#'               save_model_mean = FALSE, para_estim_mean = FALSE, para_estim_mean_nc = 0,
#'               bf_covs = c(5, 5, 5), m_covs = list(c(2, 3), c(2, 3), c(2, 3)),
#'               use_whole = FALSE, use_tri = FALSE, use_tri_constr = TRUE,
#'               use_tri_constr_weights = FALSE, np = TRUE, mp = TRUE,
#'               use_discrete_cov = FALSE,
#'               para_estim_cov = FALSE, para_estim_cov_nc = 5,
#'               var_level = 0.95, N_B = NA, N_C = NA, N_E = NA,
#'               use_famm = FALSE, use_bam_famm = TRUE,
#'               bs_int_famm = list(bs = "ps", k = 8, m = c(2, 3)),
#'               bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3)),
#'               save_model_famm = FALSE, use_discrete_famm = FALSE,
#'               para_estim_famm = FALSE, para_estim_famm_nc = 0)}
#'
#'\dontrun{
#'# whole data set with estimation in the FAMM framework
#'
#'data("acoustic")
#'acoustic_results <- sparseFLMM(curve_info = acoustic, use_RI = FALSE, use_simple = FALSE,
#'               method = "fREML", use_bam = TRUE, bs = "ps", d_grid = 100, min_grid = 0,
#'               max_grid = 1, my_grid = NULL, bf_mean = 8, bf_covariates = 8, m_mean = c(2,3),
#'               covariate = TRUE, num_covariates = 4, covariate_form = rep("by", 4),
#'               interaction = TRUE,
#'               which_interaction = matrix(c(FALSE, TRUE, TRUE, TRUE, TRUE,
#'               FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
#'               FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#'               byrow = TRUE, nrow = 4, ncol = 4),
#'               save_model_mean = FALSE, para_estim_mean = FALSE, para_estim_mean_nc = 0,
#'               bf_covs = c(5, 5, 5), m_covs = list(c(2, 3), c(2, 3), c(2, 3)),
#'               use_whole = FALSE, use_tri = FALSE, use_tri_constr = TRUE,
#'               use_tri_constr_weights = FALSE, np = TRUE, mp = TRUE,
#'               use_discrete_cov = FALSE,
#'               para_estim_cov = TRUE, para_estim_cov_nc = 5,
#'               var_level = 0.95, N_B = NA, N_C = NA, N_E = NA,
#'               use_famm = TRUE, use_bam_famm = TRUE,
#'               bs_int_famm = list(bs = "ps", k = 8, m = c(2, 3)),
#'               bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3)),
#'               save_model_famm = FALSE, use_discrete_famm = FALSE,
#'               para_estim_famm = TRUE, para_estim_famm_nc = 5)}
#'
#' @export
#' @import methods parallel mgcv MASS Matrix data.table refund
#' @importFrom grDevices terrain.colors
#' @importFrom stats approx as.formula coef coefficients fitted predict reshape var weights
#' @importFrom utils packageVersion
sparseFLMM <- function(curve_info, use_RI = FALSE, use_simple = FALSE, method = "fREML",
                               use_bam = TRUE, bs = "ps", d_grid = 100, min_grid = 0,
                               max_grid = 1, my_grid = NULL, bf_mean = 8,
                               bf_covariates = 8, m_mean = c(2,3), covariate = FALSE,
                               num_covariates, covariate_form, interaction, which_interaction = matrix(NA),
                               save_model_mean = FALSE, para_estim_mean = FALSE, para_estim_mean_nc = 0,
                               bf_covs, m_covs, use_whole = FALSE, use_tri = FALSE, use_tri_constr = TRUE,
                               use_tri_constr_weights = FALSE, np = TRUE, mp = TRUE,
                               use_discrete_cov = FALSE,
                               para_estim_cov = FALSE, para_estim_cov_nc = 0 ,
                               var_level = 0.95,
                               N_B = NA, N_C = NA, N_E = NA,
                               use_famm = FALSE, use_bam_famm = TRUE,
                               bs_int_famm = list(bs = "ps", k = 8, m = c(2, 3)),
                               bs_y_famm = list(bs = "ps", k = 8, m = c(2, 3)),
                               save_model_famm = FALSE, use_discrete_famm = FALSE,
                               para_estim_famm = FALSE, para_estim_famm_nc = 0){

  ##############################################################################
  # preparations
  ##############################################################################

  #####################
  # checks and warnings
  #####################
  if(packageVersion(pkg = "data.table") == "1.9.4")
    stop("you are using data.table version 1.9.4 which cannot be used due to internal bugs,
         please install a different version")

  if(packageVersion(pkg = "mgcv") == "1.8-12")
    warning("you are using mgcv version 1.8-12 for which errors might occur.
            In case of problems, please use a different version")

  if(!is.data.table(curve_info))
    stop("curve_info needs to be a data.table, see ?data.table")

  if(bs != "ps")
    stop("so far only P-splines are allowed")

  if(is.na(var_level) & (is.na(N_B) | is.na(N_E)))
    stop("either var_level or N_B, N_C, and N_E have to be specified")

  if(var_level > 1 & (is.na(N_B) | is.na(N_C) | is.na(N_E)))
    stop("var_level has to be smaller than 1")

  if(covariate & (length(covariate_form) != num_covariates))
    stop("covariate_form has to be of length ", num_covariates)

  if(covariate & interaction & ((nrow(which_interaction) != num_covariates) | (ncol(which_interaction) != num_covariates)))
    stop("which_interaction has to be a matrix of dimension ", num_covariates, "x", num_covariates)

  if(covariate & interaction & !isSymmetric(which_interaction))
    stop("which_interaction needs to be symmetric")

  if(covariate & (length(bf_mean) != 1))
    stop(" bf_mean needs to be a number vector of length 1")

  if(covariate & (length(bf_covariates) != 1))
    stop("so far, the same number of basis functions is used for all covariate and interaction effects,
         please set bf_covariates to a numeric vector of length 1")

  if(!use_RI & (length(bf_covs) != 3))
    stop("for use_RI == FALSE, bf_covs has to be of length 3")

  if(use_RI & (!use_simple) & ((length(bf_covs) != 2)))
    stop("for use_RI == TRUE and use_simple == FALSE, bf_covs has to be of length 2")

  if(use_simple & ((length(bf_covs) != 1)))
    stop("for use_simple == TRUE, bf_covs has to be of length 1")

  if(!is.list(m_covs))
    stop("m_covs needs to be a list")

  if(!use_RI & ((length(m_covs) != 3)))
    stop(" for use_RI == FALSE, m_covs needs to be list of length 3")

  if(use_RI & (!use_simple) & ((length(m_covs) != 2)))
    stop(" for use_RI == TRUE and use_simple == FALSE, m_covs needs to be list of length 2")

  if(use_simple & (length(m_covs) != 1))
    stop(" for use_simple == TRUE, m_covs needs to be list of length 1")

  if(!use_RI){
    if(((length(m_covs[[1]]) != 2)) | ((length(m_covs[[2]]) != 2)) | ((length(m_covs[[3]]) != 2)))
      stop("all three entries of the list m_covs need to be vectors of 2 numerical entries")
  }

  if(use_RI & (!use_simple)){
    if(((length(m_covs[[1]]) != 2)) | ((length(m_covs[[2]]) != 2)))
      stop("all two entries of the list m_covs need to be vectors of 2 numerical entries")
  }

  if(use_simple & (length(m_covs[[1]]) != 2))
    stop("m_covs need to be lists of a vector of 2 numerical entries")

  if(!use_bam)
    warning("note that when use_bam == FALSE, computation time can be large")

  if(max_grid <= min_grid)
    stop("min_grid needs to be smaller than max_grid")

  if(!use_famm)
    warning("as use_famm == FALSE, no confidence bands for mean and covariate effects are given out")

  if(use_famm & !is.list(bs_int_famm))
    stop("bs_int_famm has to be a list")

  if(use_famm & !is.list(bs_y_famm))
    stop("bs_y_famm has to be a list")

  if(use_famm & (length(bs_int_famm) != 3))
    stop("bs_int_famm has to be a list of length 3")

  if(use_famm & (length(bs_y_famm) != 3))
    stop("bs_int_famm has to be a list of length 3")

  if(use_famm & (length(bs_int_famm$k) != 1))
    stop("bs_int_famm$k has to be of length 1")

  if(use_famm & (length(bs_y_famm$k) != 1))
    stop("bs_y_famm$k has to be of length 1")

  if(use_famm & (bs_int_famm$bs == "ps") & (length(bs_int_famm$m) != 2))
    stop("when bs_int_famm$bs == 'ps' is used, bs_int_famm$m has to be of length 2")

  if(use_famm & (bs_y_famm$bs == "ps") & (length(bs_y_famm$m) != 2))
    stop("when bs_y_famm$bs == 'ps' is used, bs_y_famm$m has to be of length 2")

  if(!use_RI & use_simple)
    stop("when using use_simple == TRUE, please also set use_RI to TRUE")

  if(!is.integer(curve_info$n_long))
    stop("in curve_info, n_long needs to be an integer")

  if(!is.integer(curve_info$subject_long))
    stop("in curve_info, subject_long needs to be an integer")

  if(!use_RI & (!is.integer(curve_info$word_long)))
    stop("in curve_info, word_long needs to be an integer")

  if(use_discrete_cov  & (method != 'fREML'))
    stop("when use_discrete_cov == TRUE, please use method 'fREML'")

  if(use_famm & use_discrete_famm  & (method != 'fREML'))
    stop("when use_discrete_famm == TRUE, please use method 'fREML'")

  if(use_discrete_cov & (para_estim_cov == FALSE))
    stop("when use_discrete_cov == TRUE, please set para_estim_cov to TRUE")

  if(use_famm & use_discrete_famm & (para_estim_famm == FALSE))
    stop("when use_discrete_famm == TRUE, please set para_estim_famm to TRUE")

  if(use_discrete_cov & (para_estim_cov_nc == 0 | para_estim_cov_nc == 1))
    stop("when use_discrete_cov == TRUE, please set para_estim_cov_nc larger than 1 ")

  if(use_famm &  use_discrete_famm & (para_estim_famm_nc == 0 | para_estim_famm_nc == 1))
    stop("when use_discrete_famm == TRUE, please set para_estim_famm_nc larger than 1 ")

  if(para_estim_mean & !use_bam)
    warning("parallelization is only performed when use_bam == TRUE")

  if(para_estim_cov & !use_bam)
    warning("parallelization is only performed when use_bam == TRUE")

  if(use_famm & para_estim_famm & !use_bam)
    warning("parallelization is only performed when use_bam == TRUE")

  if(use_discrete_cov)
    warning("Note: covariance estimation is approximated using the 'discrete' option in bam{mcgv} with nthreads = para_estim_cov_nc")

  if(use_famm & use_discrete_famm)
    warning("Note: fpc-famm estimation is approximated using the 'discrete' option in bam{mcgv} with nthreads = para_estim_famm_nc")

  if(use_discrete_cov & !use_bam)
    stop("option 'discrete' can only used with bam, please set use_bam to TRUE")

  if(use_famm & use_discrete_famm & !use_bam)
    stop("option 'discrete' can only used with bam, please set use_bam_famm to TRUE")

  if(use_discrete_cov & !("discrete" %in% names(formals("bam"))))
    stop("option 'discrete' is not implemented in the loaded version of package mgcv, please install a new version")

  if(use_famm & use_discrete_famm & !("discrete" %in% names(formals("bam"))))
    stop("option 'discrete' is not implemented in the loaded version of package mgcv, please install a new version")

  if(use_famm & use_discrete_famm & covariate & any(covariate_form == "smooth"))
    warning("Caution: using the discrete option 'use_discrete_famm' together with estimating smooth effects of covariates may
            (currently) result in unexpected results which has not been clarified yet. It is thus recommended
            to avoid this combination. For more details, see https://github.com/refunders/refund/issues/70")

  ###################
  # initialize output
  ###################
  res <- list()

  ##################################
  # get levels of grouping variables
  ##################################
  I <- length(unique(curve_info$subject_long)) # number of levels of first grouping variable (e.g. speaker)
  if(!use_RI){
    J <- length(unique(curve_info$word_long))  # number of levels of second grouping variable (e.g. target words)
  }else{
    J <- NA # in case of one fRI, J is automatically set to NA
  }

  # renumberlevels of the grouping
  # variables
  ######################################

  # renumber curves
  n <- length(unique(curve_info$n_long))
  n_orig <- unique(curve_info$n_long)
  curve_info[, n_long_orig := n_long]
  for(i in seq_along(n_orig)){
    curve_info[n_long_orig == n_orig[i], n_long := (1:n)[i]]
  }
  res[["n_orig"]] <- n_orig

  # renumber levels of first grouping variable (e.g. speakers)
  subject_orig <- unique(curve_info$subject_long)
  curve_info[, subject_long_orig := subject_long]
  for(i in seq_along(subject_orig)){
    curve_info[subject_long_orig == subject_orig[i], subject_long := (1:I)[i]]
  }
  res[["subject_orig"]] <- subject_orig

  # renumber levels of second grouping variable (e.g. target words)
  if(!use_RI){
    word_orig <- unique(curve_info$word_long)
    curve_info[, word_long_orig := word_long]
    for(i in seq_along(word_orig)){
      curve_info[word_long_orig == word_orig[i], word_long := (1:J)[i]]
    }
    res[["word_orig"]] <- word_orig
  }else{
    res[["word_orig"]] <- NULL
  }

  ##########################################
  # extract measurement points of curve_info
  ##########################################
  t <- curve_info$t

  ##############
  # specify grid
  ##############
  if(is.null(my_grid)){
    # compute grid on which the covariance will later be evaluated
    res[["my_grid"]] <- seq(from = min_grid, to = max_grid, length = d_grid)
  }else{
    res[["my_grid"]] <- my_grid
    d_grid <- length(my_grid)
    if(!isTRUE(all.equal(0, var(diff(res[["my_grid"]]))))){
      warning("the evaluation grid is not equidistant")
    }
  }

  ##############################################################################
  ##############################################################################
  # estimations
  ##############################################################################
  ##############################################################################

  ##############################################################################
  # estimate mean function
  ##############################################################################
  cat("mean estimation", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
  res[["mean_hat"]] <- estimate_mean_fun(bf = bf_mean, bf_covariates = bf_covariates, method = method,
                                         covariate = covariate,num_covariates = num_covariates,
                                         covariate_form = covariate_form,
                                         save_model_mean = save_model_mean, n = n, my_grid = res[["my_grid"]], bs = bs,
                                         m = m_mean, curve_info = curve_info, interaction = interaction,
                                         which_interaction = which_interaction, use_bam = use_bam,
                                         para_estim = para_estim_mean, para_estim_nc = para_estim_mean_nc)

  ########################################
  # preparations for covariance estimation
  ########################################
  if(any(!is.na(res[["mean_hat"]][["y_tilde"]]))){
    y_tilde <- res[["mean_hat"]][["y_tilde"]]

    #########################
    # take out y_tilde
    # from res[["mean_hat]]
    # and add in curve_info
    #########################
    res[["mean_hat"]][["y_tilde"]] <- NULL
    curve_info[, y_tilde := y_tilde]

    ####################
    # get cross products
    ####################
    cat("get cross products", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
    curve_info[, id := 1:nrow(curve_info)]
    preps <- get_crossprods_fun(y_tilde = y_tilde, curve_info = curve_info,
                                my_grid = res[["my_grid"]], d_grid = d_grid, use_RI = use_RI, I = I, J = J,
                                t = t)
    set(curve_info, i = NULL, "id", NULL)
    index <- preps$index
    set(index, i = NULL, "id2", NULL)

    grid_row <- preps$grid_row
    grid_col <- preps$grid_col

    same_subject_grid <- preps$same_subject_grid
    if(!use_RI)
      same_word_grid <- preps$same_word_grid
    same_curve_grid <- preps$same_curve_grid
    same_point_grid <- preps$same_point_grid

    preps <- NULL

    ############################################################################
    # estimation of the covariances
    ############################################################################
    if(use_RI)
      same_word_grid <- NA

    #################################
    # whole: using all cross products
    # without symmetry constraint
    #################################
    if(use_whole){
      cat("covariance estimation whole", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      res[["time_cov_whole"]] <-
        system.time(res[["cov_hat_whole"]] <-
                      estimate_cov_whole_fun(index = index, bf = bf_covs, method = method,
                                             grid_col = grid_col, grid_row = grid_row, d_grid = d_grid,
                                             bs = bs, m = m_covs, use_bam = use_bam, t = t,
                                             same_subject_grid = same_subject_grid,
                                             same_word_grid = same_word_grid, same_curve_grid = same_curve_grid,
                                             same_point_grid = same_point_grid,
                                             mp = mp, para_estim = para_estim_cov, para_estim_nc = para_estim_cov_nc,
                                             use_RI = use_RI, weights = NULL, use_simple = use_simple, np = np,
                                             use_discrete = use_discrete_cov))
    }


    #######################
    # select cross products
    # on upper triangle
    #######################
    if(use_tri | use_tri_constr | use_tri_constr_weights)
      index_upperTri <- index[(n1 <= n2 & row_t_bivariate <= col_t_bivariate) |
                                (n1 > n2 & row_t_bivariate < col_t_bivariate), ]

    #############################
    # tri: using triangle only
    # without symmetry constraint
    #############################
    if(use_tri){
      cat("covariance estimation tri", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      res[["time_cov_tri"]] <-
        system.time(res[["cov_hat_tri"]] <-
                      estimate_cov_tri_fun(index_upperTri = index_upperTri, bf = bf_covs, method = method,
                                           grid_col = grid_col, grid_row = grid_row, d_grid = d_grid,
                                           bs = bs, m = m_covs, use_bam = use_bam, t = t,
                                           same_subject_grid = same_subject_grid,
                                           same_word_grid = same_word_grid, same_curve_grid = same_curve_grid,
                                           same_point_grid = same_point_grid,
                                           mp = mp, para_estim = para_estim_cov, para_estim_nc = para_estim_cov_nc,
                                           use_RI = use_RI, weights = NULL, use_simple = use_simple, np = np,
                                           use_discrete = use_discrete_cov))
    }

    #################################
    # tri constr: using triangle only
    # with our symmetry constraint
    #################################
    if(use_tri_constr){
      cat("covariance estimation tri constraint", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      res[["time_cov_tri_constr"]] <-
        system.time(res[["cov_hat_tri_constr"]] <-
                      estimate_cov_tri_constr_fun(index_upperTri = index_upperTri, bf = bf_covs, method = method,
                                                  grid_col = grid_col, grid_row = grid_row, d_grid = d_grid,
                                                  bs = bs, m = m_covs, use_bam = use_bam, t = t, same_subject_grid = same_subject_grid,
                                                  same_word_grid = same_word_grid, same_curve_grid = same_curve_grid,
                                                  same_point_grid = same_point_grid,
                                                  mp = mp, para_estim = para_estim_cov, para_estim_nc = para_estim_cov_nc,
                                                  use_RI = use_RI, weights = NULL, use_simple = use_simple, np = np,
                                                  use_discrete = use_discrete_cov))
    }

    ###################################################
    # tri constr weights: using triangle only
    # with our symmetry constraint and diagonal weights
    ###################################################
    if(use_tri_constr_weights){
      cat("covariance estimation tri constraint with weights", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      index_upperTri[, weights := 1]
      index_upperTri[row_t_bivariate == col_t_bivariate, weights := 0.5]
      res[["time_cov_tri_constr_weights"]] <-
        system.time(res[["cov_hat_tri_constr_weights"]] <-
                      estimate_cov_tri_constr_fun(index_upperTri = index_upperTri, bf = bf_covs, method = method,
                                                  grid_col = grid_col, grid_row = grid_row, d_grid = d_grid,
                                                  bs = bs, m = m_covs, use_bam = use_bam, t = t, same_subject_grid = same_subject_grid,
                                                  same_word_grid = same_word_grid, same_curve_grid = same_curve_grid,
                                                  same_point_grid = same_point_grid,
                                                  mp = mp, para_estim = para_estim_cov, para_estim_nc = para_estim_cov_nc,
                                                  use_RI = use_RI, weights = index_upperTri$weights,
                                                  use_simple = use_simple, np = np,
                                                  use_discrete = use_discrete_cov))
    }

    ##########################
    # remove unnecessary stuff
    ##########################
    same_subject_grid <- NULL
    same_word_grid <- NULL
    same_curve_grid <- NULL
    same_point_grid <- NULL
    grid_row <- NULL
    grid_col <- NULL

    set(index, i = NULL, "row_t_bivariate", NULL)
    set(index, i = NULL, "col_t_bivariate", NULL)


    ##################################
    # check if grid is equidistant
    # if not: specify equidistant grid
    ##################################
    if(!isTRUE(all.equal(0, var(diff(res[["my_grid"]]))))){
      res[["my_grid"]] <- seq(min_grid, max_grid, length = d_grid)
      warning("the pre-specified evaluation grid was adapted to be equidistant
            for eigen decomposition an following steps")
    }


    ############################################################################
    # Functional principal component analysis
    ############################################################################

    if(use_RI){
      N_C <- 0  #set to zero for ease of implementation
    }

    #################################
    # whole: using all cross products
    # without symmetry constraint
    #################################
    if(use_whole){
      cat("fpc estimation whole", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      if(!is.na(res[["cov_hat_whole"]][["sigmasq"]])){
        res[["time_fpc_whole"]] <-
          system.time(res[["fpc_hat_whole"]] <-
                        estimate_fpc_fun(cov_B = res[["cov_hat_whole"]][["grid_mat_B"]],
                                         cov_C = res[["cov_hat_whole"]][["grid_mat_C"]],
                                         cov_E = res[["cov_hat_whole"]][["grid_mat_E"]],
                                         sigmasq_int_hat = res[["cov_hat_whole"]][["sigmasq_int"]],
                                         my_grid = res[["my_grid"]],
                                         var_level = var_level, N_B = N_B, N_C = N_C, N_E = N_E,
                                         curve_info = curve_info,
                                         I = I, J = J, n = n,
                                         sigmasq_hat = res[["cov_hat_whole"]][["sigmasq"]],
                                         use_RI = use_RI))
      }else{
        res[["fpc_hat_whole"]] <- NA
      }
    }

    #############################
    # tri: using triangle only
    # without symmetry constraint
    #############################
    if(use_tri){
      cat("fpc estimation tri", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      if(!is.na(res[["cov_hat_tri"]][["sigmasq"]])){
        res[["time_fpc_tri"]] <-
          system.time(res[["fpc_hat_tri"]] <-
                        estimate_fpc_fun(cov_B = res[["cov_hat_tri"]][["grid_mat_B"]],
                                         cov_C = res[["cov_hat_tri"]][["grid_mat_C"]],
                                         cov_E = res[["cov_hat_tri"]][["grid_mat_E"]],
                                         sigmasq_int_hat = res[["cov_hat_tri"]][["sigmasq_int"]],
                                         my_grid = res[["my_grid"]],
                                         var_level = var_level, N_B = N_B, N_C = N_C, N_E = N_E,
                                         curve_info = curve_info, I = I, J = J, n = n,
                                         sigmasq_hat = res[["cov_hat_tri"]][["sigmasq"]],
                                         use_RI = use_RI))

      }else{
        res[["fpc_hat_tri"]] <- NA
      }
    }

    #################################
    # tri constr: using triangle only
    # with our symmetry constraint
    #################################
    if(use_tri_constr){
      cat("fpc estimation tri constraint", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      if(!is.na(res[["cov_hat_tri_constr"]][["sigmasq"]])){
        res[["time_fpc_tri_constr"]] <-
          system.time(res[["fpc_hat_tri_constr"]] <-
                        estimate_fpc_fun(cov_B = res[["cov_hat_tri_constr"]][["grid_mat_B"]],
                                         cov_C = res[["cov_hat_tri_constr"]][["grid_mat_C"]],
                                         cov_E = res[["cov_hat_tri_constr"]][["grid_mat_E"]],
                                         sigmasq_int_hat = res[["cov_hat_tri_constr"]][["sigmasq_int"]],
                                         my_grid = res[["my_grid"]],
                                         var_level = var_level, N_B = N_B, N_C = N_C, N_E = N_E,
                                         curve_info = curve_info, I = I, J = J, n = n,
                                         sigmasq_hat = res[["cov_hat_tri_constr"]][["sigmasq"]],
                                         use_RI = use_RI))
      }else{
        res[["fpc_hat_tri_constr"]] <- NA
      }
    }

    ###################################################
    # tri constr: using triangle only
    # with our symmetry constraint and diagonal weights
    ###################################################
    if(use_tri_constr_weights){
      cat("fpc estimation tri constraint weights", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
      if(!is.na(res[["cov_hat_tri_constr_weights"]][["sigmasq"]])){
        res[["time_fpc_tri_constr_weights"]] <-
          system.time(res[["fpc_hat_tri_constr_weights"]] <-
                        estimate_fpc_fun(cov_B = res[["cov_hat_tri_constr_weights"]][["grid_mat_B"]],
                                         cov_C = res[["cov_hat_tri_constr_weights"]][["grid_mat_C"]],
                                         cov_E = res[["cov_hat_tri_constr_weights"]][["grid_mat_E"]],
                                         sigmasq_int_hat = res[["cov_hat_tri_constr_weights"]][["sigmasq_int"]],
                                         my_grid = res[["my_grid"]],
                                         var_level = var_level, N_B = N_B, N_C = N_C, N_E = N_E,
                                         curve_info = curve_info, I = I, J = J, n = n,
                                         sigmasq_hat = res[["cov_hat_tri_constr_weights"]][["sigmasq"]],
                                         use_RI = use_RI))
      }else{
        res[["fpc_hat_tri_constr_weights"]] <- NA
      }
    }


    if(use_famm){
      ########################################################################
      # FPC-FAMM re-estimation of the mean and prediction of the basis weights
      ########################################################################
      if(use_whole){
        #################################
        # whole: using all cross products
        # without symmetry constraint
        #################################
        cat("fpc famm estimation whole", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
        if(!is.na(res[["cov_hat_whole"]][["sigmasq"]])){
          res[["time_fpc_famm_whole"]]<-
            system.time(res[["fpc_famm_hat_whole"]] <-
                          estimate_fpc_famm_fun(curve_info = curve_info, my_grid = res[["my_grid"]],
                                                phi_B_hat_grid = res[["fpc_hat_whole"]]$phi_B_hat_grid,
                                                phi_C_hat_grid = res[["fpc_hat_whole"]]$phi_C_hat_grid,
                                                phi_E_hat_grid = res[["fpc_hat_whole"]]$phi_E_hat_grid,
                                                nu_B_hat = res[["fpc_hat_whole"]]$nu_B_hat,
                                                nu_C_hat = res[["fpc_hat_whole"]]$nu_C_hat,
                                                nu_E_hat = res[["fpc_hat_whole"]]$nu_E_hat,
                                                t = t, N_B = res[["fpc_hat_whole"]]$N_B,
                                                N_C = res[["fpc_hat_whole"]]$N_C,
                                                N_E = res[["fpc_hat_whole"]]$N_E,
                                                use_bam_famm = use_bam_famm,
                                                num_covariates = num_covariates,
                                                which_interaction = which_interaction,
                                                n = n, interaction = interaction,
                                                use_RI = use_RI, method = method, bs_y_famm = bs_y_famm,
                                                bs_int_famm = bs_int_famm, sigmasq = res[["cov_hat_whole"]]$sigmasq,
                                                covariate = covariate, para_estim_famm = para_estim_famm,
                                                para_estim_famm_nc = para_estim_famm_nc, covariate_form = covariate_form,
                                                save_model_famm = save_model_famm,
                                                use_discrete = use_discrete_famm))
        }else{
          res[["fpc_famm_hat_whole"]]<-NA
        }
      }

      if(use_tri){
        #############################
        # tri: using triangle only
        # without symmetry constraint
        #############################
        cat("fpc famm estimation tri", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
        if(!is.na(res[["cov_hat_tri"]][["sigmasq"]])){
          res[["time_fpc_famm_tri"]]<-
            system.time(res[["fpc_famm_hat_tri"]] <-
                          estimate_fpc_famm_fun(curve_info = curve_info, my_grid = res[["my_grid"]],
                                                phi_B_hat_grid = res[["fpc_hat_tri"]]$phi_B_hat_grid,
                                                phi_C_hat_grid = res[["fpc_hat_tri"]]$phi_C_hat_grid,
                                                phi_E_hat_grid = res[["fpc_hat_tri"]]$phi_E_hat_grid,
                                                nu_B_hat = res[["fpc_hat_tri"]]$nu_B_hat,
                                                nu_C_hat = res[["fpc_hat_tri"]]$nu_C_hat,
                                                nu_E_hat = res[["fpc_hat_tri"]]$nu_E_hat,
                                                t = t, N_B = res[["fpc_hat_tri"]]$N_B,
                                                N_C = res[["fpc_hat_tri"]]$N_C,
                                                N_E = res[["fpc_hat_tri"]]$N_E,
                                                use_bam_famm = use_bam_famm,
                                                num_covariates = num_covariates,
                                                which_interaction = which_interaction,
                                                n = n, interaction = interaction,
                                                use_RI = use_RI, method = method, bs_y_famm = bs_y_famm,
                                                bs_int_famm = bs_int_famm, sigmasq = res[["cov_hat_tri"]]$sigmasq,
                                                covariate = covariate, para_estim_famm = para_estim_famm,
                                                para_estim_famm_nc = para_estim_famm_nc, covariate_form = covariate_form,
                                                save_model_famm = save_model_famm,
                                                use_discrete = use_discrete_famm))
        }else{
          res[["fpc_famm_hat_tri"]]<-NA
        }
      }

      if(use_tri_constr){
        ########################################
        # tri constr: using tri_constrangle only
        # with our symmetry constraint
        #########################################
        cat("fpc famm estimation tri constr", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
        if(!is.na(res[["cov_hat_tri_constr"]][["sigmasq"]])){
          res[["time_fpc_famm_tri_constr"]]<-
            system.time(res[["fpc_famm_hat_tri_constr"]] <-
                          estimate_fpc_famm_fun(curve_info = curve_info, my_grid = res[["my_grid"]],
                                                phi_B_hat_grid = res[["fpc_hat_tri_constr"]]$phi_B_hat_grid,
                                                phi_C_hat_grid = res[["fpc_hat_tri_constr"]]$phi_C_hat_grid,
                                                phi_E_hat_grid = res[["fpc_hat_tri_constr"]]$phi_E_hat_grid,
                                                nu_B_hat = res[["fpc_hat_tri_constr"]]$nu_B_hat,
                                                nu_C_hat = res[["fpc_hat_tri_constr"]]$nu_C_hat,
                                                nu_E_hat = res[["fpc_hat_tri_constr"]]$nu_E_hat,
                                                t = t, N_B = res[["fpc_hat_tri_constr"]]$N_B,
                                                N_C = res[["fpc_hat_tri_constr"]]$N_C,
                                                N_E = res[["fpc_hat_tri_constr"]]$N_E,
                                                use_bam_famm = use_bam_famm,
                                                num_covariates = num_covariates,
                                                which_interaction = which_interaction,
                                                n = n, interaction = interaction,
                                                use_RI = use_RI, method = method, bs_y_famm = bs_y_famm,
                                                bs_int_famm = bs_int_famm, sigmasq = res[["cov_hat_tri_constr"]]$sigmasq,
                                                covariate = covariate, para_estim_famm = para_estim_famm,
                                                para_estim_famm_nc = para_estim_famm_nc, covariate_form = covariate_form,
                                                save_model_famm = save_model_famm,
                                                use_discrete = use_discrete_famm))
        }else{
          res[["fpc_famm_hat_tri_constr"]]<-NA
        }
      }

      if(use_tri_constr_weights){
        #########################################
        # tri constr weights: using triangle only
        # with our symmetry constraint
        #########################################
        cat("fpc famm estimation tri constr with weights", "; time: ", format(Sys.time(), "%a %b %d %X"), "\n", sep = "")
        if(!is.na(res[["cov_hat_tri_constr_weights"]][["sigmasq"]])){
          res[["time_fpc_famm_tri_constr_weights"]]<-
            system.time(res[["fpc_famm_hat_tri_constr_weights"]] <-
                          estimate_fpc_famm_fun(curve_info = curve_info, my_grid = res[["my_grid"]],
                                                phi_B_hat_grid = res[["fpc_hat_tri_constr_weights"]]$phi_B_hat_grid,
                                                phi_C_hat_grid = res[["fpc_hat_tri_constr_weights"]]$phi_C_hat_grid,
                                                phi_E_hat_grid = res[["fpc_hat_tri_constr_weights"]]$phi_E_hat_grid,
                                                nu_B_hat = res[["fpc_hat_tri_constr_weights"]]$nu_B_hat,
                                                nu_C_hat = res[["fpc_hat_tri_constr_weights"]]$nu_C_hat,
                                                nu_E_hat = res[["fpc_hat_tri_constr_weights"]]$nu_E_hat,
                                                t = t, N_B = res[["fpc_hat_tri_constr_weights"]]$N_B,
                                                N_C = res[["fpc_hat_tri_constr_weights"]]$N_C,
                                                N_E = res[["fpc_hat_tri_constr_weights"]]$N_E,
                                                use_bam_famm = use_bam_famm,
                                                num_covariates = num_covariates,
                                                which_interaction = which_interaction,
                                                n = n, interaction = interaction,
                                                use_RI = use_RI, method = method, bs_y_famm = bs_y_famm,
                                                bs_int_famm = bs_int_famm, sigmasq = res[["cov_hat_tri_constr_weights"]]$sigmasq,
                                                covariate = covariate, para_estim_famm = para_estim_famm,
                                                para_estim_famm_nc = para_estim_famm_nc, covariate_form = covariate_form,
                                                save_model_famm = save_model_famm,
                                                use_discrete = use_discrete_famm))
        }else{
          res[["fpc_famm_hat_tri_constr_weights"]]<-NA
        }
      }
    }
  }else{
    warning("mean estimation did not succeed")
  }

  ###############
  # return output
  ###############
  return(res)
}


###############################################################################
# documentation of data sets

#' Phonetics acoustic data (complete)
#'
#' The data are part of a large study on consonant assimilation, which is
#' the phenomenon that the articulation of two consonants becomes
#' phonetically more alike when they appear subsequently in fluent speech.
#' The data set contains the audio signals of nine different speakers which
#' repeated the same sixteen German target words each five times.
#' The target words are bisyllabic noun-noun compound words which
#' contained the two abutting consonants of interest,
#' s and sh, in either order. Consonant assimilation is accompanied by a complex interplay
#' of language-specific, perceptual and articulatory factors. The aim in the study was to investigate the assimilation
#' of the two consonants as a function of their order (either first s, then sh or vice-versa),
#' syllable stress (stressed or unstressed) and vowel
#' context, i.e. which vowels are immediately adjacent to the target consonants of interest.
#' The vowels are either of the form ia or ai. For more details, see references below.
#'
#'
#' The variables are as follows:
#'
#' \itemize{
#' \item \code{subject_long}: unique identification number for each speaker.
#' \item \code{word_long}: unique identification number for each target word.
#' \item \code{combi_long}: number of the repetition of the combination of the
#' corresponding speaker and target word.
#' \item \code{y_vec}: the response values for each observation point
#' \item \code{n_long}: unique identification number for each curve.
#' \item \code{t}: the observations point locations.
#' \item \code{covariate.1}: (order of the consonants, reference category first /s/ then /sh/).
#' \item \code{covariate.2}: (stress of the final syllable of the first compound,
#' reference category 'stressed').
#' \item \code{covariate.3}: (stress of the initial syllable of the second compound,
#' reference category 'stressed').
#' \item \code{covariate.4}: (vowel context, reference category ia).
#' \item \code{word_names_long}: names of the target words.
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 24830 rows and 11 variables
#' @name acoustic
#' @references Pouplier, Marianne and Hoole, Philip (2016): Articulatory and
#' Acoustic Characteristics of German Fricative Clusters,
#' Phonetica, 73(1), 52--78.
#'
#' Cederbaum, Pouplier, Hoole, Greven (2016): Functional Linear Mixed Models
#' for Irregularly or Sparsely Sampled Data. Statistical Modelling, 16(1), 67-88.
NULL

###############################################################################################################
# documentation of data sets

#' Phonetics acoustic data (subset)
#'
#' A small subset of the phonetics acoustic
#' data set \code{\link[sparseFLMM]{acoustic}} with observations from two speakers and two items only.
#' This will not produce meaningful results but can be used as a toy
#' data set when testing the code.
#' The variables are as in the full data set, see \code{\link[sparseFLMM]{acoustic}}.
#'
#' @docType data
#' @keywords datasets
#' @format A data frame with 656 rows and 11 variables
#' @name acoustic_subset
#' @references Pouplier, Marianne and Hoole, Philip (2016): Articulatory and
#' Acoustic Characteristics of German Fricative Clusters,
#' Phonetica, 73(1), 52--78.
#'
#' Cederbaum, Pouplier, Hoole, Greven (2016): Functional Linear Mixed Models
#' for Irregularly or Sparsely Sampled Data. Statistical Modelling, 16(1), 67-88.
NULL




