#############################################################################################
# author: Jona Cederbaum
#############################################################################################
# description: useful functions that are used in get_crossprods and also for plotting results
#############################################################################################

###############################
# long data frame for bivariate
###############################
create_data_frame_bivariate_fun <- function(index){
  index[, same_subject := as.numeric(row_subject_bivariate == col_subject_bivariate)]
  index[, same_word := as.numeric(row_word_bivariate == col_word_bivariate)]
  index[, same_curve := as.numeric(row_word_bivariate == col_word_bivariate &
                        row_subject_bivariate == col_subject_bivariate & row_combi_bivariate == col_combi_bivariate)]
  index[, same_point := as.numeric(same_curve == 1 & row_t_bivariate == col_t_bivariate)]

  set(index, i = NULL, "row_subject_bivariate", NULL)
  set(index, i = NULL, "col_subject_bivariate", NULL)
  set(index, i = NULL, "row_word_bivariate", NULL)
  set(index, i = NULL, "col_word_bivariate", NULL)
  set(index, i = NULL, "row_combi_bivariate", NULL)
  set(index, i = NULL, "col_combi_bivariate", NULL)
   return(index)
}

###############################
# long data frame for bivariate
# for case of RI
###############################
create_data_frame_bivariate_RI_fun <- function(index){
  index[, same_subject := as.numeric(row_subject_bivariate == col_subject_bivariate)]
  index[, same_curve := as.numeric(row_curve_bivariate == col_curve_bivariate)]
  index[, same_point := as.numeric(same_curve == 1 & row_t_bivariate == col_t_bivariate)]
  
  set(index, i = NULL, "row_subject_bivariate", NULL)
  set(index, i = NULL, "col_subject_bivariate", NULL)
  set(index, i = NULL, "row_curve_bivariate", NULL)
  set(index, i = NULL, "col_curve_bivariate", NULL)
  
  return(index)
}


###############################
# long data frame for bivariate
# on grid
###############################
create_grid_data_fun <- function(my_grid, d_grid){
  out <- list()
  out[["grid_row"]] <- rep(my_grid, each = d_grid)
  out[["grid_col"]] <- rep(my_grid, d_grid)
  out[["same_subject_grid"]] <- rep(1, length(out[["grid_row"]]))
  out[["same_word_grid"]] <- rep(1, length(out[["grid_row"]]))
  out[["same_curve_grid"]] <- rep(1, length = length(out[["grid_row"]]))
  out[["same_point_grid"]] <- rep(0, length = length(out[["grid_row"]]))
  return(out)
}

#################################
# get nice colors for persp plots
#################################
get_color_persp_fun<-function(z){
  jet.colors <- function(n){terrain.colors(n,alpha = 1)}
  nrz <- nrow(z)
  ncz <- ncol(z)
  # generate the desired number of colors from this palette
  nbcol <- 100
  color <- jet.colors(n=nbcol)
  # compute the z-value at the facet centres
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  # recode facet z-values into color indices
  facetcol <- cut(zfacet, nbcol)
  return(color[facetcol])
}


#################################
# get covariance surface from
# eigenvalues and eigenfunctions
#################################
get_covs_fun <- function(N, phi, nu){
  cov <- phi %*% diag(nu, N, N) %*% t(phi)
  return(cov)
}
################################################################################
