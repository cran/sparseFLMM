##################################################################################################################
# author: Alexander Volkmann
##################################################################################################################
##################################################################################################################
# description: function to estimate smooth auto-covariance operators using the triangle with symmetry
# constraint for the nested model. Gives out auto-covariances evaluated on a pre-specified grid.
##################################################################################################################
estimate_cov_tri_constr_nest_fun <- function(index_upperTri,
                                             bf,
                                             method = "REML",
                                             d_grid = 100,
                                             grid_row,
                                             grid_col,
                                             same_subject_grid,
                                             same_word_grid,
                                             same_curve_grid,
                                             same_point_grid,
                                             bs,
                                             m,
                                             use_bam,
                                             t,
                                             mp = TRUE,
                                             para_estim = FALSE,
                                             para_estim_nc,
                                             weights = NULL,
                                             np = TRUE,
                                             use_discrete = TRUE){

  # local binding of nested because of data.table
  nested <- NULL

  #---------------------------
  # create new variable that is a multiplication of the indices, giving a
  # nested model structure
  index_upperTri[, nested := same_subject*same_word]


  if(use_bam == TRUE){
    #######################
    # covariance estimation
    # using bam
    #######################

    ################
    # set cluster
    # for estimation
    # if specified
    ################
    if(para_estim & !use_discrete){
      if(detectCores() > 1){
        nc_use <- min(detectCores(), para_estim_nc)
        if(.Platform$OS.type=="unix"){
          cl_estim <- makeForkCluster(nnodes = nc_use) # only runs on linux
        }else{
          cl_estim <- makeCluster(nc_use) # also runs on windows
        }
      }else{
        cl_estim <- NULL
      }
    }else{
      cl_estim <- NULL
    }


    ##################
    # for nested fRIs
    if(use_discrete){
      time_cov_estim <-
        system.time(gam1 <- try(bam(cross_vec_bivariate ~ - 1 +
                                      s(row_t_bivariate, col_t_bivariate, by = same_subject, k = bf[1],
                                        bs = "symm", m = m[[1]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                      s(row_t_bivariate, col_t_bivariate, by = nested, k = bf[2],
                                        bs = "symm", m = m[[2]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                      s(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[3],
                                        bs = "symm", m = m[[3]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                      same_point, method = method, data = index_upperTri,
                                    discrete = TRUE, nthreads = para_estim_nc, weights = weights)))
    }else{
      time_cov_estim <-
        system.time(gam1 <- try(bam(cross_vec_bivariate ~ - 1 +
                                      s(row_t_bivariate, col_t_bivariate, by = same_subject, k = bf[1],
                                        bs = "symm", m = m[[1]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                      s(row_t_bivariate, col_t_bivariate, by = nested, k = bf[2],
                                        bs = "symm", m = m[[2]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                      s(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[3],
                                        bs = "symm", m = m[[3]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                      same_point, method = method, data = index_upperTri,
                                    cluster = cl_estim, weights = weights)))
    }


    ##########################
    # stop cluster if existing
    ##########################
    if (!is.null(cl_estim)) stopCluster(cl_estim)
  }else{
    #######################
    # covariance estimation
    # using gam
    #######################

    ##################
    # for nested fRIs
    time_cov_estim <-
      system.time(gam1 <- try(gam(cross_vec_bivariate ~ - 1 +
                                    s(row_t_bivariate, col_t_bivariate, by = same_subject, k = bf[1],
                                      bs = "symm", m = m[[1]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                    s(row_t_bivariate, col_t_bivariate, by = nested, k = bf[2],
                                      bs = "symm", m = m[[2]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                    s(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[3],
                                      bs = "symm", m = m[[3]], xt = list(bsmargin = "ps", kroneckersum = mp, np = np)) +
                                    same_point, method = method, data = index_upperTri, weights = weights)))

  }


  if(class(gam1)[1] != "try-error"){

    ##############################
    # extract smoothing parameters
    ##############################
    sp <- gam1$sp

    ###################
    # extract intercept
    ###################
    intercept <- as.numeric(coefficients(gam1)[1])

    #################
    # extract sigmasq
    # and truncate
    #################
    sigmasq <- as.numeric(max(coefficients(gam1)["same_point"], 0))

    #################################
    # compute the integral of sigmasq
    #################################
    sigmasq_int <- (max(unlist(t))-min(unlist(t))) * sigmasq

    ##########################
    # construct data frame for
    # evaluation on grid
    ##########################

    ##################
    # for nested fRIs
    grid_data <- data.table(row_t_bivariate = grid_row, col_t_bivariate = grid_col, same_subject = same_subject_grid,
                            same_word = same_word_grid, same_curve = same_curve_grid, same_point = same_point_grid)
    grid_data[, nested := same_word*same_subject]

    ##########################
    # remove unnecessary stuff
    ##########################
    grid_row <- NULL
    grid_col <- NULL
    same_subject_grid <- NULL
    same_word_grid <- NULL
    same_curve_grid <- NULL
    same_point_grid <- NULL


    ####################
    # evaluation on grid
    ####################
    time_cov_pred_grid <-
      system.time(grid_smooth <- predict(gam1, newdata = grid_data[row_t_bivariate <=  col_t_bivariate, ],
                                         na.omit = TRUE, type = "terms"))

    #####################
    # extract covariances
    #####################
    ##################
    # for nested fRIs
    grid_mat_B <- matrix(NA, nrow = d_grid, ncol = d_grid)
    grid_mat_B[lower.tri(grid_mat_B, diag = TRUE)] <- grid_smooth[, 2]
    grid_mat_B <- t(grid_mat_B)
    grid_mat_B[lower.tri(grid_mat_B, diag = TRUE)] <- grid_smooth[, 2]

    grid_mat_C <- matrix(NA, nrow = d_grid, ncol = d_grid)
    grid_mat_C[lower.tri(grid_mat_C, diag = TRUE)] <- grid_smooth[, 3]
    grid_mat_C <- t(grid_mat_C)
    grid_mat_C[lower.tri(grid_mat_C, diag = TRUE)] <- grid_smooth[, 3]

    grid_mat_E <- matrix(NA, nrow = d_grid, ncol = d_grid)
    grid_mat_E[lower.tri(grid_mat_E, diag = TRUE)] <- grid_smooth[, 4]
    grid_mat_E <- t(grid_mat_E)
    grid_mat_E[lower.tri(grid_mat_E, diag = TRUE)] <- grid_smooth[, 4]

    grid_smooth <- NULL
    grid_data <- NULL

  }else{
    sigmasq <- NA
    sigmasq_int <- NA
    grid_mat_B <- NA
    grid_mat_C <- NA
    grid_mat_E <- NA
    smooth_cov_y <- NA
    sp <- NA
    time_cov_pred_grid <- NA
    time_cov_estim <- NA
    warning("covariance estimation did not succeed")
  }

  gam1 <- NULL

  ###############
  # return output
  ###############
  out <- list(sigmasq = sigmasq, sigmasq_int = sigmasq_int, grid_mat_B = grid_mat_B,
              grid_mat_C = grid_mat_C, grid_mat_E = grid_mat_E,
              sp = sp, time_cov_estim = time_cov_estim,
              time_cov_pred_grid = time_cov_pred_grid)
  return(out)
}

####################################################################################
