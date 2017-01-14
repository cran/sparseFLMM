##################################################################################################################
# author: Jona Cederbaum
##################################################################################################################
# description: function to estimate smooth auto-covariance operators using all cross products and without symmetry
# constraint. Gives out auto-covariances evaluated on a pre-specified grid.
##################################################################################################################
estimate_cov_whole_fun <- function(index, use_RI = FALSE, use_simple = FALSE, bf, method = "REML", 
                                   d_grid = 100, grid_row, grid_col, same_subject_grid, same_word_grid, 
                                   same_curve_grid, same_point_grid, bs, m, use_bam, t,
                                   mp = TRUE, para_estim = FALSE, para_estim_nc, 
                                   weights = NULL, np = TRUE, use_discrete = TRUE){
  
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
          cl_estim <- makeForkCluster(nnodes = nc_use) # only runs on unix
        }else{
          cl_estim <- makeCluster(nc_use) # also runs on windows
        }
      }else{
        cl_estim <- NULL
      }
    }else{
      cl_estim <- NULL  
    }
    
    if(!use_RI){
      ##################
      # for crossed fRIs
      if(use_discrete){
        time_cov_estim <- 
          system.time(gam1 <- try(bam(cross_vec_bivariate ~ - 1 + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_subject, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_word, k = bf[2], bs = c(bs, bs), m = m[[2]], mp = mp, np = np) + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[3], bs = c(bs, bs), m = m[[3]], mp = mp, np = np) + 
                                        same_point, method = method, data = index, 
                                      discrete = TRUE, nthreads = para_estim_nc, weights = weights)))
      }else{
        time_cov_estim <- 
          system.time(gam1 <- try(bam(cross_vec_bivariate ~ - 1 + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_subject, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_word, k = bf[2], bs = c(bs, bs), m = m[[2]], mp = mp, np = np) + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[3], bs = c(bs, bs), m = m[[3]], mp = mp, np = np) + 
                                        same_point, method = method, data = index, cluster = cl_estim, weights = weights)))
      }
      
    }else{
      if(!use_simple){
        #############
        # for one fRI
        if(use_discrete){
          time_cov_estim <- 
            system.time(gam1 <- try(bam(cross_vec_bivariate ~ 1 + 
                                          te(row_t_bivariate, col_t_bivariate, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                          te(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[2], bs = c(bs, bs), m = m[[2]], mp = mp, np = np) + 
                                          same_point, method = method, data = index, 
                                        discrete = TRUE, nthreads = para_estim_nc, weights = weights)))  
        }else{
          time_cov_estim <- 
            system.time(gam1 <- try(bam(cross_vec_bivariate ~ 1 + 
                                          te(row_t_bivariate, col_t_bivariate, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                          te(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[2], bs = c(bs, bs), m = m[[2]], mp = mp, np = np) + 
                                          same_point, method = method, data = index, cluster = cl_estim, weights = weights)))
        }
        
      }else{
        ########################
        # for independent curves
        if(use_discrete){
          time_cov_estim <- 
            system.time(gam1 <- try(bam(cross_vec_bivariate ~ 1 + 
                                          te(row_t_bivariate, col_t_bivariate, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                          same_point, method = method, data = index, 
                                        discrete = TRUE, nthreads = para_estim_nc, weights = weights)))    
        }else{
          time_cov_estim <- 
            system.time(gam1 <- try(bam(cross_vec_bivariate ~ 1 + 
                                          te(row_t_bivariate, col_t_bivariate, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                          same_point, method = method, data = index, cluster = cl_estim, weights = weights)))
        }
      
      }
    }
    
    ##########################
    # stop cluster if existing
    ##########################
    if(!is.null(cl_estim)) stopCluster(cl_estim) 
  }else{      
    #######################
    # covariance estimation
    # using gam
    #######################
    
    if(!use_RI){
      ##################
      # for crossed fRIs
      time_cov_estim <- 
        system.time(gam1 <- try(gam(cross_vec_bivariate ~ - 1 + 
                                      te(row_t_bivariate, col_t_bivariate, by = same_subject, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                      te(row_t_bivariate, col_t_bivariate, by = same_word, k = bf[2], bs = c(bs, bs), m = m[[2]], mp = mp, np = np) + 
                                      te(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[3], bs = c(bs, bs), m = m[[3]], mp = mp, np = np) + 
                                      same_point, method = method, data = index, weights = weights)))
    }else{ 
      if(!use_simple){
        #############
        # for one fRI
        time_cov_estim <- 
          system.time(gam1 <- try(gam(cross_vec_bivariate ~ 1 + 
                                        te(row_t_bivariate, col_t_bivariate, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                        te(row_t_bivariate, col_t_bivariate, by = same_curve, k = bf[2], bs = c(bs, bs), m = m[[2]], mp = mp, np = np) + 
                                        same_point, method = method, data = index, weights = weights)))
      }else{
        ########################
        # for independent curves
        time_cov_estim <- 
          system.time(gam1 <- try(gam(cross_vec_bivariate ~ 1 + 
                                        te(row_t_bivariate, col_t_bivariate, k = bf[1], bs = c(bs, bs), m = m[[1]], mp = mp, np = np) + 
                                        same_point, method = method, data = index, weights = weights)))
      }
      
    }
  }
  
  if(class(gam1)[1] != "try-error"){    
    
    ##############################
    # extract smoothing parameters
    ##############################
    sp <- gam1$sp
    
    if(use_RI){
      ###################
      # extract intercept
      ###################
      intercept <- as.numeric(coefficients(gam1)[1])
    }
    
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
    
    if(!use_RI){
      ##################
      # for crossed fRIs
      grid_data <- data.table(row_t_bivariate = grid_row, col_t_bivariate = grid_col, same_subject = same_subject_grid, 
                              same_word = same_word_grid, same_curve = same_curve_grid, same_point = same_point_grid)
      
      ##########################
      # remove unnecessary stuff
      ##########################
      grid_row <- NULL
      grid_col <- NULL
      same_subject_grid <- NULL
      same_word_grid <- NULL
      same_curve_grid <- NULL
      same_point_grid <- NULL
    }else{
      ###################################
      # for one fRI or independent curves
      grid_data <- data.table(row_t_bivariate = grid_row, col_t_bivariate = grid_col, same_subject = same_subject_grid, 
                              same_curve = same_curve_grid, same_point = same_point_grid)
      
      ##########################
      # remove unnecessary stuff
      ##########################
      grid_row <- NULL
      grid_col <- NULL
      same_subject_grid <- NULL
      same_curve_grid <- NULL
      same_point_grid <- NULL
    }
    
    ####################
    # evaluation on grid
    ####################
    time_cov_pred_grid <- 
      system.time(grid_smooth <- predict(gam1, newdata = grid_data, na.omit = TRUE, type = "terms"))
    
    grid_data <- NULL
    
    #####################
    # extract covariances
    #####################
    if(!use_RI){
      ##################
      # for crossed fRIs
      grid_mat_B <- matrix(grid_smooth[, 2], ncol = d_grid, nrow = d_grid, byrow = TRUE)  
      grid_mat_C <- matrix(grid_smooth[, 3], ncol = d_grid, nrow = d_grid, byrow = TRUE)  
      grid_mat_E <- matrix(grid_smooth[, 4], ncol = d_grid, nrow = d_grid, byrow = TRUE)
    }else{
      if(!use_simple){
        #############
        # for one fRI
        grid_mat_B <- matrix(grid_smooth[, 2], ncol = d_grid, nrow = d_grid, byrow = TRUE) + intercept
        grid_mat_C <- matrix(rep(0, length = d_grid^2), ncol = d_grid, nrow = d_grid, byrow = TRUE)  
        grid_mat_E <- matrix(grid_smooth[, 3], ncol = d_grid, nrow = d_grid, byrow = TRUE)  
      }else{
        ########################
        # for independent curves
        grid_mat_B <- matrix(grid_smooth[, 2], ncol = d_grid, nrow = d_grid, byrow = TRUE) + intercept
        grid_mat_C <- matrix(rep(0, length = d_grid^2), ncol = d_grid, nrow = d_grid, byrow = TRUE)  
        grid_mat_E <- matrix(rep(0, length = d_grid^2), ncol = d_grid, nrow = d_grid, byrow = TRUE)  
      }
      grid_smooth <- NULL    
    }
  }else{
    sigmasq <- NA
    sigmasq_int <- NA
    grid_mat_B <- NA
    grid_mat_C <- NA
    grid_mat_E <- NA
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
