#######################################################################################################################
# author: Jona Cederbaum
#######################################################################################################################
# description: estimates smooth mean function and smooth covariate and interaction effects and gives out centered data.
# NOTE: so far all covariates need to enter the mean in the same way.
#######################################################################################################################
estimate_mean_fun <- function(bf, bf_covariates, method, save_model_mean,
                        n, my_grid, bs, m, use_bam, curve_info, num_covariates, 
                        covariate_form, interaction, which_interaction, covariate,
                        para_estim, para_estim_nc){
  
  results <- list() 
  dat_help <- copy(curve_info)
  
  ###############
  # estimate mean
  ###############
  if(covariate){
    names <- vector()
    
    for(i in 1:num_covariates){
      if(covariate_form[i] == "by"){
        name_help <- paste0("s(t, k = bf_covariates, bs = bs, m = m, 
                            by = covariate.", i, ")")
      }
      if(covariate_form[i] == "smooth"){
        if(all(dat_help[[paste0("covariate.", i)]] %in% c(0, 1))){
          stop("no smooth effects for dummy covariates allowed, 
               please use covariate_form = 'by' for dummy covariates")
        }
        name_help <- paste0("ti(t, covariate.", i, ", k = bf_covariates, bs = bs, 
                            m = m, mc = c(0, 1)", ")")
      }
      names <- cbind(names, name_help)
    }
    
    if(interaction == FALSE){
      listofbys <- as.vector(names)
      pred <- as.formula(paste("y_vec ~ ", paste0(listofbys, collapse = "+"), 
                               " + s(t, k = bf, bs = bs, m = m)"))
    }else{
      inter_names <- vector(mode = "character")
      inter_by <- numeric()
      for(i in 1:num_covariates){
        for(k in 1:num_covariates){
          if(which_interaction[i, k] & (i < k)){
            
            if(!all(dat_help[[paste0("covariate.", i)]] %in% c(0, 1))|!all(dat_help[[paste0("covariate.", k)]] %in% c(0, 1))){
              stop("interaction effects are only implemented between dummy covariates")
            }
            
            prod_help <- curve_info[[paste0("covariate.", i)]] * curve_info[[paste0("covariate.", k)]]
            dat_help[, paste0("inter_", i, "_", k) := prod_help]
            if(covariate_form[i] == "by" & covariate_form[k] == "by"){
              inter_names <- cbind(inter_names, paste0("s(t, k = bf_covariates, 
                                                       bs = bs, m = m, by = inter_", i, "_", k, ")"))
            }else{
              warning("interaction effects are only implemented between dummy covariates acting as varying-coefficients")
            }
          }
        }
      }
      listofbys <- c(as.vector(names), c(inter_names))
      pred <- as.formula(paste("y_vec ~ ", paste(listofbys, collapse = "+"), 
                               " + s(t, k = bf, bs = bs, m = m)", sep = ""))
    }
  }else{
    ys <- curve_info$y_vec
    t <- curve_info$t
    pred <- ys ~ s(t, k = bf, bs = bs, m = m)
  }
  
  
  ################
  # set cluster
  # for estimation
  # if specified
  ################
  if(para_estim){
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
  
  ############
  # estimation
  ############
  if(use_bam == TRUE){
    gam1 <- try(bam(pred, data = dat_help, method = method))
  }else{
    gam1 <- try(gam(pred, data = dat_help, method = method))
  } 
  
  ##########################
  # stop cluster if existing
  ##########################
  if (!is.null(cl_estim)) stopCluster(cl_estim) 
  dat_help <- NULL
  
  
  ########################
  # estimation successfull
  ########################
  
  if(class(gam1)[1] != "try-error"){
    
    ###################
    ##extract intercept
    ###################
    intercept <- coefficients(gam1)[1]
    
    #######################
    # Evaluate mean on grid
    #######################
    # make data frame for prediction
    if(covariate){
      newdata <- data.table(t = my_grid)
      for(i in 1:num_covariates){
        if(covariate_form[i] == "by"){
          newdata[, paste0("covariate.", i) := rep(1, length(my_grid))]
        }else{
          range_mean <- range(curve_info[[paste0("covariate.", i)]])
          newdata[, paste0("covariate.", i) := seq(from = range_mean[1], 
                                                   to = range_mean[2], length = length(my_grid))]
        }
        
        if(interaction){
          for(k in 1:num_covariates){
            if(which_interaction[i, k] & (i < k)){
              if(all(curve_info[[paste0("covariate.", i)]] %in% c(0, 1)) & all(curve_info[[paste0("covariate.", k)]] %in% c(0, 1))){
                newdata[, paste0("inter_", i, "_", k) := rep(1, length(my_grid))]
              }else{
                warning("interaction effects are only implemented between dummy covariates")
              }
            }    
          }
        }
      }
      
      
      # predict all components at once with type = iterms
      mean_pred <- predict(gam1, newdata = newdata, type = "iterms")
      
      if(any(covariate_form == "smooth")){
        use_grid <- seq(min(my_grid), max(my_grid), length = length(my_grid))
        newdata_smooth <- data.table(t = expand.grid(use_grid, use_grid)[, 1])
        newdata_smooth_mean <- data.table(t = use_grid)
        for(i in 1:num_covariates){
          if(covariate_form[i] == "by"){
            newdata_smooth[, paste0("covariate.", i) := rep(1, nrow(newdata_smooth))]
            mean_use <- mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", i)]])
            newdata_smooth_mean[, paste0("covariate.", i) := rep(mean_use, nrow(newdata_smooth_mean))]
            if(interaction){
              for(k in 1:num_covariates){
                if(which_interaction[i, k] & (i < k)){
                  if(all(curve_info[[paste0("covariate.", i)]] %in% c(0, 1)) & all(curve_info[[paste0("covariate.", k)]] %in% c(0, 1))){
                    newdata_smooth[, paste0("inter_", i, "_", k) := rep(1, nrow(newdata_smooth))]
                    mean_use <- mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", i)]]) * mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", k)]])
                    newdata_smooth_mean[, paste0("inter_", i, "_", k) := rep(mean_use, nrow(newdata_smooth_mean))]
                  }else{
                    warning("interaction effects are only implemented between dummy covariates")
                  }
                }    
              }
            }
          }
          if(covariate_form[i] == "smooth"){
            use_cov <- seq(min(newdata[[paste0("covariate.", i)]]), max(newdata[[paste0("covariate.", i)]]), length = length(my_grid))
            newdata_smooth[, paste0("covariate.", i) := expand.grid(use_grid, use_cov)[, 2]]
            mean_use <- mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", i)]])
            newdata_smooth_mean[, paste0("covariate.", i) := rep(mean_use, nrow(newdata_smooth_mean))]
          }
        }
        
        mean_pred_smooth <- list()
        mean_pred_smooth$predict <- predict(gam1, newdata = newdata_smooth, type = "iterms")
        mean_pred_smooth$predict_mean <- predict(gam1, newdata = newdata_smooth_mean, type = "iterms")
        mean_pred_smooth$newdata <- newdata_smooth
      }else{
        mean_pred_smooth <- NA
      }
      
    }else{
      newdat <- data.frame(t = my_grid)
      mean_pred <- predict(gam1, newdata = newdat)
      mean_pred_smooth <- NA
    }
    
    
    # construct estimated mean on original data points 
    eta_hat <- fitted(gam1)
    
    ##########
    # center y
    ##########
    y_tilde <- curve_info$y_vec-eta_hat
    
  }else{
    y_tilde <- rep(NA, length = nrow(curve_info))
    intercept <- NA
    mean_pred <- NA
  }
  
  ########
  # Output
  ########
  results[["y_tilde"]] <-  y_tilde
  results[["intercept"]] <- intercept
  results[["mean_pred"]] <- mean_pred
  results[["mean_pred_smooth"]] <- mean_pred_smooth
  
  ###################
  # save model object
  # if specified
  ###################
  if(save_model_mean == TRUE){
    results[["gam_object"]] <- gam1
  }
  
  gam1 <- NULL
  
  return(results)
}

################################################################################
