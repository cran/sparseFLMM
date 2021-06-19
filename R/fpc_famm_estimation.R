##################################################################################################
# author: Jona Cederbaum
##################################################################################################
# description: FPCA-FAMM which predicts FPC weights and B, C, E and re-estimates the mean function
# including covariates as functional additive mixed model.
##################################################################################################
estimate_fpc_famm_fun <- function(curve_info, my_grid, phi_B_hat_grid, phi_C_hat_grid, phi_E_hat_grid, 
                                  nu_B_hat, nu_C_hat, nu_E_hat, t, N_B, N_C, N_E, use_bam_famm, 
                                  num_covariates, interaction, which_interaction, n, use_RI, 
                                  method, bs_y_famm, bs_int_famm, sigmasq, covariate, para_estim_famm, 
                                  para_estim_famm_nc, covariate_form, save_model_famm,
                                  use_discrete){
  
  ###################
  # initialize output
  ###################
  results <- list() 
  
  ######################
  # extract curve values
  ######################
  y_vec <- curve_info$y_vec
  
  ##########################
  # construct data and ydata 
  # for function pffr
  ##########################
  if(!use_RI){
    data <- data.table(id_subject = as.factor(curve_info$subject_long), 
                       id_word = as.factor(curve_info$word_long), 
                       id_n = as.factor(curve_info$n_long))
  }else{
    data <- data.table(id_subject = as.factor(curve_info$subject_long), 
                       id_n = as.factor(curve_info$n_long))
  }
  
  if(covariate){
    for(i in 1:num_covariates){
      data[, paste0("covariate.", i) := curve_info[[paste0("covariate.", i)]]]      
      if(interaction){
        for(k in 1:num_covariates){
          if(which_interaction[i, k] & (i < k)){
            prod_help <- curve_info[[paste0("covariate.", i)]] * curve_info[[paste0("covariate.", k)]]
            data[, paste0("inter_", i, "_", k) := prod_help]              
          }
        }  
      }
    }
  }
  
  
  data <- data[!duplicated(data$id_n), ]  # only need one row per curve not per observation in data
  ydata <- data.frame(.obs = curve_info$n_long, .index = t, .value = y_vec)
  rownames(data) <- 1:n
  
  ################
  # construct pcre
  # terms for famm
  ################
  names <- vector()
  N_vec <- c(N_B, N_C, N_E)
  
  if(!is.null(N_vec)){
    funs_names <- c("_B", "_C", "_E")
    for(i in seq_along(funs_names)){
      if(i == 1){
        id <- "id_subject"
      }
      if(i == 2){
        id <- "id_word"
      }
      if(i == 3){
        id <- "id_n"
      }
      if(N_vec[i] > 0){
        names <- cbind(names, paste("pcre(id = ", id, ", efunctions = phi", funs_names[i], 
                                    "_hat_grid, evalues = nu", funs_names[i], "_hat, yind = my_grid)", 
                                    sep = ""))
      }
    }
    listofbys_pc <- c(as.vector(names))
    
    if(covariate){
      names <- vector()
      
      for(i in 1:num_covariates){
        if(covariate_form[i] == "by"){
          name_help <- paste("covariate.", i, "", sep = "")
        }
        if(covariate_form[i] == "smooth"){
          if(all(data[[paste0("covariate.", i)]] %in% c(0, 1))){
            stop("no smooth effects for dummy covariates allowed,  
               please use covariate_form = 'by' for dummy covariates")
          }
          name_help <- paste0("s(covariate.", i, ")")
        }
        names <- cbind(names, name_help)
      }
      
      if(interaction == FALSE){
        listofbys_cov <- as.vector(names)
      }else{
        inter_names <- vector(mode = "character")
        inter_by <- numeric()
        for(i in 1:num_covariates){
          for(k in 1:num_covariates){
            if(which_interaction[i, k] & (i < k)){
              if(!all(data[[paste0("covariate.", i)]] %in% c(0, 1))|!all(data[[paste0("covariate.", k)]] %in% c(0, 1))){
                stop("interaction effects are only implemented between dummy covariates")
              }
              if(covariate_form[i] == "by" & covariate_form[k] == "by"){
                inter_names <- cbind(inter_names, paste("inter_", i, "_", k, "", sep = ""))
              }else{
                warning("interaction effects are only implemented between dummy covariates acting as varying-coefficients")
              }
            }
          }
        }
        listofbys_cov <- c(as.vector(names), c(inter_names))
      }
      
      listofbys <- c(listofbys_pc, listofbys_cov)
    }else{
      listofbys <- listofbys_pc
    }
    
    pred <- as.formula(paste("y_vec ~ 1 + ", paste(listofbys, collapse = " + "), sep = ""))
    
    ####################################
    # get design matrices to get S.scale
    ####################################
    if(use_bam_famm){
      #######################
      # FPC famm estimation/
      # prediction using bam
      #######################
      famm_setup_get_scale <- pffr(pred, yind = my_grid, data = data, ydata = ydata, 
                                   algorithm = "bam", control = gam.control(trace = TRUE), 
                                   method = method, bs.yindex = bs_y_famm, 
                                   bs.int = bs_int_famm, fit = FALSE)
    }else{
      #######################
      # FPC famm estimation/
      # prediction using gam
      #######################
      famm_setup_get_scale <- pffr(pred, yind = my_grid, data = data, ydata = ydata, 
                                   algorithm = "gam", control = gam.control(trace = TRUE), 
                                   method = method, bs.yindex = bs_y_famm, 
                                   bs.int = bs_int_famm, fit = FALSE)
    } 
    
    #################
    # extract S.scale
    #################
    if(N_B > 0){
      scale_B <- famm_setup_get_scale$smooth[[2]]$S.scale    
      if(N_C > 0){
        scale_C <- famm_setup_get_scale$smooth[[3]]$S.scale    
        if(N_E > 0){
          scale_E <- famm_setup_get_scale$smooth[[4]]$S.scale    
          scale_res <- c(scale_B, scale_C, scale_E)
        }else{
          scale_res <- c(scale_B, scale_C)
        }
      }else{
        if(N_E > 0){
          scale_E <- famm_setup_get_scale$smooth[[3]]$S.scale    
          scale_res <- c(scale_B, scale_E)
        }else{
          scale_res <- c(scale_B)
        }
      }  
    }else{
      if(N_C > 0){
        scale_C <- famm_setup_get_scale$smooth[[2]]$S.scale  
        if(N_E > 0){
          scale_E <- famm_setup_get_scale$smooth[[3]]$S.scale  
          scale_res <- c(scale_C, scale_E)
        }else{
          scale_res <- c(scale_C)
        }
      }else{
        if(N_E > 0){
          scale_E <- famm_setup_get_scale$smooth[[2]]$S.scale  
          scale_res <- c(scale_E)  
        }else{
          scale_res <- NULL
          warning("no FPCs chosen at all")
        }
      }
    }
    
    #################
    # specify sp_fix
    #################
    if(covariate){
      num_smooth <- sum(covariate_form == "smooth") # add -1 for smooth effects
      sp_fix <- c(-1, scale_res * sigmasq, rep(-1, (length(listofbys_cov) + num_smooth)))  
    }else{
      sp_fix <- c(-1, scale_res * sigmasq)  
    }
    
    ##############
    # set cluster
    # if specified
    ##############
    if(para_estim_famm & !use_discrete){
      if(detectCores() > 1){
        nc_use <- min(detectCores(), para_estim_famm_nc)
        if(.Platform$OS.type == "unix"){
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
    
    
    ###############
    # estimate famm
    ###############
    if(use_bam_famm){
      if(use_discrete){
        famm_estim <- pffr(pred, sp = sp_fix, yind = my_grid, data = data, ydata = ydata, 
                           algorithm = "bam", control = gam.control(trace = TRUE), method = method, 
                           bs.yindex = bs_y_famm, bs.int = bs_int_famm, discrete = TRUE, nthreads = para_estim_famm_nc)
      }else{
        famm_estim <- pffr(pred, sp = sp_fix, yind = my_grid, data = data, ydata = ydata, 
                           algorithm = "bam", control = gam.control(trace = TRUE), method = method, 
                           bs.yindex = bs_y_famm, bs.int = bs_int_famm, cluster = cl_estim)  
      }
      
      
    }else{
      famm_estim <- pffr(pred, sp = sp_fix, yind = my_grid, data = data, ydata = ydata, 
                         algorithm = "gam", control = gam.control(trace = TRUE), method = method, 
                         bs.yindex = bs_y_famm, bs.int = bs_int_famm)
    }
    
    
    ##########################
    # stop cluster if existing
    ##########################
    if(!is.null(cl_estim)) stopCluster(cl_estim) 
    
    results[["intercept"]] <- famm_estim$coefficients[1]
    results[["sigmasq"]] <- famm_estim$sig2
    
    ###############
    # get residuals
    ###############
    results[["residuals"]] <- famm_estim$residuals
    
    #####################
    # predict in parts
    # to avoid many
    # doubles
    #####################
    if(N_B > 0){
      pred1 <- coef(famm_estim, n2 = length(my_grid))$smterms[[2]]$coef
      pred1_sorted <- cbind(pred1[order(pred1$id_subject.vec), c("id_subject.vec", "value")],
                            use = rep(1:length(my_grid), times = length(unique(curve_info$subject_long))))
      
      pred1_sorted_reshaped <- as.matrix(subset(reshape(pred1_sorted, idvar = "id_subject.vec", direction = "wide", 
                                                        timevar = "use", times = "value", v.names = "value"), 
                                                select = - c(id_subject.vec)),
                                         nrow = length(unique(curve_info$subject_long)), ncol = length(my_grid))
      
      dimnames(pred1_sorted_reshaped)[[1]] <- unique(curve_info$subject_long_orig)
      dimnames(pred1_sorted_reshaped)[[2]] <- NULL
      results[["famm_predict_B"]] <- pred1_sorted_reshaped
    }else{
      results[["famm_predict_B"]] <- NA
    }
    if(N_C > 0){
      if(N_B > 0){
        pred2 <- coef(famm_estim, n2 = length(my_grid))$smterms[[3]]$coef
      }else{
        pred2 <- coef(famm_estim, n2 = length(my_grid))$smterms[[2]]$coef
      }
      pred2_sorted <- cbind(pred2[order(pred2$id_word.vec), c("id_word.vec", "value")],
                            use = rep(1:length(my_grid), times = length(unique(curve_info$word_long))))
      
      pred2_sorted_reshaped <- as.matrix(subset(reshape(pred2_sorted, idvar = "id_word.vec", direction = "wide", 
                                                        timevar = "use", times = "value", v.names = "value"), 
                                                select = - c(id_word.vec)),
                                         nrow = length(unique(curve_info$word_long)), ncol = length(my_grid))
      
      dimnames(pred2_sorted_reshaped)[[1]] <- unique(curve_info$word_long_orig)
      dimnames(pred2_sorted_reshaped)[[2]] <- NULL
      
      results[["famm_predict_C"]] <- pred2_sorted_reshaped
    }else{
      results[["famm_predict_C"]] <- NA
    }
    if(N_E > 0){
      if(N_B > 0){
        if(N_C > 0){
          pred3 <- coef(famm_estim, n2 = length(my_grid))$smterms[[4]]$coef
        }else{
          pred3 <- coef(famm_estim, n2 = length(my_grid))$smterms[[3]]$coef
        }
      }else{
        pred3 <- coef(famm_estim, n2 = length(my_grid))$smterms[[3]]$coef
      }
      pred3_sorted <- cbind(pred3[order(pred3$id_n.vec), c("id_n.vec", "value")],
                            use = rep(1:length(my_grid), times = length(unique(curve_info$n_long))))
      
      pred3_sorted_reshaped <- as.matrix(subset(reshape(pred3_sorted, idvar = "id_n.vec", direction = "wide", 
                                                        timevar = "use", times = "value", v.names = "value"), 
                                                select = - c(id_n.vec)),
                                         nrow = length(unique(curve_info$n_long)), ncol = length(my_grid))
      
      dimnames(pred3_sorted_reshaped)[[1]] <- unique(curve_info$n_long_orig)
      dimnames(pred3_sorted_reshaped)[[2]] <- NULL
      
      results[["famm_predict_E"]] <- pred3_sorted_reshaped
      
    }else{
      results[["famm_predict_E"]] <- NA
    }
    
    ##################
    # get covariate
    # effects and
    # confidence bands
    ##################
    coef_use <- coef(famm_estim, n1 = length(my_grid), n2 = length(my_grid))  
    if(any(covariate_form == "smooth")){
      newdata_smooth_mean <- data.table(id_subject = rep(1, length = length(my_grid)), 
                                        id_word = rep(1, length = length(my_grid)), 
                                        id_n = rep(1, length = length(my_grid)))
      if(covariate){
        for(i in 1:num_covariates){
          if(covariate_form[i] == "by"){
            mean_use <- mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", i)]])
            newdata_smooth_mean[, paste0("covariate.", i) := rep(mean_use, length(my_grid))]    
          }else{
            mean_use <- mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", i)]])
            newdata_smooth_mean[, paste0("covariate.", i) := rep(mean_use, length(my_grid))]  
          }
          
          if(interaction){
            for(k in 1:num_covariates){
              if(which_interaction[i, k] & (i < k)){
                if(all(data[[paste0("covariate.", i)]] %in% c(0, 1)) & all(data[[paste0("covariate.", k)]] %in% c(0, 1))){
                  mean_use <- mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", i)]]) * mean(curve_info[!duplicated(n_long), ][[paste0("covariate.", k)]])
                  newdata_smooth_mean[, paste0("inter_", i, "_", k) := rep(mean_use, length(my_grid))]
                }else{
                  warning("interaction effects are only implemented between dummy covariates")
                }
              }
            }  
          } 
        }
      }
      
      pred_smooth_mean <- predict(famm_estim, type = "iterms", newdata = newdata_smooth_mean)
      results[["pred_smooth_mean"]] <- pred_smooth_mean
    }
    
    coef_use_sm <- coef_use$smterms
    results[["famm_cb_mean"]] <- coef_use_sm[[1]]$coef
    
    if(covariate){
      for(i in 1:num_covariates){
        if(covariate_form[i] == "by"){
          results[[paste0("famm_cb_covariate.", i)]] <- coef_use_sm[[paste0("covariate.", i, "(yindex)")]][["coef"]]
        }
        if(covariate_form[i] == "smooth"){
          results[[paste0("famm_cb_covariate.", i)]] <- coef_use_sm[[paste0("s(covariate.", i, ")")]][["coef"]]
          results[[paste0("famm_cb_covariate.", i, "_mean")]] <- pred_smooth_mean[[paste0("ti(covariate.", i, ", yindex.vec)")]]
        }
        if(interaction){
          for(k in 1:num_covariates){
            if(which_interaction[i, k] & (i < k)){
              if(all(curve_info[[paste0("covariate.", i)]] %in% c(0, 1)) & all(curve_info[[paste0("covariate.", k)]] %in% c(0, 1))){
                results[[paste0("famm_cb_inter_", i, "_", k)]] <- coef_use_sm[[paste0("inter_", i, "_", k, "(yindex)")]][["coef"]]
              }else{
                warning("interaction effects are only implemented between dummy covariates")
              }
            }
          }
        }
      }
    }
    
    ###################################################
    # get original basis weights of FPC-FAMM estimation
    ###################################################
    #######
    # for B
    if(N_B > 0){
      xid <- diag(length(unique(data$id_subject)))
      xef <- phi_B_hat_grid
      xi_B_hat_famm_help <- qr.coef(qr(xid %x% xef), as.vector(t(results[["famm_predict_B"]])))
      xi_B_hat_famm <- t(matrix(xi_B_hat_famm_help, ncol = length(unique(data$id_subject)), nrow = N_B))
      row.names(xi_B_hat_famm) <- unique(curve_info$subject_long_orig)
    }else{
      xi_B_hat_famm <- NA
    }
    results[["xi_B_hat_famm"]] <- xi_B_hat_famm
    
    #######
    # for C
    if(N_C > 0){
      xid <- diag(length(unique(data$id_word)))
      xef <- phi_C_hat_grid
      xi_C_hat_famm_help <- qr.coef(qr(xid %x% xef), as.vector(t(results[["famm_predict_C"]])))
      xi_C_hat_famm <- t(matrix(xi_C_hat_famm_help, ncol = length(unique(data$id_word)), nrow = N_C))
      row.names(xi_C_hat_famm) <- unique(curve_info$word_long_orig)
    }else{
      xi_C_hat_famm <- NA
    }
    results[["xi_C_hat_famm"]] <- xi_C_hat_famm
    
    #######
    # for E
    if(N_E > 0){
      xid <- diag(length(unique(data$id_n)))
      xef <- phi_E_hat_grid
      xi_E_hat_famm_help <- qr.coef(qr(xid %x% xef), as.vector(t(results[["famm_predict_E"]])))
      xi_E_hat_famm <- t(matrix(xi_E_hat_famm_help, ncol = length(unique(data$id_n)), nrow = N_E))
      row.names(xi_E_hat_famm) <- unique(curve_info$n_long_orig)
    }else{
      xi_E_hat_famm <- NA
    }
    results[["xi_E_hat_famm"]] <- xi_E_hat_famm
    
    
    ##################
    # save famm object
    # if specified
    ##################
    if(save_model_famm){
      results[["famm_estim"]] <- famm_estim
    }
    famm_estim <- NULL
  }else{
    warning(" no FPCs chosen at all")
  }
  
  ###############
  # return output
  ###############
  return(results)
}

###################################################################################################
