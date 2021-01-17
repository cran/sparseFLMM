##################################################################################################
# author : Jona Cederbaum
##################################################################################################
# description : eigen decomposition of the covariance estimates evaluated on a pre-specified grid
# and prediction of the FPC weights.
##################################################################################################
estimate_fpc_fun <- function(cov_B, cov_C, cov_E, sigmasq_int_hat, my_grid, var_level,
                               N_B, N_C, N_E, curve_info, I, J, n, sigmasq_hat, use_RI,
                             nested){

  ##################################################
  # compute interval length of Riemann approximation
  ##################################################
  interv <- my_grid[2] - my_grid[1]

  ############################
  # extract measurement points
  ############################
  t <- curve_info$t

  ##############################
  # eigen decomposition of cov_B
  ##############################
  cov_B <- symmpart(cov_B)
  my_eigen_B <- eigen(cov_B, symmetric = TRUE)
  cov_B <- NULL

  #####################
  # rescale eigenvalues
  #####################
  nu_B_hat <- my_eigen_B$values * interv

  ######################
  # truncate eigenvalues
  ######################
  neg_nu_B <- which(nu_B_hat < 10^(-10) * max(nu_B_hat))
  if(length(neg_nu_B) > 0)
    warning(paste0(length(neg_nu_B), " negative eigenvalues of B are truncated"))
  nu_B_hat[neg_nu_B] <- 0


  ##############################
  # eigen decomposition of cov_C
  ##############################
  if(!use_RI){
    cov_C <- symmpart(cov_C)
    my_eigen_C <- eigen(cov_C, symmetric = TRUE)
    cov_C <- NULL

    #####################
    # rescale eigenvalues
    #####################
    nu_C_hat <- my_eigen_C$values * interv

    ######################
    # truncate eigenvalues
    ######################
    neg_nu_C <- which(nu_C_hat < 10^(-10) * max(nu_C_hat))
    if(length(neg_nu_C) > 0)
      warning(paste0(length(neg_nu_C), " negative eigenvalues of C are truncated"))
    nu_C_hat[neg_nu_C] <- 0

  }else{
    my_eigen_C <- eigen(cov_C)
    nu_C_hat <- 0
  }

  ##############################
  # eigen decomposition of cov_E
  ##############################
  cov_E <- symmpart(cov_E)
  my_eigen_E <- eigen(cov_E, symmetric = TRUE)
  cov_E <- NULL

  #####################
  # rescale eigenvalues
  #####################
  nu_E_hat <- my_eigen_E$values * interv

  ######################
  # truncate eigenvalues
  ######################
  neg_nu_E <- which(nu_E_hat < 10^(-10) * max(nu_E_hat))
  if(length(neg_nu_E) > 0)
    warning(paste0(length(neg_nu_E), " negative eigenvalues of E are truncated"))
  nu_E_hat[neg_nu_E] <- 0


  ###############################################
  # compute total variance and variance explained
  ###############################################
  total_var <- sum(nu_B_hat * (nu_B_hat > 0)) + sum(nu_C_hat * (nu_C_hat > 0)) +
    sum(nu_E_hat * (nu_E_hat > 0)) + sigmasq_int_hat

  ##############################################
  # specify number of components to keep
  # (depends on var_level and N if it is not NA)
  ##############################################
  if(is.na(N_B)|is.na(N_C)|is.na(N_E)){
    prop <- N_B <- N_C <- N_E <- 0
    while(prop<var_level){
      if(!use_RI){
        nu_all <- c(nu_B_hat[N_B + 1], nu_C_hat[N_C + 1], nu_E_hat[N_E + 1])
        N_all <- c(N_B, N_C, N_E)
        maxi <- which.max(nu_all)
        N_all[maxi] <-  N_all[maxi] + 1
        N_B <- N_all[1]
        N_C <- N_all[2]
        N_E <- N_all[3]
        prop <- (sum(nu_B_hat[seq(len = N_B)]) + sum(nu_C_hat[seq(len = N_C)]) +
                   sum(nu_E_hat[seq(len = N_E)]) + sigmasq_int_hat) / total_var
      }else{
        nu_all <- c(nu_B_hat[N_B + 1], nu_E_hat[N_E + 1])
        N_all <- c(N_B, N_E)
        maxi <- which.max(nu_all)
        N_all[maxi] <-  N_all[maxi] + 1
        N_B <- N_all[1]
        N_E <- N_all[2]
        prop <- (sum(nu_B_hat[seq(len = N_B)]) + sum(nu_E_hat[seq(len = N_E)]) +
                   sigmasq_int_hat) / total_var

      }
    }
  }

  if(N_B != 0|N_C != 0|N_E != 0){
    ######################################################
    # truncate eigen values to level of explained variance
    ######################################################
    nu_B_hat <- nu_B_hat[seq(len = N_B), drop = FALSE]
    nu_C_hat <- nu_C_hat[seq(len = N_C), drop = FALSE]
    nu_E_hat <- nu_E_hat[seq(len = N_E), drop = FALSE]
    var_explained <- (sum(nu_B_hat) + sum(nu_C_hat) + sum(nu_E_hat) + sigmasq_int_hat) / total_var

    ###############################
    # truncate and recalculate NPC
    # if one chosen eigenvalue is 0
    ###############################

    #######
    # for B
    if(N_B > 0){
      while(nu_B_hat[N_B]<10^(-8)){
        N_B <- N_B-1
        if(N_B > 0){
          nu_B_hat <- nu_B_hat[seq(len = N_B), drop = FALSE]
        }else{
          nu_B_hat <- 0
        }
        if(N_B == 0) break  # if this was the only eigenvalue than stop
      }
    }

    #######
    # for C
    if(N_C > 0){
      while(nu_C_hat[N_C]<10^(-8)){
        N_C <- N_C-1
        if(N_C > 0){
          nu_C_hat <- nu_C_hat[seq(len = N_C), drop = FALSE]
        }else{
          nu_C_hat <- 0
        }
        if(N_C == 0) break  # if this was the only eigenvalue than stop
      }
    }

    #######
    # for E
    if(N_E > 0){
      while(nu_E_hat[N_E]<10^(-8)){
        N_E <- N_E-1
        if(N_E > 0){
          nu_E_hat <- nu_E_hat[seq(len = N_E), drop = FALSE]
        }else{
          nu_E_hat <- 0
        }
        if(N_E == 0) break  # if this was the only eigenvalue than stop
      }
    }


    if(N_B != 0){
      ########################
      # rescale eigenfunctions
      ########################
      phi_B_hat_grid <- (1 / sqrt(interv)) * my_eigen_B$vectors[, seq(len = N_B), drop = FALSE]

      ############################
      # interpolate eigenfunctions
      # and evaluate on original
      # measurement points
      ############################
      phi_B_hat_orig <- matrix(NA, ncol = N_B, nrow = length(unlist(t)))
      for(k in 1 : N_B){
        phi_B_hat_orig[, k] <- approx(x = my_grid, y = phi_B_hat_grid[, k],
                                      xout = unlist(t), method = "linear")$y
      }
    }else{
      phi_B_hat_grid <- matrix(0, 0, 0)
      phi_B_hat_orig <- matrix(0, 0, 0)
    }
    my_eigen_B <- NULL

    if(N_C != 0){
      ########################
      # rescale eigenfunctions
      ########################
      phi_C_hat_grid <- (1 / sqrt(interv)) * my_eigen_C$vectors[, seq(len = N_C), drop = FALSE]

      ############################
      # interpolate eigenfunctions
      # and evaluate on original
      # measurement points
      ############################
      phi_C_hat_orig <- matrix(NA, ncol = N_C, nrow = length(unlist(t)))
      for(k in 1 : N_C){
        phi_C_hat_orig[, k] <- approx(x = my_grid, y = phi_C_hat_grid[, k],
                                      xout = unlist(t), method = "linear")$y
      }
    }else{
      phi_C_hat_grid <- matrix(0, 0, 0)
      phi_C_hat_orig <- matrix(0, 0, 0)
    }
    my_eigen_C <- NULL

    if(N_E != 0){
      ########################
      # rescale eigenfunctions
      ########################
      phi_E_hat_grid <- (1 / sqrt(interv)) * my_eigen_E$vectors[, seq(len = N_E), drop = FALSE]

      ############################
      # interpolate eigenfunctions
      # and evaluate on original
      # measurement points
      ############################
      phi_E_hat_orig <- matrix(NA, ncol = N_E, nrow = length(unlist(t)))
      for(k in 1 : N_E){
        phi_E_hat_orig[, k] <- approx(x = my_grid, y = phi_E_hat_grid[, k],
                                      xout = unlist(t), method = "linear")$y
      }
    }else{
      phi_E_hat_grid <- matrix(0, 0, 0)
      phi_E_hat_orig <- matrix(0, 0, 0)
    }
    my_eigen_E <- NULL

    ####################
    # esimate covariance
    # of basis weights
    ####################
    if(!use_RI){
      if(!nested) {
        N <- I * N_B + J * N_C + n * N_E
      } else {
        N <- I * N_B + I * J * N_C + n * N_E
      }

    }else{
      N <- I * N_B + n * N_E
    }
    if(N_B > 0){
      G_B <- diag(rep(nu_B_hat, times = I))
    }else{
      G_B <- matrix(NA, ncol = 0, nrow = 0)
    }
    if(N_C > 0){
      if(!nested) {
        G_C <- diag(rep(nu_C_hat, times = J))
      } else {
        G_C <- diag(rep(nu_C_hat, times = J*I))
      }
    }else{
      G_C <- matrix(NA, ncol = 0, nrow = 0)
    }
    if(N_E > 0){
      G_E <- diag(rep(nu_E_hat, times = n))
    }else{
      G_E <- matrix(NA, ncol = 0, nrow = 0)
    }
    G <- bdiag(G_B, G_C, G_E)

    ###################
    # invert covariance
    # of basis weights
    ###################
    G_inverse <- try(solve(G, sparse = TRUE))
    G <- NULL

    if(class(G_inverse)[1] != "try-error"){
      ##############################################
      # construct design matrix for EBLUP prediction
      ##############################################

      #######
      # for B
      if(N_B > 0){
        if(!use_RI){
          help_blocks <- data.table(subject_long = curve_info$subject_long,
                                    word_long = curve_info$word_long, phi_B_hat_orig)
          setorder(help_blocks, subject_long, word_long)
        }else{
          help_blocks <- data.table(subject_long = curve_info$subject_long, phi_B_hat_orig)
          setorder(help_blocks, subject_long)
        }

        blocks <- list()
        for(i in 1 : I){
          if(!use_RI){
            blocks[[i]] <- as.matrix(subset(help_blocks, subset = subject_long == i,
                                            select = -c(subject_long, word_long)), ncol = N_B)
          }else{
            blocks[[i]] <- as.matrix(subset(help_blocks, subset = subject_long == i,
                                            select = -c(subject_long)), ncol = N_B)
          }
        }
        phi_B_block <- bdiag(blocks)
      }else{
        phi_B_block <- matrix(0, 0, 0)
      }

      #######
      # for C
      if(N_C > 0){
        help_blocks <- data.table(subject_long = curve_info$subject_long,
                                  word_long = curve_info$word_long, phi_C_hat_orig)
        setorder(help_blocks, subject_long, word_long)
        blocks <- list()
        for(i in 1 : I){
          blocks[[i]] <- list()
          for(j in 1 : J){
            blocks[[i]][[j]] <- as.matrix(subset(help_blocks, subset = subject_long == i & word_long == j,
                                                 select = -c(subject_long, word_long)), ncol = N_C)
          }
          blocks[[i]] <- bdiag(blocks[[i]])
        }
        if(!nested){
          phi_C_block <- do.call("rbind", blocks)
        } else {
          phi_C_block <- bdiag(blocks)
        }

      }else{
        phi_C_block <- matrix(0, 0, 0)
      }

      ######
      #for E
      if(N_E > 0){
        if(!use_RI){
          help_blocks <- data.table(subject_long = curve_info$subject_long, word_long = curve_info$word_long,
                                    n_long = curve_info$n_long, phi_E_hat_orig)
          setorder(help_blocks, subject_long, word_long)
        }else{
          help_blocks <- data.table(subject_long = curve_info$subject_long,
                                    n_long = curve_info$n_long, phi_E_hat_orig)
          setorder(help_blocks, subject_long)
        }
        blocks <- list()
        for(i in seq_along(unique(help_blocks$n_long))){
          if(!use_RI){
            blocks[[i]] <- as.matrix(subset(help_blocks, subset = n_long == unique(help_blocks$n_long)[i],
                                            select = -c(n_long, subject_long, word_long)),
                                     ncol = N_E)
          }else{
            blocks[[i]] <- as.matrix(subset(help_blocks, subset = n_long == unique(help_blocks$n_long)[i],
                                            select = -c(n_long, subject_long)),
                                     ncol = N_E)
          }
        }
        phi_E_block <- bdiag(blocks)
      }else{
        phi_E_block <- matrix(0, 0, 0)
      }

      help_blocks <- NULL
      blocks <- NULL

      ################
      # combine blocks
      if(N_B > 0){
        if(N_C > 0){
          if(N_E > 0){
            phi_all <- cbind(phi_B_block, phi_C_block, phi_E_block)
          }else{
            phi_all <- cbind(phi_B_block, phi_C_block)
          }
        }else{
          if(N_E > 0){
            phi_all <- cbind(phi_B_block, phi_E_block)
          }else{
            phi_all <- phi_B_block
          }
        }
      }else{
        if(N_C > 0){
          if(N_E > 0){
            phi_all <- cbind(phi_C_block, phi_E_block)
          }else{
            phi_all <- phi_C_block
          }
        }else{
          phi_all <- phi_E_block
        }
      }

      ######################
      # compute bracket
      # for EBLUP prediction
      ######################
      # bracket  <-  sigmasq_hat * G_inverse + t(phi_all) %*% phi_all
      bracket <- sigmasq_hat * G_inverse + crossprod(phi_all)

      #########################
      # extract singular values
      #########################
      svd <- svd(bracket, nu = 0, nv = 0)

      ##########################
      # compute condition number
      ##########################
      cond <- abs(max(svd$d)) / abs(min(svd$d))

      ############################
      # get y_tilde in right order
      ############################
      curve_info_sort <- copy(curve_info)
      if(!use_RI){
        setorder(curve_info_sort, subject_long, word_long)
      }else{
        setorder(curve_info_sort, subject_long)
      }

      if(cond <= 1e+10){
        xi_all_hat <- solve(bracket, t(phi_all) %*% curve_info_sort$y_tilde)
      }else{
        bracket_inverse <- ginv(matrix(bracket, nrow = nrow(bracket), ncol = ncol(bracket)))
        xi_all_hat <- bracket_inverse %*% t(phi_all) %*% curve_info_sort$y_tilde
        warning("ginv is used in prediction of basis weights as bracket not invertible")
      }

      ########################################
      # determine which basis weights belong
      # to what level based on original levels
      # in curve_info
      ########################################
      if(N_B > 0){
        xi_B_hat <- matrix(xi_all_hat[1 : (N_B * I)], ncol = N_B, byrow = T)
        row.names(xi_B_hat) <- unique(curve_info_sort$subject_long_orig)
        if(N_C > 0){
          if(!nested){
            xi_C_hat <- matrix(xi_all_hat[(N_B * I + 1) : (N_B * I + N_C * J)], ncol = N_C, byrow = T)
            row.names(xi_C_hat) <- unique(curve_info_sort$word_long_orig)
          } else {
            xi_C_hat <- matrix(xi_all_hat[(N_B * I + 1) : (N_B * I + N_C * J * I)], ncol = N_C, byrow = T)
            row.names(xi_C_hat) <- unique(paste(curve_info_sort$subject_long_orig,
                                                curve_info_sort$word_long_orig,
                                                sep = "."))
          }
          if(N_E > 0){
            if(!nested){
              xi_E_hat <- matrix(xi_all_hat[(N_B * I + N_C * J + 1) : N], ncol = N_E, byrow = T)
              row.names(xi_E_hat) <- unique(curve_info_sort$n_long_orig)
            } else {
              xi_E_hat <- matrix(xi_all_hat[(N_B * I + N_C * J * I + 1) : N], ncol = N_E, byrow = T)
              row.names(xi_E_hat) <- unique(curve_info_sort$n_long_orig)
            }
          }else{
            xi_E_hat <- rep(NA, n)
          }
        }else{
          if(use_RI){
            xi_C_hat <- NA
          }else{
            xi_C_hat <- rep(NA, J)
          }

          if(N_E > 0){
            xi_E_hat <- matrix(xi_all_hat[(N_B * I + 1) : N], ncol = N_E, byrow = T)
            row.names(xi_E_hat) <- unique(curve_info_sort$n_long_orig)
          }else{
            xi_E_hat <- rep(NA, n)
          }
        }
      }else{
        xi_B_hat <- rep(NA, I)
        if(N_C > 0){
          xi_C_hat <- matrix(xi_all_hat[1 : (N_C * J)], ncol = N_C, byrow = T)
          row.names(xi_C_hat) <- unique(curve_info_sort$word_long_orig)
          if(N_E > 0){
            xi_E_hat <- matrix(xi_all_hat[(N_C * J + 1) : N], ncol = N_E, byrow = T)
            row.names(xi_E_hat) <- unique(curve_info_sort$n_long_orig)
          }else{
            xi_E_hat <- rep(NA, n)
          }
        }else{
          if(use_RI){
            xi_C_hat <- NA
          }else{
            if(!nested) {
              xi_C_hat <- rep(NA, J)
            } else {
              xi_C_hat <- rep(NA, J*I)
            }
          }
          xi_E_hat <- matrix(xi_all_hat[1 : N], ncol = N_E, byrow = T)
          row.names(xi_E_hat) <- unique(curve_info_sort$n_long_orig)
        }
      }
      phi_B_hat_orig <- NULL
      phi_C_hat_orig <- NULL
      phi_E_hat_orig <- NULL
    }else{
      warning("basis weights cannot be computed due to inversion of their covariance")
      xi_B_hat <- NA
      xi_C_hat <- NA
      xi_E_hat <- NA
    }

  }else{
    xi_B_hat <- NA
    xi_C_hat <- NA
    xi_E_hat <- NA
    phi_B_hat_grid <- NA
    phi_C_hat_grid <- NA
    phi_E_hat_grid <- NA
    print(warning("no FPCs chosen at all"))
  }


  ##############
  # return ouput
  ##############
  results <- list(phi_B_hat_grid = phi_B_hat_grid, phi_C_hat_grid = phi_C_hat_grid, phi_E_hat_grid = phi_E_hat_grid,
                  nu_B_hat = nu_B_hat, nu_C_hat = nu_C_hat, nu_E_hat = nu_E_hat, N_B= N_B, N_C = N_C, N_E = N_E,
                  total_var = total_var, var_explained  =  var_explained,
                  xi_B_hat = xi_B_hat, xi_C_hat = xi_C_hat, xi_E_hat = xi_E_hat)
  return(results)
}
########################################################################################################################
