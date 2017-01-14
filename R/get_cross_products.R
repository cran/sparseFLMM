#####################################################################################################
# author: Jona Cederbaum and Fabian Scheipl
#####################################################################################################
# description: preparations for the covariance where only cross products of interest are constructed.
# uses functions in useful_functions.R.
#####################################################################################################

get_crossprods_fun <- function(y_tilde, curve_info, t, my_grid, d_grid, use_RI, I, J){
  
  ###################
  # initialize output
  ###################
  output <- list()  
  
  if(!use_RI){
    ##################
    # for crossed fRIs
    res <- vector(mode = "list", length = I + J)  
  }else{
    ###################################
    # for one fRI or independent curves
    res <- vector(mode = "list", length = I )  
  }
  
  ##################################################
  # loop over subjects (first groupin_variable)
  # get all combinations on the same subjects
  # (same words or different words)
  # calls function make_crossprod_dt (defined below)
  ##################################################
  for(i in seq_len(I)){
    res[[i]] <- make_crossprod_dt(curve_info[subject_long == i, ], use_RI = use_RI)
  }
  
  if(!use_RI){  
    ##################################################
    # loop over words (second grouping variable)
    # get all combinations on the same words
    # (same words or different subjects)
    # calls function make_crossprod_dt (defined below)
    ##################################################
    for(i in seq_len(J)){
      res[[I + i]] <- make_crossprod_dt(curve_info[word_long == i, ], 
                                        preselection = "word", use_RI = use_RI)
    }
  }
  
  ret <- do.call(rbind, res)
  setkey(ret, id1, id2)
  
  #####################
  # take out id, y1, y2
  #####################
  set(ret, i = NULL, "id1", NULL)
  set(ret, i = NULL, "id2", NULL)
  
  set(ret, i = NULL, "y1", NULL)
  set(ret, i = NULL, "y2", NULL)
  
  ##################
  # rename t1 and t2
  ##################
  setnames(ret, old = c("t1", "t2"), new = c("row_t_bivariate", "col_t_bivariate"))
  
  #####################################
  # create indicators
  # same_word, same_subject, same_point
  #####################################
  if(!use_RI){
    ##################
    # for crossed fRIs
    output[["index"]] <- create_data_frame_bivariate_fun(index = ret)   
  }else{
    ###################################
    # for one fRI or independent curves
    output[["index"]] <- create_data_frame_bivariate_RI_fun(index = ret)   
  }
  
  #####################
  # construct grid data
  #####################
  grid_help <- create_grid_data_fun(my_grid = my_grid, d_grid = d_grid)
  
  output[["grid_row"]] <- grid_help$grid_row
  output[["grid_col"]] <- grid_help$grid_col
  
  output[["same_subject_grid"]] <- grid_help$same_subject
  if(!use_RI)
    output[["same_word_grid"]] <- grid_help$same_word
  output[["same_curve_grid"]] <- grid_help$same_curve_grid
  output[["same_point_grid"]] <- grid_help$same_point_grid
  
  rm(grid_help)
  
  ###############
  # return output
  ###############
  output
}  

############################################################################################

make_crossprod_dt <- function(curve_info, preselection = c("none", "subject", "word"), use_RI){
  
  preselection <- match.arg(preselection)
  setkey(curve_info, id)
  
  ####################
  # take  combinations
  ####################
  combinations <- with(curve_info, CJ(id = id, id2 = id))
  
  if(!use_RI){
    tmp1 <- curve_info[combinations, list(id1 = id, subj1 = subject_long, word1 = word_long, 
                                          rep1 = combi_long)]  
  }else{
    tmp1 <- curve_info[combinations, list(id1 = id, subj1 = subject_long, n1 = n_long)]  
  }
  
  if(!use_RI){
    tmp2 <- curve_info[combinations[, list(id = id2)], 
                       list(id2 = id, subj2 = subject_long, word2 = word_long, rep2 = combi_long)]
  }else{
    tmp2 <- curve_info[combinations[, list(id = id2)], 
                       list(id2 = id, subj2 = subject_long, n2 = n_long)]
  } 
  
  if(!use_RI){
    crosstable <- tmp1[, `:=`(id2 = tmp2$id2, subj2 = tmp2$subj2, word2 = tmp2$word2, rep2 = tmp2$rep2)]
  }else{
    crosstable <- tmp1[, `:=`(id2 = tmp2$id2, subj2 = tmp2$subj2, n2 = tmp2$n2)]
  }
  
  ################################
  # remove irrelevant combinations
  ################################
  # NOTE: if preselection is "none" all subjects and words are used
  ## if preselection is "subject" only the subset of crosstable is used for which word1 != word2
  ## if preselection is "word" only  the subset of crosstable is used for which subj1 != subj2
  
  if(preselection == "none") {
    
  }
  if(preselection == "subject") {
    crosstable <- crosstable[word1 != word2, ]
  }    
  if(preselection == "word") {
    crosstable <- crosstable[subj1 != subj2, ]  
  }
  
  ###############################
  # add y and t to the data.table
  ###############################
  # once for id1 ordering and once for id2 ordering
  # leading to t1, t2, y1, y2 (when use_tri = TRUE also an indicator for t is added for both ordering)
  # NOTE: using the id key, we can directly assign the right values
  
  # add t1, y1, and (t_ind1)
  crosstable[, id := id1]
  setkey(crosstable, "id")
  crosstable <- crosstable[curve_info[, list(id, y_tilde, t)], ]
  setnames(crosstable, old = c("id1", "y_tilde", "t"), 
           new = c("id1", "y1", "t1"))
  
  # add t2, y2, and (t_ind2)
  crosstable[, id := id2]
  setkey(crosstable, "id")
  crosstable <- crosstable[curve_info[, list(id, y_tilde, t)], ]
  if(!use_RI){
    setnames(crosstable, old = c("id", "y_tilde", "t", "subj1", "subj2", "word1", "word2", "rep1", "rep2"), 
             new = c("id2", "y2", "t2", "row_subject_bivariate", 
                   "col_subject_bivariate", "row_word_bivariate", "col_word_bivariate", "row_combi_bivariate", "col_combi_bivariate"))
  }else{
    setnames(crosstable, old = c("id", "y_tilde", "t", "subj1", "subj2", "n1", "n2"), 
             new = c("id2", "y2", "t2", "row_subject_bivariate", 
                   "col_subject_bivariate", "row_curve_bivariate", 
                   "col_curve_bivariate"))
  }
  
  #################################
  # redo sort (again sorted by id1)
  #################################
  setkey(crosstable, "id1") 
  
  (crosstable[, cross_vec_bivariate := y1 * y2])
  
  ##############################
  # sort data table first by id1
  # and then by id2
  ##############################
  setkey(crosstable, id1, id2)
  
  ###############
  # return output
  ###############
  crosstable
}

