.onLoad <- function(libname = find.package("sparseFLMM"), pkgname = "sparseFLMM"){

  # CRAN note avoidance
  utils::globalVariables(c("col_combi_bivariate", "col_curve_bivariate", "col_subject_bivariate",
                                "col_t_bivariate", "col_word_bivariate", "combi_long", "cross_vec_bivariate", "id",
                                "id1", "id2", "id_n.vec", "id_subject.vec", "id_word.vec", "n_long", "n_long_orig",
                                "row_combi_bivariate", "row_curve_bivariate", "row_subject_bivariate",
                                "row_t_bivariate", "row_word_bivariate", "same_curve", "same_point", "same_subject",
                                "same_word", "subj1", "subj2", "subject_long", "subject_long_orig", "word1", "word2",
                                "word_long", "word_long_orig", "y1", "y2", "y_tilde", "n1", "n2"))
  invisible()
}

