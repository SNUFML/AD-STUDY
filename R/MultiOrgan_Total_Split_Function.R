group_na_summary  <- function(data) {
  grouped_data <- split(data, data$Group)
  mark_na_columns <- function(df) {
    na_ratio <- sapply(df, function(col) mean(is.na(col)))
    na_ratio <= 0.2
  }
  na_marked <- lapply(grouped_data, mark_na_columns)
  na_matrix <- as.data.frame(na_marked)
  trans_data <- as.data.frame(t(na_matrix))
  trans_data
}


impu <- function(data) {
  # 데이터의 열 개수
  Data_col <- ncol(data)
  
  # 3번째 열부터 마지막 열까지의 데이터 선택
  Data <- data[, 3:Data_col]
  
  # 0 값을 확인
  if (any(Data == 0, na.rm = TRUE)) {
    # 0을 NA로 대체
    Data[Data == 0] <- NA
    
    # 전체 최소값 계산 (NA 제외)
    overall_minimum <- min(Data, na.rm = TRUE)
    
    # NA 값을 전체 최소값의 1/5로 대체
    Data[is.na(Data)] <- overall_minimum / 5
    
    # 조합된 데이터프레임 생성
    impu_Data <- cbind(data[, 1, drop = FALSE], data[, 2, drop = FALSE], Data)
    
    # 결과를 CSV 파일로 저장
    write.csv(impu_Data, "Imputed_data.csv", row.names = FALSE)
    
    # 처리된 데이터 반환
    return(impu_Data)
  } else {
    # 0 값이 없는 경우 메시지 출력
    message("0값이 없습니다.")
    return(data)
  }
}

YMSTAT<-function (Data,group_combinations=group_combinations,Adjust_p_value = TRUE, Adjust_method = "BH", 
                  parallel = FALSE) {  
  colnames(Data)[1:2] <- c("Sample", "Group")
  Data$Group <- as.character(Data$Group)
  Data <- dplyr::arrange(Data, Data$Group)
  Data_ori <- Data
  nmet <- ncol(Data) - 2
  nmet_seq <- seq_len(nmet)
  cols_to_convert <- 3:ncol(Data)
  Data[, cols_to_convert] <- lapply(Data[, cols_to_convert], 
                                    as.numeric)
  Data_final <- data.table::as.data.table(Data)
  Data_final_raw <- Data_final[, setdiff(names(Data_final), 
                                         c("Sample", "Group")), with = FALSE]
  data.table::set(Data_final, j = "Group", value = as.factor(as.character(Data_final$Group)))
  groups_split <- split(Data_final_raw, Data_final$Group)
  group_nottwo <- length(unique(Data_final$Group)) > 2
  Data_final <- as.data.frame(Data_final)
  colnames(Data_final) <- c("Sample", "Group", paste0("V", 
                                                      nmet_seq))
  if (parallel) {
    num_core <- parallel::detectCores() - 1
    num_col <- ncol(Data_final)
    if (num_col < 1000) {
      parallel <- FALSE
    }
    else if (num_col < 5000) {
      cl <- parallel::makeCluster(3)
    }
    else {
      cl <- parallel::makeCluster(min(5, num_core))
    }
    if (parallel) {
      on.exit(parallel::stopCluster(cl))
    }
  }
  group_matrices <- lapply(groups_split, as.matrix)
  # 추가된 Cohen's d 함수
  vectorized_cohens_d <- function(mat1, mat2) {
    mean_diff <- colMeans(mat1) - colMeans(mat2)
    pooled_sd <- sqrt(((nrow(mat1) - 1) * apply(mat1, 2, var) + 
                         (nrow(mat2) - 1) * apply(mat2, 2, var)) / 
                        (nrow(mat1) + nrow(mat2) - 2))
    d_values <- mean_diff / pooled_sd
    return(d_values)
  }
  # 추가된 Hedge's g 함수
  vectorized_hedges_g <- function(mat1, mat2) {
    cohens_d <- vectorized_cohens_d(mat1, mat2)
    correction_factor <- 1 - (3 / (4 * (nrow(mat1) + nrow(mat2)) - 9))
    hedges_g <- cohens_d * correction_factor
    return(hedges_g)
  }
  
  # 추가된 Fold Change 계산 함수
  calculate_fold_change <- function(mat1, mat2) {
    mean1 <- colMeans(mat1)
    mean2 <- colMeans(mat2)
    return(mean1 / mean2)
  }
  
  # 추가된 Odds Ratio 계산 함수
  
  calculate_cliffs_delta <- function(mat1, mat2) {
    sapply(1:ncol(mat1), function(i) {
      effsize::cliff.delta(mat1[, i], mat2[, i])$estimate
    })
  }
  
  vectorized_t_test <- function(mat1, mat2) {
    n1 <- nrow(mat1)
    n2 <- nrow(mat2)
    means1 <- colMeans(mat1)
    means2 <- colMeans(mat2)
    vars1 <- apply(mat1, 2, var)
    vars2 <- apply(mat2, 2, var)
    t_stats <- (means1 - means2)/sqrt(vars1/n1 + vars2/n2)
    df <- (((vars1/n1) + (vars2/n2))^2)/((vars1/n1)^2/(n1 - 
                                                         1) + (vars2/n2)^2/(n2 - 1))
    p_values <- 2 * (1 - pt(abs(t_stats), df = df))
    return(p_values)
  }

  if (parallel) {
    result_list_t <- parallel::parLapply(cl, group_combinations, 
                                         function(combo) {
                                           mat1 <- group_matrices[[combo[1]]]
                                           mat2 <- group_matrices[[combo[2]]]
                                           list(
                                             t_test = vectorized_t_test(mat1, mat2),
                                             cohens_d = vectorized_cohens_d(mat1, mat2),
                                             hedges_g = vectorized_hedges_g(mat1, mat2),
                                             fold_change = calculate_fold_change(mat1, mat2),
                                             cliffs_delta = calculate_cliffs_delta(mat1, mat2))
                                        
                                         })
  }  else {
    result_list_t <- lapply(group_combinations, function(combo) {
      mat1 <- group_matrices[[combo[1]]]
      mat2 <- group_matrices[[combo[2]]]
      list(
        t_test = vectorized_t_test(mat1, mat2),
        cohens_d = vectorized_cohens_d(mat1, mat2),
        hedges_g = vectorized_hedges_g(mat1, mat2),
        fold_change = calculate_fold_change(mat1, mat2),
        cliffs_delta = calculate_cliffs_delta(mat1, mat2)
      )
      #vectorized_t_test(mat1, mat2)
    })
  }
  # T-test, Cohen's d, Fold Change, Odds Ratio 결과 데이터 프레임 생성
  df_ttest <- (data.frame(sapply(result_list_t, `[[`, "t_test")))
  df_cohens_d <- (data.frame(sapply(result_list_t, `[[`, "cohens_d")))
  df_hedges_g <- (data.frame(sapply(result_list_t, `[[`, "hedges_g")))
  df_cliffs_delta <- (data.frame(sapply(result_list_t, `[[`, "cliffs_delta")))
  df_fold_change <- (data.frame(sapply(result_list_t, `[[`, "fold_change")))
  print("T-test has finished")
  vectorized_u_test <- function(mat1, mat2) {
    result <- matrixTests::col_wilcoxon_twosample(mat1, mat2)
    return(result$pvalue)
  }
  if (parallel) {
    result_list_u <- parallel::parLapply(cl, group_combinations, 
                                         function(combo) {
                                           mat1 <- group_matrices[[combo[1]]]
                                           mat2 <- group_matrices[[combo[2]]]
                                           list(
                                             U_test = vectorized_u_test(mat1, mat2),
                                             cohens_d = vectorized_cohens_d(mat1, mat2),
                                             hedges_g = vectorized_hedges_g(mat1, mat2),
                                             fold_change = calculate_fold_change(mat1, mat2),
                                             cliffs_delta = calculate_cliffs_delta(mat1, mat2))
                                         })
  } else {
    result_list_u <- lapply(group_combinations, function(combo) {
      mat1 <- group_matrices[[combo[1]]]
      mat2 <- group_matrices[[combo[2]]]
      list(
        u_test = vectorized_u_test(mat1, mat2),
        cohens_d = vectorized_cohens_d(mat1, mat2),
        hedges_g = vectorized_hedges_g(mat1, mat2),
        fold_change = calculate_fold_change(mat1, mat2),
        cliffs_delta = calculate_cliffs_delta(mat1, mat2))
    })
  }
  df_utest <-  (data.frame(sapply(result_list_u, `[[`, "u_test")))
  df_cohens_d <- (data.frame(sapply(result_list_u, `[[`, "cohens_d")))
  df_hedges_g <- (data.frame(sapply(result_list_t, `[[`, "hedges_g")))
  df_cliffs_delta <- (data.frame(sapply(result_list_u, `[[`, "cliffs_delta")))
  df_fold_change <- (data.frame(sapply(result_list_u, `[[`, "fold_change")))
  print("U-test has finished")
  if (group_nottwo) {
    metabolite_names <- colnames(Data_final)[nmet_seq + 2]
    perform_anova_tests <- function(data_subset) {
      colname <- names(data_subset)[2]
      formula_str <- as.formula(paste0(colname, " ~ Group"))
      model <- aov(formula_str, data = data_subset)
      anova_res <- summary(model)[[1]]$"Pr(>F)"[1]
      posthoc_res <- DescTools::PostHocTest(model, method = "scheffe")$Group[, 
                                                                             4]
      list(anova_res = anova_res, posthoc_res = posthoc_res)
    }
    if (parallel) {
      results <- parallel::parLapply(cl, metabolite_names, 
                                     function(col) {
                                       subset_data <- Data_final[, c("Group", col)]
                                       perform_anova_tests(subset_data)
                                     })
    }
    else {
      results <- lapply(metabolite_names, function(col) {
        subset_data <- Data_final[, c("Group", col)]
        perform_anova_tests(subset_data)
      })
    }
    p_anova <- sapply(results, function(x) {
      x$anova_res
    })
    df_anova <- data.frame(p_anova)
    df_anova_post <- do.call(rbind, lapply(results, function(x) {
      x$posthoc_res
    }))
    print("Anova & PostHoc has finished")
    kw_results <- matrixTests::col_kruskalwallis(as.matrix(Data_final[metabolite_names]), 
                                                 Data_final$Group)
    perform_posthoc_tests <- function(data_subset) {
      colname <- names(data_subset)[2]
      formula <- as.formula(paste0(colname, " ~ Group"))
      dunn_res <- FSA::dunnTest(formula, data = data_subset, 
                                method = "none")
      post_pvals <- dunn_res[["res"]][["P.unadj"]]
      adjusted_pvals <- p.adjust(post_pvals, method = "BH")
      list(adjusted_pvals = adjusted_pvals, dunn_res = dunn_res)
    }
    if (parallel) {
      results_kw_posthoc <- parallel::parLapply(cl, metabolite_names, 
                                                function(col) {
                                                  subset_data <- Data_final[, c("Group", col)]
                                                  perform_posthoc_tests(subset_data)
                                                })
    }
    else {
      results_kw_posthoc <- lapply(metabolite_names, function(col) {
        subset_data <- Data_final[, c("Group", col)]
        perform_posthoc_tests(subset_data)
      })
    }
    df_kw <- data.frame(p_kw = kw_results$pvalue)
    df_kw_post <- do.call(rbind, lapply(results_kw_posthoc, 
                                        function(x) {
                                          x$adjusted_pvals
                                        }))
    print("Kruskal Wallis & PostHoc has finished")
  }
  Names <- NULL
  Names <- matrix(sapply(group_combinations, function(x) paste(x, collapse = "-")))
  rownamechange <- colnames(Data)[nmet_seq + 2]
  rownames <- rownamechange
  rownames(df_cohens_d) <- rownamechange
  rownames(df_fold_change) <- rownamechange
  rownames(df_cliffs_delta) <- rownamechange
  rownames(df_hedges_g) <- rownamechange
  rownames(df_ttest) <- rownamechange
  rownames(df_utest) <- rownamechange
  colnames(df_utest) <- paste(Names[, 1], "u-test", sep = "___")
  colnames(df_ttest) <- paste(Names[, 1], "t-test", sep = "___")
  colnames(df_cohens_d) <- paste(Names[, 1], "cohens_d", sep = "___")
  colnames(df_fold_change) <- paste(Names[, 1], "fold_change", sep = "___")
  colnames(df_cliffs_delta) <- paste(Names[, 1], "cliffs_delta", sep = "___")
  colnames(df_hedges_g) <- paste(Names[, 1], "hedges_g", sep = "___")
  
  if (group_nottwo) {
    AN_post_names <- colnames(df_anova_post)
    DU_post_names <- results_kw_posthoc[[1]]$dunn_res$res$Comparison
    rownames(df_anova) <- rownamechange
    rownames(df_anova_post) <- rownamechange
    rownames(df_kw) <- rownamechange
    rownames(df_kw_post) <- rownamechange
    colnames(df_anova) <- "Anova"
    colnames(df_anova_post) <- paste(AN_post_names, "ANO_posthoc", 
                                     sep = "___")
    colnames(df_kw) <- "Kruskal_Wallis"
    colnames(df_kw_post) <- paste(DU_post_names, "Kru_posthoc(Dunn)", 
                                  sep = "___")
  }
  if (Adjust_p_value == TRUE) {
    print("###########################################")
    print(paste0("adjusted according to the ", Adjust_method, 
                 " method"))
    print("###########################################")
    adj_func <- function(x) {
      p.adjust(x, method = Adjust_method)
    }
    df_ttest <- apply(df_ttest, 2, adj_func)
    df_utest <- apply(df_utest, 2, adj_func)
    if (group_nottwo) {
      df_anova <- apply(df_anova, 2, adj_func)
      df_kw <- apply(df_kw, 2, adj_func)
    }
  } else {
    print("###########################################")
    print("p_value not adjusted")
    print("###########################################")
  }
  print("statistical test has finished")
  Result <- cbind(df_ttest, df_utest, df_cohens_d, df_hedges_g, df_cliffs_delta,df_fold_change)
  if (group_nottwo) {
    Result <- cbind(Result, df_anova, df_anova_post, df_kw, 
                    df_kw_post)
  }
  Result <- as.matrix(Result)
  if (group_nottwo) {
    df_anova <- as.data.frame(df_anova)
    df_kw <- as.data.frame(df_kw)
    df_kw_post <- as.matrix(df_kw_post)
  }
  Final <- list()
  Final$Data <- Data_ori
  Final$Data_renamed <- Data_final
  Final$Result <- Result
  if (group_nottwo) {
    Final$Anova <- df_anova
    Final$Anova_PostHoc <- df_anova_post
    Final$KW <- df_kw
    Final$Dunn <- df_kw_post
    Final$t_test <- as.data.frame(df_ttest)
    Final$u_test <- as.data.frame(df_utest)
    Final$cohens_d <- as.data.frame(df_cohens_d)
    Final$cliffs_delta <- as.data.frame(df_cliffs_delta)
    Final$hedges_g<- as.data.frame(df_hedges_g)
    Final$fold_change <- as.data.frame(df_fold_change)
  } else {
    Final$t_test <- df_ttest[, 1]
    Final$u_test <- df_utest[, 1]
    Final$cohens_d<-df_cohens_d
    Final$hedges_g<-df_hedges_g
    Final$fold_change<-df_fold_change
    
  }
  stat<-Final
}

