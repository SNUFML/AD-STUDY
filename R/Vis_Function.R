PCA_YM_fviz <- function(Data, 
                        color = c("#FF3300", "#660099", "#FFCC00", "#99CC00", "#0066CC", "#FF6600"), 
                        legend_position = "none", 
                        fig_width = 24, 
                        fig_height = 20, 
                        components = c(1, 2), 
                        names = F, 
                        dotsize = 3, 
                        labsize = 3, 
                        ellipse = T, 
                        order = NULL, 
                        title = "PCA", 
                        legend_title = "Group",
                        title_size = 30, 
                        legend_size = 15, 
                        legend_title_size = 18, 
                        axis_text_size = 12, 
                        axis_title_size = 14, 
                        border_color = "black", 
                        border_size = 1) {
  # Check and create directory
  if (!dir.exists(file.path(getwd(), "PCA"))) {
    dir.create(file.path(getwd(), "PCA"))
  }
  
  # Data preprocessing
  Data <- Data %>% dplyr::arrange(Group)
  if (!is.null(order)) {
    Data$Group <- factor(Data$Group, levels = order)
  }
  
  rownames(Data) <- Data[, 1]
  Data <- Data[, -1]
  Data[, 2:ncol(Data)] <- sapply(Data[, 2:ncol(Data)], as.numeric)
  
  
  # Modify Group levels to include superscripts
  group_levels <- levels(Data$Group)
  group_levels_parsed <- lapply(group_levels, function(x) parse(text = x))  # Parse group labels
  Data$Group <- factor(Data$Group, levels = group_levels)
  
  PCA_Data_input <- Data
  PCA_Data_input[, 1] <- as.factor(PCA_Data_input[, 1])
  
  # PCA analysis
  Bef_merge <- prcomp(PCA_Data_input[, 2:ncol(PCA_Data_input)], scale = TRUE)
  explained_variance <- Bef_merge$sdev^2 / sum(Bef_merge$sdev^2)
  
  # PCA plot
  plot <- factoextra::fviz_pca_ind(Bef_merge,
                                   geom.ind = "point",
                                   mean.point = FALSE,
                                   pointshape = 21,
                                   pointsize = dotsize,
                                   fill.ind = Data$Group,
                                   col.ind = "black",
                                   palette = color,
                                   addEllipses = ellipse,
                                   ellipse.level = 0.68,
                                   legend.title = "Group",
                                   axes.labs = c(paste("PC", components[1]), paste("PC", components[2]))) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(color = border_color, linewidth = border_size, fill = NA),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold"),
                   legend.position = legend_position,
                   legend.text = ggplot2::element_text(size = legend_size),
                   legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = axis_text_size),
                   axis.text.y = ggplot2::element_text(size = axis_text_size),
                   axis.title.x = ggplot2::element_text(size = axis_title_size),
                   axis.title.y = ggplot2::element_text(size = axis_title_size)) 
  # Save the plot
  ggplot2::ggsave(paste(components[1], components[2], "label", names, "PCA.pdf", sep = "_"), 
                  path = file.path(getwd(), "PCA"), 
                  height = fig_height, 
                  width = fig_width, 
                  units = "cm")
  
  # Results
  Result <- list()
  Result$explained_variance <- explained_variance
  Result$coordinates <- as.data.frame(Bef_merge$x)
  Result$plot <- plot
  
  return(Result)
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
  # t_test, Cohen's d, Fold Change, Odds Ratio 결과 데이터 프레임 생성
  df_ttest <- (data.frame(sapply(result_list_t, `[[`, "t_test")))
  df_cohens_d <- (data.frame(sapply(result_list_t, `[[`, "cohens_d")))
  df_hedges_g <- (data.frame(sapply(result_list_t, `[[`, "hedges_g")))
  df_cliffs_delta <- (data.frame(sapply(result_list_t, `[[`, "cliffs_delta")))
  df_fold_change <- (data.frame(sapply(result_list_t, `[[`, "fold_change")))
  print("t_test has finished")
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
  colnames(df_ttest) <- paste(Names[, 1], "t_test", sep = "___")
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






Volcano_YM_S <- function(data, asterisk = "t_test", comb, reverse = FALSE, fig_width = 7.5, 
                         FC_log = 2, pval_log = 10, fig_height = 5, dotsize = 2.5, 
                         fixed_limit = FALSE, y_limit = c(-2, 6), pval_intercept = 0.05, 
                         max_overlap = 9, FC_range = c(-1.5, 1.5), sig_label = TRUE, 
                         color = c("blue", "black", "red")) {
  
  # Volcano 폴더 생성
  if (!dir.exists(file.path(getwd(), "Volcano"))) {
    dir.create(file.path(getwd(), "Volcano"))
  }
  
  ## 1. ASTERISK 값에 따른 p-value 데이터 선택
  if (asterisk == "Scheffe") {
    p_val_data <- data[["Anova_PostHoc"]]
  } else if (asterisk == "t_test") {
    p_val_data <- data[["t_test"]]
  } else if (asterisk == "u_test") {
    p_val_data <- data[["u_test"]]
  } else if (asterisk == "Dunn") {
    p_val_data <- data[["Dunn"]]
  } else {
    stop("asterisk 인자는 (Scheffe, t_test, u-test, Dunn) 중 하나여야 합니다.")
  }
  
  ## 2. LOG 변환된 FC_range 값 설정
  log_FC_range_1 <- log(-1/FC_range[1], base = FC_log)  # DOWN 판별 기준 (음수)
  log_FC_range_2 <- log(FC_range[2], base = FC_log)         # UP 판별 기준
  
  ## 3. 각 비교에 대해 Volcano plot 작성 (comb는 리스트 형태)
  for (i in seq_along(comb)) {
    # 그룹 순서 결정 (reverse 인자에 따라)
    if (reverse) {
      G1 <- comb[[i]][2]
      G2 <- comb[[i]][1]
    } else {
      G1 <- comb[[i]][1]
      G2 <- comb[[i]][2]
    }
    
    ## 4. YMSTAT에서 미리 계산된 fold_change 사용  
    colname1 <- paste0(G1, "-", G2, "___fold_change")
    colname2 <- paste0(G2, "-", G1, "___fold_change")
    
    if (colname1 %in% colnames(data[["fold_change"]])) {
      fold_change_vec <- data[["fold_change"]][, colname1]
    } else if (colname2 %in% colnames(data[["fold_change"]])) {
      fold_change_vec <- data[["fold_change"]][, colname2]
    } else {
      stop("해당 비교에 대한 fold_change 데이터를 찾을 수 없습니다.")
    }
    
    # Log 변환된 Fold Change 값 계산
    logFC <- log(fold_change_vec, base = FC_log)
    
    ## 5. P-value 데이터 선택 (YMSTAT 결과)
    if (length(unique(data[["Data"]]$Group)) > 2) {
      colname1_p <- paste0(G1, "-", G2, "___", asterisk)
      colname2_p <- paste0(G2, "-", G1, "___", asterisk)
      
      if (colname1_p %in% colnames(p_val_data)) {
        part_p_val <- p_val_data[, colname1_p]
      } else {
        part_p_val <- p_val_data[, colname2_p]
      }
    } else {
      part_p_val <- p_val_data
    }
    
    ## 6. Volcano 데이터 프레임 구성  
    Volc_Dat <- data.frame(
      log2FoldChange = logFC,
      pvalue = part_p_val,
      stringsAsFactors = FALSE
    )
    ## 8. 대사체 이름 적용 (최대 25자까지 표시)
    metabolite_names <- substr(colnames(data[["Data"]][, 3:ncol(data[["Data"]])]), 1, 200)
    rownames(Volc_Dat) <- metabolite_names[1:nrow(Volc_Dat)]
    
    # 1) p-value threshold를 만족하지 않는 데이터 삭제
    Volc_Dat <- Volc_Dat[Volc_Dat$pvalue < pval_intercept, ]
    
    ## 7. 유의미성 분류: LOG 변환된 FC_range 적용
    Volc_Dat$Direction <- "Not Significant"
    Volc_Dat$Direction[Volc_Dat$log2FoldChange > log_FC_range_2] <- "Up"
    Volc_Dat$Direction[Volc_Dat$log2FoldChange < log_FC_range_1] <- "Down"
    
    ## -log 변환 p-value
    Volc_Dat$pval_trans <- -log(Volc_Dat$pvalue, base = pval_log)
    
    
    
    ## 9. Volcano plot 작성
    # x축 최대 절대값 계산 후, 2의 배수로 올림
    max_abs <- max(abs(Volc_Dat$log2FoldChange), na.rm = TRUE)
    
    sem_res <- ggplot2::ggplot(data = Volc_Dat, ggplot2::aes(x = log2FoldChange, y = pval_trans, col = Direction)) +
      # 3) 닷의 크기 증가 및 4) 투명도 0.2 적용
      ggplot2::geom_point(size = dotsize, alpha = 0.5) +
      ggplot2::labs(title = paste0(G1, " divided by ", G2, " Volcano plot"),
                    x = eval(parse(text = paste0("expression(Log[", FC_log, "]~(Fold~Change))"))),
                    y = eval(parse(text = paste0("expression(Log[", pval_log, "]~(P-value))")))) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 0, color = "Black", hjust = 0.5, lineheight = 1.2),
                     axis.title.x = ggplot2::element_text(size = 15),
                     axis.title.y = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.text.y = ggplot2::element_text(size = 10),
                     legend.text = ggplot2::element_text(size = 10),
                     legend.title = ggplot2::element_text(size = 12),
                     legend.position = "bottom") +
      ggplot2::scale_color_manual(values = c("Down" = "blue", "Not Significant" = "grey", "Up" = "red")) +
      ggplot2::geom_vline(xintercept = 0, col = "black") +
      ggplot2::geom_vline(xintercept = log_FC_range_1, linetype = "dashed", col = "grey") +
      ggplot2::geom_vline(xintercept = log_FC_range_2, linetype = "dashed", col = "grey") +
      ggplot2::geom_hline(yintercept = -log(pval_intercept, base = pval_log), col = "black") 
    sem_res <- sem_res + ggplot2::xlim(-max_abs, max_abs)
    
    if (fixed_limit && !is.na(y_limit)) {
      sem_res <- sem_res + ggplot2::ylim(y_limit)
    }
    
    ## 10. 대사체 이름 레이블 추가 (최대 20자; 글씨색 검정)
    
    if (sig_label) {
      # UP 레이블 데이터
      UP_labels <- subset(Volc_Dat, Direction == "Up")
      # DOWN 레이블 데이터
      DOWN_labels <- subset(Volc_Dat, Direction == "Down")
      
      sem_res <- sem_res +
        ggrepel::geom_text_repel(
          data = UP_labels,
          aes(x = log2FoldChange, y = pval_trans, 
              label = substr(rownames(UP_labels), 1, 200)),
          color = "darkred", size = 3, show.legend = FALSE, 
          max.overlaps = max_overlap
        ) +
        ggrepel::geom_text_repel(
          data = DOWN_labels,
          aes(x = log2FoldChange, y = pval_trans, 
              label = substr(rownames(DOWN_labels), 1, 200)),
          color = "navy", size = 3, show.legend = FALSE, 
          max.overlaps = max_overlap
        ) +ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    }
    ## 11. 저장 (5) 피겨의 너비증가, 높이 감소
    ggplot2::ggsave(filename = paste(G1, "divided by", G2, "Volcano plot.pdf", collapse = ""),
                    path = file.path(getwd(), "Volcano"), width = fig_width, height = fig_height,dpi=1000)
  }
}

Boxplot_ym_pc<-function (data, asterisk = "t_test", significant_variable_only = F, 
                         color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", 
                                   "#660099"), legend_position = "none", order = NULL, tip_length = 0.01, 
                         label_size = 2.88, step_increase = 0.05, width = 0.3, size = 0.5, 
                         fig_width = NA, fig_height = NA, X_text = 10, Y_text = 0, 
                         Y_lab = 10, T_size = 15, sig_int = c(0.05, 0.01),x_text_bold = "plain", x_text_angle = 45,y_text_bold="plain") 
{
  {
    ifelse(!dir.exists(file.path(getwd(), "boxplot_pc")), dir.create(file.path(getwd(), 
                                                                            "boxplot_pc")), FALSE)
    data[["Data_renamed"]] <- data[["Data_renamed"]] %>% 
      plyr::mutate(ZZZZ = data[["Data_renamed"]][, 2])
    data[["Data_renamed"]] <- data[["Data_renamed"]][, c(-1, 
                                                         -2)]
    data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 
                                  1)] <- sapply(data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 
                                                                              1)], function(x) {
                                                                                as.numeric(x)
                                                                              })
    NAMES <- colnames(data[["Data"]][3:ncol(data[["Data"]])])
    data[["Data_renamed"]]$ZZZZ <- as.factor(data[["Data_renamed"]]$ZZZZ)
    colnames(data[["Data_renamed"]])[ncol(data[["Data_renamed"]])] <- "Group"
    Comb <- gtools::combinations(length(unique(data[["Data"]]$Group)), 
                                 2, unique(data[["Data"]]$Group))
    ckey <- color[1:length(unique(data[["Data"]]$Group))]
    if (asterisk == "Scheffe") {
      p_val_data <- data[["Anova_PostHoc"]]
    }
    else if (asterisk == "t_test") {
      p_val_data <- data[["t_test"]]
    }
    else if (asterisk == "u_test") {
      p_val_data <- data[["u_test"]]
    }
    else if (asterisk == "Dunn") {
      p_val_data <- data[["Dunn"]]
    }
    else {
      print("Wrong asterisk input must be one of (Dunn,Scheffe,u_test,t_test)")
    }
    if (length(unique(data[["Data"]]$Group)) != 2) {
      for (a in 1:nrow(p_val_data)) {
        for (b in 1:ncol(p_val_data)) { 
          if (is.nan(p_val_data[a, b]) == T) {
            p_val_data[a, b] <- 1
          }
        }
      }
    }
    else if (length(unique(data[["Data"]]$Group)) == 2) {
      for (b in 1:length(p_val_data)) {
        if (is.nan(p_val_data[b]) == T) {
          p_val_data[b] <- 1
        }
      }
    }
  }
  group_not_two <- length(unique(data[["Data_renamed"]][["Group"]])) != 
    2
  tryCatch(expr = {
    num_core <- parallel::detectCores()
    final_cores <- min(5, num_core)
    cl <- parallel::makeCluster(final_cores)
    doSNOW::registerDoSNOW(cl)
    if (group_not_two) {
      pb <- txtProgressBar(max = nrow(p_val_data), style = 3)
    }
    else {
      pb <- txtProgressBar(max = length(p_val_data), style = 3)
    }
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    suppressWarnings(if (group_not_two) {
      foreach::foreach(number = 1:nrow(p_val_data), .options.snow = opts, 
                       .packages = c("LMSstat", "dplyr", "plyr")) %dopar% 
        {
          stat.test <- as.data.frame(matrix(data = NA, 
                                            nrow = nrow(Comb), ncol = 2))
          stat.test[, 1:2] <- Comb
          if (asterisk == "Dunn") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], " - ", stat.test[, 2])
          }
          else if (asterisk == "Scheffe") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    2], "-", stat.test[, 1])
          }
          else {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], "-", stat.test[, 2])
          }
          colnames(p_val_data) <- as.data.frame(strsplit(colnames(p_val_data), 
                                                         "___"))[1, ]
          t1 <- as.data.frame(t(as.data.frame(p_val_data)[number, 
          ]))
          t2 <- stat.test
          t3 <- merge(t2, t1, by = 0)
          stat.test <- t3[, 2:4]
          colnames(stat.test) <- c("group1", "group2", 
                                   "p")
          stat.test <- stat.test %>% plyr::mutate(p.adj.signif = case_when(p > 
                                                                             sig_int[1] ~ "NS", p <= sig_int[1] & p > 
                                                                             sig_int[2] ~ "*", p <= sig_int[2] ~ "**"))
          stat.test <- stat.test[stat.test$p.adj.signif != 
                                   "NS", ]
          if (length(stat.test > 4)) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group", fill = "Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) +
              ggplot2::scale_y_continuous(label = scales::label_scientific()) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55)) +
              ggplot2::labs(title = NAMES[number], x = NULL, y = NULL) +
              ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, lineheight = 1.2), 
                             plot.subtitle = ggplot2::element_text(size = 15, hjust = 0.5),
                             plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, size = 15),
                             axis.title.y = ggplot2::element_text(size = Y_text, face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust = 1, face = x_text_bold, angle = x_text_angle,hjust = 1),
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position) + 
              ggpubr::stat_pvalue_manual(
                stat.test,
                y.position = seq(
                  from = 1.05 * max(data[["Data_renamed"]][, number]),
                  by = step_increase * max(data[["Data_renamed"]][, number]),
                  length.out = nrow(stat.test)
                ),
                label = "p.adj.signif",
                size = 3.5,
                vjust = 0.025,
                label.size = label_size,
                tip.length = tip_length
              ) +
              ggplot2::scale_x_discrete(labels = function(x) parse(text = x))
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot_pc"), width = fig_width, height = fig_height)
          }
          else if (significant_variable_only == F) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group", fill = "Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55)) + 
              ggplot2::scale_y_continuous(label = scales::label_scientific()) + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = NULL) + 
              ggplot2::theme_classic() + 
              ggplot2::theme(
                plot.title = ggplot2::element_text(size = T_size, face = "bold", color = "Black", hjust = 0.5, lineheight = 1.2), 
                plot.subtitle = ggplot2::element_text(size = 15, hjust = 0.5), 
                plot.caption = ggplot2::element_text(size = 15), 
                axis.title.x = ggplot2::element_text(vjust = 10, size = 15),
                axis.title.y = ggplot2::element_text(size = Y_text, face = y_text_bold), 
                axis.text.x = ggplot2::element_text(size = X_text, vjust = 1, face = x_text_bold, angle = x_text_angle,hjust = 1),
                axis.text.y = ggplot2::element_text(size = Y_lab), 
                legend.position = legend_position) + 
              ggplot2::scale_x_discrete(labels = function(x) parse(text = x))
            
            ggplot2::ggsave(filename = paste(NAMES[number], "boxplot.pdf", collapse = ""), 
                            path = paste0(getwd(), "/boxplot_pc"), 
                            width = fig_width, height = fig_height)
          }
        }
    }
    else {
      p_val_data <- as.data.frame(p_val_data)
      if (asterisk == "t_test") {
        colnames(p_val_data) <- colnames(data[["Result"]])[1]
      }
      else if (asterisk == "u_test") {
        colnames(p_val_data) <- colnames(data[["Result"]])[2]
      }
      foreach::foreach(number = 1:nrow(p_val_data), .options.snow = opts, 
                       .packages = c("LMSstat", "dplyr", "plyr")) %dopar% 
        {
          df <- data_summary(data[["Data_renamed"]], 
                             varname = colnames(data[["Data_renamed"]])[number], 
                             groupnames = c("Group"))
          df <- df %>% plyr::mutate(ymax = len + SEM)
          colnames(df)[1] <- "Group"
          stat.test <- as.data.frame(matrix(data = NA, 
                                            nrow = nrow(Comb), ncol = 2))
          stat.test[, 1:2] <- Comb
          if (asterisk == "Dunn") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], " - ", stat.test[, 2])
          }
          else if (asterisk == "Scheffe") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    2], "-", stat.test[, 1])
          }
          else {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], "-", stat.test[, 2])
          }
          colnames(p_val_data) <- as.data.frame(strsplit(colnames(p_val_data), 
                                                         "___"))[1, ]
          t1 <- as.data.frame(t(as.data.frame(p_val_data[number, 
          ])))
          rownames(t1) <- colnames(p_val_data)
          t2 <- stat.test
          t3 <- merge(t2, t1, by = 0)
          stat.test <- t3[, 2:4]
          colnames(stat.test) <- c("group1", "group2", 
                                   "p")
          stat.test <- stat.test %>% plyr::mutate(p.adj.signif = case_when(p > 
                                                                             sig_int[1] ~ "NS", p <= sig_int[1] & p > 
                                                                             sig_int[2] ~ "*", p <= sig_int[2] ~ "**"))
          stat.test <- stat.test[stat.test$p.adj.signif != 
                                   "NS", ]
          if (length(stat.test > 4)) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group",fill="Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) + 
              ggplot2::scale_y_continuous(label =  scales::label_scientific()) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55))    + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = NULL) + ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, 
                                                                lineheight = 1.2), plot.subtitle = ggplot2::element_text(size = 15, 
                                                                                                                         hjust = 0.5), plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, 
                                                                  size = 15), 
                             axis.title.y = ggplot2::element_text(size = Y_text,face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust =1,face = x_text_bold, angle = x_text_angle),
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position) +
              ggpubr::stat_pvalue_manual(
                stat.test,
                y.position = seq(
                  from = 1.05 * max(data[["Data_renamed"]][, number]),
                  by = step_increase * max(data[["Data_renamed"]][, number]),
                  length.out = nrow(stat.test)
                ),
                label = "p.adj.signif",
                size = 3.5,
                vjust = 0.025,
                label.size = label_size,
                tip.length = tip_length
              )
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot"), width = fig_width, height = fig_height)
          }
          else if (significant_variable_only == F) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group",fill="Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55))    + 
              ggplot2::scale_y_continuous(label =  scales::label_scientific()) + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = NULL) + ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, 
                                                                lineheight = 1.2), plot.subtitle = ggplot2::element_text(size = 15, 
                                                                                                                         hjust = 0.5), plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, 
                                                                  size = 15), 
                             axis.title.y = ggplot2::element_text(size = Y_text,face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust =1,face = x_text_bold, angle = x_text_angle),
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position)+ 
              ggpubr::stat_pvalue_manual(stat.test, y.position = 1.05 * 
                                           max(data[["Data_renamed"]][, number]), 
                                         label = "p.adj.signif", size = 3.5, vjust = 0.05, 
                                         label.size = label_size, tip.length = tip_length)
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot"), width = fig_width, height = fig_height)
          }
        }
    })
    close(pb)
  }, finally = {
    parallel::stopCluster(cl)
  })
}




Boxplot_ym<-function (data, asterisk = "t_test", significant_variable_only = F, 
                      color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", 
                                "#660099"), legend_position = "none", order = NULL, tip_length = 0.01, 
                      label_size = 2.88, step_increase = 0.05, width = 0.3, size = 0.5, 
                      fig_width = NA, fig_height = NA, X_text = 10, Y_text = 12, 
                      Y_lab = 10, T_size = 15, sig_int = c(0.05, 0.01),x_text_bold = "plain", x_text_angle = 45,y_text_bold="plain") 
{
  {
    ifelse(!dir.exists(file.path(getwd(), "boxplot")), dir.create(file.path(getwd(), 
                                                                            "boxplot")), FALSE)
    data[["Data_renamed"]] <- data[["Data_renamed"]] %>% 
      plyr::mutate(ZZZZ = data[["Data_renamed"]][, 2])
    data[["Data_renamed"]] <- data[["Data_renamed"]][, c(-1, 
                                                         -2)]
    data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 
                                  1)] <- sapply(data[["Data_renamed"]][, 1:(ncol(data[["Data_renamed"]]) - 
                                                                              1)], function(x) {
                                                                                as.numeric(x)
                                                                              })
    NAMES <- colnames(data[["Data"]][3:ncol(data[["Data"]])])
    data[["Data_renamed"]]$ZZZZ <- as.factor(data[["Data_renamed"]]$ZZZZ)
    colnames(data[["Data_renamed"]])[ncol(data[["Data_renamed"]])] <- "Group"
    Comb <- gtools::combinations(length(unique(data[["Data"]]$Group)), 
                                 2, unique(data[["Data"]]$Group))
    ckey <- color[1:length(unique(data[["Data"]]$Group))]
    if (asterisk == "Scheffe") {
      p_val_data <- data[["Anova_PostHoc"]]
    }
    else if (asterisk == "t_test") {
      p_val_data <- data[["t_test"]]
    }
    else if (asterisk == "u_test") {
      p_val_data <- data[["u_test"]]
    }
    else if (asterisk == "Dunn") {
      p_val_data <- data[["Dunn"]]
    }
    else {
      print("Wrong asterisk input must be one of (Dunn,Scheffe,u_test,t_test)")
    }
    if (length(unique(data[["Data"]]$Group)) != 2) {
      for (a in 1:nrow(p_val_data)) {
        for (b in 1:ncol(p_val_data)) { 
          if (is.nan(p_val_data[a, b]) == T) {
            p_val_data[a, b] <- 1
          }
        }
      }
    }
    else if (length(unique(data[["Data"]]$Group)) == 2) {
      for (b in 1:length(p_val_data)) {
        if (is.nan(p_val_data[b]) == T) {
          p_val_data[b] <- 1
        }
      }
    }
  }
  group_not_two <- length(unique(data[["Data_renamed"]][["Group"]])) != 
    2
  tryCatch(expr = {
    num_core <- parallel::detectCores()
    final_cores <- min(5, num_core)
    cl <- parallel::makeCluster(final_cores)
    doSNOW::registerDoSNOW(cl)
    if (group_not_two) {
      pb <- txtProgressBar(max = nrow(p_val_data), style = 3)
    }
    else {
      pb <- txtProgressBar(max = length(p_val_data), style = 3)
    }
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    suppressWarnings(if (group_not_two) {
      foreach::foreach(number = 1:nrow(p_val_data), .options.snow = opts, 
                       .packages = c("LMSstat", "dplyr", "plyr")) %dopar% 
        {
          stat.test <- as.data.frame(matrix(data = NA, 
                                            nrow = nrow(Comb), ncol = 2))
          stat.test[, 1:2] <- Comb
          if (asterisk == "Dunn") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], " - ", stat.test[, 2])
          }
          else if (asterisk == "Scheffe") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    2], "-", stat.test[, 1])
          }
          else {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], "-", stat.test[, 2])
          }
          colnames(p_val_data) <- as.data.frame(strsplit(colnames(p_val_data), 
                                                         "___"))[1, ]
          t1 <- as.data.frame(t(as.data.frame(p_val_data)[number, 
          ]))
          t2 <- stat.test
          t3 <- merge(t2, t1, by = 0)
          stat.test <- t3[, 2:4]
          colnames(stat.test) <- c("group1", "group2", 
                                   "p")
          stat.test <- stat.test %>% plyr::mutate(p.adj.signif = case_when(p > 
                                                                             sig_int[1] ~ "NS", p <= sig_int[1] & p > 
                                                                             sig_int[2] ~ "*", p <= sig_int[2] ~ "**"))
          stat.test <- stat.test[stat.test$p.adj.signif != 
                                   "NS", ]
          if (length(stat.test > 4)) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group",fill="Group", palette = ckey, add = "jitter", 
                              size = size, order = order, width = width) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55))    + 
              ggplot2::scale_y_continuous(label =  scales::label_scientific()) + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = "Intensity") + ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, 
                                                                lineheight = 1.2), plot.subtitle = ggplot2::element_text(size = 15, 
                                                                                                                         hjust = 0.5), plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, 
                                                                  size = 15), 
                             axis.title.y = ggplot2::element_text(size = Y_text,face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust=1,face = x_text_bold, angle = x_text_angle,hjust = 1), 
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position) + 
              ggpubr::stat_pvalue_manual(stat.test, y.position = 1.05 * 
                                           max(data[["Data_renamed"]][, number]), 
                                         step.increase = step_increase, label.size = label_size, 
                                         tip.length = tip_length, label = "p.adj.signif", 
                                         size = 3.5, vjust = 0.05)+ 
              ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) 
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot"), width = fig_width, height = fig_height)
          }
          else if (significant_variable_only == F) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group",fill="Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55))    + 
              ggplot2::scale_y_continuous(label =  scales::label_scientific()) + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = "Intensity") + ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, 
                                                                lineheight = 1.2), plot.subtitle = ggplot2::element_text(size = 15, 
                                                                                                                         hjust = 0.5), plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, size = 15),
                             axis.title.y = ggplot2::element_text(size = Y_text,face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust=1,face = x_text_bold, angle = x_text_angle,hjust = 1),
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position)+
              ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) 
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot"), width = fig_width, height = fig_height)
          }
        }
    }
    else {
      p_val_data <- as.data.frame(p_val_data)
      if (asterisk == "t_test") {
        colnames(p_val_data) <- colnames(data[["Result"]])[1]
      }
      else if (asterisk == "u_test") {
        colnames(p_val_data) <- colnames(data[["Result"]])[2]
      }
      foreach::foreach(number = 1:nrow(p_val_data), .options.snow = opts, 
                       .packages = c("LMSstat", "dplyr", "plyr")) %dopar% 
        {
          df <- data_summary(data[["Data_renamed"]], 
                             varname = colnames(data[["Data_renamed"]])[number], 
                             groupnames = c("Group"))
          df <- df %>% plyr::mutate(ymax = len + SEM)
          colnames(df)[1] <- "Group"
          stat.test <- as.data.frame(matrix(data = NA, 
                                            nrow = nrow(Comb), ncol = 2))
          stat.test[, 1:2] <- Comb
          if (asterisk == "Dunn") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], " - ", stat.test[, 2])
          }
          else if (asterisk == "Scheffe") {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    2], "-", stat.test[, 1])
          }
          else {
            rownames(stat.test) <- paste0(stat.test[, 
                                                    1], "-", stat.test[, 2])
          }
          colnames(p_val_data) <- as.data.frame(strsplit(colnames(p_val_data), 
                                                         "___"))[1, ]
          t1 <- as.data.frame(t(as.data.frame(p_val_data[number, 
          ])))
          rownames(t1) <- colnames(p_val_data)
          t2 <- stat.test
          t3 <- merge(t2, t1, by = 0)
          stat.test <- t3[, 2:4]
          colnames(stat.test) <- c("group1", "group2", 
                                   "p")
          stat.test <- stat.test %>% plyr::mutate(p.adj.signif = case_when(p > 
                                                                             sig_int[1] ~ "NS", p <= sig_int[1] & p > 
                                                                             sig_int[2] ~ "*", p <= sig_int[2] ~ "**"))
          stat.test <- stat.test[stat.test$p.adj.signif != 
                                   "NS", ]
          if (length(stat.test > 4)) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group",fill="Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) + 
              ggplot2::scale_y_continuous(label =  scales::label_scientific()) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55))    + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = "Intensity") + ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, 
                                                                lineheight = 1.2), plot.subtitle = ggplot2::element_text(size = 15, 
                                                                                                                         hjust = 0.5), plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, 
                                                                  size = 15), 
                             axis.title.y = ggplot2::element_text(size = Y_text,face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust=1,face = x_text_bold, angle = x_text_angle,hjust = 1),
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position) + 
              ggpubr::stat_pvalue_manual(stat.test, y.position = 1.05 * 
                                           max(data[["Data_renamed"]][, number]), 
                                         label = "p.adj.signif", size = 3.5, vjust = 0.05, 
                                         label.size = label_size, tip.length = tip_length)
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot"),  width = fig_width, height = fig_height)
          }
          else if (significant_variable_only == F) {
            ggpubr::ggboxplot(data[["Data_renamed"]], 
                              x = "Group", y = colnames(data[["Data_renamed"]])[number], 
                              color = "Group",fill="Group", palette = ckey, size = size, 
                              add = "jitter", order = order, width = width) + 
              ggplot2::scale_fill_manual(values = scales::alpha(ckey, 0.55))    + 
              ggplot2::scale_y_continuous(label =  scales::label_scientific()) + 
              ggplot2::labs(title = NAMES[number], x = NULL, 
                            y = "Intensity") + ggplot2::theme_classic() + 
              ggplot2::theme(plot.title = ggplot2::element_text(size = T_size, 
                                                                face = "bold", color = "Black", hjust = 0.5, 
                                                                lineheight = 1.2), plot.subtitle = ggplot2::element_text(size = 15, 
                                                                                                                         hjust = 0.5), plot.caption = ggplot2::element_text(size = 15), 
                             axis.title.x = ggplot2::element_text(vjust = 10, 
                                                                  size = 15), 
                             axis.title.y = ggplot2::element_text(size = Y_text,face = y_text_bold), 
                             axis.text.x = ggplot2::element_text(size = X_text, vjust=1,face = x_text_bold, angle = x_text_angle,hjust = 1),
                             axis.text.y = ggplot2::element_text(size = Y_lab), 
                             legend.position = legend_position)+ 
              ggplot2::scale_x_discrete(labels = function(x) parse(text = x)) 
            ggplot2::ggsave(filename = paste(NAMES[number], 
                                             "boxplot.pdf", collapse = ""), path = paste0(getwd(), 
                                                                                          "/boxplot"), width = fig_width, height = fig_height)
          }
        }
    })
    close(pb)
  }, finally = {
    parallel::stopCluster(cl)
  })
}


