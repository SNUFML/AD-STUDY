###############################################################################
# Renamed Functions (same behavior/IO as originals)
###############################################################################

# ---------- small helpers ----------
.melt_mat <- function(mat) {
  mat <- as.matrix(mat)
  data.frame(
    Var1  = rownames(mat)[row(mat)],
    Var2  = colnames(mat)[col(mat)],
    value = c(mat),
    row.names = NULL,
    check.names = FALSE
  )
}
.safe_name <- function(x) paste0(deparse(substitute(x)), collapse = "")

# ---------- Normality test (NT -> run_normality_tests) ----------
run_normality_tests <- function(data) {
  stopifnot(is.data.frame(data) || is.matrix(data))
  dn <- .safe_name(data)
  
  Wp <- sapply(seq_len(ncol(data)), function(i) {
    st <- shapiro.test(as.numeric(data[, i]))
    c(W = unname(st$statistic), p = unname(st$p.value))
  })
  outputs <- t(Wp)
  colnames(outputs) <- c("shapiro test W", "shapiro test p value")
  rownames(outputs) <- colnames(data)
  
  write.csv(outputs, paste0("normality.test_", dn, ".csv"))
  
  pdf(paste0("normality.test_", dn, ".pdf"))
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
  par(mfrow = c(3, 2))
  for (j in seq_len(ncol(data))) {
    x <- as.numeric(data[, j])
    hist(x, xlab = colnames(data)[j],
         main = paste("Histogram of", colnames(data)[j]))
    qqnorm(x, main = paste("Normal Q-Q plot of", colnames(data)[j]))
    qqline(x, col = "red")
  }
  dev.off(invisible = TRUE)
  invisible(outputs)
}

# ---------- Spearman cross-correlation (SC -> spearman_crosscorr) ----------
spearman_crosscorr <- function(data1, data2) {
  data1 <- as.data.frame(data1, check.names = FALSE)
  data2 <- as.data.frame(data2, check.names = FALSE)
  
  coef_mat <- suppressWarnings(cor(data1, data2, method = "spearman",
                                   use = "pairwise.complete.obs"))
  
  i_idx <- rep(seq_len(ncol(data1)), times = ncol(data2))
  j_idx <- rep(seq_len(ncol(data2)), each  = ncol(data1))
  pv <- mapply(function(i, j) {
    suppressWarnings(cor.test(data1[[i]], data2[[j]],
                              method = "spearman", exact = FALSE)$p.value)
  }, i_idx, j_idx)
  
  p_mat <- matrix(pv, nrow = ncol(data1), ncol = ncol(data2),
                  dimnames = list(colnames(data1), colnames(data2)))
  rownames(coef_mat) <- colnames(data1); colnames(coef_mat) <- colnames(data2)
  list(coef_mat, p_mat)
}

# ---------- Conditional logistic OR (OR -> clogit_or_scan) ----------
clogit_or_scan <- function(form, meta.d, clinic.d, outcome.d) {
  stopifnot(nrow(meta.d) == nrow(clinic.d), nrow(meta.d) == nrow(outcome.d))
  meta.d   <- as.data.frame(meta.d, check.names = FALSE)
  clinic.d <- as.data.frame(clinic.d, check.names = FALSE)
  outcome.d<- as.data.frame(outcome.d, check.names = FALSE)
  
  rows <- lapply(seq_len(ncol(meta.d)), function(i) {
    tmp <- data.frame(
      meta = meta.d[, i] / sd(meta.d[, i], na.rm = TRUE),
      clinic.d, outcome.d, check.names = FALSE
    )
    fit <- try(survival::clogit(form, data = tmp), silent = TRUE)
    if (inherits(fit, "try-error")) return(c(0, 0, 0, 0))
    
    ci <- try(exp(confint(fit))[1, ], silent = TRUE)
    if (inherits(ci, "try-error")) ci <- c(NA_real_, NA_real_)
    or <- try(exp(coef(fit))[1], silent = TRUE)
    if (inherits(or, "try-error")) or <- NA_real_
    p  <- try(as.data.frame(summary(fit)$coefficient)$`Pr(>|z|)`[1],
              silent = TRUE)
    if (inherits(p, "try-error")) p <- NA_real_
    c(ci[1], ci[2], or, p)
  })
  
  outputs1 <- do.call(rbind, rows)
  colnames(outputs1) <- c("2.5% CI", "97.5% CI", "OR", "p")
  rownames(outputs1) <- colnames(meta.d)
  
  fdr <- p.adjust(outputs1[, "p"], "fdr")
  outputs <- data.frame(outputs1, fdr)
  colnames(outputs) <- c("2.5% CI", "97.5% CI", "OR", "p", "fdr")
  outputs
}

# ---------- Dataset split (split_dataset -> split_by_match) ----------
split_by_match <- function(n.seed, proportion.subdataset1) {
  ids <- names(table(outcome.d$match))
  set.seed(n.seed)
  n_train <- round(length(ids) * proportion.subdataset1)
  train_ID <- sample(ids, n_train)
  
  idx_train <- which(outcome.d$match %in% train_ID)
  idx_valid <- setdiff(seq_len(nrow(outcome.d)), idx_train)
  
  t_meta.d    <- meta.d[idx_train, , drop = FALSE]
  t_clinic.d  <- clinic.d[idx_train, , drop = FALSE]
  t_outcome.d <- outcome.d[idx_train, , drop = FALSE]
  
  v_meta.d    <- meta.d[idx_valid, , drop = FALSE]
  v_clinic.d  <- clinic.d[idx_valid, , drop = FALSE]
  v_outcome.d <- outcome.d[idx_valid, , drop = FALSE]
  
  list(t_meta.d, t_clinic.d, t_outcome.d, v_meta.d, v_clinic.d, v_outcome.d)
}

# ---------- Random subset (random_dataset -> sample_by_match) ----------
sample_by_match <- function(n.seed, percentage) {
  ids <- names(table(outcome.d$match))
  set.seed(n.seed)
  sel_ID <- sample(ids, round(length(ids) * percentage))
  idx <- which(outcome.d$match %in% sel_ID)
  
  s_meta.d    <- meta.d[idx, , drop = FALSE]
  s_clinic.d  <- clinic.d[idx, , drop = FALSE]
  s_outcome.d <- outcome.d[idx, , drop = FALSE]
  list(s_meta.d, s_clinic.d, s_outcome.d)
}

# ---------- Random OR counts (random_OR -> random_or_counts) ----------
random_or_counts <- function(times.random, percentage) {
  OR_list <- replicate(times.random, {
    rd <- sample_by_match(n.seed = sample.int(1e6, 1), percentage = percentage)
    clogit_or_scan(form1, rd[[1]], rd[[2]], rd[[3]])
  }, simplify = FALSE)
  
  OR_mat  <- do.call(cbind, lapply(OR_list, `[[`, "OR"))
  p_mat   <- do.call(cbind, lapply(OR_list, `[[`, "p"))
  fdr_mat <- do.call(cbind, lapply(OR_list, `[[`, "fdr"))
  
  gt1_p   <- rowSums(OR_mat > 1 & p_mat < 0.05,  na.rm = TRUE)
  lt1_p   <- rowSums(OR_mat < 1 & p_mat < 0.05,  na.rm = TRUE)
  gt1_fdr <- rowSums(OR_mat > 1 & fdr_mat < 0.05, na.rm = TRUE)
  lt1_fdr <- rowSums(OR_mat < 1 & fdr_mat < 0.05, na.rm = TRUE)
  
  out <- data.frame(
    `OR>1,p<0.05`   = gt1_p,
    `OR>1,fdr<0.05` = gt1_fdr,
    `OR<1,p<0.05`   = lt1_p,
    `OR<1,fdr<0.05` = lt1_fdr,
    check.names = FALSE
  )
  rownames(out) <- colnames(meta.d)
  out
}

# ---------- Mediation grid (CM -> causal_mediation_grid) ----------
causal_mediation_grid <- function(form.a, form.b, form.c) {
  out_rows <- list(); k <- 1
  for (i in seq_len(ncol(meta.d))) {
    for (j in seq_len(ncol(clinic.d))) {
      tmp <- data.frame(
        case   = outcome.d$case,
        meta   = (meta.d[, i] - mean(meta.d[, i], na.rm = TRUE)) /
          sd(meta.d[, i], na.rm = TRUE),
        clinic = clinic.d[, j],
        check.names = FALSE
      )
      cm.a <- lm(form.a, data = tmp)
      cm.b <- glm(form.b, data = tmp, family = binomial())
      cm.c <- glm(form.c, data = tmp, family = binomial())
      
      set.seed(1)
      md <- mediation::mediate(cm.a, cm.b, treat = "meta",
                               mediator = "clinic", boot = TRUE)
      
      a <- c(unname(coef(summary(cm.a))[2, 1]),
             unname(coef(summary(cm.a))[2, 4]))
      b <- c(unname(coef(summary(cm.b))[3, 1]),
             unname(coef(summary(cm.b))[3, 4]))
      c_ <- c(unname(coef(summary(cm.c))[2, 1]),
              unname(coef(summary(cm.c))[2, 4]))
      
      out_rows[[k]] <- c(
        metabolites = colnames(meta.d)[i],
        `clinical parameters` = colnames(clinic.d)[j],
        `beta a` = a[1], `p a` = a[2],
        `beta b` = b[1], `p b` = b[2],
        `beta c` = c_[1], `p c` = c_[2],
        acme = md$d0, `2.5%CI` = md$d0.ci[1], `97.5%CI` = md$d0.ci[2], `p acme` = md$d0.p,
        ade  = md$z0, `2.5%CI.1` = md$z0.ci[1], `97.5%CI.1` = md$z0.ci[2], `p ade`  = md$z0.p,
        te   = md$tau.coef, `2.5%CI.2` = md$tau.ci[1], `97.5%CI.2` = md$tau.ci[2], `p te` = md$tau.p,
        prop = md$n0, `2.5%CI.3` = md$n0.ci[1], `97.5%CI.3` = md$n0.ci[2], `p prop` = md$n0.p
      )
      k <- k + 1
    }
  }
  outputs <- as.data.frame(do.call(rbind, out_rows), check.names = FALSE)
  colnames(outputs) <- c("metabolites","clinical parameters",
                         "beta a","p a","beta b","p b","beta c","p c",
                         "acme","2.5%CI","97.5%CI","p acme",
                         "ade","2.5%CI","97.5%CI","p ade",
                         "te","2.5%CI","97.5%CI","p te",
                         "prop","2.5%CI","97.5%CI","p prop")
  outputs
}

# ---------- Variance partition (VP -> variance_partition_scan) ----------
# ---- Normality diagnostics ----
run_normality_tests <- function(data) {
  stopifnot(is.data.frame(data) || is.matrix(data))
  dn <- deparse(substitute(data))
  
  Wp <- sapply(seq_len(ncol(data)), function(i) {
    st <- shapiro.test(as.numeric(data[, i]))
    c(W = unname(st$statistic), p = unname(st$p.value))
  })
  outputs <- t(Wp)
  colnames(outputs) <- c("shapiro test W", "shapiro test p value")
  rownames(outputs) <- colnames(data)
  
  write.csv(outputs, paste0("normality.test_", dn, ".csv"))
  
  pdf(paste0("normality.test_", dn, ".pdf"))
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
  par(mfrow = c(3, 2))
  for (j in seq_len(ncol(data))) {
    x <- as.numeric(data[, j])
    hist(x, xlab = colnames(data)[j],
         main = paste("Histogram of", colnames(data)[j]))
    qqnorm(x, main = paste("Normal Q-Q plot of", colnames(data)[j]))
    qqline(x, col = "red")
  }
  dev.off()   
  invisible(outputs)
}