# ============================================================================
# 11_clinical_group_x_age_interaction_MI.R
# 临床分组 x 年龄 交互分析 —— 多重插补版 (mice, m=10)
# 与 10 相同的问题: depression/anxiety 缺失 -> 插补
# 每个插补数据集内按 cutoff 重新二分, 拟合交互模型, Rubin's rule 合并
#
# 分组: depression (PHQ-9) >= 10 -> clinical;  anxiety (GAD-7) >= 10 -> clinical
# 模型: ROI ~ age * group + gender + mri_equip + EMO_mean_FD
# 仅线性, 不做可视化
#
# 输出: _tables/depression_group_x_age_MI.xlsx, anxiety_group_x_age_MI.xlsx
# ============================================================================
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(broom)
library(writexl)
library(mice)
library(emmeans)

# ---- 0. 读取数据 ----
df <- read.csv("_data/simulated_mahalanobis_data.csv", stringsAsFactors = FALSE)
df$mri_equip <- factor(df$mri_equip)
df$gender    <- factor(df$gender, levels = c(1, 2), labels = c("male", "female"))
df$age_z     <- as.numeric(scale(df$age))

roi_cols <- grep("^ROI[0-9]+$", colnames(df), value = TRUE)
stopifnot(length(roi_cols) == 432)
df$mean_md <- rowMeans(df[, roi_cols], na.rm = TRUE)

PHQ9_CUT <- 10
GAD7_CUT <- 10

cat("数据维度:", nrow(df), "x", ncol(df), "\n")
cat("ROI数量:", length(roi_cols), "\n")
cat("anxiety 缺失:", sum(is.na(df$anxiety)), "\n")
cat("depression 缺失:", sum(is.na(df$depression)), "\n")

# ---- 1. 多重插补 (m=10) ----
cat("\n========== 多重插补 (mice, m=10) ==========\n")
m_mi <- 10
set.seed(2024)

imp_df <- df %>%
  select(anxiety, depression, age, gender, mri_equip, EMO_mean_FD, mean_md)

imp <- mice(imp_df, m = m_mi, method = "pmm", maxit = 5, seed = 123, print = FALSE)

comp_list <- lapply(1:m_mi, function(k) {
  d <- complete(imp, k)
  d$mri_equip <- factor(d$mri_equip)
  d$gender    <- factor(d$gender, levels = c("male", "female"))
  # 补回ROI列 (插补模型未包含)
  d[, roi_cols] <- df[, roi_cols]
  # 补回标准化年龄 (age无缺失, 与总体一致)
  d$age_z <- as.numeric(scale(d$age))
  d
})
cat("插补完成:", m_mi, "个数据集 (分析用 N =", nrow(comp_list[[1]]), ")\n")

# ---- 2. Rubin's rule 合并 (mice::pool.scalar, Barnard-Rubin df) ----
# n = 分析样本量, k = 模型参数个数 (用于完整数据自由度 n-k)
pool_est <- function(betas, ses, n, k) {
  ok <- is.finite(betas) & is.finite(ses)
  if (sum(ok) < 2) {
    return(c(beta = NA, se = NA, t_value = NA, p_value = NA,
             fmi = NA, df = NA))
  }
  betas <- betas[ok]; ses <- ses[ok]
  ps   <- mice::pool.scalar(Q = betas, U = ses^2, n = n, k = k)
  qbar <- ps$qbar; tt <- ps$t; dfv <- ps$df
  tval <- qbar / sqrt(tt)
  pval <- 2 * pt(-abs(tval), df = dfv)
  c(beta = qbar, se = sqrt(tt), t_value = tval, p_value = pval,
    fmi = ps$fmi, df = dfv)
}

# ---- 3. emtrends 稳健提取 ----
emtrends_df <- function(fit, by_var, var = "age") {
  f <- as.formula(sprintf("~ %s", by_var))
  d <- as.data.frame(summary(emtrends(fit, f, var = var), infer = c(TRUE, TRUE)))
  trend_col <- grep("\\.trend$", names(d), value = TRUE)[1]
  if (!is.na(trend_col) && trend_col != "age.trend") {
    names(d)[names(d) == trend_col] <- "age.trend"
  }
  stat_cands <- intersect(c("t.ratio", "z.ratio", "t.statistic", "statistic"), names(d))
  if (length(stat_cands) > 0 && stat_cands[1] != "t.ratio") {
    names(d)[names(d) == stat_cands[1]] <- "t.ratio"
  }
  if (by_var %in% names(d) && by_var != "group") {
    names(d)[names(d) == by_var] <- "group"
  }
  need <- c("age.trend", "SE", "df", "lower.CL", "upper.CL", "t.ratio", "p.value")
  for (cn in need) if (!cn %in% names(d)) d[[cn]] <- NA
  d
}

# 合并多个插补的简单斜率 (按 group 分别 pool)
pool_slopes <- function(sl_all, n, k) {
  groups <- unique(sl_all$group)
  bind_rows(lapply(groups, function(g) {
    sub <- sl_all[sl_all$group == g, ]
    pv  <- pool_est(sub$age.trend, sub$SE, n, k)
    data.frame(group = g, t(pv), stringsAsFactors = FALSE)
  }))
}

# ---- 4. 单组分组的 MI 分析 ----
analyze_group_mi <- function(score_var, cutoff, measure_label, out_xlsx) {

  cat(sprintf("\n\n############ %s (MI) ############\n", measure_label))
  dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

  # 每个插补数据集派生组别
  comps <- lapply(comp_list, function(d) {
    d$grp <- factor(ifelse(d[[score_var]] >= cutoff, "clinical", "healthy"),
                    levels = c("healthy", "clinical"))
    d
  })

  # --- 观察样本描述 (去缺失后) ---
  obs_group <- ifelse(is.na(df[[score_var]]), NA,
                      ifelse(df[[score_var]] >= cutoff, "clinical", "healthy"))
  d_obs <- df[!is.na(obs_group), ]
  d_obs$.grp <- obs_group[!is.na(obs_group)]
  make_row <- function(sub, label) {
    data.frame(
      measure    = measure_label,
      cutoff     = cutoff,
      group      = label,
      n          = nrow(sub),
      n_male     = sum(sub$gender == "male"),
      n_female   = sum(sub$gender == "female"),
      pct_female = round(100 * mean(sub$gender == "female"), 1),
      age_mean   = round(mean(sub$age), 2),
      age_sd     = round(sd(sub$age), 2),
      age_min    = round(min(sub$age), 1),
      age_max    = round(max(sub$age), 1),
      score_mean = round(mean(sub[[score_var]]), 2),
      score_sd   = round(sd(sub[[score_var]]), 2),
      stringsAsFactors = FALSE
    )
  }
  sample_summary <- rbind(
    make_row(d_obs, "Total (observed, excl. missing)"),
    make_row(d_obs[d_obs$.grp == "healthy", ], "Healthy (observed)"),
    make_row(d_obs[d_obs$.grp == "clinical", ], "Clinical (observed)")
  )
  cat(sprintf("\n%s 观察缺失: %d | 去缺失 n = %d | MI分析 N = %d (m=%d)\n",
              measure_label, sum(is.na(obs_group)), nrow(d_obs),
              nrow(comps[[1]]), m_mi))
  print(sample_summary)

  # --- 全局交互 (pooled) ---
  fml <- as.formula("mean_md ~ age * grp + gender + mri_equip + EMO_mean_FD")
  fml_z <- as.formula("mean_md ~ age_z * grp + gender + mri_equip + EMO_mean_FD")
  n_global <- nrow(comps[[1]])
  k_global <- length(coef(lm(fml, data = comps[[1]])))

  gb <- gs <- numeric(m_mi); gbz <- gsz <- numeric(m_mi)
  sl_global <- vector("list", m_mi)
  coef_list <- vector("list", m_mi)
  for (k in 1:m_mi) {
    d <- comps[[k]]
    fit  <- lm(fml, data = d)
    fitz <- lm(fml_z, data = d)
    sm  <- summary(fit)$coefficients
    smz <- summary(fitz)$coefficients
    coef_list[[k]] <- sm
    i1 <- grep("^age:grp", rownames(sm));  i2 <- grep("^age_z:grp", rownames(smz))
    gb[k]  <- sm[i1, "Estimate"];  gs[k]  <- sm[i1, "Std. Error"]
    gbz[k] <- smz[i2, "Estimate"]; gsz[k] <- smz[i2, "Std. Error"]
    sl_global[[k]] <- emtrends_df(fit, by_var = "grp", var = "age")
  }
  global_int  <- pool_est(gb, gs, n_global, k_global)
  global_intz <- pool_est(gbz, gsz, n_global, k_global)
  global_slopes <- pool_slopes(bind_rows(sl_global), n_global, k_global)

  # 全局模型所有系数 pooled
  terms_all <- rownames(coef_list[[1]])
  global_model_tbl <- bind_rows(lapply(terms_all, function(tt) {
    b <- sapply(coef_list, function(m) m[tt, "Estimate"])
    s <- sapply(coef_list, function(m) m[tt, "Std. Error"])
    pv <- pool_est(b, s, n_global, k_global)
    data.frame(term = tt, t(pv), stringsAsFactors = FALSE)
  })) %>%
    mutate(across(where(is.numeric), ~ round(., 6)))

  cat(sprintf("\n全局交互(MI): beta = %.4f, se = %.4f, p = %.2e, FMI = %.3f\n",
              global_int["beta"], global_int["se"], global_int["p_value"],
              global_int["fmi"]))
  cat(sprintf("全局交互(MI, std): beta_std = %.4f, p = %.2e\n",
              global_intz["beta"], global_intz["p_value"]))
  cat("全局简单斜率(MI):\n"); print(global_slopes)

  # --- ROI水平交互 (pooled) ---
  n_roi <- nrow(comps[[1]])
  k_roi <- length(coef(lm(as.formula(sprintf(
    "`%s` ~ age * grp + gender + mri_equip + EMO_mean_FD", roi_cols[1])),
    data = comps[[1]])))

  run_roi_pooled <- function(age_var) {
    pat <- sprintf("^%s:grp", age_var)
    beta_mat <- se_mat <- matrix(NA, nrow = length(roi_cols), ncol = m_mi)
    for (k in 1:m_mi) {
      d <- comps[[k]]
      for (i in seq_along(roi_cols)) {
        f <- as.formula(sprintf(
          "`%s` ~ %s * grp + gender + mri_equip + EMO_mean_FD",
          roi_cols[i], age_var))
        fit <- tryCatch(lm(f, data = d), error = function(e) NULL)
        if (is.null(fit)) next
        sm <- summary(fit)$coefficients
        idx <- grep(pat, rownames(sm))
        if (!length(idx)) next
        beta_mat[i, k] <- sm[idx[1], "Estimate"]
        se_mat[i, k]   <- sm[idx[1], "Std. Error"]
      }
    }
    pooled <- t(sapply(seq_along(roi_cols), function(i) {
      pool_est(beta_mat[i, ], se_mat[i, ], n_roi, k_roi)
    }))
    out <- data.frame(roi = roi_cols, index = seq_along(roi_cols), pooled,
                      stringsAsFactors = FALSE)
    for (cn in c("beta", "se", "t_value", "p_value", "fmi", "df"))
      out[[cn]] <- as.numeric(out[[cn]])
    out
  }

  cat("\nROI水平交互 (raw) ...\n")
  roi_raw <- run_roi_pooled("age")
  roi_raw$p_fdr <- p.adjust(roi_raw$p_value, method = "fdr")
  roi_raw$significant <- roi_raw$p_fdr < 0.05
  roi_raw$system <- ifelse(roi_raw$index <= 400, "Cortical", "Subcortical")
  cat(sprintf("  FDR<0.05: %d / 432\n", sum(roi_raw$significant, na.rm = TRUE)))

  cat("ROI水平交互 (std) ...\n")
  roi_std <- run_roi_pooled("age_z")
  roi_std$p_fdr <- p.adjust(roi_std$p_value, method = "fdr")
  roi_std$significant <- roi_std$p_fdr < 0.05
  cat(sprintf("  FDR<0.05: %d / 432\n", sum(roi_std$significant, na.rm = TRUE)))

  roi_merged <- roi_raw %>%
    select(roi, index, system, beta, se, t_value, p_value, p_fdr, significant, fmi) %>%
    rename(beta_raw = beta, se_raw = se, t_raw = t_value, p_raw = p_value,
           p_fdr_raw = p_fdr, significant_raw = significant, fmi_raw = fmi) %>%
    left_join(
      roi_std %>%
        select(roi, beta, se, t_value, p_value, p_fdr, significant, fmi) %>%
        rename(beta_std = beta, se_std = se, t_std = t_value, p_std = p_value,
               p_fdr_std = p_fdr, significant_std = significant, fmi_std = fmi),
      by = "roi") %>%
    arrange(p_fdr_raw)

  # --- 显著ROI的简单斜率 (pooled) ---
  sig_rois <- roi_merged %>% filter(significant_raw) %>% pull(roi)
  cat(sprintf("显著交互ROI数: %d\n", length(sig_rois)))

  roi_slopes <- NULL
  if (length(sig_rois) > 0) {
    roi_slopes <- bind_rows(lapply(sig_rois, function(rn) {
      sl <- bind_rows(lapply(comps, function(d) {
        fit <- lm(as.formula(sprintf(
          "`%s` ~ age * grp + gender + mri_equip + EMO_mean_FD", rn)), data = d)
        emtrends_df(fit, by_var = "grp", var = "age")
      }))
      pooled <- pool_slopes(sl, n_roi, k_roi)
      pooled$roi <- rn
      pooled
    })) %>%
      left_join(roi_merged %>% select(roi, p_fdr_raw), by = "roi") %>%
      select(roi, group, p_fdr_raw, beta, se, t_value, p_value, fmi)
    print(head(roi_slopes))
  }

  # --- 保存 ---
  gi_tbl <- data.frame(
    measure   = measure_label,
    term      = "age:group",
    beta_raw  = global_int["beta"],  se_raw  = global_int["se"],
    t_raw     = global_int["t_value"], p_raw = global_int["p_value"], fmi_raw = global_int["fmi"],
    beta_std  = global_intz["beta"], se_std  = global_intz["se"],
    t_std     = global_intz["t_value"], p_std = global_intz["p_value"], fmi_std = global_intz["fmi"],
    stringsAsFactors = FALSE
  ) %>% mutate(across(where(is.numeric), ~ round(., 6)))

  global_model_tbl <- global_model_tbl %>%
    mutate(measure = measure_label, .before = 1)

  sheet_list <- list(
    sample_summary       = sample_summary,
    global_model         = global_model_tbl,
    global_interaction   = gi_tbl,
    global_simple_slopes = global_slopes %>%
      mutate(across(where(is.numeric), ~ round(., 6))),
    roi_interaction      = roi_merged %>%
      select(roi, index, system, beta_raw, se_raw, t_raw, p_raw, p_fdr_raw,
             significant_raw, fmi_raw, beta_std, se_std, t_std, p_std,
             p_fdr_std, significant_std, fmi_std) %>%
      mutate(across(c(beta_raw, se_raw, t_raw, beta_std, se_std, t_std, fmi_raw, fmi_std),
                    ~ round(., 6)),
             across(c(p_raw, p_fdr_raw, p_std, p_fdr_std), ~ signif(., 4))),
    roi_interaction_sig  = roi_merged %>%
      filter(significant_raw) %>%
      select(roi, index, system, beta_raw, se_raw, t_raw, p_raw, p_fdr_raw,
             beta_std, p_fdr_std, fmi_raw) %>%
      mutate(across(c(beta_raw, se_raw, t_raw, beta_std, fmi_raw), ~ round(., 6)),
             across(c(p_raw, p_fdr_raw, p_fdr_std), ~ signif(., 4)))
  )
  if (!is.null(roi_slopes)) {
    sheet_list[["roi_simple_slopes"]] <- roi_slopes %>%
      mutate(across(c(beta, se, t_value, fmi), ~ round(., 6)),
             across(c(p_value, p_fdr_raw), ~ signif(., 4)))
  }

  write_xlsx(sheet_list, out_xlsx)
  cat(sprintf("\n已保存: %s\n", out_xlsx))
}

# ---- 5. 运行 ----
analyze_group_mi("depression", PHQ9_CUT, "Depression",
                 "_tables/depression_group_x_age_MI.xlsx")
analyze_group_mi("anxiety", GAD7_CUT, "Anxiety",
                 "_tables/anxiety_group_x_age_MI.xlsx")

cat("\n========== 全部完成 ==========\n")
