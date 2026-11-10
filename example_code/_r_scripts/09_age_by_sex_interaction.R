# ============================================================================
# 09_age_by_sex_interaction.R
# 年龄对马氏距离的效应是否存在性别差异 (方案B: 交互模型)
#
# 核心模型: ROI ~ age * gender + mri_equip + EMO_mean_FD
#   - age:gender 交互项 = 两性年龄斜率的差异 (性别调节的正式检验)
#   - 交互显著(FDR)的ROI, 用 emmeans::emtrends 提取男女各自的简单斜率
# 仅线性, 不做可视化
#
# 输出: _tables/age_by_sex_interaction.xlsx
# ============================================================================
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(broom)
library(writexl)
library(emmeans)

# ---- 0. 读取数据 ----
df <- read.csv("_data/simulated_mahalanobis_data.csv", stringsAsFactors = FALSE)
df$mri_equip <- factor(df$mri_equip)
# gender: 1 = 男, 2 = 女 (reference = 1)
df$gender    <- factor(df$gender, levels = c(1, 2), labels = c("male", "female"))
df$age_z     <- as.numeric(scale(df$age))

roi_cols <- grep("^ROI[0-9]+$", colnames(df), value = TRUE)
stopifnot(length(roi_cols) == 432)

df$mean_md <- rowMeans(df[, roi_cols], na.rm = TRUE)

cat("数据维度:", nrow(df), "x", ncol(df), "\n")
cat("ROI数量:", length(roi_cols), "\n")
cat("性别分布:\n"); print(table(df$gender, useNA = "ifany"))
cat("性别 x 设备 (检查混杂):\n"); print(table(df$gender, df$mri_equip))

# ---- 0.5 稳健提取 emtrends 结果的辅助函数 ----
# 不同 emmeans 版本统计列名可能为 t.ratio / z.ratio / statistic 等,
# 这里统一列名, 并保证所需列都存在 (缺失补 NA)
emtrends_df <- function(fit, var = "age") {
  d <- as.data.frame(
    summary(emtrends(fit, ~ gender, var = var), infer = c(TRUE, TRUE))
  )

  # 统一斜率列名 -> age.trend
  trend_col <- grep("\\.trend$", names(d), value = TRUE)[1]
  if (!is.na(trend_col) && trend_col != "age.trend") {
    names(d)[names(d) == trend_col] <- "age.trend"
  }
  # 统一统计量列名 -> t.ratio
  stat_cands <- intersect(c("t.ratio", "z.ratio", "t.statistic", "statistic"),
                          names(d))
  if (length(stat_cands) > 0 && stat_cands[1] != "t.ratio") {
    names(d)[names(d) == stat_cands[1]] <- "t.ratio"
  }
  # 保证所需列存在
  need <- c("age.trend", "SE", "df", "lower.CL", "upper.CL", "t.ratio", "p.value")
  for (cn in need) if (!cn %in% names(d)) d[[cn]] <- NA
  d
}

# ---- 1. 全局: 平均神经距离的年龄 x 性别交互 ----
cat("\n========== 全局: mean_md ~ age * gender ==========\n")
m_global <- lm(mean_md ~ age * gender + mri_equip + EMO_mean_FD, data = df)
print(summary(m_global))

# 提取交互项
sm_g <- summary(m_global)$coefficients
idx_g <- grep("^age:gender", rownames(sm_g))
global_int <- data.frame(
  term      = rownames(sm_g)[idx_g],
  beta      = sm_g[idx_g, "Estimate"],
  se        = sm_g[idx_g, "Std. Error"],
  t_value   = sm_g[idx_g, "t value"],
  p_value   = sm_g[idx_g, "Pr(>|t|)"],
  stringsAsFactors = FALSE
)

# 标准化年龄版本 (beta尺度可解释, p与raw一致)
m_global_z <- lm(mean_md ~ age_z * gender + mri_equip + EMO_mean_FD, data = df)
sm_gz <- summary(m_global_z)$coefficients
idx_gz <- grep("^age_z:gender", rownames(sm_gz))
global_int$beta_std <- sm_gz[idx_gz, "Estimate"]
global_int$se_std   <- sm_gz[idx_gz, "Std. Error"]

cat(sprintf("\n全局交互: beta = %.4f, p = %.2e | beta_std = %.4f\n",
            global_int$beta, global_int$p_value, global_int$beta_std))

# 全局简单斜率 (男女各自的 age 斜率, 基于交互模型, 全样本)
global_slopes <- emtrends_df(m_global, var = "age")
global_slopes$roi <- "mean_md"
cat("全局简单斜率:\n"); print(global_slopes)

# ---- 2. 每个ROI: 年龄 x 性别交互 ----
cat("\n========== ROI水平: 432个交互模型 ==========\n")
cat("进度: ")

# 提取交互项的小函数
extract_int <- function(fit, pattern) {
  sm <- summary(fit)$coefficients
  idx <- grep(pattern, rownames(sm))
  if (length(idx) == 0) {
    return(c(beta = NA, se = NA, t_value = NA, p_value = NA))
  }
  c(beta    = sm[idx[1], "Estimate"],
    se      = sm[idx[1], "Std. Error"],
    t_value = sm[idx[1], "t value"],
    p_value = sm[idx[1], "Pr(>|t|)"])
}

run_interaction <- function(age_var) {
  pat <- sprintf("^%s:gender", age_var)
  lapply(seq_along(roi_cols), function(i) {
    if (i %% 50 == 0) cat(i, ".. ")
    roi_name <- roi_cols[i]
    fml <- as.formula(sprintf("`%s` ~ %s * gender + mri_equip + EMO_mean_FD",
                              roi_name, age_var))
    fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
    if (is.null(fit)) {
      return(data.frame(roi = roi_name, index = i,
                        beta = NA, se = NA, t_value = NA, p_value = NA,
                        stringsAsFactors = FALSE))
    }
    v <- extract_int(fit, pat)
    data.frame(roi = roi_name, index = i,
               beta = v["beta"], se = v["se"],
               t_value = v["t_value"], p_value = v["p_value"],
               stringsAsFactors = FALSE)
  })
}

roi_raw <- bind_rows(run_interaction("age"))
cat("\n")
roi_raw <- as.data.frame(roi_raw)
roi_raw$p_fdr       <- p.adjust(roi_raw$p_value, method = "fdr")
roi_raw$significant <- roi_raw$p_fdr < 0.05
roi_raw$system      <- ifelse(roi_raw$index <= 400, "Cortical", "Subcortical")
cat(sprintf("Raw age 交互 FDR<0.05: %d / %d\n",
            sum(roi_raw$significant, na.rm = TRUE), nrow(roi_raw)))

roi_std <- bind_rows(run_interaction("age_z"))
cat("\n")
roi_std <- as.data.frame(roi_std)
roi_std$p_fdr       <- p.adjust(roi_std$p_value, method = "fdr")
roi_std$significant <- roi_std$p_fdr < 0.05
cat(sprintf("Std age 交互 FDR<0.05: %d / %d\n",
            sum(roi_std$significant, na.rm = TRUE), nrow(roi_std)))

# 合并 raw + std
roi_merged <- roi_raw %>%
  select(roi, index, system, beta, se, t_value, p_value, p_fdr, significant) %>%
  rename(beta_raw = beta, se_raw = se, t_raw = t_value,
         p_raw = p_value, p_fdr_raw = p_fdr, significant_raw = significant) %>%
  left_join(
    roi_std %>%
      select(roi, beta, se, t_value, p_value, p_fdr, significant) %>%
      rename(beta_std = beta, se_std = se, t_std = t_value,
             p_std = p_value, p_fdr_std = p_fdr, significant_std = significant),
    by = "roi"
  ) %>%
  arrange(p_fdr_raw)

# ---- 3. 交互显著ROI的简单斜率 ----
cat("\n========== 显著的ROI计算简单斜率 ==========\n")
sig_rois <- roi_merged %>% filter(significant_raw) %>% pull(roi)
cat(sprintf("显著交互ROI数: %d\n", length(sig_rois)))

roi_slopes <- NULL
if (length(sig_rois) > 0) {
  roi_slopes <- bind_rows(lapply(sig_rois, function(rn) {
    fml <- as.formula(sprintf("`%s` ~ age * gender + mri_equip + EMO_mean_FD", rn))
    fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    et <- emtrends_df(fit, var = "age")
    et$roi <- rn
    et
  }))
  # 附加完整模型里交互的p_fdr
  roi_slopes <- roi_slopes %>%
    left_join(roi_merged %>% select(roi, p_fdr_raw), by = "roi") %>%
    select(roi, gender, p_fdr_raw,
           any_of(c("age.trend", "SE", "df", "lower.CL", "upper.CL",
                    "t.ratio", "p.value")))
  print(head(roi_slopes))
} else {
  cat("无交互显著ROI, 跳过简单斜率。\n")
}

# ---- 4. 保存xlsx ----
dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

round_cols <- function(d, cols, digits = 6) {
  for (cn in cols) if (cn %in% names(d)) d[[cn]] <- round(as.numeric(d[[cn]]), digits)
  d
}

global_model_tbl <- tidy(m_global, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(., 6)))

sheet_list <- list(
  global_model          = global_model_tbl,
  global_interaction    = global_int %>% mutate(across(where(is.numeric), ~ round(., 6))),
  global_simple_slopes  = global_slopes %>%
    select(roi, gender,
           any_of(c("age.trend", "SE", "df", "lower.CL", "upper.CL",
                    "t.ratio", "p.value"))) %>%
    mutate(across(where(is.numeric), ~ round(., 6))),
  roi_interaction       = roi_merged %>%
    select(roi, index, system, beta_raw, se_raw, t_raw, p_raw, p_fdr_raw,
           significant_raw, beta_std, se_std, t_std, p_std, p_fdr_std,
           significant_std) %>%
    mutate(across(c(beta_raw, se_raw, t_raw, beta_std, se_std, t_std), ~ round(., 6)),
           across(c(p_raw, p_fdr_raw, p_std, p_fdr_std), ~ signif(., 4))),
  roi_interaction_sig   = roi_merged %>%
    filter(significant_raw) %>%
    select(roi, index, system, beta_raw, se_raw, t_raw, p_raw, p_fdr_raw,
           beta_std, p_fdr_std) %>%
    mutate(across(c(beta_raw, se_raw, t_raw, beta_std), ~ round(., 6)),
           across(c(p_raw, p_fdr_raw, p_fdr_std), ~ signif(., 4)))
)

if (!is.null(roi_slopes)) {
  sheet_list[["roi_simple_slopes"]] <- roi_slopes %>%
    mutate(across(c(age.trend, SE, df, lower.CL, upper.CL, t.ratio), ~ round(., 6)),
           across(c(p.value, p_fdr_raw), ~ signif(., 4)))
}

write_xlsx(sheet_list, "_tables/age_by_sex_interaction.xlsx")
cat("\n已保存: _tables/age_by_sex_interaction.xlsx\n")

# ---- 5. 摘要 ----
cat("\n========== 摘要 ==========\n")
cat(sprintf("全局交互: beta = %.4f, p = %.2e\n",
            global_int$beta, global_int$p_value))
cat(sprintf("ROI交互显著 (FDR<0.05): %d / 432\n",
            sum(roi_merged$significant_raw, na.rm = TRUE)))
cat(sprintf("  正向(女斜率>男): %d, 负向(女斜率<男): %d\n",
            sum(roi_merged$significant_raw & roi_merged$beta_raw > 0, na.rm = TRUE),
            sum(roi_merged$significant_raw & roi_merged$beta_raw < 0, na.rm = TRUE)))
cat("  注: gender reference = male, 交互beta = female斜率 - male斜率\n")
cat("==================================\n")
