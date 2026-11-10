# ============================================================================
# 02_age_effect_analysis.R
# 年龄对马氏距离的线性效应分析
#
# 分析1: 全局 - 年龄 vs 全脑平均马氏距离 (控制 gender + mri_equip + FD)
# 分析2: ROI水平 - 每个ROI单独做 mass-univariate + FDR校正
# ============================================================================
rm(list = ls())

gc()


library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(writexl)

# ---- 0. 读取数据 ----
df <- read.csv("_data/simulated_mahalanobis_data.csv", stringsAsFactors = FALSE)
df$mri_equip <- factor(df$mri_equip)
df$gender    <- factor(df$gender)

# ROI列索引
roi_cols <- grep("^ROI[0-9]+$", colnames(df), value = TRUE)
stopifnot(length(roi_cols) == 432)

cat("数据维度:", nrow(df), "x", ncol(df), "\n")
cat("ROI数量:", length(roi_cols), "\n")

# ---- 1. 全局分析: 年龄 vs 全脑平均马氏距离 ----
df$mean_md <- rowMeans(df[, roi_cols], na.rm = TRUE)
df$age_z <- as.numeric(scale(df$age))

m_global <- lm(mean_md ~ age + gender + mri_equip + EMO_mean_FD, data = df)
cat("\n========== 全局模型 ==========\n")
print(summary(m_global))

# ---- 2. 全局可视化 ----
dir.create("_plots", showWarnings = FALSE, recursive = TRUE)

# 2a. Partial regression plot: 先回归掉协变量, 再画年龄与残差
m_cov <- lm(mean_md ~ gender + mri_equip + EMO_mean_FD, data = df)
df$mean_md_resid <- residuals(m_cov)
m_age_only <- lm(mean_md_resid ~ age, data = df)

b_global <- coef(summary(m_age_only))["age", ]
anno_text <- sprintf("beta = %.4f, p = %.2e",
                     b_global["Estimate"], b_global["Pr(>|t|)"])

p_global <- ggplot(df, aes(x = age, y = mean_md_resid)) +
  geom_point(alpha = 0.35, size = 0.8, color = "#3366B0") +
  geom_smooth(method = "lm", se = TRUE, color = "#D7263D", linewidth = 1.2,
              fill = "#D7263D", alpha = 0.12) +
  annotate("text", x = -Inf, y = Inf, label = anno_text,
           hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic") +
  labs(x = "Age (years)", y = "Mean Mahalanobis Distance (residualized)",
       title = "Age Effect on Whole-Brain Pattern Dissimilarity") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("_plots/age_vs_global_meanMD.jpg", p_global,
       width = 8, height = 6, dpi = 800)
cat("已保存: _plots/age_vs_global_meanMD.jpg\n")

# 2b. Standardized age 版本
m_age_z <- lm(mean_md_resid ~ age_z, data = df)
b_global_z <- coef(summary(m_age_z))["age_z", ]
anno_text_z <- sprintf("beta_std = %.4f, p = %.2e",
                       b_global_z["Estimate"], b_global_z["Pr(>|t|)"])

p_global_z <- ggplot(df, aes(x = age_z, y = mean_md_resid)) +
  geom_point(alpha = 0.35, size = 0.8, color = "#3366B0") +
  geom_smooth(method = "lm", se = TRUE, color = "#D7263D", linewidth = 1.2,
              fill = "#D7263D", alpha = 0.12) +
  annotate("text", x = -Inf, y = Inf, label = anno_text_z,
           hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic") +
  labs(x = "Age (z-scored)", y = "Mean Mahalanobis Distance (residualized)",
       title = "Age Effect on Whole-Brain Pattern Dissimilarity (Standardized)") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("_plots/age_vs_global_meanMD_std.jpg", p_global_z,
       width = 8, height = 6, dpi = 800)
cat("已保存: _plots/age_vs_global_meanMD_std.jpg\n")

# ---- 3. ROI水平 mass-univariate (raw age) ----
cat("\n========== ROI水平分析 (432个模型, raw age) ==========\n")
cat("进度: ")

run_roi_models <- function(age_var, label) {
  lapply(seq_along(roi_cols), function(i) {
    if (i %% 50 == 0) cat(i, ".. ")
    roi_name <- roi_cols[i]
    fml <- as.formula(sprintf("`%s` ~ %s + gender + mri_equip + EMO_mean_FD",
                              roi_name, age_var))
    fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
    if (is.null(fit)) {
      return(data.frame(roi = roi_name, index = i,
                        beta = NA, se = NA, t_value = NA, p_value = NA,
                        stringsAsFactors = FALSE))
    }
    sm <- summary(fit)$coefficients
    idx <- which(rownames(sm) == age_var)
    if (length(idx) == 0) idx <- grep(age_var, rownames(sm), fixed = TRUE)[1]
    data.frame(
      roi     = roi_name,
      index   = i,
      beta    = sm[idx, "Estimate"],
      se      = sm[idx, "Std. Error"],
      t_value = sm[idx, "t value"],
      p_value = sm[idx, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  })
}

roi_raw <- bind_rows(run_roi_models("age", "raw"))
cat("\n")
roi_raw <- as.data.frame(roi_raw)
roi_raw$p_fdr <- p.adjust(roi_raw$p_value, method = "fdr")

cat(sprintf("Raw age FDR < 0.05: %d / %d ROIs\n",
            sum(roi_raw$p_fdr < 0.05), nrow(roi_raw)))

# ---- 4. ROI水平 mass-univariate (standardized age) ----
cat("\n========== ROI水平分析 (432个模型, standardized age) ==========\n")
cat("进度: ")

roi_std <- bind_rows(run_roi_models("age_z", "std"))
cat("\n")
roi_std <- as.data.frame(roi_std)
roi_std$p_fdr <- p.adjust(roi_std$p_value, method = "fdr")

cat(sprintf("Std age FDR < 0.05: %d / %d ROIs\n",
            sum(roi_std$p_fdr < 0.05), nrow(roi_std)))

# ---- 5. 合并标注 ----
prepare_roi_df <- function(dat, label) {
  dat$significant <- dat$p_fdr < 0.05
  dat$direction   <- ifelse(dat$beta > 0, "positive", "negative")
  dat$system      <- ifelse(dat$index <= 400, "Cortical", "Subcortical")
  dat$label       <- label
  # 颜色分类: n.s. 统一灰色, 显著用红/蓝; 形状区分皮层/皮层下
  dat <- dat %>%
    mutate(
      color_group = case_when(
        significant & beta > 0 ~ "FDR < 0.05 (+)",
        significant & beta < 0 ~ "FDR < 0.05 (-)",
        TRUE                     ~ "n.s."
      ),
      shape_group = system
    )
  dat
}

roi_raw <- prepare_roi_df(roi_raw, "raw")
roi_std <- prepare_roi_df(roi_std, "std")

# ---- 6. 可视化函数 ----
sig_colors <- c(
  "FDR < 0.05 (+)" = "#D7263D",
  "FDR < 0.05 (-)" = "#1B5299",
  "n.s."            = "grey70"
)
sig_shapes <- c("Cortical" = 16, "Subcortical" = 17)

make_manhattan <- function(dat, title_suffix, filepath) {
  fdr_thresh <- max(dat$p_value[dat$significant], na.rm = TRUE)
  if (!is.finite(fdr_thresh)) fdr_thresh <- 0.05 / 432
  p_line <- -log10(fdr_thresh)
  ymax <- max(-log10(dat$p_value), na.rm = TRUE)
  ggplot(dat, aes(x = index, y = -log10(p_value))) +
    geom_vline(xintercept = 400.5, color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = color_group, shape = shape_group),
               size = 0.9, alpha = 0.85) +
    geom_hline(yintercept = p_line, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    scale_color_manual(values = sig_colors) +
    scale_shape_manual(values = sig_shapes) +
    annotate("text", x = 200, y = ymax * 0.96,
             label = "Cortical (Schaefer400)", size = 3.2, color = "grey40") +
    annotate("text", x = 416, y = ymax * 0.96,
             label = "Subcortical",           size = 3.2, color = "grey40") +
    labs(x = "ROI Index", y = expression(-log[10](italic(p))),
         title = paste("Manhattan Plot:", title_suffix),
         color = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = c(0.82, 0.85),
          legend.background = element_rect(fill = "white", color = "grey80"),
          legend.key.size = unit(0.35, "cm"))
}

make_volcano <- function(dat, title_suffix, filepath) {
  n_pos <- sum(dat$significant & dat$beta > 0, na.rm = TRUE)
  n_neg <- sum(dat$significant & dat$beta < 0, na.rm = TRUE)
  fdr_thresh <- max(dat$p_value[dat$significant], na.rm = TRUE)
  if (!is.finite(fdr_thresh)) fdr_thresh <- 0.05 / 432
  p_line <- -log10(fdr_thresh)
  ggplot(dat, aes(x = beta, y = -log10(p_value))) +
    geom_point(aes(color = color_group, shape = shape_group),
               size = 1.1, alpha = 0.85) +
    geom_hline(yintercept = p_line, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "grey70", linewidth = 0.5) +
    scale_color_manual(values = sig_colors) +
    scale_shape_manual(values = sig_shapes) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("FDR (+): %d  |  FDR (-): %d", n_pos, n_neg),
             hjust = -0.05, vjust = 1.5, size = 4) +
    labs(x = expression(beta), y = expression(-log[10](italic(p))),
         title = paste("Volcano Plot:", title_suffix),
         color = NULL, shape = NULL) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = c(0.82, 0.85),
          legend.background = element_rect(fill = "white", color = "grey80"),
          legend.key.size = unit(0.35, "cm"))
}

make_beta_hist <- function(dat, title_suffix, filepath) {
  ggplot(dat, aes(x = beta)) +
    geom_histogram(bins = 60, fill = "#4682B4", color = "white", alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#D7263D", linewidth = 0.8) +
    geom_vline(xintercept = mean(dat$beta, na.rm = TRUE),
               linetype = "dotted", color = "grey30", linewidth = 0.8) +
    annotate("text", x = mean(dat$beta, na.rm = TRUE),
             y = Inf, vjust = 2, size = 3.8, color = "grey30",
             label = sprintf("Mean = %.4f", mean(dat$beta, na.rm = TRUE))) +
    labs(x = expression(beta), y = "Number of ROIs",
         title = paste("Beta Distribution:", title_suffix)) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
}

# ---- 7. 生成全部图表 ----
# Manhattan & Volcano: p值不变, 只画一套
p_man <- make_manhattan(roi_raw, "", "_plots/manhattan_age.jpg")
ggsave("_plots/manhattan_age.jpg", p_man, width = 14, height = 6, dpi = 800)
cat("已保存: _plots/manhattan_age.jpg\n")

p_vol <- make_volcano(roi_raw, "", "_plots/volcano_age.jpg")
ggsave("_plots/volcano_age.jpg", p_vol, width = 9, height = 7, dpi = 800)
cat("已保存: _plots/volcano_age.jpg\n")

p_vol_std <- make_volcano(roi_std, "(Standardized)", "_plots/volcano_age_std.jpg")
ggsave("_plots/volcano_age_std.jpg", p_vol_std, width = 9, height = 7, dpi = 800)
cat("已保存: _plots/volcano_age_std.jpg\n")

# Beta分布: raw 和 std 各一套
p_hist_raw <- make_beta_hist(roi_raw, "Raw Age",
                             "_plots/beta_distribution_raw_age.jpg")
ggsave("_plots/beta_distribution_raw_age.jpg", p_hist_raw, width = 8, height = 6, dpi = 800)
cat("已保存: _plots/beta_distribution_raw_age.jpg\n")

p_hist_std <- make_beta_hist(roi_std, "Standardized Age",
                             "_plots/beta_distribution_std_age.jpg")
ggsave("_plots/beta_distribution_std_age.jpg", p_hist_std, width = 8, height = 6, dpi = 800)
cat("已保存: _plots/beta_distribution_std_age.jpg\n")

# ---- 8. 保存结果到xlsx ----
dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

# 8a. 全局模型
global_tbl <- tidy(m_global, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(., 6)))

# 8b. Raw age 每个ROI完整结果
roi_raw_out <- roi_raw %>%
  select(roi, index, system, beta, se, t_value, p_value, p_fdr, significant, direction) %>%
  mutate(across(c(beta, se, t_value), ~ round(., 6)),
         across(c(p_value, p_fdr), ~ signif(., 4))) %>%
  rename(beta_raw = beta, se_raw = se, t_raw = t_value,
         p_raw = p_value, p_fdr_raw = p_fdr)

# 8c. Standardized age 每个ROI完整结果
roi_std_out <- roi_std %>%
  select(roi, index, system, beta, se, t_value, p_value, p_fdr, significant, direction) %>%
  mutate(across(c(beta, se, t_value), ~ round(., 6)),
         across(c(p_value, p_fdr), ~ signif(., 4))) %>%
  rename(beta_std = beta, se_std = se, t_std = t_value,
         p_std = p_value, p_fdr_std = p_fdr)

# 8d. 合并两张ROI表
roi_merged <- roi_raw_out %>%
  left_join(
    roi_std_out %>% select(roi, beta_std, se_std, t_std, p_std, p_fdr_std),
    by = "roi"
  ) %>%
  arrange(p_fdr_raw)

# 8e. 显著ROI
sig_raw <- roi_merged %>% filter(significant == TRUE)
sig_std <- roi_merged %>% filter(p_fdr_std < 0.05)

write_xlsx(
  list(
    global_model      = global_tbl,
    roi_raw_results   = roi_raw_out %>% arrange(p_fdr_raw),
    roi_std_results   = roi_std_out %>% arrange(p_fdr_std),
    roi_merged        = roi_merged,
    sig_raw           = sig_raw,
    sig_std           = sig_std
  ),
  "_tables/age_effect_results.xlsx"
)
cat("已保存: _tables/age_effect_results.xlsx\n")

# ---- 9. 输出摘要 ----
cat("\n========== 分析摘要 ==========\n")
cat(sprintf("全局模型 age beta = %.4f, p = %.2e\n",
            coef(summary(m_global))["age", "Estimate"],
            coef(summary(m_global))["age", "Pr(>|t|)"]))
cat(sprintf("\n--- Raw Age ---\n"))
cat(sprintf("显著ROI (FDR<0.05): %d / %d\n",
            sum(roi_raw$significant), nrow(roi_raw)))
cat(sprintf("  正向: %d, 负向: %d\n",
            sum(roi_raw$significant & roi_raw$beta > 0, na.rm = TRUE),
            sum(roi_raw$significant & roi_raw$beta < 0, na.rm = TRUE)))
cat(sprintf("Beta范围: [%.4f, %.4f]\n",
            min(roi_raw$beta, na.rm = TRUE), max(roi_raw$beta, na.rm = TRUE)))
cat(sprintf("\n--- Standardized Age ---\n"))
cat(sprintf("显著ROI (FDR<0.05): %d / %d\n",
            sum(roi_std$significant), nrow(roi_std)))
cat(sprintf("  正向: %d, 负向: %d\n",
            sum(roi_std$significant & roi_std$beta > 0, na.rm = TRUE),
            sum(roi_std$significant & roi_std$beta < 0, na.rm = TRUE)))
cat(sprintf("Beta范围: [%.4f, %.4f]\n",
            min(roi_std$beta, na.rm = TRUE), max(roi_std$beta, na.rm = TRUE)))
cat("==================================\n")
