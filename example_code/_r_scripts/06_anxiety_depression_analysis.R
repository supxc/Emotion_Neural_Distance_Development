# ============================================================================
# 06_anxiety_depression_analysis.R
# 神经距离 (ROI马氏距离) 对焦虑/抑郁得分的效应分析
# 主分析 = 多重插补 (mice, m=10) + Rubin's rule pooled 估计
# 所有可视化与表格均使用 pooled 统计量
#
# 分析1 (全局): 全脑平均马氏距离 -> anxiety / depression
# 分析2 (ROI):  每个ROI马氏距离 -> anxiety / depression
# 控制变量: gender + mri_equip + EMO_mean_FD
# ============================================================================
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(writexl)
library(mice)

# ---- 0. 读取数据 ----
df <- read.csv("_data/simulated_mahalanobis_data.csv", stringsAsFactors = FALSE)
df$mri_equip <- factor(df$mri_equip)
df$gender    <- factor(df$gender)

roi_cols <- grep("^ROI[0-9]+$", colnames(df), value = TRUE)
stopifnot(length(roi_cols) == 432)

df$mean_md <- rowMeans(df[, roi_cols], na.rm = TRUE)

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
  d$gender    <- factor(d$gender)
  # 补回ROI列 (插补模型未包含, 但分析需要)
  d[, roi_cols] <- df[, roi_cols]
  d
})
cat("插补完成:", m_mi, "个数据集\n")

# ---- 2. Rubin's rule 合并函数 ----
rubin_pool <- function(betas, ses, m) {
  Qbar <- mean(betas)
  Ubar <- mean(ses^2)
  B    <- sum((betas - Qbar)^2) / (m - 1)
  T    <- Ubar + (1 + 1 / m) * B
  r    <- (1 + 1 / m) * B / Ubar
  v    <- (m - 1) * (1 + 1 / r)^2
  t    <- Qbar / sqrt(T)
  p    <- 2 * pt(-abs(t), df = v)
  fmi  <- (r + 2 / (v + 3)) / (r + 1)
  c(beta = Qbar, se = sqrt(T), t_value = t, p_value = p, fmi = fmi)
}

# ---- 3. 主分析函数 (基于MI pooled) ----
analyze_behavior <- function(behavior_var) {

  cat(sprintf("\n========== 分析: %s (MI pooled) ==========\n", behavior_var))

  dir.create("_plots", showWarnings = FALSE, recursive = TRUE)
  dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

  # 3a. 全局模型: 每个插补数据集跑, 然后pooled
  fml_global <- as.formula(sprintf("%s ~ mean_md + gender + mri_equip + EMO_mean_FD",
                                   behavior_var))
  global_fits <- lapply(comp_list, function(d) {
    tryCatch(lm(fml_global, data = d), error = function(e) NULL)
  })
  global_beta <- sapply(global_fits, function(f)
    coef(summary(f))["mean_md", "Estimate"])
  global_se   <- sapply(global_fits, function(f)
    coef(summary(f))["mean_md", "Std. Error"])
  global_pooled <- rubin_pool(global_beta, global_se, m_mi)

  cat(sprintf("全局模型: beta = %.4f, se = %.4f, t = %.3f, p = %.2e, FMI = %.3f\n",
              global_pooled["beta"], global_pooled["se"], global_pooled["t_value"],
              global_pooled["p_value"], global_pooled["fmi"]))

  # 3b. 全局 partial regression 图
  # 散点用第一个插补数据集, 回归线用pooled斜率
  d_plot <- comp_list[[1]]
  m_cov <- lm(as.formula(sprintf("%s ~ gender + mri_equip + EMO_mean_FD",
                                 behavior_var)), data = d_plot)
  d_plot$behav_resid <- residuals(m_cov)
  d_plot$mean_md_z <- as.numeric(scale(d_plot$mean_md))

  # pooled 回归线: 截距取残差均值, 斜率取 pooled beta
  intercept_pool <- mean(d_plot$behav_resid) - global_pooled["beta"] * mean(d_plot$mean_md)

  anno_text <- sprintf("pooled beta = %.4f, p = %.2e",
                       global_pooled["beta"], global_pooled["p_value"])

  p_global <- ggplot(d_plot, aes(x = mean_md, y = behav_resid)) +
    geom_point(alpha = 0.35, size = 0.8, color = "#3366B0") +
    geom_abline(intercept = intercept_pool, slope = global_pooled["beta"],
                color = "#D7263D", linewidth = 1.2) +
    annotate("text", x = -Inf, y = Inf, label = anno_text,
             hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic") +
    labs(x = "Mean Mahalanobis Distance", y = sprintf("%s (residualized)", behavior_var),
         title = sprintf("Mean MD Effect on %s (MI pooled)", behavior_var)) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))

  fname <- sprintf("_plots/%s_vs_global_meanMD.jpg", behavior_var)
  ggsave(fname, p_global, width = 8, height = 6, dpi = 800)
  cat(sprintf("已保存: %s\n", fname))

  # 3c. ROI水平 mass-univariate: 每个插补数据集跑, pooled + FDR
  for (scale_type in c("raw", "std")) {
    cat(sprintf("\n--- ROI水平 (%s) ---\n", scale_type))
    cat("进度: ")

    beta_mat <- se_mat <- matrix(NA, nrow = length(roi_cols), ncol = m_mi)

    for (k in 1:m_mi) {
      d <- comp_list[[k]]
      for (i in seq_along(roi_cols)) {
        if (i %% 100 == 0 && k == 1) cat(i, ".. ")
        roi_raw <- d[[roi_cols[i]]]
        roi_sd  <- sd(roi_raw, na.rm = TRUE)
        val     <- if (scale_type == "raw") roi_raw else {
          if (roi_sd > 1e-10) as.numeric(scale(roi_raw)) else rep(0, nrow(d))
        }
        fml <- as.formula(sprintf("%s ~ roi_val + gender + mri_equip + EMO_mean_FD",
                                  behavior_var))
        fit <- tryCatch(lm(fml, data = data.frame(d, roi_val = val)),
                        error = function(e) NULL)
        if (is.null(fit)) next
        sm <- summary(fit)$coefficients
        if (!"roi_val" %in% rownames(sm)) next
        beta_mat[i, k] <- sm["roi_val", "Estimate"]
        se_mat[i, k]   <- sm["roi_val", "Std. Error"]
      }
    }
    cat("\n")

    pooled <- t(sapply(seq_along(roi_cols), function(i) {
      if (all(is.na(beta_mat[i, ]))) {
        c(beta = NA, se = NA, t_value = NA, p_value = NA, fmi = NA)
      } else {
        rubin_pool(beta_mat[i, ], se_mat[i, ], m_mi)
      }
    }))

    roi_df <- data.frame(roi = roi_cols, index = 1:length(roi_cols), pooled,
                         stringsAsFactors = FALSE)
    # 确保数值列为 numeric (sapply合并可能强转为字符)
    num_cols <- c("beta", "se", "t_value", "p_value", "fmi")
    for (cn in num_cols) roi_df[[cn]] <- as.numeric(roi_df[[cn]])
    roi_df$p_fdr       <- p.adjust(roi_df$p_value, method = "fdr")
    roi_df$significant <- roi_df$p_fdr < 0.05
    roi_df$direction   <- ifelse(roi_df$beta > 0, "positive", "negative")
    roi_df$system      <- ifelse(roi_df$index <= 400, "Cortical", "Subcortical")

    roi_df <- roi_df %>%
      mutate(
        color_group = case_when(
          significant & beta > 0 ~ "FDR < 0.05 (+)",
          significant & beta < 0 ~ "FDR < 0.05 (-)",
          TRUE                     ~ "n.s."
        ),
        shape_group = system
      )

    cat(sprintf("FDR < 0.05: %d / %d\n", sum(roi_df$significant, na.rm = TRUE),
                nrow(roi_df)))

    # 可视化
    sig_colors <- c("FDR < 0.05 (+)" = "#D7263D", "FDR < 0.05 (-)" = "#1B5299",
                    "n.s."            = "grey70")
    sig_shapes <- c("Cortical" = 16, "Subcortical" = 17)

    # Manhattan (pooled p)
    fdr_thresh <- max(roi_df$p_value[roi_df$significant], na.rm = TRUE)
    if (!is.finite(fdr_thresh)) fdr_thresh <- 0.05 / 432
    p_line <- -log10(fdr_thresh)
    ymax   <- max(-log10(roi_df$p_value), na.rm = TRUE)
    p_man  <- ggplot(roi_df, aes(x = index, y = -log10(p_value))) +
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
               label = "Subcortical", size = 3.2, color = "grey40") +
      labs(x = "ROI Index", y = expression(-log[10](italic(p))),
           title = sprintf("Manhattan: ROI MD -> %s (MI pooled)", behavior_var),
           color = NULL, shape = NULL) +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = c(0.82, 0.85),
            legend.background = element_rect(fill = "white", color = "grey80"),
            legend.key.size = unit(0.35, "cm"))

    f_man <- sprintf("_plots/manhattan_%s_%s.jpg", behavior_var, scale_type)
    ggsave(f_man, p_man, width = 14, height = 6, dpi = 800)
    cat(sprintf("已保存: %s\n", f_man))

    # Volcano (pooled beta/p)
    n_pos <- sum(roi_df$significant & roi_df$beta > 0, na.rm = TRUE)
    n_neg <- sum(roi_df$significant & roi_df$beta < 0, na.rm = TRUE)
    p_vol <- ggplot(roi_df, aes(x = beta, y = -log10(p_value))) +
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
           title = sprintf("Volcano: ROI MD -> %s (MI pooled)", behavior_var),
           color = NULL, shape = NULL) +
      theme_bw(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = c(0.82, 0.85),
            legend.background = element_rect(fill = "white", color = "grey80"),
            legend.key.size = unit(0.35, "cm"))

    f_vol <- sprintf("_plots/volcano_%s_%s.jpg", behavior_var, scale_type)
    ggsave(f_vol, p_vol, width = 9, height = 7, dpi = 800)
    cat(sprintf("已保存: %s\n", f_vol))

    # Beta 分布 (pooled beta)
    p_hist <- ggplot(roi_df, aes(x = beta)) +
      geom_histogram(bins = 60, fill = "#4682B4", color = "white", alpha = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#D7263D", linewidth = 0.8) +
      geom_vline(xintercept = mean(roi_df$beta, na.rm = TRUE),
                 linetype = "dotted", color = "grey30", linewidth = 0.8) +
      annotate("text", x = mean(roi_df$beta, na.rm = TRUE),
               y = Inf, vjust = 2, size = 3.8, color = "grey30",
               label = sprintf("Mean = %.4f", mean(roi_df$beta, na.rm = TRUE))) +
      labs(x = expression(beta), y = "Number of ROIs",
           title = sprintf("Beta Distribution: ROI MD -> %s (MI pooled)", behavior_var)) +
      theme_bw(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5))

    f_hist <- sprintf("_plots/beta_distribution_%s_%s.jpg", behavior_var, scale_type)
    ggsave(f_hist, p_hist, width = 8, height = 6, dpi = 800)
    cat(sprintf("已保存: %s\n", f_hist))

    # 汇总打印
    cat(sprintf("  正向: %d,  负向: %d\n  beta范围: [%.4f, %.4f]\n",
                n_pos, n_neg,
                min(roi_df$beta, na.rm = TRUE), max(roi_df$beta, na.rm = TRUE)))

    # 保存结果
    assign(sprintf("roi_%s_%s", behavior_var, scale_type), roi_df, envir = .GlobalEnv)
    assign(sprintf("global_pooled_%s", behavior_var), global_pooled, envir = .GlobalEnv)
  }
}

# ---- 4. 运行两个临床变量 ----
analyze_behavior("anxiety")
analyze_behavior("depression")

# ---- 5. 保存xlsx (全部为pooled统计量) ----
prepare_out <- function(dat) {
  dat %>%
    select(roi, index, system, beta, se, t_value, p_value, p_fdr,
           significant, direction, fmi) %>%
    mutate(across(c(beta, se, t_value, fmi), ~ round(., 6)),
           across(c(p_value, p_fdr), ~ signif(., 4)))
}

sheet_list <- list()

# 全局模型 (pooled)
global_list <- list()
for (behav in c("anxiety", "depression")) {
  gp <- get(sprintf("global_pooled_%s", behav))
  global_list[[behav]] <- data.frame(
    term = "mean_md", t(gp)
  ) %>% mutate(across(where(is.numeric), ~ round(., 6)))
}

for (behav in c("anxiety", "depression")) {
  for (scale_type in c("raw", "std")) {
    obj_name <- sprintf("roi_%s_%s", behav, scale_type)
    dat <- get(obj_name)
    out <- prepare_out(dat) %>%
      rename(t = t_value, p = p_value)
    sheet_list[[sprintf("%s_%s", behav, scale_type)]] <- out %>% arrange(p_fdr)
  }
}

write_xlsx(c(global_list, sheet_list), "_tables/anxiety_depression_effect_results.xlsx")
cat("\n已保存: _tables/anxiety_depression_effect_results.xlsx\n")

cat("\n========== 全部完成 ==========\n")
