# ============================================================================
# 04_behavior_effect_analysis.R
# 神经距离 (ROI马氏距离) 对情绪行为变量的效应分析
#
# 分析1 (全局): 全脑平均马氏距离 -> regulate_effect / generation_effect
# 分析2 (ROI):  每个ROI马氏距离 -> regulate_effect / generation_effect
# 控制变量: gender + mri_equip + EMO_mean_FD
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

roi_cols <- grep("^ROI[0-9]+$", colnames(df), value = TRUE)
stopifnot(length(roi_cols) == 432)

df$mean_md <- rowMeans(df[, roi_cols], na.rm = TRUE)

cat("数据维度:", nrow(df), "x", ncol(df), "\n")
cat("ROI数量:", length(roi_cols), "\n")

# ---- 1. 主分析函数 ----
analyze_behavior <- function(behavior_var, df, roi_cols) {

  cat(sprintf("\n========== 分析: %s ==========\n", behavior_var))

  dir.create("_plots", showWarnings = FALSE, recursive = TRUE)

  # 1a. 全局模型 (raw scale)
  fml_global <- as.formula(sprintf("%s ~ mean_md + gender + mri_equip + EMO_mean_FD",
                                   behavior_var))
  m_global <- lm(fml_global, data = df)
  cat("\n--- 全局模型 ---\n")
  print(summary(m_global))

  # 1b. 全局 partial regression 图
  m_cov <- lm(as.formula(sprintf("%s ~ gender + mri_equip + EMO_mean_FD",
                                 behavior_var)), data = df)
  df$behav_resid <- residuals(m_cov)

  m_md <- lm(behav_resid ~ mean_md, data = df)
  b_global <- coef(summary(m_md))["mean_md", ]
  anno_text <- sprintf("beta = %.4f, p = %.2e",
                       b_global["Estimate"], b_global["Pr(>|t|)"])

  p_global <- ggplot(df, aes(x = mean_md, y = behav_resid)) +
    geom_point(alpha = 0.35, size = 0.8, color = "#3366B0") +
    geom_smooth(method = "lm", se = TRUE, color = "#D7263D", linewidth = 1.2,
                fill = "#D7263D", alpha = 0.12) +
    annotate("text", x = -Inf, y = Inf, label = anno_text,
             hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic") +
    labs(x = "Mean Mahalanobis Distance", y = sprintf("%s (residualized)", behavior_var),
         title = sprintf("Mean MD Effect on %s", behavior_var)) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))

  fname <- sprintf("_plots/%s_vs_global_meanMD.jpg", behavior_var)
  ggsave(fname, p_global, width = 8, height = 6, dpi = 800)
  cat(sprintf("已保存: %s\n", fname))

  # 1c. 全局模型 std 版本
  df$mean_md_z <- as.numeric(scale(df$mean_md))
  m_md_z <- lm(behav_resid ~ mean_md_z, data = df)
  b_global_z <- coef(summary(m_md_z))["mean_md_z", ]
  anno_text_z <- sprintf("beta_std = %.4f, p = %.2e",
                         b_global_z["Estimate"], b_global_z["Pr(>|t|)"])

  p_global_z <- ggplot(df, aes(x = mean_md_z, y = behav_resid)) +
    geom_point(alpha = 0.35, size = 0.8, color = "#3366B0") +
    geom_smooth(method = "lm", se = TRUE, color = "#D7263D", linewidth = 1.2,
                fill = "#D7263D", alpha = 0.12) +
    annotate("text", x = -Inf, y = Inf, label = anno_text_z,
             hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic") +
    labs(x = "Mean Mahalanobis Distance (z-scored)",
         y = sprintf("%s (residualized)", behavior_var),
         title = sprintf("Mean MD Effect on %s (Standardized)", behavior_var)) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))

  fname_z <- sprintf("_plots/%s_vs_global_meanMD_std.jpg", behavior_var)
  ggsave(fname_z, p_global_z, width = 8, height = 6, dpi = 800)
  cat(sprintf("已保存: %s\n", fname_z))

  # 1d. ROI水平 mass-univariate: raw + std
  for (scale_type in c("raw", "std")) {
    predictor <- if (scale_type == "raw") "roi_val" else "roi_val_z"
    cat(sprintf("\n--- ROI水平 (%s) ---\n", scale_type))
    cat("进度: ")

    roi_results <- lapply(seq_along(roi_cols), function(i) {
      if (i %% 50 == 0) cat(i, ".. ")
      roi_raw <- df[[roi_cols[i]]]
      roi_sd  <- sd(roi_raw, na.rm = TRUE)
      roi_z   <- if (roi_sd > 1e-10) as.numeric(scale(roi_raw)) else rep(0, length(roi_raw))
      val     <- if (scale_type == "raw") roi_raw else roi_z

      fml <- as.formula(sprintf("%s ~ roi_val + gender + mri_equip + EMO_mean_FD",
                                behavior_var))
      fit <- tryCatch(
        lm(fml, data = data.frame(df, roi_val = val)),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(data.frame(roi = roi_cols[i], index = i,
                          beta = NA, se = NA, t_value = NA, p_value = NA,
                          stringsAsFactors = FALSE))
      }
      sm <- summary(fit)$coefficients
      if (!"roi_val" %in% rownames(sm)) {
        return(data.frame(roi = roi_cols[i], index = i,
                          beta = NA, se = NA, t_value = NA, p_value = NA,
                          stringsAsFactors = FALSE))
      }
      data.frame(
        roi     = roi_cols[i],
        index   = i,
        beta    = sm["roi_val", "Estimate"],
        se      = sm["roi_val", "Std. Error"],
        t_value = sm["roi_val", "t value"],
        p_value = sm["roi_val", "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    })
    cat("\n")

    roi_df <- bind_rows(roi_results)
    roi_df <- as.data.frame(roi_df)
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

    cat(sprintf("FDR < 0.05: %d / %d\n", sum(roi_df$significant), nrow(roi_df)))

    # 可视化
    sig_colors <- c("FDR < 0.05 (+)" = "#D7263D", "FDR < 0.05 (-)" = "#1B5299",
                    "n.s."            = "grey70")
    sig_shapes <- c("Cortical" = 16, "Subcortical" = 17)

    # Manhattan
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
           title = sprintf("Manhattan: ROI MD -> %s", behavior_var),
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

    # Volcano
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
           title = sprintf("Volcano: ROI MD -> %s", behavior_var),
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

    # Beta 分布
    p_hist <- ggplot(roi_df, aes(x = beta)) +
      geom_histogram(bins = 60, fill = "#4682B4", color = "white", alpha = 0.9) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#D7263D", linewidth = 0.8) +
      geom_vline(xintercept = mean(roi_df$beta, na.rm = TRUE),
                 linetype = "dotted", color = "grey30", linewidth = 0.8) +
      annotate("text", x = mean(roi_df$beta, na.rm = TRUE),
               y = Inf, vjust = 2, size = 3.8, color = "grey30",
               label = sprintf("Mean = %.4f", mean(roi_df$beta, na.rm = TRUE))) +
      labs(x = expression(beta), y = "Number of ROIs",
           title = sprintf("Beta Distribution: ROI MD -> %s", behavior_var)) +
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

    # 保存结果供外部合并
    assign(sprintf("roi_%s_%s", behavior_var, scale_type), roi_df, envir = .GlobalEnv)
    assign(sprintf("m_global_%s", behavior_var), m_global, envir = .GlobalEnv)
  }
}

# ---- 2. 运行两个行为变量 ----
analyze_behavior("regulate_effect", df, roi_cols)
analyze_behavior("generation_effect", df, roi_cols)

# ---- 3. 保存xlsx ----
dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

prepare_out <- function(dat) {
  dat %>%
    select(roi, index, system, beta, se, t_value, p_value, p_fdr, significant, direction) %>%
    mutate(across(c(beta, se, t_value), ~ round(., 6)),
           across(c(p_value, p_fdr), ~ signif(., 4)))
}

# 各结果列表
sheet_list <- list()

for (behav in c("regulate_effect", "generation_effect")) {
  for (scale_type in c("raw", "std")) {
    obj_name <- sprintf("roi_%s_%s", behav, scale_type)
    if (exists(obj_name)) {
      dat <- get(obj_name)
      out <- prepare_out(dat) %>%
        rename(beta = beta, se = se, t = t_value, p = p_value)
      suffix <- if (scale_type == "raw") "raw" else "std"
      sheet_list[[sprintf("%s_%s", behav, suffix)]] <- out %>% arrange(p)
    }
  }
}

# 全局模型
global_list <- list()
for (behav in c("regulate_effect", "generation_effect")) {
  obj_name <- sprintf("m_global_%s", behav)
  if (exists(obj_name)) {
    tbl <- tidy(get(obj_name), conf.int = TRUE) %>%
      mutate(across(where(is.numeric), ~ round(., 6)))
    global_list[[behav]] <- tbl
  }
}

write_xlsx(c(global_list, sheet_list), "_tables/behavior_effect_results.xlsx")
cat("已保存: _tables/behavior_effect_results.xlsx\n")

cat("\n========== 全部完成 ==========\n")
