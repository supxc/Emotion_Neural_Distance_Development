# ============================================================================
# 05_behavior_sig_visualization.R
# 提取 regulate_effect / generation_effect 模型中的显著ROI
# 流程复刻 03: xlsx + 皮层7网络图 + 皮层下bar plot
# ============================================================================
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)

# ---- 0. 读取ROI标签 ----
labels <- read.csv("_outputs/400P7N_Subcortexlabel.csv", header = FALSE,
                   col.names = c("index", "label"), stringsAsFactors = FALSE)
labels$roi <- sprintf("ROI%d", labels$index)

# 7网络映射
network_abbrev <- c(
  "Vis" = "Visual", "SomMot" = "Somatomotor",
  "DorsAttn" = "Dorsal\nAttention", "SalVentAttn" = "Ventral\nAttention",
  "Limbic" = "Limbic", "Cont" = "Frontoparietal", "Default" = "Default"
)

dir.create("_plots", showWarnings = FALSE, recursive = TRUE)
dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

# ---- 1. 处理函数 ----
process_behavior <- function(sheet_name, prefix, title_label) {

  cat(sprintf("\n========== %s ==========\n", title_label))

  # 读取std结果
  dat <- read_excel("_tables/behavior_effect_results.xlsx", sheet = sheet_name)

  # 合并标签
  dat <- dat %>% left_join(labels, by = c("index","roi"))
  dat$significant_fdr <- dat$p_fdr < 0.05

  cat(sprintf("ROI总数: %d, FDR < 0.05: %d\n",
              nrow(dat), sum(dat$significant_fdr)))

  # 分离皮层/皮层下
  sig_all <- dat %>%
    filter(significant_fdr) %>%
    arrange(p_fdr)

  sig_cortical <- sig_all %>% filter(index <= 400)
  sig_subcort  <- sig_all %>% filter(index > 400)

  cat(sprintf("皮层显著: %d, 皮层下显著: %d\n",
              nrow(sig_cortical), nrow(sig_subcort)))

  # ---- 保存xlsx ----
  sig_out <- sig_all %>%
    select(roi, index, label, beta, se, t, p_fdr) %>%
    mutate(across(c(beta, se, t), ~ round(., 6)),
           across(p_fdr, ~ signif(., 4)))

  fname_xlsx <- sprintf("_tables/%s_sig_ROI_brain_labels.xlsx", prefix)
  write_xlsx(
    list(significant    = sig_out,
       cortical       = sig_cortical %>% select(roi, label, beta, p_fdr),
       subcortical    = sig_subcort  %>% select(roi, label, beta, p_fdr)),
    fname_xlsx
  )
  cat(sprintf("已保存: %s\n", fname_xlsx))

  # ---- 皮层7网络统计图 ----
  sig_cortical <- sig_cortical %>%
    mutate(
      network = gsub("^7Networks_(LH|RH)_([A-Za-z]+)_.*", "\\2", label),
      direction_label = ifelse(beta > 0, "Positive", "Negative")
    )
  sig_cortical$network_full <- network_abbrev[sig_cortical$network]

  net_count <- sig_cortical %>%
    count(network_full, direction_label, .drop = FALSE) %>%
    tidyr::complete(network_full, direction_label, fill = list(n = 0)) %>%
    filter(!is.na(network_full))

  p_net <- ggplot(net_count, aes(x = network_full, y = n, fill = direction_label)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_dodge(width = 0.8), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("Positive" = "#D7263D", "Negative" = "#1B5299"),
                      name = expression(beta[std])) +
    labs(x = NULL, y = "Number of Significant ROIs",
         title = sprintf("Significant Cortical ROIs by Yeo 7-Network (%s)", title_label),
         subtitle = sprintf("FDR < 0.05, %d cortical ROIs", nrow(sig_cortical))) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          legend.position = c(0.92, 0.88),
          legend.background = element_rect(fill = "white", color = "grey80"))

  fname_net <- sprintf("_plots/%s_cortical_network.jpg", prefix)
  ggsave(fname_net, p_net, width = 10, height = 5, dpi = 800)
  cat(sprintf("已保存: %s\n", fname_net))

  # ---- 皮层下 bar plot ----
  if (nrow(sig_subcort) > 0) {
    sig_subcort <- sig_subcort %>%
      mutate(label_clean = gsub("-lh$", " (L)", gsub("-rh$", " (R)", label)),
             direction_label = ifelse(beta > 0, "Positive", "Negative"))

    ylim_max <- max(abs(sig_subcort$beta), na.rm = TRUE) * 1.2

    p_sub <- ggplot(sig_subcort,
                    aes(x = reorder(label_clean, beta), y = beta,
                        fill = direction_label)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      scale_fill_manual(values = c("Positive" = "#D7263D", "Negative" = "#1B5299"),
                        guide = "none") +
      ylim(-ylim_max, ylim_max) +
      coord_flip() +
      labs(x = NULL, y = expression(beta[std]),
           title = sprintf("Significant Subcortical ROIs (%s, FDR < 0.05)", title_label),
           subtitle = sprintf("%d subcortical ROIs", nrow(sig_subcort))) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

    fname_sub <- sprintf("_plots/%s_subcortical_bar.jpg", prefix)
    ggsave(fname_sub, p_sub,
           width = 8, height = max(4, nrow(sig_subcort) * 0.45), dpi = 800)
    cat(sprintf("已保存: %s\n", fname_sub))
  } else {
    cat("无显著皮层下ROI, 跳过bar plot。\n")
  }
}

# ---- 2. 运行 ----
process_behavior("regulate_effect_std", "regulate", "Regulate Effect")
process_behavior("generation_effect_std", "generation", "Generation Effect")

cat("\nDone.\n")
