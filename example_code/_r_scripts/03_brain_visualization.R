# ============================================================================
# 03_brain_visualization.R
# 将显著ROI的标准化beta可视化到脑图上
# 皮层 (ROI1-400): ggseg Schaefer7 400
# 皮层下 (ROI401-432): bar plot
# 直接从 age_effect_results.xlsx 读取已有结果
# ============================================================================
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)

# ---- 0. 读取已有结果 ----
# roi_std_results sheet列名: roi, index, system, beta_std, se_std, t_std,
#                             p_std, p_fdr_std, significant, direction
roi_std <- read_excel("_tables/age_effect_results.xlsx", sheet = "roi_std_results")

# 读取ROI标签
labels <- read.csv("_outputs/400P7N_Subcortexlabel.csv", header = FALSE,
                   col.names = c("index", "label"), stringsAsFactors = FALSE)
labels$roi <- sprintf("ROI%d", labels$index)

# 合并标签
roi_std <- roi_std %>% left_join(labels, by = c("index","roi"))
roi_std$significant_fdr <- roi_std$p_fdr_std < 0.05

cat(sprintf("ROI总数: %d, FDR < 0.05: %d\n",
            nrow(roi_std), sum(roi_std$significant_fdr)))

# ---- 2. 分离皮层和皮层下 ----
sig_all <- roi_std %>%
  filter(significant_fdr) %>%
  arrange(p_fdr_std)

sig_cortical <- sig_all %>% filter(index <= 400)
sig_subcort  <- sig_all %>% filter(index > 400)

cat(sprintf("皮层显著: %d, 皮层下显著: %d\n",
            nrow(sig_cortical), nrow(sig_subcort)))

# ---- 3. 保存显著ROI列表 ----
dir.create("_tables", showWarnings = FALSE, recursive = TRUE)

sig_out <- sig_all %>%
  select(roi, index, label, beta_std, se_std, t_std, p_std, p_fdr_std) %>%
  rename(beta = beta_std, se = se_std, t = t_std, p = p_std, p_fdr = p_fdr_std) %>%
  mutate(across(c(beta, se, t), ~ round(., 6)),
         across(c(p, p_fdr), ~ signif(., 4)))

write_xlsx(
  list(significant_std = sig_out,
       cortical        = sig_cortical %>% select(roi, label, beta_std, p_fdr_std),
       subcortical     = sig_subcort  %>% select(roi, label, beta_std, p_fdr_std)),
  "_tables/age_sig_ROI_brain_labels.xlsx"
)
cat("已保存: _tables/age_sig_ROI_brain_labels.xlsx\n")

# ---- 4. 皮层7网络统计图 ----
# 从标签提取网络名: "7Networks_LH_Vis_1" -> "Vis"
sig_cortical <- sig_cortical %>%
  mutate(
    network = gsub("^7Networks_(LH|RH)_([A-Za-z]+)_.*", "\\2", label),
    direction_label = ifelse(beta_std > 0, "Positive", "Negative")
  )

# 7网络标准缩写映射
network_abbrev <- c(
  "Vis" = "Visual", "SomMot" = "Somatomotor",
  "DorsAttn" = "Dorsal\nAttention", "SalVentAttn" = "Ventral\nAttention",
  "Limbic" = "Limbic", "Cont" = "Frontoparietal", "Default" = "Default"
)
sig_cortical$network_full <- network_abbrev[sig_cortical$network]

# 统计各网络正负ROI数
net_count <- sig_cortical %>%
  count(network_full, direction_label, .drop = FALSE) %>%
  tidyr::complete(network_full, direction_label, fill = list(n = 0)) %>%
  filter(!is.na(network_full))

dir.create("_plots", showWarnings = FALSE, recursive = TRUE)

p_net <- ggplot(net_count, aes(x = network_full, y = n, fill = direction_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Positive" = "#D7263D", "Negative" = "#1B5299"),
                    name = expression(beta[std])) +
  labs(x = NULL, y = "Number of Significant ROIs",
       title = "Significant Cortical ROIs by Yeo 7-Network",
       subtitle = sprintf("FDR < 0.05, %d cortical ROIs", nrow(sig_cortical))) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position = c(0.92, 0.88),
        legend.background = element_rect(fill = "white", color = "grey80"))

ggsave("_plots/cortical_network_sig_count.jpg", p_net,
       width = 10, height = 5, dpi = 800)
cat("已保存: _plots/cortical_network_sig_count.jpg\n")

# ---- 5. 皮层下 bar plot ----
if (nrow(sig_subcort) > 0) {
  sig_subcort <- sig_subcort %>%
    mutate(label_clean = gsub("-lh$", " (L)", gsub("-rh$", " (R)", label)),
           direction_label = ifelse(beta_std > 0, "Positive", "Negative"))

  ylim_max <- max(abs(sig_subcort$beta_std), na.rm = TRUE) * 1.2

  p_sub <- ggplot(sig_subcort,
                  aes(x = reorder(label_clean, beta_std), y = beta_std,
                      fill = direction_label)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("Positive" = "#D7263D", "Negative" = "#1B5299"),
                      guide = "none") +
    ylim(-ylim_max, ylim_max) +
    coord_flip() +
    labs(x = NULL, y = expression(beta[std]),
         title = "Significant Subcortical ROIs (FDR < 0.05)",
         subtitle = sprintf("%d subcortical ROIs", nrow(sig_subcort))) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave("_plots/subcortical_sig_bar.jpg", p_sub,
         width = 8, height = max(4, nrow(sig_subcort) * 0.45), dpi = 800)
  cat("已保存: _plots/subcortical_sig_bar.jpg\n")
} else {
  cat("无显著皮层下ROI, 跳过bar plot。\n")
}

cat("\nDone.\n")
