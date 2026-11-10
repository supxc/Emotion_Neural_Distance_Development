# ============================================================================
# 08_scatter_roi_behaviors.R
# 单独绘制指定 ROI 马氏距离与情绪行为变量的偏回归散点图
# (先回归掉 gender + age + mri_equip + EMO_mean_FD, 再画残差散点)
#
# 图1: ROI206 ~ generation_effect   (颜色 #362350)
# 图2: ROI208 ~ generation_effect   (颜色 #457AC6)
# 图3: ROI206 ~ regulate_effect     (颜色 #362350)
# ============================================================================
rm(list = ls())
gc()

library(dplyr)
library(ggplot2)

# ---- 0. 读取数据 ----
df <- read.csv("_data/simulated_mahalanobis_data.csv", stringsAsFactors = FALSE)
df$mri_equip <- factor(df$mri_equip)
df$gender    <- factor(df$gender)

# 确认列存在
stopifnot(all(c("ROI206", "ROI208", "generation_effect", "regulate_effect",
                "gender", "age", "mri_equip", "EMO_mean_FD") %in% colnames(df)))

cat("数据维度:", nrow(df), "x", ncol(df), "\n")

# ---- 1. 偏回归绘图函数 ----
# 把 y 和 x 都先对协变量回归取残差, 再画残差-残差散点
# (残差回归的斜率 = 完整模型中该预测子的 beta)
make_partial_scatter <- function(data, x_var, y_var, color_hex, out_file,
                                 x_lab, y_lab) {

  cov_fml <- "~ gender + mri_equip + EMO_mean_FD + age"

  # x 对协变量回归取残差
  d_tmp <- data.frame(
    y = data[[y_var]],
    x = data[[x_var]],
    gender    = data$gender,
    mri_equip = data$mri_equip,
    EMO_mean_FD = data$EMO_mean_FD,
    age       = data$age
  )
  d_tmp <- d_tmp[complete.cases(d_tmp), ]

  resid_x <- residuals(lm(as.formula(paste0("x", cov_fml)), data = d_tmp))
  resid_y <- residuals(lm(as.formula(paste0("y", cov_fml)), data = d_tmp))

  d_plot <- data.frame(x = resid_x, y = resid_y)

  # 残差回归的 beta 与完整模型一致, 用它做标注
  fit <- lm(y ~ x, data = d_plot)
  sm  <- summary(fit)$coefficients["x", ]
  r   <- cor(d_plot$x, d_plot$y)
  anno <- sprintf("beta = %.4f, p = %.2e", sm["Estimate"], sm["Pr(>|t|)"])

  p <- ggplot(d_plot, aes(x = x, y = y)) +
    geom_point(alpha = 0.35, size = 1, color = color_hex) +
    geom_smooth(method = "lm", se = TRUE, color = color_hex, linewidth = 1.3,
                fill = color_hex, alpha = 0.15) +
    annotate("text", x = -Inf, y = Inf, label = anno,
             hjust = -0.1, vjust = 1.5, size = 4.2, fontface = "italic") +
    labs(x = paste0(x_lab, " (residualized)"),
         y = paste0(y_lab, " (residualized)")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank())

  ggsave(out_file, p, width = 7, height = 6, dpi = 800)
  cat("已保存:", out_file, "\n")
}

# ---- 2. 生成三张偏回归散点图 ----
dir.create("_plots", showWarnings = FALSE, recursive = TRUE)

# 图1: ROI206 vs generation_effect
make_partial_scatter(
  data    = df,
  x_var   = "generation_effect",
  y_var   = "ROI206",
  color_hex = "#362350",
  out_file = "_plots/scatter_ROI206_generation_effect.jpg",
  x_lab   = "Generation Effect",
  y_lab   = "ROI206 SMD"
)

# 图2: ROI208 vs generation_effect
make_partial_scatter(
  data    = df,
  x_var   = "generation_effect",
  y_var   = "ROI208",
  color_hex = "#457AC6",
  out_file = "_plots/scatter_ROI208_generation_effect.jpg",
  x_lab   = "Generation Effect",
  y_lab   = "ROI208 SMD"
)

# 图3: ROI206 vs regulate_effect
make_partial_scatter(
  data    = df,
  x_var   = "regulate_effect",
  y_var   = "ROI206",
  color_hex = "#362350",
  out_file = "_plots/scatter_ROI206_regulate_effect.jpg",
  x_lab   = "Regulation Effect",
  y_lab   = "ROI206 SMD"
)

cat("\nDone.\n")
