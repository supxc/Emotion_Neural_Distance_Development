# ============================================================================
# 01_simulate_data.R
# 模拟儿童情绪调节fMRI马氏距离数据
# 
# 实验设计:
#   - 1464名6-19岁儿童, 三种条件: 中性面孔 / 负性面孔 / 带调节策略看负性面孔
#   - 情绪产生 (generation)  = 负性面孔神经活动 - 中性面孔神经活动
#   - 情绪调节 (regulation)  = 带调节策略看负性面孔神经活动 - 负性面孔神经活动
#   - 脑模板: Schaefer400 (400 ROIs) + 皮层下结构 (32 ROIs) = 432 ROIs
#   - 每个被试每个ROI计算两个条件间的马氏距离
#
# 数据列:
#   [1] subject_id      被试编号
#   [2] mri_equip       MRI设备型号
#   [3] age             年龄 (小数)
#   [4] gender          性别 (1=男, 2=女)
#   [5] EMO_mean_FD     平均头动
#   [6] regulate_effect 情绪调节效果 (负性打分 - 调节后打分, 越高越好)
#   [7] generation_effect 负性情绪产生 (负性打分 - 中性打分)
#   [8-407] ROI1 ~ ROI400  Schaefer400各ROI的马氏距离
#   [408-439] ROI401 ~ ROI432  皮层下各ROI的马氏距离
# ============================================================================

library(MASS)       # mvrnorm
library(tidyverse)

set.seed(42)

# ---- 1. 基础参数 ----
n_subjects  <- 1464
n_schaefer  <- 400
n_subcort   <- 32
n_roi       <- n_schaefer + n_subcort

# ---- 2. 生成人口学和行为变量 ----
sim <- tibble(
  subject_id      = sprintf("sub-%04d", 1:n_subjects),
  mri_equip       = sample(c("Siemens_Prisma", "Siemens_Trio", "GE_Discovery"),
                           n_subjects, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
  age             = runif(n_subjects, min = 6, max = 19),
  gender          = sample(1:2, n_subjects, replace = TRUE)
)

# 头动: 年幼儿童更大, 加入噪声
sim$EMO_mean_FD <- 0.35 - 0.012 * sim$age +
  rnorm(n_subjects, 0, 0.06)
sim$EMO_mean_FD <- pmax(sim$EMO_mean_FD, 0.02)

# 行为变量: generation_effect 随年龄递减(前额叶成熟, 负性情绪产生减弱)
sim$generation_effect <- 3.2 - 0.08 * sim$age +
  rnorm(n_subjects, 0, 1.2)
sim$generation_effect <- pmax(sim$generation_effect, 0.2)  # 正值

# regulate_effect 随年龄递增(调控能力增强), 并与generation_effect轻微正相关
sim$regulate_effect <- 0.5 + 0.12 * sim$age +
  0.3 * sim$generation_effect +
  rnorm(n_subjects, 0, 1.0)
sim$regulate_effect <- pmax(sim$regulate_effect, 0)

# ---- 3. 生成ROI马氏距离 ----

# 3a. 构建空间相关结构 (邻近ROI更相似)
# 使用指数衰减相关矩阵: cor(i,j) = exp(-|i-j| / decay)
build_spatial_cor <- function(n, decay = 20, noise = 0.3) {
  d <- abs(outer(1:n, 1:n, "-"))
  R <- exp(-d / decay)
  # 加入噪声使结构不完全平滑
  R <- (1 - noise) * R + noise * diag(n)
  R
}

# Schaefer 400: 较大衰减(皮层梯度)
R_schaefer <- build_spatial_cor(n_schaefer, decay = 50, noise = 0.2)
# 皮层下 32: 较小衰减
R_subcort  <- build_spatial_cor(n_subcort, decay = 8, noise = 0.15)

# 皮层下与皮层之间的跨系统相关(弱)
cross_cor <- matrix(0.05, nrow = n_subcort, ncol = n_schaefer)

# 拼成完整 432x432 相关矩阵
R_full <- rbind(
  cbind(R_schaefer,           t(cross_cor)),
  cbind(cross_cor,            R_subcort)
)
# 确保对称正定
R_full <- (R_full + t(R_full)) / 2
# 用近邻法修正非正定
eig <- eigen(R_full)
eig$values <- pmax(eig$values, 1e-6)
R_full <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
R_full <- cov2cor(R_full)

# 3b. 定义均值和协方差尺度
# 每个ROI有基础马氏距离均值(不同脑区基础激活不同)
base_md <- c(
  # Schaefer 400: 视觉/默认/注意/额顶网络等, 模拟一个合理的分布
  0.3 + 0.5 * sin(seq(0, 4*pi, length.out = n_schaefer)) +
    rnorm(n_schaefer, 0, 0.1),
  # 皮层下: 杏仁核/海马等情绪相关结构基础值稍高
  0.5 + 0.3 * sin(seq(0, 2*pi, length.out = n_subcort)) +
    rnorm(n_subcort, 0, 0.08)
)
# 各ROI标准差(乘性噪声)
roi_sd <- pmax(0.15, base_md * 0.4 + rnorm(n_roi, 0, 0.05))
Sigma <- diag(roi_sd) %*% R_full %*% diag(roi_sd)

# 3c. 为每个被试生成432维马氏距离
roi_matrix <- mvrnorm(n_subjects, mu = base_md, Sigma = Sigma)
colnames(roi_matrix) <- c(
  paste0("ROI", 1:n_schaefer),
  paste0("ROI", (n_schaefer + 1):n_roi)
)

# ---- 4. 加入年龄效应 ----
# 某些ROI的马氏距离与年龄有关(发展效应), 某些无关
# 选择约40%的ROI有年龄效应(发展敏感的脑区)

# Schaefer中年龄效应ROI: 模拟额顶网络随年龄变化更强
age_effect_rois <- sort(c(
  sample(1:n_schaefer, size = 150),          # 150个皮层ROI
  n_schaefer + sample(1:n_subcort, size = 12) # 12个皮层下ROI
))

for (j in age_effect_rois) {
  # 每个ROI的年龄效应量不同, 正负均有
  beta_j <- rnorm(1, 0, 0.02)
  roi_matrix[, j] <- roi_matrix[, j] +
    beta_j * (sim$age - mean(sim$age))
}

# 确保马氏距离为正值
roi_matrix <- pmax(roi_matrix, 0.01)

# ---- 5. 加入MRI设备效应(批次效应) ----
for (j in 1:n_roi) {
  # 大约一半ROI有微弱设备效应
  if (runif(1) < 0.5) {
    shift <- rnorm(1, 0, 0.03)
    idx2 <- sim$mri_equip == "Siemens_Trio"
    idx3 <- sim$mri_equip == "GE_Discovery"
    roi_matrix[idx2, j] <- roi_matrix[idx2, j] + shift
    roi_matrix[idx3, j] <- roi_matrix[idx3, j] + rnorm(1, 0, 0.04)
  }
}
roi_matrix <- pmax(roi_matrix, 0.01)

# ---- 6. 合并最终数据 ----
sim_data <- cbind(sim, as.data.frame(roi_matrix))

# ---- 7. 保存 ----
out_dir <- "_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



# ---- 8. 输出数据描述 ----
cat("\n========== 数据生成完毕 ==========\n")
cat("样本量:", nrow(sim_data), "\n")
cat("总列数:", ncol(sim_data), "\n")
cat("文件保存至:", file.path(out_dir, "simulated_mahalanobis_data.csv"), "\n\n")

cat("前7列 (人口学+行为):\n")
print(head(sim_data[, 1:7]))
cat("\nROI列 (前5个 + 后5个):\n")
roi_start <- which(colnames(sim_data) == "ROI1")
print(colnames(sim_data)[roi_start:(roi_start+4)])
print(tail(colnames(sim_data), 5))

cat("\n描述统计:\n")
cat("年龄范围:", round(range(sim$age), 2), "\n")
cat("generate_effect:", round(mean(sim$generation_effect), 2),
    "(SD =", round(sd(sim$generation_effect), 2), ")\n")
cat("regulate_effect:", round(mean(sim$regulate_effect), 2),
    "(SD =", round(sd(sim$regulate_effect), 2), ")\n")
cat("FD:", round(mean(sim$EMO_mean_FD), 3),
    "(SD =", round(sd(sim$EMO_mean_FD), 3), ")\n")
cat("mri_equip分布:\n")
print(table(sim$mri_equip))
cat("\nanxiety 缺失数:", sum(is.na(sim_data$anxiety)), "\n")
cat("depression 缺失数:", sum(is.na(sim_data$depression)), "\n")
cat("\n==================================\n")


sim_data <- sim_data %>% 
  select(1:5,7,6,8:439)

# ---- 7.5 生成焦虑和抑郁得分 ----
# 列440: anxiety, 列441: depression
# 约100名被试缺失 (NA)
n_missing <- 100
set.seed(123)

# 焦虑: 与 regulate_effect 负相关 (调节差->焦虑高), 与 generation_effect 正相关
sim_data$anxiety <- 12 +
  1.5 * sim_data$generation_effect -
  2.0 * sim_data$regulate_effect +
  rnorm(nrow(sim_data), 0, 3)
sim_data$anxiety <- pmax(sim_data$anxiety, 0)

# 抑郁: 与焦虑正相关, 与年龄负相关(儿童期年龄大抑郁低)
sim_data$depression <- 10 +
  0.6 * sim_data$anxiety -
  0.4 * sim_data$age +
  rnorm(nrow(sim_data), 0, 3)
sim_data$depression <- pmax(sim_data$depression, 0)

# 随机选100名被试缺失两个变量
missing_idx <- sample(1:nrow(sim_data), n_missing)
sim_data$anxiety[missing_idx]    <- NA
sim_data$depression[missing_idx] <- NA

# 列顺序: anxiety = 440, depression = 441
sim_data <- sim_data %>%
  select(1:439, anxiety, depression)

write.csv(sim_data,
          file.path(out_dir, "simulated_mahalanobis_data.csv"),
          row.names = FALSE)
