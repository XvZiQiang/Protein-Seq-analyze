# 查看当前有效库路径（新库在前，原库在后）
.libPaths()
# 设置新的包路径，同时保留原来的路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")  # 原来的库路径
.libPaths(c(new_lib, old_libs))  # 以后先搜索 new_lib，再搜索原来的库


# 检查2个R包路径是否存在
.libPaths()
#以后所有R安装包时会默认安装到新的包路径，以后加载包时，R 会从 new_lib 和原来的库里查找





# 0. 加载必要的包
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(ggrepel)

# 1. 参数设置 (请在此处修改文件路径和列名)

# 1.1 文件路径 (注意路径中的斜杠方向 / )
file_path <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/proteome_discover3.2/P20200802342_F.xlsx"
save_csv_path <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/proteome_discover3.2/WT_SP_vs_WT_NT_Diff_Proteins.csv"
save_img_path <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/proteome_discover3.2/WT_SP_vs_WT_NT_Volcano.png"

# 1.2 定义列名 (必须与 Excel 表头完全一致)
# 对照组 (Control)
cols_nt <- c("WT_NT_1", "WT_NT_2", "WT_NT_3")
# 处理组 (Treatment)
cols_st <- c("WT_ST_1", "WT_ST_2", "WT_ST_3")

# 1.3 统计阈值 (严格参考文章)
fc_cutoff <- 1.2        # Fold Change > 1.2 或 < 1/1.2 (0.833)
p_cutoff  <- 0.05       # P-value < 0.05

# 2. 读取数据与预处理
raw_data <- read_excel(file_path) # 默认读取第一个 Sheet

# 2.1 提取简短 ID (Accession)为ShortID
raw_data$ShortID <- raw_data$Accession


# 2.2 过滤无效数据
# 确保用于计算的列都是数值型，允许3个重复有1个缺失值，生成calc_data
# 定义一个辅助函数计算非NA数量
count_valid <- function(x) sum(!is.na(x))

calc_data <- raw_data %>%
  rowwise() %>%
  filter(count_valid(c_across(all_of(cols_nt))) >= 2 & 
           count_valid(c_across(all_of(cols_st))) >= 2) %>%
  ungroup()

message(paste("放宽过滤后的行数:", nrow(calc_data)))


# 2.3 数据标准化 (复现文章算法)
# 提取纯数值矩阵
tmt_matrix <- as.matrix(calc_data[, c(cols_nt, cols_st)])


# --- Step A: 列校正 (消除上样误差) ---
# 计算每列中位数
col_medians <- apply(tmt_matrix, 2, median, na.rm = TRUE)
# 计算参考值 (所有中位数的平均)
ref_value <- mean(col_medians)
# 执行校正
norm_factors <- ref_value / col_medians
sln_matrix <- sweep(tmt_matrix, 2, norm_factors, "*")

# --- Step B: 行归一化 (获得你想要的 1.几 的数值) ---
# 计算每一行(每个蛋白)在所有样品中的平均强度
row_means <- rowMeans(sln_matrix, na.rm = TRUE)

# 让每个数值除以该行的平均值
# 结果：如果一个样品是平均水平，它的值就是 1.0；如果是平均的1.2倍，就是 1.2
ratio_matrix <- sln_matrix / row_means

# 检查一下：看看前几行是不是变成了 1.0 左右的数
print(head(ratio_matrix, 3))

# 将这些 "1.x" 的数值保存回数据框，用于后续 T-test
# 注意：我们将原始的绝对强度替换为了相对比值
calc_data[, c(cols_nt, cols_st)] <- as.data.frame(ratio_matrix)




# 3. 统计计算 (均值, FC, T-test)
# 初始化结果列
calc_data$Mean_NT <- NA
calc_data$Mean_ST <- NA
calc_data$FC      <- NA
calc_data$Log2FC  <- NA
calc_data$P_value <- NA

# 循环计算 (稳健且易于调试)
for(i in 1:nrow(calc_data)) {
  # 提取数值
  vals_nt <- as.numeric(calc_data[i, cols_nt])
  vals_st <- as.numeric(calc_data[i, cols_st])
  
  # 1. 计算均值
  mean_nt <- mean(vals_nt, na.rm = TRUE)
  mean_st <- mean(vals_st, na.rm = TRUE)
  
  # 2. 计算 FC (处理/对照)
  # 防止分母为0
  if(mean_nt == 0) {
    fc <- NA
  } else {
    fc <- mean_st / mean_nt
  }
  
  # 3. T 检验 (双尾，等方差，符合文章 "using standard deviation")
  # tryCatch 防止因为数据完全一致(方差为0)导致的报错
  log_nt <- log2(vals_nt)
  log_st <- log2(vals_st)
  
  t_res <- tryCatch({
    t.test(log_nt, log_st, var.equal = TRUE)
  }, error = function(e) return(NULL))
  
  # 写入结果
  calc_data$Mean_NT[i] <- mean_nt
  calc_data$Mean_ST[i] <- mean_st
  calc_data$FC[i]      <- fc
  calc_data$Log2FC[i]  <- ifelse(is.na(fc) | fc <= 0, NA, log2(fc))
  calc_data$P_value[i] <- if(!is.null(t_res)) t_res$p.value else NA
}

# 4. 分组与打标签
# 去除计算失败的行
final_df <- calc_data %>% filter(!is.na(P_value) & !is.na(Log2FC))

# 分类: Up, Down, NoDiff
final_df <- final_df %>%
  mutate(
    Group = case_when(
      P_value < p_cutoff & FC > fc_cutoff ~ "Up",       # 显著上调
      P_value < p_cutoff & FC < (1/fc_cutoff) ~ "Down", # 显著下调
      TRUE ~ "NoDiff"                                   # 不显著
    ),
    MinusLog10P = -log10(P_value)
  )

# 统计数量
n_up <- sum(final_df$Group == "Up")
n_down <- sum(final_df$Group == "Down")
n_nodiff <- sum(final_df$Group == "NoDiff")

message(paste0("统计完成: 上调 ", n_up, " 个, 下调 ", n_down, " 个"))

# 生成图例标签
label_up <- paste0("Up (", n_up, ")")
label_down <- paste0("Down (", n_down, ")")
label_nodiff <- paste0("NoDiff (", n_nodiff, ")")

# 提取用于在图中显示名字的子集 (仅显著差异蛋白)
label_data <- final_df %>% filter(Group != "NoDiff")

# 5. 绘制火山图
p <- ggplot(final_df, aes(x = Log2FC, y = MinusLog10P, color = Group)) +
  # 1. 散点
  geom_point(alpha = 0.6, size = 2) +
  
  # 2. 颜色与图例
  scale_color_manual(
    values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NoDiff" = "grey70"), # 红蓝灰配色
    labels = c("Up" = label_up, "Down" = label_down, "NoDiff" = label_nodiff)
  ) +
  
  # 3. 辅助阈值线
  geom_vline(xintercept = c(log2(fc_cutoff), log2(1/fc_cutoff)), 
             linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(p_cutoff), 
             linetype = "dashed", color = "black", linewidth = 0.4) +
  
  # 4. 标注蛋白 ID (使用 ggrepel 避免重叠)
  geom_text_repel(
    data = label_data,
    aes(label = ShortID),
    size = 3,
    box.padding = 0.5,
    max.overlaps = 20, # 如果标签太多导致报错，可调大此数值或设为 Inf
    show.legend = FALSE
  ) +
  
  # 5. 坐标轴与标题
  labs(
    title = "Volcano Plot: WT-SP vs WT-NSP",
    x = expression(log[2] ~ "(Fold Change)"),
    y = expression(-log[10] ~ "(P-value)"),
    color = "Expression"
  ) +
  
  # 6. 主题美化
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.text = element_text(size = 11)
  )

# 显示图片
print(p)

# 6. 保存结果

# 保存差异蛋白列表
write.csv(label_data, save_csv_path, row.names = FALSE)

# 保存图片 (可选，取消注释以保存)
ggsave(save_img_path, plot = p, width = 8, height = 6, dpi = 300)
