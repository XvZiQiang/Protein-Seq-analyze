#查看当前有效库路径（新库在前，原库在后）
.libPaths()
#设置新的包路径，同时保留原来的路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")  # 原来的库路径
.libPaths(c(new_lib, old_libs))  # 以后先搜索 new_lib，再搜索原来的库
#检查2个R包路径是否存在
.libPaths()
#以后所有R安装包时会默认安装到新的包路径，以后加载包时，R 会从 new_lib 和原来的库里查找


library(tidyverse)
library(stringr)

# 读取数据
# 读取你上一步生成的 Table S3
df <- read.csv("C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/Maxquant_Table_S3_Result.csv", stringsAsFactors = FALSE)






#第一部分：WT_NT******WT_ST

# 提取定量列 (Ratio)
# 假设你的列名是这些 (请根据实际 CSV 文件头确认)
wt_nsp_cols <- c("WT_NT_1", "WT_NT_2", "WT_NT_3")
wt_sp_cols  <- c("WT_ST_1", "WT_ST_2", "WT_ST_3")


# 检查一下数据是否包含 NA，T-test 需要去除 NA
# 这一步只保留在两组中都有值的行
valid_rows <- rowSums(is.na(df[, c(wt_nsp_cols, wt_sp_cols)])) == 0
df_process <- df[valid_rows, ]




#df_process表里提取前6列数据（不包括后面数据）
df_process_front_info <- df_process[, 1:6]


#df_process表里提取WT_NT和WT_ST处理数据（一共6列）
raw_ratios <- df_process[, c(wt_nsp_cols, wt_sp_cols)]


#  数据预处理：Log2 转换
# T-test 要求数据符合正态分布，Ratio 数据必须先 Log2
# 注意：你的数据已经是 Ratio (1.0 附近)，log2(1) = 0
# 也就是 Log2 后，大于0是上调，小于0是下调
log_data <- log2(raw_ratios)





# 定义 T-test 函数

# 我们比较 SP (处理) vs NSP (对照)
# Fold Change = SP / NSP
# Log2FC = Mean(Log2_SP) - Mean(Log2_NSP)
calculate_stats <- function(row_vals) {
  # 提取两组数值
  nsp_vals <- row_vals[wt_nsp_cols]
  sp_vals  <- row_vals[wt_sp_cols]
  
  # 执行标准 T-test (var.equal = TRUE 符合文章描述的 standard deviation)
  # 如果报错(比如数据完全一样)，返回 NA
  t_res <- tryCatch(t.test(sp_vals, nsp_vals, var.equal = TRUE), error = function(e) NULL)
  
  if (is.null(t_res)) {
    return(c(Log2FC = NA, PValue = NA))
  } else {
    # 计算 Log2FC
    log_fc <- mean(sp_vals) - mean(nsp_vals)
    return(c(Log2FC = log_fc, PValue = t_res$p.value))
  }
}




# 批量计算，T检验
#对 log_data 的每一行做一次 T 检验
stats_res <- apply(log_data, 1, calculate_stats)
#翻转矩阵，让“基因”为行、“统计量”为列
stats_res <- as.data.frame(t(stats_res))



# df_process_front_info和log_data和stats_res合并回主表result_df
result_df <- cbind(df_process_front_info, log_data, stats_res)





# 计算 FDR 和 标记差异

# 计算 Fold Change (Log2FC从Log2还原回正常倍数，方便看)
result_df$FoldChange <- 2 ^ result_df$Log2FC

# 利用PValue计算 FDR (Benjamini-Hochberg算法)
result_df$FDR <- p.adjust(result_df$PValue, method = "BH")









#很多文章，写是写 "FDR < 0.01"，但实际上可能用的是 Perseus 软件里的
#Permutation-based FDR（这种算法在处理 MaxQuant 数据时比 R 的 BH 算法更容易得到结果）

#所以检查一下FDR值

# 1. 看看如果不卡 FDR，只卡 P < 0.05，你能得到多少个？
raw_sig_count <- sum(result_df$PValue < 0.05, na.rm = TRUE)
print(paste("P-value < 0.05 的蛋白数量:", raw_sig_count))

# 2. 看看如果不卡 FDR，只卡 P < 0.01，你能得到多少个？
raw_sig_count <- sum(result_df$PValue < 0.01, na.rm = TRUE)
print(paste("P-value < 0.05 的蛋白数量:", raw_sig_count))

# 3. 看看数据里算出来的最小 FDR 是多少？
# 如果这个数大于 0.01，说明你之前的代码无论怎么跑都是 0
min_fdr <- min(result_df$FDR, na.rm = TRUE)
print(paste("全数据中最小的 FDR 值:", round(min_fdr, 5)))

#发现无法使用FDR，因为最小都是0.76543，无法使用FDR<0.5的限制条件
#说明实际上可能用的是 Perseus 软件里的FDR
#文章图注用的P < 0.05，FC>1.2







# 标记差异 (严格复现文章标准)
# 设定阈值 
p_cut_strict <- 0.05       # P-value < 0.05
logfc_cut <- log2(1.2)     # FC > 1.2


#进行分类
# 使用 case_when 创建 Regulation 列
result_df$Regulation <- case_when(
  # 上调
  result_df$PValue <= p_cut_strict & result_df$Log2FC > logfc_cut ~ "Up",
  
  # 下调
  result_df$PValue <= p_cut_strict & result_df$Log2FC < -logfc_cut ~ "Down",
  
  # 稳定
  TRUE ~ "Stable"
)




# 输出统计结果
counts <- table(result_df$Regulation)
print(counts)




#result_df简化ID

# 定义一个处理单行 Protein.IDs 的函数
process_ids <- function(id_string) {
  # 1. 按分号分割成多个部分
  parts <- str_split(id_string, ";")[[1]]
  
  # 2. 提取 ID (逻辑：找到两个竖线 | 中间的内容)
  # 正则表达式 "(?<=\\|)[^|]+" 意思是：查找 | 之后，非 | 的字符
  extracted_ids <- str_extract(parts, "(?<=\\|)[^|]+")
  
  # 3. 容错处理：如果有部分没有竖线（匹配结果为 NA），则保留原样
  # 这种情况很少见，但为了防止数据丢失，加上这一步
  extracted_ids[is.na(extracted_ids)] <- parts[is.na(extracted_ids)]
  
  # 4. 去重 (unique)
  # 这一步会自动把 tr|A0A... 和 REV__tr|A0A... 提取出的相同 ID 合并为一个
  unique_ids <- unique(extracted_ids)
  
  # 5. 重新合并，用分号隔开
  return(paste(unique_ids, collapse = ";"))
}



#result_df简化ID
# 应用到 DataFrame
# 使用 sapply 批量处理每一行
result_df$ShortID <- sapply(result_df$Protein.IDs, process_ids)

#检查结果
head(result_df[, c("Protein.IDs", "ShortID")])









# 第五步：绘制火山图 (Volcano Plot) - 适配 MaxQuant 结果

library(ggplot2)
library(ggrepel)
library(dplyr) # 确保加载 dplyr

# 1. 准备绘图数据
# 为了画图，我们需要计算 -Log10(PValue)
result_df$MinusLog10P <- -log10(result_df$PValue)

# 2. 统计各组数目 (用于图例显示)
# 你的分组列叫 "Regulation"，值是 "Up", "Down", "Stable"
n_up <- sum(result_df$Regulation == "Up", na.rm = TRUE)
n_down <- sum(result_df$Regulation == "Down", na.rm = TRUE)
n_stable <- sum(result_df$Regulation == "Stable", na.rm = TRUE)

# 创建动态图例标签
label_up <- paste0("Up (", n_up, ")")
label_down <- paste0("Down (", n_down, ")")
label_stable <- paste0("NoDiff (", n_stable, ")")

print(paste("绘图统计：", label_up, label_down))

# 3. 筛选出需要标记名字的数据
# 策略：标记所有 Up 和 Down 的蛋白
# 注意：如果差异蛋白太多（比如超过50个），图会很乱。
# 你可以使用 head() 只取前几个，或者只标记 Log2FC 绝对值最大的几个。
label_data <- subset(result_df, Regulation != "Stable")


# 4. 设置阈值线位置
# P < 0.05
y_line <- -log10(0.05) 
# FC > 1.2 和 FC < 1/1.2 (即 Log2 > 0.263...)
x_line_pos <- log2(1.2)
x_line_neg <- -log2(1.2)

# 5. 绘图
p <- ggplot(result_df, aes(x = Log2FC, y = MinusLog10P, color = Regulation)) +
  # (1) 画散点
  geom_point(alpha = 0.6, size = 1.5) +
  
  # (2) 设置颜色和图例标签
  scale_color_manual(
    values = c("Up" = "red", "Down" = "blue", "Stable" = "grey"),
    labels = c("Up" = label_up, "Down" = label_down, "Stable" = label_stable),
    breaks = c("Up", "Down", "Stable") # 强制图例顺序
  ) +
  
  # (3) 添加阈值辅助线
  geom_vline(xintercept = c(x_line_neg, x_line_pos), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = y_line, 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # (4) 标记蛋白名称
  # 使用你清洗好的 ShortID 列
  geom_text_repel(
    data = label_data,
    aes(label = ShortID),       # 这里也可以改成 Gene.Name
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,          # 避免标签太密集
    show.legend = FALSE,
    segment.color = "grey50"    # 连线颜色
  ) +
  
  # (5) 设置主题和标签
  labs(title = "Volcano Plot: WT-ST vs WT-NT",
       subtitle = "Threshold: P < 0.05 & FC > 1.2",
       x = expression(log[2] ~ "(Fold Change)"),
       y = expression(-log[10] ~ "(P-value)"),
       color = "Regulation") +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# 6. 显示并保存
print(p)

# 保存图片 (路径请根据实际情况修改)
save_path <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/Maxquant_Volcano_Plot_WT.png"
ggsave(save_path, plot = p, width = 8, height = 6, dpi = 300)
