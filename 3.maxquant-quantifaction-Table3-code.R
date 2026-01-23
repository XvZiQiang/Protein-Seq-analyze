library(tidyverse)


# 读取与清洗数据

# 请修改为你的 proteinGroups.txt 路径
raw_data <- read.delim("C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/proteinGroups.txt", stringsAsFactors = FALSE)

# 过滤掉 污染物(+)、反库(+)、仅由位点鉴定(+)
df_clean <- raw_data %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+")

print(paste("清洗后剩余蛋白数量:", nrow(df_clean)))


# 提取 TMT 强度并重命名


# 自动寻找 "Reporter intensity corrected" 列
# 无论后缀是 .first.run 还是其他，这里都能抓取到
# 注意：MaxQuant 输出的是 "Reporter intensity corrected 1, 2..."

tmt_cols_found <- grep("Reporter\\.intensity\\.corrected\\.\\d+", colnames(df_clean), value = TRUE)

# 确保只取前 12 列 (你的 12 个样本)
tmt_data <- df_clean[, tmt_cols_found[1:12]]



# 将 0 转换为 NA (关键！否则中位数计算会偏低)
tmt_data[tmt_data == 0] <- NA

# 重命名列 (请务必核对顺序)
my_sample_names <- c(
  "WT_NT_1", "WT_NT_2", "WT_NT_3",
  "WT_ST_1", "WT_ST_2", "WT_ST_3",
  "ANK_NT_1", "ANK_NT_2", "ANK_NT_3",
  "ANK_ST_1", "ANK_ST_2", "ANK_ST_3"
)
colnames(tmt_data) <- my_sample_names


# 提取df_clean中元数据 (Gene Name, ID 等)
meta_data <- df_clean %>%
  select(
    Protein.IDs, 
    Fasta.headers,
    Peptides, 
    Unique.peptides, 
    Mol..weight..kDa.
  )


# 尝试从 Fasta header 提取 Gene Name (例如 GN=TaHSP70)
meta_data$Gene.Name <- str_extract(meta_data$Fasta.headers, "GN=[^\\s]+")
meta_data$Gene.Name <- gsub("GN=", "", meta_data$Gene.Name)


# 核心计算：复现文章的标准化逻辑

# --- 步骤 A: 列校正 (消除上样量误差) ---
# 计算每列的中位数
col_medians <- apply(tmt_data, 2, median, na.rm = TRUE)

# 计算所有中位数的平均值 (Reference)
ref_value <- mean(col_medians)

# 计算校正因子
norm_factors <- ref_value / col_medians

# 执行校正 (每一列乘以对应的因子)
data_sln <- sweep(tmt_data, 2, norm_factors, "*")



# --- 步骤 B: 行归一化 (生成 Table S3 风格的 Ratio) ---

# 计算每一行（每个蛋白）的平均强度
row_means <- rowMeans(data_sln, na.rm = TRUE)

# 每个数值除以该行的平均值 -> 得到以 1.0 为中心的比值
data_ratio <- data_sln / row_means


# 5. 合并与保存
final_table_s3 <- cbind(meta_data, data_ratio)

# 预览前几行
print(head(final_table_s3))

# 保存
write.csv(final_table_s3, "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/Maxquant_Table_S3_Result.csv", row.names = FALSE)


# 6. 画图验证
# 画个箱线图看看，中心线应该都在 0 (因为 log2(1)=0)
boxplot(log2(data_ratio), main="Table S3 Data Distribution (Log2 Ratio)", 
        ylab="Log2 (Ratio to Average)", las=2, cex.axis=0.7)
abline(h=0, col="red")





#检验一下
# 三组数据折线图对比代码 (Line Charts)

# 1. 准备数据
# --------------------------------------------------------
# Group 1: Wheatwin-1 (PR4A) - 完美一致
art_1 <- c(0.832, 0.766, 0.839, 0.996, 0.864, 1.002, 0.981, 1.212, 1.166, 1.280, 1.092, 0.938)
my_1  <- c(0.837325966, 0.752424339, 0.83996871, 0.980049025, 0.872321922, 1.0032255, 0.98482781, 1.23030758, 1.164618458, 1.301528854, 1.085475783, 0.947926055)

# Group 2: Uncharacterized (A0A3B6PRC9) - 趋势一致
art_2 <- c(0.95152, 0.90739, 0.95800, 0.92269, 0.92330, 0.96098, 1.03656, 1.07650, 1.08732, 1.06666, 1.05558, 1.06192)
my_2  <- c(0.965009699, 0.883640896, 0.990194593, 0.973326834, 0.965668738, 0.978520922, 1.035342552, 1.010468943, 1.076427282, 1.041487103, 1.056004672, 1.023907765)

# Group 3: Beta-amylase (A0A3B6J0P1) - 压缩效应 (趋势一致但幅度不同)
art_3 <- c(0.95148, 0.90316, 0.89859, 0.98258, 0.96064, 0.79175, 1.17774, 0.92821, 1.17145, 1.01168, 1.15240, 1.07172)
my_3  <- c(1.023543, 0.972046, 0.949516, 0.974406, 0.975091, 0.869259, 1.059976, 1.050022, 1.045526, 1.040293, 1.037501, 1.002820)

# 2. 设置绘图布局 (3行1列)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))

# 3. 定义一个绘图函数 (方便复用)
plot_trend <- function(art_data, my_data, title_text) {
  # 计算 Y 轴范围 (保证两条线都能显示全)
  y_limits <- range(c(art_data, my_data)) * c(0.95, 1.05)
  # 计算相关系数
  r_val <- cor(art_data, my_data)
  
  # 画文章的数据 (蓝色线)
  plot(art_data, type = "b", pch = 19, col = "blue", lwd = 2, 
       ylim = y_limits, xaxt = "n", xlab = "", ylab = "Ratio",
       main = paste0(title_text, " (R = ", round(r_val, 4), ")"))
  
  # 画你的数据 (红色线)
  lines(my_data, type = "b", pch = 19, col = "red", lwd = 2)
  
  # 添加 X 轴坐标 (1-12)
  axis(1, at = 1:12, labels = 1:12)
  
  # 添加辅助线 (Ratio = 1)
  abline(h = 1, col = "gray", lty = 2)
  
  # 添加图例
  legend("topleft", legend = c("Article (PD)", "My Result (MQ)"),
         col = c("blue", "red"), lwd = 2, pch = 19, bty = "n", cex = 0.8)
}

# 4. 开始绘图
plot_trend(art_1, my_1, "1.  Wheatwin-1 (PR4A) ")
plot_trend(art_2, my_2, "2. Uncharacterized (A0A3B6PRC9)")
plot_trend(art_3, my_3, "3. Beta-amylase (A0A3B6J0P1) ")

# 5. 恢复默认布局
par(mfrow = c(1, 1))
