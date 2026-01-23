# 查看当前有效库路径（新库在前，原库在后）
.libPaths()
# 设置新的包路径，同时保留原来的路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")  # 原来的库路径
.libPaths(c(new_lib, old_libs))  # 以后先搜索 new_lib，再搜索原来的库


# 检查2个R包路径是否存在
.libPaths()
#以后所有R安装包时会默认安装到新的包路径，以后加载包时，R 会从 new_lib 和原来的库里查找

library(tidyverse)

# 读取数据 (关键修正)
# ⚠️ 注意：Table S2 是肽段表，必须读取 'peptides.txt'
# ⚠️ 注意：我把变量名统一改为 pep_raw，以防后面报错
pep_raw <- read.delim("C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/peptides.txt", stringsAsFactors = FALSE)


# 数据清洗
pep_clean <- pep_raw %>%
  filter(Potential.contaminant != "+") %>% # 去除污染物
  filter(Reverse != "+")                   # 去除反向库

# 提取 Table S2 所需的列
# 提取基础信息
base_info <- pep_clean %>%
  select(
    Sequence,           # 序列
    PEP,                # 后验误差概率
    Score,              # Andromeda Score
    Mass,               # 质量
    Proteins            # 对应蛋白ID
  )

##################################
#每个人maxquant里experiment命名不同，我的命名是first run
# 1. 打印所有的列名，看看 TMT 强度到底叫什么
all_cols <- colnames(pep_clean)
print(paste("总列数:", length(all_cols)))


# 寻找所有包含 "Reporter" 字眼的列
possible_tmt_cols <- grep("Reporter", all_cols, value = TRUE)
print(possible_tmt_cols)





# 精确提取 TMT 列 (解决顺序和后缀问题)
# 注意：MaxQuant 输出的是 "Reporter intensity corrected 1, 2..."


# 此实验是12个样品，4个处理*3个重复，利用TMT 16plex标记，但16个没有全用完，只用了前12个
#只用了前12个（如果不确定，等会单独，对比实验结果看一下相关性）
# 我们手动构建这 12 个列名
#每个人maxquant里experiment命名不同，我的命名是first run
target_cols <- paste0("Reporter.intensity.corrected.", 1:12, ".first.run")

# 检查一下这些列是否存在于你的数据中
missing_cols <- setdiff(target_cols, colnames(pep_clean))
if(length(missing_cols) > 0) {
  stop(paste("错误：以下列名未找到", paste(missing_cols, collapse=", ")))
}


# 提取这 12 列数据
tmt_data <- pep_clean[, target_cols]




# 将 0 转换为 NA (以便后续统计)
tmt_data[tmt_data == 0] <- NA


# 重命名列 (请根据你的实验设计修改这里!)
# 假设你的12个通道顺序如下：
colnames(tmt_data) <- c(
  "WT_NT_1", "WT_NT_2", "WT_NT_3",
  "WT_ST_1", "WT_ST_2", "WT_ST_3",
  "ANK_NT_1", "ANK_NT-2", "ANK_NT_3",
  "ANK_ST_1", "ANK_ST_2", "ANK_ST_3"
)


# 4. 合并base_info和_data，并保存结果
final_table_s2 <- cbind(base_info, tmt_data)

# 预览
head(final_table_s2)

# 保存表格
write.csv(final_table_s2, "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/Maxquant_Table_S2_Result.csv", row.names = FALSE)





#检验一下maxquant输出结果跟文章的结果差别

#EQFPGANDFGSEVIPGATSTGMR
# 输入你的数据 (MaxQuant)
mq_vals <- c(230110, 239160, 280010, 246390, 320870, 317020, 296000, 309480, 331290, 385950, 325950, 322180)

# 输入文章的数据  (Proteome Discoverer)
pd_vals <- c(6611127, 6975094, 8192432, 7411078, 9471128, 9128349, 8583192, 8649082, 9608482, 11094547, 9386626, 9140692)
# 输入文章的数据  (Proteome Discoverer)
pd_vals_2 <- c(2333880.531, 2592940.813, 3242404, 2544598.125, 3572238.875, 3366642.5, 3285263.25, 3218294.875, 3656611.563, 3987416.563, 3374531.75, 3343558.188
)

# 计算相关系数
cor_val <- cor(mq_vals, pd_vals)
print(paste("12个点的相关系数 R =", round(cor_val, 4)))

# 归一化画图
plot(mq_vals/mq_vals[1], type="b", col="red", lwd=2, ylim=c(0.8, 1.8), 
     ylab="Ratio to Sample 1", xlab="Sample Index", main="MaxQuant (Red) vs Paper (Blue) Trend")
lines(pd_vals/pd_vals[1], type="b", col="blue", lwd=2)
legend("topleft", legend=c("Your Result", "Paper Result"), col=c("red", "blue"), lwd=2)






# 计算相关系数
cor_val <- cor(mq_vals, pd_vals_2)
print(paste("12个点的相关系数 R =", round(cor_val, 4)))

# 归一化画图
plot(mq_vals/mq_vals[1], type="b", col="red", lwd=2, ylim=c(0.8, 1.8), 
     ylab="Ratio to Sample 1", xlab="Sample Index", main="MaxQuant (Red) vs Paper (Blue) Trend")
lines(pd_vals_2/pd_vals_2[1], type="b", col="blue", lwd=2)
legend("topleft", legend=c("Your Result", "Paper Result"), col=c("red", "blue"), lwd=2)



#SSAEEDDELSKLPAAQR 

# 你的 MaxQuant 数据 (Raw Intensity)
mq2_vals <- c(15143, 15611, 22311, 16882, 23317, 20583, 19521, 17898, 20947, 22973, 21074, 19367)

# 文章的 Proteome Discoverer 数据 (Raw Intensity)
pd2_vals <- c(300036.53, 268266.19, 372751.78, 309528.06, 405145.88, 341841.00, 
              342423.31, 306035.84, 394627.03, 455559.97, 415051.38, 314415.13)

# 计算相关系数
cor_val <- cor(mq2_vals, pd2_vals)
print(paste("12个点的相关系数 R =", round(cor_val, 4)))

# 归一化画图
plot(mq2_vals/mq2_vals[1], type="b", col="red", lwd=2, ylim=c(0.8, 1.8), 
     ylab="Ratio to Sample 1", xlab="Sample Index", main="MaxQuant (Red) vs Paper (Blue) Trend")
lines(pd2_vals/pd2_vals[1], type="b", col="blue", lwd=2)
legend("topleft", legend=c("Your Result", "Paper Result"), col=c("red", "blue"), lwd=2)



