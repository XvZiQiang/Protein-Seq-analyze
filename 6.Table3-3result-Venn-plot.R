# 查看当前有效库路径（新库在前，原库在后）
.libPaths()
# 设置新的包路径，同时保留原来的路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")  # 原来的库路径
.libPaths(c(new_lib, old_libs))  # 以后先搜索 new_lib，再搜索原来的库
.libPaths() 

# 检查2个R包路径是否存在
.libPaths()
#以后所有R安装包时会默认安装到新的包路径，以后加载包时，R 会从 new_lib 和原来的库里查找

install.packages("ggvenn")

# ==============================================================================
# 1. 加载必要的包
# ==============================================================================
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggvenn) # 用于画韦恩图

# ==============================================================================
# 2. 定义文件路径 (请根据你的实际路径修改，注意Windows路径要用双反斜杠 \\ 或反斜杠 /)
# ==============================================================================

# A. 文章结果 (PD 1.4)
file_paper <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/ScienceDirect_files_01Dec2025_15-27-11.695/上下调基因17+4ID.xlsx"

# B. 你的 PD 3.2 结果
file_pd32 <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/proteome_discover3.2/WT_SP_vs_WT_NT_Diff_Proteins.csv"

# C. 你的 MaxQuant 结果
file_mq <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/quantifaction/Maxquant_Table_S3_Result_lable_data.csv"

# ==============================================================================
# 3. 数据读取与清洗
# ==============================================================================

# --- A. 处理文章数据 (Paper) ---
# 假设第一列是 Accession ID
data_paper <- read_excel(file_paper)
# 提取第一列作为 ID 列表，并去除空格
list_paper <- data_paper %>%
  pull(1) %>%             # 提取第一列
  str_trim() %>%          # 去除前后空格
  unique()                # 去重

print(paste("文章 (Paper) 差异蛋白数量:", length(list_paper)))


# --- B. 处理 PD 3.2 数据 ---
# 读取 CSV
data_pd32 <- read.csv(file_pd32)

# 提取 Group 为 Up 或 Down 的 Accession
list_pd32 <- data_pd32 %>%
  filter(Group %in% c("Up", "Down")) %>% # 只选差异蛋白
  pull(Accession) %>%                    # 提取 Accession 列
  str_trim() %>%
  unique()

print(paste("PD 3.2 复现差异蛋白数量:", length(list_pd32)))


# --- C. 处理 MaxQuant 数据 (难点：一行多个ID) ---
data_mq <- read.csv(file_mq)

# 提取 Regulation 为 Up 或 Down 的行，并拆分 ShortID
list_mq <- data_mq %>%
  filter(Regulation %in% c("Up", "Down")) %>% # 筛选差异蛋白
  select(ShortID) %>%                         # 只选 ID 列
  # 关键步骤：将 "P09863;A0A3B6GLF7" 这种格式拆分成多行
  separate_rows(ShortID, sep = ";") %>%
  pull(ShortID) %>%                           # 提取为向量
  str_trim() %>%
  unique()

print(paste("MaxQuant 复现差异蛋白数量 (拆分Group后):", length(list_mq)))


# 4. 绘制韦恩图 (Venn Diagram)

# 将三个列表放入一个 list 中
venn_data <- list(
  Paper_PD1.4 = list_paper,
  My_PD3.2 = list_pd32,
  My_MaxQuant = list_mq
)

# 绘图
p <- ggvenn(
  venn_data, 
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"), # 红、蓝、绿
  stroke_size = 0.5, 
  set_name_size = 4,
  text_size = 4
) + 
  ggtitle("Comparison of DEPs: Paper vs. Re-analysis") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=14))

print(p)

# 保存图片
ggsave("Venn_Comparison.png", plot = p, width = 8, height = 6, dpi = 300)



# 5.  查看具体的重叠蛋白

# 找出三个方法都鉴定到的蛋白 (核心交集)
common_ids <- intersect(intersect(list_paper, list_pd32), list_mq)
print("三个方法共有的蛋白 ID:")
print(common_ids)

