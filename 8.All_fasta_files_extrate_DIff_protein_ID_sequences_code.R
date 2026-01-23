#环境设置与包加载
# 设置库路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")
.libPaths(c(new_lib, old_libs))
.libPaths()

# CELLO 预测数据准备：提取序列 FASTA + 提取 ID/Group 映射表

# 1. 加载必要的包
library(dplyr)
library(Biostrings) 

# 2. 设置路径
input_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/paper-Volcano/"
output_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/subcellular-location/"
db_path <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/download-data/uniprot_Triticum_aestivum_143178_20200217.fasta"

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 定义 4 个对比组名称
comp_names <- c(
  "WT-SP_vs_WT-NSP", 
  "ANK-SP_vs_ANK-NSP", 
  "ANK-NSP_vs_WT-NSP", 
  "ANK-SP_vs_WT-SP"
)

# 3. 加载 FASTA 数据库
fasta_db <- readAAStringSet(db_path)


# 清洗数据库的 ID (从 >tr|ID|... 中只提取 ID)
names(fasta_db) <- sapply(strsplit(names(fasta_db), "\\|"), function(x) {
  if(length(x) >= 2) return(x[2]) 
  return(x[1])
})
message("数据库加载完成！")






# 4. 循环处理 4 个文件
for (comp in comp_names) {
  
  # 输入文件路径
  infile <- paste0(input_dir, comp, "_Sig_DEPs.csv")
  
  if (file.exists(infile)) {
    message(paste("正在处理:", comp, "..."))
    
    # 读取差异蛋白表
    df <- read.csv(infile)
    
    # --- 步骤 A: 提取 ID 和 Group (用于后续画图) ---
    map_df <- df %>% select(Accession, Group)
    map_file <- paste0(output_dir, comp, "_ID_Group_Map.csv")
    write.csv(map_df, map_file, row.names = FALSE)
    
    # --- 步骤 B: 提取序列 (用于 CELLO 网站预测) ---
    target_ids <- df$Accession
    # 在数据库中找到存在的 ID
    valid_ids <- target_ids[target_ids %in% names(fasta_db)]
    
    if (length(valid_ids) > 0) {
      sub_fasta <- fasta_db[valid_ids]
      fasta_out_file <- paste0(output_dir, comp, "_for_CELLO.fasta")
      writeXStringSet(sub_fasta, fasta_out_file)
      message(paste("  -> 成功生成:", fasta_out_file))
    } else {
      warning("  -> 未找到匹配的序列，请检查 ID 格式！")
    }
    
  } else {
    warning(paste("找不到文件:", infile))
  }
}






