#环境设置与包加载
# 设置库路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")
.libPaths(c(new_lib, old_libs))
.libPaths()

# 批量提取差异蛋白的表达量和注释信息 (用于绘制热图)
library(readxl)
library(dplyr)
library(stringr)


# 定义输入和输出路径
# 1. 差异蛋白 ID 列表所在的文件夹 (生成的 Map 文件，****_ID_Group_Map.csv)
id_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/subcellular-location/"

# 2. 总数据表 Excel 路径 (Table S3)
master_file <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/ScienceDirect_files_01Dec2025_15-27-11.695/1-s2.0-S0168945222002539-mmc4.xlsx"

# 3. 结果保存文件夹
out_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/Diff-Anonation-Heatmap/"

# 如果输出目录不存在，自动创建
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}



# 2. 读取总数据表 (Master Table)
# skip = 1 跳过第一行标题说明，直接读取表头
raw_data <- read_excel(master_file, sheet = "Table S3", skip = 1)

# 确保列名正确 (防止 Excel 第一列名字不对)
if(!"Accession" %in% colnames(raw_data)) {
  colnames(raw_data)[1] <- "Accession"
}

# 定义你需要提取的表达量列名
# 注意：请确保这些名字与 Excel 表头完全一致
sample_cols <- c(
  "WT-NT-1", "WT-NT-2", "WT-NT-3", 
  "WT-ST-1", "WT-ST-2", "WT-ST-3",
  "ANK-NT-1", "ANK-NT-2", "ANK-NT-3", 
  "ANK-ST-1", "ANK-ST-2", "ANK-ST-3"
)

# 检查列名是否存在
missing_cols <- setdiff(sample_cols, colnames(raw_data))
if(length(missing_cols) > 0) {
  stop(paste("错误：在 Excel 中找不到以下列名，请检查拼写：", paste(missing_cols, collapse=", ")))
}



# 3. 批量处理 4 个文件
# 获取目录下所有的 ID_Group_Map.csv 文件
map_files <- list.files(id_dir, pattern = "_ID_Group_Map.csv", full.names = TRUE)

for (f in map_files) {
  
  # 获取文件名 (不带路径)
  fname <- basename(f)
  # 生成新的文件名 (例如 WT-SP_vs_WT-NSP_Diff_Expression.csv)
  out_name <- sub("_ID_Group_Map.csv", "_Diff_Expression_Data.csv", fname)
  
  message(paste("正在处理:", fname))
  
  # A. 读取差异蛋白 ID 列表
  id_list <- read.csv(f)
  
  # B. 从总表中提取数据
  # 使用 inner_join 确保只提取列表里有的蛋白
  extracted_data <- id_list %>%
    select(Accession, Group) %>% # 保留 Accession 和 Group (Sig_Up/Down)
    inner_join(raw_data, by = "Accession") %>%
    select(
      Accession, 
      Group, 
      Description, # 保留描述，画图时可能需要看功能
      all_of(sample_cols) # 提取 12 个样品的表达量
    )
  
  # C. 保存结果
  out_path <- paste0(out_dir, out_name)
  write.csv(extracted_data, out_path, row.names = FALSE)
  
  message(paste("  -> 提取成功，共", nrow(extracted_data), "行，已保存至:", out_name))
}





# 4. 针对图片 Figure 4 的三类特定蛋白提取
# 文章图4的三个热图里的 ID 

# 定义文章图4中的 ID 列表
fig4_ids <- c(
  # (a) Redox homeostasis
  "W5ET85", "W5GW81", "Q41579", "A0A3B6LKF5", "V5NW73", "A0A3B6JMF7", "D0PRB4", "W5FUB4", 
  "A0A3B6ISK3", "A0A3B6PUE9", "Q6W8Q2", "A0A3B6IMK6", "Q6PLQ7", "H2DPU1", "A0A3B6JPF1", 
  "A0A2X0U184", "A0A1D5ULD4", "A0A3B6N2Z1", "W5B8D6", "A0A3B6TEG7", "A0A3B6PRG2", "A0T1D6", 
  "A0A3B6PLT1", "A0A3B6PSL8", "A0A3B6I4M8", "A0A3B6R8P8", "A0A3B6RLR1", "A0A1B5GEP8", 
  "A0A3B6ESH8", "A0A3B6REP2", "A0A3B6GVN0", "D4P8R8", "P04464",
  
  # (b) Photosynthesis
  "A0A3B6AWA4", "A0A3B6DNQ7", "A0A3B6TZK7", "A0A3B6MZJ7", "A0A3B6RQ08", "A0A077S241", "A0A3B6I2T5",
  
  # (c) Carbohydrate
  "A0A077RZ28", "A0A3B6QPD1", "A0A3B6MII9", "W5G100", "A0A3B6C8Z7", "A0A3B6HUV0", 
  "A0A3B6B0U0", "A0A3B6C1I7", "W5D0E3", "W5C4N4", "A0A3B6KHT1", "A0A3B6PFH8", "Q8L5C6", 
  "Q4W6G2", "A0A3B6QI71", "A0A3B6D5T4", "A0A3B6ASR5", "A0A3B6SLM3", "A0A1D5V0H4", 
  "A0A3B5XW30", "A0A1D5V0T8", "A7BJ78", "A0A3B6H0W1", "A0A3B6T820", "A0A3B5XU98", 
  "A0A3B6KG51", "A7BJ77", "A0A3B6TLI0", "A0A023W638", "A0A3B6GNS8", "A0A3B6ES99"
)

# 刚刚得到的4个表里面，提取这部分特定蛋白ALL，一起保存
fig4_data <- raw_data %>%
  filter(Accession %in% fig4_ids) %>%
  select(Accession, Description, all_of(sample_cols))

# 保存图4专用数据
write.csv(fig4_data, paste0(out_dir, "Figure4_All_Specific_Proteins_to_heatmap.csv"), row.names = FALSE)





# 定义三个列表，分别对应 Figure 4 的 a, b, c，分别提取到3个表格里
# 列表 1: Redox homeostasis (Figure 4a)
ids_redox <- c(
  "W5ET85", "W5GW81", "Q41579", "A0A3B6LKF5", "V5NW73", "A0A3B6JMF7", "D0PRB4", 
  "W5FUB4", "A0A3B6ISK3", "A0A3B6PUE9", "Q6W8Q2", "A0A3B6IMK6", "Q6PLQ7", "H2DPU1", 
  "A0A3B6JPF1", "A0A2X0U184", "A0A1D5ULD4", "A0A3B6N2Z1", "W5B8D6", "A0A3B6TEG7", 
  "A0A3B6PRG2", "A0T1D6", "A0A3B6PLT1", "A0A3B6PSL8", "A0A3B6I4M8", "A0A3B6R8P8", 
  "A0A3B6RLR1", "A0A1B5GEP8", "A0A3B6ESH8", "A0A3B6REP2", "A0A3B6GVN0", "D4P8R8", 
  "P04464"
)

# 列表 2: Photosynthesis and electron transfer (Figure 4b)
ids_photo <- c(
  "A0A3B6AWA4", "A0A3B6DNQ7", "A0A3B6TZK7", "A0A3B6MZJ7", "A0A3B6RQ08", 
  "A0A077S241", "A0A3B6I2T5"
)

# 列表 3: Carbohydrate metabolism (Figure 4c)
ids_carbo <- c(
  "A0A077RZ28", "A0A3B6QPD1", "A0A3B6PSL8", "A0A3B6MII9", "W5G100", "A0A3B6C8Z7", 
  "A0A3B6HUV0", "A0A3B6B0U0", "A0A3B6C1I7", "W5D0E3", "W5C4N4", "A0A3B6KHT1", 
  "A0A3B6PFH8", "Q8L5C6", "Q4W6G2", "A0A3B6QI71", "A0A3B6D5T4", "A0A3B6ASR5", 
  "A0A3B6SLM3", "A0A1D5V0H4", "A0A3B5XW30", "A0A1D5V0T8", "A7BJ78", "A0A3B6H0W1", 
  "A0A3B6T820", "A0A3B5XU98", "A0A3B6KG51", "A7BJ77", "A0A3B6TLI0", "A0A023W638", 
  "A0A3B6GNS8", "A0A3B6ES99"
)

# 定义提取函数
extract_and_save <- function(id_list, category_name, raw_df, output_path) {
  
  # 提取
  extracted <- raw_df %>%
    filter(Accession %in% id_list) %>%
    select(Accession, Description, all_of(sample_cols))
  
  # 检查提取数量
  message(paste0("[", category_name, "] 目标ID数: ", length(id_list), " -> 实际提取数: ", nrow(extracted)))
  
  # 保存
  file_name <- paste0(output_path, "Fig4_", category_name, ".csv")
  write.csv(extracted, file_name, row.names = FALSE)
  message(paste("  -> 已保存至:", file_name))
}



# 执行提取
# 1. 提取 Redox
extract_and_save(ids_redox, "Redox_homeostasis", raw_data, out_dir)

# 2. 提取 Photosynthesis
extract_and_save(ids_photo, "Photosynthesis", raw_data, out_dir)

# 3. 提取 Carbohydrate
extract_and_save(ids_carbo, "Carbohydrate_metabolism", raw_data, out_dir)









