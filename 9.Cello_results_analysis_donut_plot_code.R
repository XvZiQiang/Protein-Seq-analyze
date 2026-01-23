#环境设置与包加载
# 设置库路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")
.libPaths(c(new_lib, old_libs))
.libPaths()


library(dplyr)
library(stringr)
library(ggplot2)


# 1. 设置路径
# 存放 CELLO 结果 txt 的文件夹
cello_result_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/subcellular-location/cello_results/"
# 结果保存路径
output_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/subcellular-location/cello_plot/"

# 获取所有 txt 文件
txt_files <- list.files(cello_result_dir, pattern = "\\.txt$", full.names = TRUE)

if(length(txt_files) == 0) stop("未找到 txt 文件！请检查路径和文件名。")


# 2. 定义配色 (参考文章 Figure 3e-h)
# ------------------------------------------------------------------------------
cello_colors <- c(
  "Nuclear" = "#748DC5",          # 蓝紫色
  "Cytoplasmic" = "#94D483",      # 浅绿色
  "Extracellular" = "#C995C7",    # 紫粉色
  "Chloroplast" = "#D9F09E",      # 黄绿色 (亮)
  "Mitochondrial" = "#4CA8A2",    # 青色
  "PlasmaMembrane" = "#3B7A9E",   # 深蓝色 (CELLO输出是连写的PlasmaMembrane)
  "Plasma Membrane" = "#3B7A9E",  # 兼容带空格的写法
  "Lysosomal" = "#D66F6D",        # 红色
  "Vacuole" = "#6D5A8C",          # 深紫
  "Golgi" = "#A65C88",            # 暗红/紫
  "Cytoskeletal" = "#445588",     # 深蓝灰
  "ER" = "#553344",               # 暗色
  "Peroxisomal" = "#224422"       # 深绿
)



# 3. 解析parse_cello_txt_final函数：读取单个 txt 并提取最佳定位
# 3. 修复后的解析函数：按 ID 分块处理，防止越界
# ------------------------------------------------------------------------------
parse_cello_txt_final <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  id_indices <- grep("^SeqID:", lines)
  
  if (length(id_indices) == 0) return(data.frame())
  
  results <- data.frame(Accession = character(), Location = character(), stringsAsFactors = FALSE)
  
  for (k in 1:length(id_indices)) {
    # 确定当前蛋白块范围
    start_line <- id_indices[k]
    if (k < length(id_indices)) {
      end_line <- id_indices[k+1] - 1
    } else {
      end_line <- length(lines)
    }
    chunk <- lines[start_line:end_line]
    
    # 提取 ID
    current_id <- str_trim(sub("SeqID:", "", chunk[1]))
    
    # 查找预测结果起始点
    pred_header_idx <- grep("CELLO Prediction:", chunk)
    
    if (length(pred_header_idx) > 0) {
      idx <- pred_header_idx[1]
      check_limit <- min(length(chunk), idx + 20)
      
      if (idx + 1 <= check_limit) {
        for (j in (idx + 1):check_limit) {
          line_content <- chunk[j]
          clean_line <- str_trim(line_content)
          
          # 【核心修复 1】如果行是空的，跳过
          if (nchar(clean_line) == 0) next
          
          # 【核心修复 2】如果行以星号开头 (例如 *******)，说明是分隔线，直接结束当前蛋白
          if (grepl("^\\*", clean_line)) break
          
          # 【核心修复 3】正常的预测行必须包含星号，且不以星号开头
          if (grepl("\\*", clean_line)) {
            # 提取第一个单词
            loc <- str_split(clean_line, "\\s+")[[1]][1]
            
            # 双重保险：确保提取出来的 loc 不是一串星号
            if (!grepl("\\*", loc)) {
              results <- rbind(results, data.frame(Accession = current_id, Location = loc))
            }
          }
        }
      }
    }
  }
  return(results)
}

# 4. 循环处理并绘图
# ------------------------------------------------------------------------------
for (file in txt_files) {
  
  comp_name <- sub("_CELLO\\.txt$", "", basename(file))
  comp_name <- sub("_for_CELLO", "", comp_name)
  
  message(paste("正在处理:", comp_name))
  
  # A. 解析
  df <- parse_cello_txt_final(file)
  
  if (nrow(df) == 0) {
    warning(paste("  -> 跳过空文件:", comp_name))
    next
  }
  
  # 保存干净的解析表格
  write.csv(df, paste0(output_dir, comp_name, "_Location_Parsed.csv"), row.names = FALSE)
  
  # B. 统计数量
  plot_data <- df %>%
    group_by(Location) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count)) %>%
    mutate(
      fraction = Count / sum(Count),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      label_pos = (ymax + ymin) / 2
    )
  
  total_count <- sum(plot_data$Count)
  
  # C. 绘图
  p <- ggplot(plot_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = Location)) +
    geom_rect(color = "black", size = 0.2) + 
    scale_fill_manual(values = cello_colors) +
    coord_polar(theta = "y") +
    xlim(c(0.5, 4)) +
    
    annotate("text", x = 0.5, y = 0, label = paste0("Total=", total_count), 
             size = 5, fontface = "bold") +
    
    scale_fill_manual(
      values = cello_colors, 
      labels = paste0(plot_data$Count, " ", plot_data$Location),
      breaks = plot_data$Location
    ) +
    
    labs(title = comp_name, fill = NULL) +
    
    theme_void() + 
    theme(
      plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 10)
    )
  
  # D. 保存图片
  ggsave(paste0(output_dir, comp_name, "_Donut.png"), width = 7, height = 5, dpi = 300)
  print(paste("  -> 图片已保存:", comp_name, "_Donut.png"))
}







