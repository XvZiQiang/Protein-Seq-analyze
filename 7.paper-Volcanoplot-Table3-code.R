# 1. 环境设置与包加载
# 设置库路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")
.libPaths(c(new_lib, old_libs))
.libPaths()
# 加载包
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

# 2. 数据读取与预处理

# 定义工作目录和文件路径
work_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/ScienceDirect_files_01Dec2025_15-27-11.695/"
file_name <- "1-s2.0-S0168945222002539-mmc4.xlsx"
full_path <- paste0(work_dir, file_name)
save_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/paper-Volcano/"


# 读取数据
raw_data <- read_excel(full_path, sheet = "Table S3", skip = 1)










# 3. 定义通用分析绘图函数

analyze_and_plot <- function(data, group_ctrl_cols, group_treat_cols, comparison_name,
                             x_breaks = waiver(), x_limits = NULL) {
  
  message(paste("正在处理对比组:", comparison_name, "..."))
  
  # --- A. 计算 FC 和 P-value ---
  results <- data.frame(
    Accession = data$Accession,
    Mean_Ctrl = NA,
    Mean_Treat = NA,
    FC = NA,
    Log2FC = NA,
    P_value = NA
  )
  
  # 循环计算 (虽然慢但逻辑清晰)
  for(i in 1:nrow(data)) {
    vals_ctrl <- as.numeric(data[i, group_ctrl_cols])
    vals_treat <- as.numeric(data[i, group_treat_cols])
    
    # 均值计算
    mean_c <- mean(vals_ctrl, na.rm = TRUE)
    mean_t <- mean(vals_treat, na.rm = TRUE)
    
    # FC 计算
    fc <- mean_t / mean_c
    
    # T检验 (使用 tryCatch 防止数据全一样时报错)
    t_res <- tryCatch({
      t.test(vals_ctrl, vals_treat, var.equal = TRUE)
    }, error = function(e) return(NULL))
    
    results$Mean_Ctrl[i] <- mean_c
    results$Mean_Treat[i] <- mean_t
    results$FC[i] <- fc
    results$Log2FC[i] <- log2(fc)
    results$P_value[i] <- if(!is.null(t_res)) t_res$p.value else NA
  }
  
  # --- B. 6种分类设定 (核心步骤) ---
  fc_cut_up <- 1.2
  fc_cut_down <- 0.833  # 即 1/1.2
  p_cut <- 0.05
  
  results <- results %>%
    mutate(
      Group = case_when(
        # 1. Sig_Up: FC > 1.2 且 P < 0.05
        FC > fc_cut_up & P_value < p_cut ~ "Sig_Up",
        
        # 2. FC_Up_Only: FC > 1.2 且 P >= 0.05
        FC > fc_cut_up & P_value >= p_cut ~ "FC_Up_Only",
        
        # 3. Sig_Down: FC < 0.833 且 P < 0.05
        FC < fc_cut_down & P_value < p_cut ~ "Sig_Down",
        
        # 4. FC_Down_Only: FC < 0.833 且 P >= 0.05
        FC < fc_cut_down & P_value >= p_cut ~ "FC_Down_Only",
        
        # 5. pVal_Only: FC 在 0.833~1.2 之间 且 P < 0.05
        FC >= fc_cut_down & FC <= fc_cut_up & P_value < p_cut ~ "pVal_Only",
        
        # 6. NoDiff: 其他情况
        TRUE ~ "NoDiff"
      ),
      MinusLog10P = -log10(P_value)
    )
  
  # --- C. 统计数量用于图例 ---
  # 定义因子顺序，保证图例颜色和标签一一对应
  results$Group <- factor(results$Group, levels = c("Sig_Up", "FC_Up_Only", "pVal_Only", "NoDiff", "FC_Down_Only", "Sig_Down"))
  
  counts <- table(results$Group)
  
  # 动态标签
  my_labels <- c(
    "Sig_Up" = paste0("Sig_Up (", counts["Sig_Up"], ")"),
    "FC_Up_Only" = paste0("FC_Up_Only (", counts["FC_Up_Only"], ")"),
    "pVal_Only" = paste0("pVal_Only (", counts["pVal_Only"], ")"),
    "NoDiff" = paste0("NoDiff (", counts["NoDiff"], ")"),
    "FC_Down_Only" = paste0("FC_Down_Only (", counts["FC_Down_Only"], ")"),
    "Sig_Down" = paste0("Sig_Down (", counts["Sig_Down"], ")")
  )
  my_labels <- gsub("\\(NA\\)", "(0)", my_labels)
  
  # --- D. 保存 CSV 数据 ---
  csv_name <- paste0(save_dir, comparison_name, "_Analysis_Data.csv")
  write.csv(results, csv_name, row.names = FALSE)
  print(paste("数据表已保存:", csv_name))
  
  # --- E. 绘制火山图 ---
  
  my_colors <- c(
    "Sig_Up" = "#D62728",       # 红色
    "FC_Up_Only" = "#FF9896",   # 粉色
    "pVal_Only" = "#E377C2",    # 紫粉色
    "NoDiff" = "#C7C7C7",       # 灰色
    "FC_Down_Only" = "#AEC7E8", # 浅蓝
    "Sig_Down" = "#1F77B4"      # 深蓝
  )
  
  # 提取要标记的基因 (只标显著的，防止太乱)
  label_data <- subset(results, Group %in% c("Sig_Up", "Sig_Down"))
  # 仅标记 P值最小的前 30 个 (可注释掉这行以标记所有)
  label_data <- label_data %>% arrange(P_value) %>% head(30)
  
  p <- ggplot(results, aes(x = Log2FC, y = MinusLog10P, color = Group)) +
    geom_point(alpha = 0.7, size = 1.5) +
    
    scale_color_manual(values = my_colors, labels = my_labels) +
    
    # 使用函数传入的参数设置 X 轴
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    
    geom_vline(xintercept = c(log2(fc_cut_up), log2(fc_cut_down)), 
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cut), 
               linetype = "dashed", color = "black", linewidth = 0.4) +
    
    geom_text_repel(
      data = label_data,
      aes(label = Accession),
      size = 2.5,
      box.padding = 0.4,
      point.padding = 0.2,
      max.overlaps = 30,
      show.legend = FALSE,
      color = "black"
    ) +
    
    labs(title = comparison_name,
         x = expression(log[2] ~ "(Fold Change)"),
         y = expression(-log[10] ~ "(P-value)"),
         color = "Significant") +
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )
  
  # --- F. 保存图片 ---
  plot_name <- paste0(save_dir, comparison_name, "_Volcano.png")
  ggsave(plot_name, plot = p, width = 7, height = 6, dpi = 300)
  print(paste("火山图已保存:", plot_name))
  
  return(p)
}









# 4. 执行四组对比分析

# 1. WT-SP vs WT-NSP
analyze_and_plot(
  data = raw_data,
  group_ctrl_cols = c("WT-NT-1", "WT-NT-2", "WT-NT-3"),
  group_treat_cols = c("WT-ST-1", "WT-ST-2", "WT-ST-3"),
  comparison_name = "WT-SP_vs_WT-NSP",
  x_breaks = c(-2, 0, 2),
  x_limits = c(-3.5, 3.5)
)

# 2. ANK-SP vs ANK-NSP
analyze_and_plot(
  data = raw_data,
  group_ctrl_cols = c("ANK-NT-1", "ANK-NT-2", "ANK-NT-3"),
  group_treat_cols = c("ANK-ST-1", "ANK-ST-2", "ANK-ST-3"),
  comparison_name = "ANK-SP_vs_ANK-NSP",
  x_breaks = c(-2, 0, 2),
  x_limits = c(-3.5, 3.5)
)

# 3. ANK-NSP vs WT-NSP
analyze_and_plot(
  data = raw_data,
  group_ctrl_cols = c("WT-NT-1", "WT-NT-2", "WT-NT-3"),
  group_treat_cols = c("ANK-NT-1", "ANK-NT-2", "ANK-NT-3"),
  comparison_name = "ANK-NSP_vs_WT-NSP",
  x_breaks = seq(-4, 4, 2),
  x_limits = c(-5, 5)
)

# 4. ANK-SP vs WT-SP
analyze_and_plot(
  data = raw_data,
  group_ctrl_cols = c("WT-ST-1", "WT-ST-2", "WT-ST-3"),
  group_treat_cols = c("ANK-ST-1", "ANK-ST-2", "ANK-ST-3"),
  comparison_name = "ANK-SP_vs_WT-SP",
  x_breaks = seq(-5, 5, 2.5),
  x_limits = c(-6, 6)
)


# 5. 批量提取并保存差异蛋白总表 (Sig_Up + Sig_Down)

# 定义 4 个对比组名称
comp_names <- c(
  "WT-SP_vs_WT-NSP", 
  "ANK-SP_vs_ANK-NSP", 
  "ANK-NSP_vs_WT-NSP", 
  "ANK-SP_vs_WT-SP"
)

# 确保保存目录存在
if (!exists("save_dir")) {
  save_dir <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/paper-Volcano/"
}

# 循环处理
for (comp in comp_names) {
  
  # 1. 读取全量数据
  full_file_path <- paste0(save_dir, comp, "_Analysis_Data.csv")
  #检查是否存在表格
  if (file.exists(full_file_path)) {
    
    full_data <- read.csv(full_file_path)
    
    
    # 2. 筛选 Sig_Up 和 Sig_Down
    # 使用 %in% 同时匹配两种情况
    dep_data <- full_data %>% 
      filter(Group %in% c("Sig_Up", "Sig_Down")) %>%
      arrange(Group, P_value) #按分组和P值排序，方便查看
    
    # 3. 保存文件
    if (nrow(dep_data) > 0) {
      # 命名为 _Sig_DEPs.csv (Differentially Expressed Proteins)
      out_file <- paste0(save_dir, comp, "_Sig_DEPs.csv")
      write.csv(dep_data, out_file, row.names = FALSE)
      
      message(paste("成功保存:", comp, "差异蛋白表 ->", nrow(dep_data), "行"))
    } else {
      warning(paste(comp, "没有检测到显著差异蛋白，未保存文件。"))
    }
    
  } else {
    warning(paste("找不到文件:", full_file_path))
  }
}
