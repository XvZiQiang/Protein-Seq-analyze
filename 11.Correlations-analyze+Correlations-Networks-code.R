# 0. 环境设置与包加载
# 设置库路径
new_lib <- "C:/R packages"
old_libs <- .libPaths("C:/生信软件/R/R-4.5.1/library")
.libPaths(c(new_lib, old_libs))
.libPaths()

# 加载必要的包
library(psych)
library(dplyr)
library(reshape2)
library(readr)
library(tibble)
library(scales) # 用于缩放线条宽度和点大小
library(ggraph)
library(igraph)
library(tidygraph)
library(ggplot2)
library(ggnewscale)




# 1. 设置路径与读取数据

# 定义路径
physio_wt_file  <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/Physiological-index-chart-cropping/生理指标-WT.csv"
physio_ank_file <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/Physiological-index-chart-cropping/生理指标-ANK.csv"
protein_dir     <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/Correlations/"
out_dir         <- "C:/Users/XZQ/Desktop/NCBI/pep-seq-project/plant-science-paper/Correlations/"

# 确保输出目录存在
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 定义样本名
wt_samples  <- c("WT-NT-1", "WT-NT-2", "WT-NT-3", "WT-ST-1", "WT-ST-2", "WT-ST-3")
ank_samples <- c("ANK-NT-1", "ANK-NT-2", "ANK-NT-3", "ANK-ST-1", "ANK-ST-2", "ANK-ST-3")


# --- A. 读取生理数据 ---
phy_wt <- read_csv(physio_wt_file, show_col_types = FALSE) %>% 
  column_to_rownames("Accession") %>% t()

phy_ank <- read_csv(physio_ank_file, show_col_types = FALSE) %>% 
  column_to_rownames("Accession") %>% t()

# 强制校正行名
if(nrow(phy_wt) == 6) rownames(phy_wt) <- wt_samples
if(nrow(phy_ank) == 6) rownames(phy_ank) <- ank_samples

# 获取所有生理指标名称
all_physio_names <- colnames(phy_wt)




# --- B. 读取蛋白数据 ---
# 注意：pattern 使用 "Fig4_.*\\.csv" 匹配
files <- list.files(protein_dir, pattern = "Fig4_.*\\.csv", full.names = TRUE)
if(length(files) == 0) stop("未找到蛋白文件，请检查路径！")

#合并三个EXCEL表得到pro_data_all
pro_data_all <- bind_rows(lapply(files, read_csv, show_col_types = FALSE)) %>% 
  distinct(Accession, .keep_all = TRUE)

# 构建矩阵
pro_wt <- pro_data_all %>% select(Accession, all_of(wt_samples)) %>% 
  column_to_rownames("Accession") %>% t()

pro_ank <- pro_data_all %>% select(Accession, all_of(ank_samples)) %>% 
  column_to_rownames("Accession") %>% t()



# 2. 第一步：执行 Pearson 相关性分析并导出表格
#定义相关性分析函数
export_correlation_table <- function(pro_mat, phy_mat, prefix) {
  
  message(paste0(">>> 正在计算 [", prefix, "] 组 Pearson 相关性..."))
  
  # 1. 数据清洗 (去除零方差)
  pro_mat <- pro_mat[, apply(pro_mat, 2, var, na.rm=TRUE) != 0]
  phy_mat <- phy_mat[, apply(phy_mat, 2, var, na.rm=TRUE) != 0]
  
  # 2. 计算相关性 (adjust = "none", 只看原始P值)
  cor_res <- corr.test(pro_mat, phy_mat, method = "pearson", adjust = "none")
  
  # 3. 整理数据
  edges <- melt(cor_res$r) %>% rename(Data1 = Var1, Data2 = Var2, rho = value)
  pvals <- melt(cor_res$p) %>% rename(Data1 = Var1, Data2 = Var2, pvalue = value)
  
  full_table <- left_join(edges, pvals, by = c("Data1", "Data2"))
  
  # 4. 筛选 (P < 0.05 且 |rho| > 0.8)
  final_table <- full_table %>%
    filter(abs(rho) > 0.8 & pvalue < 0.05) %>%
    mutate(
      relation = ifelse(rho > 0, "positive", "negative")
    ) %>%
    arrange(desc(abs(rho))) %>%
    select(Data1, Data2, rho, pvalue, relation)
  
  # 5. 导出
  file_name <- paste0(out_dir, prefix, "_Correlation_Result.csv")
  write.csv(final_table, file_name, row.names = FALSE)
  
  message(paste("    筛选出", nrow(final_table), "条连线 ->", file_name))
}

# --- 执行导出 ---
export_correlation_table(pro_wt, phy_wt, "WT")
export_correlation_table(pro_ank, phy_ank, "ANK")










# 3. 第二步：生成 EDGE 和 NODE 文件 (详细注释版)
# 定义生成文件的函数
generate_graph_files_custom <- function(prefix, all_physio_names) {
  
  # 1. 读取上一步生成的 Pearson 相关性分析结果
  corr_file <- paste0(out_dir, prefix, "_Correlation_Result.csv")
  
  if (!file.exists(corr_file)) {
    message(paste("错误：找不到文件:", corr_file))
    return(NULL)
  }
  
  # 读取数据
  corr_data <- read.csv(corr_file)
  message(paste0(">>> 正在为 [", prefix, "] 组生成 Edge 和 Node 文件..."))
  
  # ----------------------------------------------------------------------------
  # A. 生成 Edge 文件 (全部保留)
  # ----------------------------------------------------------------------------
  # 要求：from, to, relation, edge_color, edge_width, edge_type
  
  edge_df <- corr_data %>%
    rename(from = Data1, to = Data2) %>% # 重命名以匹配要求
    mutate(
      # edge_color: positive为#DFC27D, negative为#9F9F9F
      edge_color = ifelse(rho > 0, "#DFC27D", "#9F9F9F"),
      
      # edge_width: 跟rho内容一样，全部保留 (取绝对值，否则ggplot无法正确渲染宽度)
      edge_width = abs(rho),
      
      # edge_type: positive为solid, negative为dashed
      edge_type = ifelse(rho > 0, "solid", "dashed")
    ) %>%
    # 只保留指定的列
    select(from, to, relation, edge_color, edge_width, edge_type)
  
  # 保存 Edge 文件
  edge_out <- paste0(out_dir, prefix, "_Edge.csv")
  write.csv(edge_df, edge_out, row.names = FALSE)
  
  # ----------------------------------------------------------------------------
  # B. 生成 Node 文件 (筛选 Top 35)
  # ----------------------------------------------------------------------------
  # 要求表头：name, group, node_type, label, label_size, label_colour, 
  #           label_family, label_face, degree, colour_value, node_color, node_size, node_label
  
  # 1. 筛选 Top 35 蛋白 (基于 Data1 列出现次数)
  # 统计每个蛋白在相关性表格中出现的次数 (即度 Degree)
  protein_counts <- table(corr_data$Data1)
  
  # 排序并取前 35 个 (如果不足 35 个则取全部)
  sorted_proteins <- sort(protein_counts, decreasing = TRUE)
  top_n <- min(35, length(sorted_proteins))
  top_35_proteins <- names(sorted_proteins)[1:top_n]
  
  # 2. 获取所有生理指标 (Data2 列去重)
  # 这里我们取 Data2 中出现的所有生理指标
  physio_in_data <- unique(corr_data$Data2)
  all_physios_final <- unique(c(physio_in_data))
  
  # 3. 合并生成最终节点列表 Name
  final_node_names <- c(top_35_proteins, all_physios_final)
  
  # 4. 计算 Degree (用于 Node_Size 计算)
  # 这里的 Degree 定义为：Data1列出现次数 + Data2列出现次数
  # 注意：我们需要在完整的 corr_data 中计算，反映其真实连接数
  all_nodes_in_corr <- c(corr_data$Data1, corr_data$Data2)
  degree_lookup <- table(all_nodes_in_corr)
  
  # 5. 构建 Node 数据框
  node_df <- data.frame(name = final_node_names) %>%
    mutate(
      # group: 蛋白ID为Data1，生理指标为Data2
      group = ifelse(name %in% all_physios_final, "Data2", "Data1"),
      
      # node_type: 均为 16
      node_type = 16,
      
      # label: 同 name 列
      label = name,
      
      # label 样式固定值
      label_size = 3,
      label_colour = "#000000",
      label_family = "Arial",
      label_face = "bold",
      
      # degree: 从全局统计中查找，没找到的设为0
      degree = as.integer(degree_lookup[name]),
      degree = ifelse(is.na(degree), 0, degree),
      
      # colour_value: 同 degree
      colour_value = degree,
      
      # node_color: 蛋白#FB8861，生理#92C5DE
      node_color = ifelse(group == "Data1", "#FB8861", "#92C5DE"),
      
      # node_label: 同 group
      node_label = group
    )
  
  # 6. 计算 node_size (严格线性映射)
  # 公式: 5 + (Degree - Min) * (15 - 5) / (Max - Min)
  min_degree <- min(node_df$degree)
  max_degree <- max(node_df$degree)
  
  # 防止 Max == Min 导致除以零错误
  if (max_degree == min_degree) {
    node_df$node_size <- 10 # 如果所有点度数一样，给个中间值
  } else {
    node_df$node_size <- 5 + (node_df$degree - min_degree) * (10) / (max_degree - min_degree)
  }
  
  # 保存 Node 文件
  node_out <- paste0(out_dir, prefix, "_Node.csv")
  write.csv(node_df, node_out, row.names = FALSE)
  
  message(paste("    Edge文件已保存 ->", edge_out))
  message(paste("    Node文件已保存 ->", node_out, "(包含", nrow(node_df), "个节点)"))
}



# 执行生成函数
generate_graph_files_custom("WT", all_physio_list)
generate_graph_files_custom("ANK", all_physio_list)








# 4. 第三步：利用 Edge/Node 文件绘图 (完整注释版)
# 定义绘图函数
# 终极绘图函数：自定义图例 + 样式精修 (Degree黑圈/Rho分级/布局调整)

plot_network_final_custom_legend <- function(prefix, title_text) {
  
  message(paste0(">>> 正在绘制 [", prefix, "] 网络图..."))
  
  # 1. 读取文件
  # ----------------------------------------------------------------------------
  node_path <- paste0(out_dir, prefix, "_Node.csv")
  edge_path <- paste0(out_dir, prefix, "_Edge.csv")
  
  if (!file.exists(node_path) || !file.exists(edge_path)) {
    message("错误：找不到CSV文件，请先运行生成步骤。")
    return()
  }
  
  # 读取数据
  nodes <- read.csv(node_path)
  edges <- read.csv(edge_path)
  
  # 2. 数据过滤与准备
  # ----------------------------------------------------------------------------
  # 确保 Edge 的两端都在 Node 表里
  edges_filtered <- edges %>%
    filter(from %in% nodes$name & to %in% nodes$name)
  
  if (nrow(edges_filtered) == 0) {
    message("警告：没有有效的连线，跳过。")
    return()
  }
  
  # 3. 构建图对象 & 排序
  # ----------------------------------------------------------------------------
  # 转换为 tidygraph 对象
  g <- graph_from_data_frame(d = edges_filtered, vertices = nodes, directed = FALSE)
  g_tbl <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    # 关键排序：先按 Group 分类(左右分屏)，再按 Degree 降序(大点在中间)
    mutate(group = factor(group, levels = c("Data1", "Data2"))) %>% 
    arrange(group, desc(degree))
  
  # 4. 开始绘图 (ggraph)
  # ----------------------------------------------------------------------------
  
  # 定义颜色方案
  cols_node <- c("Data1" = "#FB8861", "Data2" = "#92C5DE") # 橙色和蓝色
  labels_node <- c("Data1" = "Protein of seeds", "Data2" = "Parameters of seedlings")
  
  cols_edge <- c("positive" = "#E5C494", "negative" = "#9F9F9F", "negtive" = "#9F9F9F")
  types_edge <- c("positive" = "solid", "negative" = "dashed", "negtive" = "dashed")
  
  p <- ggraph(g_tbl, layout = 'linear', circular = TRUE) + 
    
    # --- A. 绘制连线 (Edge) ---
    # edge_width 映射到 CSV 中的 rho 值 (绝对值)
    geom_edge_link(aes(edge_width = edge_width,   
                       color = relation,          
                       linetype = relation),      
                   alpha = 0.8) + # 设置透明度，防止太乱
    
    # [关键] 设置 Rho 图例刻度：0.85, 0.90, 0.95
    scale_edge_width_continuous(
      range = c(0.2, 1.5),             # 线条粗细范围
      breaks = c(0.85, 0.90, 0.95),    # 图例上显示的刻度
      limits = c(0.8, 1.0),            # 数据的范围
      name = "abs_rho"                 # 图例标题
    ) +
    
    # 设置连线颜色 (正=金, 负=灰)
    scale_edge_color_manual(values = cols_edge, name = "relation") +
    
    # 设置连线虚实 (正=实, 负=虚)
    scale_edge_linetype_manual(values = types_edge, name = "relation") +
    
    # --- B. 绘制节点 (Node) ---
    # shape=21 表示带边框的圆，fill控制填充色，color控制边框色
    geom_node_point(aes(size = degree,            
                        fill = group),            
                    shape = 21, color = "white", stroke = 0.5) +
    
    # [关键] 设置 Degree 图例刻度：1, 6, 18
    scale_size_continuous(
      range = c(3, 10),                # 点的大小范围
      breaks = c(1, 7, 14),            # 图例上显示的刻度
      name = "degree"                  # 图例标题
    ) +
    
    # 设置节点颜色
    scale_fill_manual(
      values = cols_node, 
      labels = labels_node,
      name = "group"
    ) +
    
    # --- C. 绘制标签 (Label) ---
    geom_node_text(aes(x = x * 1.06, y = y * 1.06, 
                       label = label, 
                       angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                   hjust = 'outward', size = 3, fontface = "bold", family = "sans") +
    
    # 人为扩展画布范围
    # 圆的半径通常是1，标签大概延伸到1.2。
    # 这里我们把 x 轴范围强制设为 -1.5 到 1.8 (向右多留空)
    # y 轴也稍微扩一点防止标题太挤
    # 如果调整图例位置，调整这行代码
    expand_limits(x = c(-1.5, 1.8), y = c(-1.2, 1.2)) +
    
    # --- D. 主题设置 (Theme) ---
    coord_fixed(clip = "off") + 
    theme_void() + # 去除背景
    
    labs(title = title_text) +
    
    theme(
      # 背景纯白
      plot.background = element_rect(fill = "white", color = NA),
      
      # 标题位置：hjust=0.5居中，margin(b=30) 增加标题到底部的距离(向上移)
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 30)),
      
      # 边距设置：t=top, r=right, b=bottom, l=left
      # 增加 r (右边距) 以便容纳图例，增加 t (上边距) 让标题不顶格
      plot.margin = margin(t = 50, r = 20, b = 50, l = 50),
      
      # 图例整体设置
      legend.position = "right",        # 图例在右侧
      legend.box = "vertical",          # 图例垂直排列
      legend.spacing.y = unit(0.5, "cm"), # 图例之间的间距
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    
    # --- E. 深度自定义图例 (Guides) ---
    guides(
      # 1. Group 图例：点变大，保持原来的颜色
      fill = guide_legend(
        order = 1, 
        override.aes = list(size = 5)
      ),
      
      # 2. Degree 图例：【关键】强制显示为黑色圆圈
      size = guide_legend(
        order = 2,
        override.aes = list(
          fill = "black",  # 填充变黑
          color = "black", # 边框变黑
          shape = 21       # 保持圆圈形状
        )
      ),
      
      # 3. abs_rho 图例：显示线条粗细，颜色统一为黑色或灰色以便观察
      edge_width = guide_legend(
        order = 3,
        override.aes = list(edge_color = "black")
      ),
      
      # 4. Relation 图例：合并 color 和 linetype
      edge_color = guide_legend(order = 4),
      edge_linetype = guide_legend(order = 4)
    )
  
  # 保存图片
  out_png <- paste0(out_dir, prefix, "_Network_Final_LegendFixed.png")
  ggsave(out_png, plot = p, width = 13, height = 10, dpi = 300) # 宽度增加一点给图例
  message(paste("    图片已保存 ->", out_png))
}

# --- 执行绘图 ---
plot_network_final_custom_legend("WT", "Figure 7a: WT Correlation Network")
plot_network_final_custom_legend("ANK", "Figure 7b: ANK Correlation Network")

