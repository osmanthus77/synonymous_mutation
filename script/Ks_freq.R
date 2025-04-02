# 加载R包
library(pheatmap)
library(readr)
library(dplyr)

# define function: 将 Ks值[0,1.5]划分为不同区间，统计每个区间的Ks值的频率
bins_ks_freq <- function(column_data, breaks) {
    # 过滤NA
    column_data <- column_data[!is.na(column_data)]
    # 统计总个数
    total_count <- length(column_data)
    # 依据 Ks值的大小分组
    bins <- cut(column_data, breaks = breaks, include.lowest = TRUE)
    # 每个分组中元素的频数
    counts <- table(bins)
    counts_df <- as.data.frame(counts)
    colnames(counts_df) <- c("ks_bin", "count")
    # 每个分组中元素的频率 [0,100]
    counts_df$count_normolized <- (counts_df$count / total_count)*100
    return(counts_df)
}

# 导入输入数据(即 Ks 值的tsv文件)
in_df <- read_tsv("Ks_Cdomain/Ks_Cdomain_all.tsv")
## 将n/c值全部转变为NA
in_df[] <- lapply(in_df, function(x) as.numeric(as.character(x)))
## 删除全为NA的行
in_df <- in_df[rowSums(is.na(in_df)) < ncol(in_df), ]

# 打印 Ks值 的最小值、最大值和 >1.5的值
print(max(in_df, na.rm = TRUE))    # 1.79848
print(min(in_df, na.rm = TRUE))    #  0.00796442
print(in_df[in_df > 1.5], na.rm = TRUE)

# >1.5的值定义为1.5, 处理前过滤 NA
in_df[] <- lapply(in_df, function(x) {
  x[!is.na(x) & x > 1.5] <- 1.5
  return(x)
})

# 检查ks值分布 —— 是否存在小于0或大于1.5的值
all_in_range <- all(in_df >= 0 & in_df <= 1.5, na.rm = TRUE)
if (all_in_range) {
  print("YES")
} else {
  print("NO")
}

# breaks
breaks <- seq(0,1.5, length.out = 61)
ks_freq_list <- lapply(in_df, bins_ks_freq, breaks = breaks)

# combine
ks_bins <- ks_freq_list[[1]]$ks_bin
ks_freq_df <- data.frame(ks_bins = ks_bins)
for (i in seq_along(ks_freq_list)) {
    ks_freq_df[[names(in_df)[i]]] <- ks_freq_list[[i]]$count_normolized
}
ks_matrix <- as.matrix(ks_freq_df[, -1])
rownames(ks_matrix) <- ks_freq_df$ks_bin

# plot
plot <- pheatmap(
    t(ks_matrix),
    color = colorRampPalette(c("#ffffff", "#e62d2d"))(60),
    scale = "none",  # 不进行标准化处理
    cluster_rows = FALSE, cluster_cols = FALSE,
    border_color = NA,  # 设置单元格边框颜色
    fontsize_row = 7, fontsize_col = 7,
    main = "Ks Heatmap",
    # display_numbers = TRUE, fontsize_number = 5 # 显示单元格数值
    cellwidth = 12, cellheight = 18,  # 设置单元格大小
)