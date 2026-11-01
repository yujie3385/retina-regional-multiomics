# 最终整理脚本 01 · 胚胎+成年合并 QC 箱线图与 cellcount 柱状图 (2026-10-27)
# 来源: 00use_scripts/18_2026-10-20_merge_dev_adult_qc_cellcount.R (原脚本即最终版, 内容原样收录)
# 产出: Plot/2026-10-20_dev_adult_qc.pdf                (20x12.4, 4 指标箱线, p3/p4 高度压缩)
#       Plot/2026-10-20_dev_adult_cellcount_barplot.pdf (20x8, 18 个年龄细胞数)
#
# 合并胚胎与成年 ATAC QC / cellcount 图（年龄统一排序，参考 RNA 项目 2026-10-19_merge_dev_adult_qc_cellcount.R 的布局）
# 数据：2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds（胚胎，12 个年龄）
#       2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds（成年，6 个年龄）
# 输出：Plot/2026-10-20_dev_adult_qc.pdf
#       Plot/2026-10-20_dev_adult_cellcount_barplot.pdf
# QC 面板与 04/08 summary 脚本一致（nCount_ATAC / TSS.enrichment / blacklist_ratio / nucleosome_signal）
# 2026-10-25 调整: p3/p4 Y 轴收窄 (blacklist_ratio ≤ 0.01, nucleosome_signal ≤ 3, 均覆盖全部箱须)
#                 + p3/p4 面板高度压缩 (0.65×), 总高 15→12.4 (p1/p2 尺寸不变, 宽 20 不变)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

date <- "2026-10-20_dev_adult"
data_dir <- "/project/work/2026-05-08_ATAC_signac_all_sample/Data/"
plot_dir <- "/project/work/2026-05-08_ATAC_signac_all_sample/Plot/"

cols.age <- c(
  `70days`   = "#66C2A5",
  `79days`   = "#FC8D62",
  `87days`   = "#8DA0CB",
  `91days`   = "#E78AC3",
  `100days`  = "#A6D854",
  `103days`  = "#FFD92F",
  `116days`  = "#E5C494",
  `137days`  = "#B3B3B3",
  `141days`  = "#7FC97F",
  `142days`  = "#BEAED4",
  `162days`  = "#FDC086",
  `165days`  = "#FFFF99",
  `21years`  = "#FFFF33",
  `25years`  = "#E41A1C",
  `30years`  = "#377EB8",
  `32years`  = "#4DAF4A",
  `67years`  = "#984EA3",
  `70years`  = "#FF7F00",
  `75years`  = "#F781BF",
  `86years`  = "#E6AB02"
)

# ATAC 数据实际存在的 18 个年龄（30years / 32years 仅在 RNA 数据中，不参与绘图）
age_order <- c(
  "70days", "79days", "87days", "91days", "100days", "103days",
  "116days", "137days", "141days", "142days", "162days", "165days",
  "21years", "25years", "67years", "70years", "75years", "86years"
)

dev_meta <- readRDS(file.path(data_dir, "2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds"))
adult_meta <- readRDS(file.path(data_dir, "2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds"))

metadata <- rbind(
  data.frame(
    donor_age = as.character(dev_meta$donor_age),
    nCount_ATAC = as.numeric(dev_meta$nCount_ATAC),
    TSS.enrichment = as.numeric(dev_meta$TSS.enrichment),
    blacklist_ratio = as.numeric(dev_meta$blacklist_ratio),
    nucleosome_signal = as.numeric(dev_meta$nucleosome_signal),
    stringsAsFactors = FALSE
  ),
  data.frame(
    donor_age = as.character(adult_meta$donor_age),
    nCount_ATAC = as.numeric(adult_meta$nCount_ATAC),
    TSS.enrichment = as.numeric(adult_meta$TSS.enrichment),
    blacklist_ratio = as.numeric(adult_meta$blacklist_ratio),
    nucleosome_signal = as.numeric(adult_meta$nucleosome_signal),
    stringsAsFactors = FALSE
  )
)
metadata$donor_age <- factor(metadata$donor_age, levels = age_order)
stopifnot(!any(is.na(metadata$donor_age)))

## ---- 1. QC: 4 个 boxplot 纵向排列（nCount_ATAC / TSS.enrichment / blacklist_ratio / nucleosome_signal） ----
p1 <- ggplot(metadata, aes(fill = donor_age, y = nCount_ATAC, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    coord_cartesian(ylim = c(NA, 30000)) +
    labs(y = "nCount_ATAC") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 15))

p2 <- ggplot(metadata, aes(fill = donor_age, y = TSS.enrichment, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    labs(y = "TSS.enrichment") +
    coord_cartesian(ylim = c(NA, 12)) +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 15))

p3 <- ggplot(metadata, aes(fill = donor_age, y = blacklist_ratio, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    coord_cartesian(ylim = c(NA, 0.01)) +
    labs(y = "blacklist_ratio") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 15))

p4 <- ggplot(metadata, aes(fill = donor_age, y = nucleosome_signal, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    coord_cartesian(ylim = c(NA, 3)) +
    labs(y = "nucleosome_signal") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
          axis.title.x = element_blank(),
          text = element_text(size = 15))

qc_pdf <- file.path(plot_dir, paste0(date, "_qc", ".pdf"))
pdf(qc_pdf, width = 20, height = 12.4)
  wrap_plots(list(p1, p2, p3, p4), ncol = 1, byrow = F, heights = c(1, 1, 0.65, 0.65))
dev.off()

## ---- 2. cellcount barplot ----
plot_data <- metadata %>%
  count(donor_age)
plot_data$donor_age <- factor(plot_data$donor_age, levels = age_order)

p <- ggplot(plot_data, aes(x = donor_age, y = n, fill = donor_age)) +
  geom_bar(stat = "identity", width = 0.85) +
  geom_text(aes(label = n),
            vjust = -0.5,
            size = 5) +
  scale_fill_manual(values = cols.age) +
  labs(y = "number", x = "age") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid = element_line(color = "grey92")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

cc_pdf <- file.path(plot_dir, paste0(date, "_cellcount_barplot", ".pdf"))
pdf(cc_pdf, width = 20, height = 8)
  print(p)
dev.off()

message("Done: ", qc_pdf)
message("Done: ", cc_pdf)
message("dev cells = ", nrow(dev_meta), "; adult cells = ", nrow(adult_meta))
message("ages = ", nlevels(metadata$donor_age))
