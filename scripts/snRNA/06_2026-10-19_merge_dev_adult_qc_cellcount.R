# 合并胚胎与成年 QC / cellcount 图（年龄统一排序，参考 2026-08-28_cellratio_by_tissue_age_central 的排序方式）
# 数据：2026-06-05_develop_all_human_harmony_anno_meta.rds（胚胎，12 个年龄）
#       2026-06-05_new_adult_human_harmony_anno_meta.rds（成年，8 个年龄）
# 输出：Plot/2026-10-19_dev_adult_qc.pdf
#       Plot/2026-10-19_dev_adult_cellcount_barplot.pdf
# 样式与 02/05 summary 脚本一致，配色统一引用 00_colors.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

date <- "2026-10-19_dev_adult"
data_dir <- "/project/work/2026-05-25_create_summary/Data/"
plot_dir <- "/project/work/2026-05-25_create_summary/Plot/"

source("/project/work/00use_scripts/00_colors.R")

age_order <- c(
  "70days", "79days", "87days", "91days", "100days", "103days",
  "116days", "137days", "141days", "142days", "162days", "165days",
  "21years", "25years", "30years", "32years",
  "67years", "70years", "75years", "86years"
)

dev_meta <- readRDS(file.path(data_dir, "2026-06-05_develop_all_human_harmony_anno_meta.rds"))
adult_meta <- readRDS(file.path(data_dir, "2026-06-05_new_adult_human_harmony_anno_meta.rds"))

metadata <- rbind(
  data.frame(
    donor_age = as.character(dev_meta$donor_age),
    percent.mt = as.numeric(dev_meta$percent.mt),
    nCount_RNA = as.numeric(dev_meta$nCount_RNA),
    nFeature_RNA = as.numeric(dev_meta$nFeature_RNA),
    stringsAsFactors = FALSE
  ),
  data.frame(
    donor_age = as.character(adult_meta$donor_age),
    percent.mt = as.numeric(adult_meta$percent.mt),
    nCount_RNA = as.numeric(adult_meta$nCount_RNA),
    nFeature_RNA = as.numeric(adult_meta$nFeature_RNA),
    stringsAsFactors = FALSE
  )
)
metadata$donor_age <- factor(metadata$donor_age, levels = age_order)
stopifnot(!any(is.na(metadata$donor_age)))

## ---- 1. QC: 4 个 boxplot 纵向排列（percent.mt / nCount_RNA / nFeature_RNA / nFeature_nCount） ----
p1 <- ggplot(metadata, aes(fill = donor_age, y = percent.mt, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    scale_y_continuous()+
    labs(y = "percent.mt") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
         text = element_text(size = 15))

p2 <- ggplot(metadata, aes(fill = donor_age, y = nCount_RNA, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    scale_y_log10()+
    labs(y = "nCount_RNA") +
    coord_cartesian(ylim = c(NA, 15000)) +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
         text = element_text(size = 15))

p3 <- ggplot(metadata, aes(fill = donor_age, y = nFeature_RNA, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    scale_y_log10()+
    labs(y = "nFeature_RNA") +
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
         text = element_text(size = 15))

p4 <- ggplot(metadata, aes(fill = donor_age, y = nFeature_RNA / nCount_RNA, x = donor_age)) +
    geom_boxplot(outlier.shape = NA, show.legend = T)+
    scale_fill_manual(values = cols.age)+
    labs(y = "nFeature / nCount") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
          axis.title.x = element_blank(),
         text = element_text(size = 15))

qc_pdf <- file.path(plot_dir, paste0(date, "_qc", ".pdf"))
pdf(qc_pdf, width = 20, height = 15)
  wrap_plots(list(p1, p2, p3, p4), ncol = 1, byrow = F)
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
