# 最终整理脚本 03 · snATAC celltype 占比堆叠柱状图 (central/peripheral 分组, 按分期 / 按年龄) (2026-10-27)
# 来源: 00use_scripts/21_2026-09-08_snATAC_cellratio_central.ipynb
# 产出: Plot/2026-08-28_cellratio_by_tissue_stage_central.pdf (10x8,  x = central/peripheral × Early~Old 五期)
#       Plot/2026-08-28_cellratio_by_tissue_age_central.pdf   (20x10, x = central/peripheral × 18 个年龄)
# 说明: fovea/macula 统一为 central; 占比按各分组列内归一 (position = "fill")
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

date     <- "2026-08-28"
plot_dir <- "/project/work/2026-05-08_ATAC_signac_all_sample/Plot/"

cols.palette2 <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"))
## celltype colors (palette2[1:8] + palette2[13:14])
cols.celltype <- c(cols.palette2[1:8], cols.palette2[13:14])
names(cols.celltype) <- c('Rod', 'Cone', 'AC', 'BC', 'HC', 'MG', 'RGC', 'RPC', 'Micro', 'Astro')

develop_meta <- readRDS("/project/work/2026-05-08_ATAC_signac_all_sample/Data/2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds")
adult_meta   <- readRDS("/project/work/2026-05-08_ATAC_signac_all_sample/Data/2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds")

cat("develop cells:", nrow(develop_meta), "\n")
cat("adult cells:  ", nrow(adult_meta),   "\n")
cat("develop celltypes:", paste(sort(unique(as.character(develop_meta$celltype))), collapse=" | "), "\n")
cat("adult   celltypes:", paste(sort(unique(as.character(adult_meta$celltype))),   collapse=" | "), "\n")

## snATAC 实际存在的 age 顺序（12 个发育 days + 6 个成年 years）
age_order <- c('70days','79days','87days','91days','100days','103days',
               '116days','137days','141days','142days','162days','165days',
               '21years','25years','67years','70years','75years','86years')
stage_order          <- c('Early', 'Mid', 'Late', 'Young', 'Old')
tissue_merged_order  <- c('central', 'peripheral')
tissue_age_merged_order   <- as.vector(outer(tissue_merged_order, age_order,   function(t, a) paste(t, a, sep = "_")))
tissue_stage_merged_order <- as.vector(outer(tissue_merged_order, stage_order, function(t, s) paste(t, s, sep = "_")))

meta.merged <- bind_rows(develop_meta, adult_meta) |>
  dplyr::mutate(
    donor_age     = as.character(donor_age),
    tissue        = as.character(tissue),
    tissue_merged = ifelse(tissue %in% c("fovea", "macula"), "central", "peripheral"),
    stage = dplyr::case_when(
      donor_age %in% c("70days",  "79days",  "87days",  "91days")              ~ "Early",
      donor_age %in% c("100days", "103days", "116days")                        ~ "Mid",
      donor_age %in% c("137days", "141days", "142days", "162days", "165days") ~ "Late",
      donor_age %in% c("21years", "25years")                                   ~ "Young",
      donor_age %in% c("67years", "70years", "75years", "86years")             ~ "Old"
    )
  )

## ---- 图 1: by tissue × stage (10x8) ----
plot.data.stage.merged <- meta.merged |>
  dplyr::mutate(tissue_stage_merged = paste(tissue_merged, stage, sep = "_")) |>
  dplyr::select(celltype, tissue_stage_merged) |>
  dplyr::filter(!is.na(celltype), !is.na(tissue_stage_merged)) |>
  dplyr::mutate(tissue_stage_merged = factor(tissue_stage_merged, levels = tissue_stage_merged_order)) |>
  dplyr::count(tissue_stage_merged, celltype, name = "value") |>
  dplyr::mutate(
    tissue_stage_merged = droplevels(tissue_stage_merged),
    celltype = factor(celltype, levels = names(cols.celltype))
  )

p_stage <- ggplot(plot.data.stage.merged, aes(x = tissue_stage_merged, y = value, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill", show.legend = TRUE) +
  scale_fill_manual(name = NULL, values = cols.celltype) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),
    axis.text.y  = element_text(color = "black"),
    text         = element_text(size = 15),
    axis.line.y  = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)
  ) +
  labs(x = NULL, y = "Cell proportion") +
  ggtitle("celltype (snATAC)")

pdf(paste0(plot_dir, date, "_cellratio_by_tissue_stage_central", ".pdf"), width = 10, height = 8)
  print(p_stage)
dev.off()

## ---- 图 2: by tissue × age (20x10) ----
plot.data.age.merged <- meta.merged |>
  dplyr::mutate(tissue_age_merged = paste(tissue_merged, donor_age, sep = "_")) |>
  dplyr::select(celltype, tissue_age_merged) |>
  dplyr::filter(!is.na(celltype), !is.na(tissue_age_merged)) |>
  dplyr::mutate(tissue_age_merged = factor(tissue_age_merged, levels = tissue_age_merged_order)) |>
  dplyr::count(tissue_age_merged, celltype, name = "value") |>
  dplyr::mutate(
    tissue_age_merged = droplevels(tissue_age_merged),
    celltype = factor(celltype, levels = names(cols.celltype))
  )

p_age <- ggplot(plot.data.age.merged, aes(x = tissue_age_merged, y = value, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill", show.legend = TRUE) +
  scale_fill_manual(name = NULL, values = cols.celltype) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y  = element_text(color = "black"),
    text         = element_text(size = 15),
    axis.line.y  = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    plot.margin  = margin(t = 5, r = 5, b = 20, l = 5)
  ) +
  labs(x = NULL, y = "Cell proportion") +
  ggtitle("celltype (snATAC)")

pdf(paste0(plot_dir, date, "_cellratio_by_tissue_age_central", ".pdf"), width = 20, height = 10)
  print(p_age)
dev.off()

cat("done: 2 cellratio PDFs\n")
