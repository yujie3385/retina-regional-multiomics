#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-21_11  AUCell violin byCT（development / adult 两分面版, 按用户要求）
# 在 2026-10-20_5b 基础上修改:
#   1. 已整合的成年数据(GSE277326 + GSE155288)合并为一个 "adult" 分面,
#      不再按数据集分面; 胚胎(GSE142526)为 "development" 分面
#      (strip 标题与主项目口径一致: development / adult)
#   2. 去掉 21y 供体(Macula_21y/Retina_21y, 此前质控确认数据质量异常)
#   3. 输出为新文件(_devadult 后缀), 不覆盖/删除 2026-10-20_5 版本
#   - 成年分面 x 轴: 27y -> 43y -> 78y -> 90y (各 central/peripheral)
#   - 分数沿用 2026-10-20_4 各数据集独立打分(与现版图完全同源)
#   - 其余风格不变: celltype 分页、共Y、45°黑x标签、居中标题、仅PDF
# 输出: Plot/2026-10-21_11_AUCell_violin_byCT_{DM,DAL,AAM}_{Central,Peripheral}_devadult.pdf
# =====================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
plot_dir <- file.path(proj_dir, "Plot")
sig_file <- "/project/work/2026-09-01_DEG_summary/Data/2026-09-01_2_AUCell_signature_definition.csv"

## ---- 设定 ----
period_lab <- c("development" = "development", "adult" = "adult")
age_order  <- c("FD125", "27y", "43y", "78y", "90y")   # 21y 已去除
celltype_order <- c("Rod", "Cone", "AC", "BC", "HC", "MG")
subclass_order <- c("Developmental maintained", "Developmental aging-lost",
                    "Adult acquired and maintained")
subclass_short <- c("Developmental maintained" = "DM",
                    "Developmental aging-lost" = "DAL",
                    "Adult acquired and maintained" = "AAM")
direction_order <- c("Macular direction", "Peripheral direction")
cols.tissue <- c(central = "#6A3D9A", peripheral = "#E69F00")

signature_key <- read.csv(sig_file, stringsAsFactors = FALSE) %>%
  select(signature_id, celltype, subclass, direction_class)

## ---- 读入三个数据集 meta + score ----
meta_all <- bind_rows(lapply(c("GSE142526", "GSE277326", "GSE155288"), function(ds) {
  md <- readRDS(file.path(data_dir, paste0("2026-10-20_4_AUCell_", ds, "_meta_data.rds")))
  md %>% select(cell_id, dataset, sample_id, tissue, donor, age_group, celltype,
                starts_with("AUCell_"))
}))
cat("cells loaded:", nrow(meta_all), "\n")

## ---- 年龄标签 + 去除21y + 分面标签 ----
plot_long <- meta_all %>%
  mutate(age_lab = case_when(
    dataset == "GSE142526" ~ "FD125",
    dataset == "GSE277326" ~ paste0(sub(".*_(\\d+)y$", "\\1", sample_id), "y"),
    dataset == "GSE155288" ~ paste0(sub("donor(\\d+)y$", "\\1", donor), "y")
  )) %>%
  filter(age_lab != "21y") %>%
  mutate(period_lab = ifelse(dataset == "GSE142526", "development", "adult")) %>%
  pivot_longer(starts_with("AUCell_"), names_to = "signature_id",
               values_to = "AUCell_score") %>%
  mutate(signature_id = sub("^AUCell_", "", signature_id)) %>%
  filter(!is.na(AUCell_score)) %>%
  left_join(signature_key, by = "signature_id") %>%
  ## 每个 signature 只画对应 celltype 的细胞; 胚胎 PR 同时进入 Rod/Cone 行
  filter(celltype.x == celltype.y |
           (celltype.x == "PR" & celltype.y %in% c("Rod", "Cone"))) %>%
  rename(celltype = celltype.x, sig_celltype = celltype.y) %>%
  mutate(
    stage_tissue = paste(age_lab, tissue, sep = " "),
    subclass     = factor(subclass, levels = subclass_order),
    direction_class = factor(direction_class, levels = direction_order)
  )

st_levels <- paste(rep(age_order, each = 2),
                   rep(c("central", "peripheral"), times = length(age_order)), sep = " ")
plot_long$stage_tissue <- factor(plot_long$stage_tissue, levels = st_levels)
plot_long$period_lab <- factor(plot_long$period_lab, levels = c("development", "adult"))

cat("plot rows:", nrow(plot_long), "\n")
cat("x slots:", nlevels(droplevels(plot_long$stage_tissue)), "\n")
cat("各分面细胞数:\n"); print(table(plot_long$period_lab, plot_long$age_lab))

## ---- 绘图: 每个 subclass × direction 一个 PDF, celltype 分页 ----
for (one_subclass in subclass_order) {
  for (one_direction in direction_order) {
    dir_tag   <- if (one_direction == "Macular direction") "Central" else "Peripheral"
    dir_label <- if (one_direction == "Macular direction") "Central direction" else "Peripheral direction"
    sub_df <- plot_long %>%
      filter(subclass == one_subclass, direction_class == one_direction)

    pdf_file <- file.path(plot_dir, paste0(
      "2026-10-21_11_AUCell_violin_byCT_",
      subclass_short[[one_subclass]], "_", dir_tag, "_devadult.pdf"))
    stopifnot(!file.exists(pdf_file))   # 不覆盖保护

    n_slots <- nlevels(droplevels(sub_df$stage_tissue))
    pdf_w <- n_slots * 0.44 + 3.4   # 0.44 inch/槽 + 轴/strip/图例边距

    n_pages <- 0
    pdf(pdf_file, width = pdf_w, height = 2.6, useDingbats = FALSE)
    for (one_ct in celltype_order) {
      ct_df <- sub_df %>% filter(sig_celltype == one_ct)
      if (nrow(ct_df) == 0) next

      p <- ggplot(ct_df, aes(x = stage_tissue, y = AUCell_score, fill = tissue)) +
        geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
        facet_grid(~ period_lab, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = cols.tissue, drop = FALSE) +
        scale_x_discrete(drop = TRUE) +
        labs(
          title = paste0(subclass_short[[one_subclass]], " - ",
                         dir_label, " - ", one_ct),
          x = NULL, y = "AUCell score", fill = "Tissue"
        ) +
        theme_bw(base_size = 9) +
        theme(
          axis.text.x        = element_text(angle = 45, hjust = 1,
                                            size = 7.5, color = "black"),
          axis.text.y        = element_text(size = 7.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          strip.background   = element_rect(fill = "gray90", color = "gray70"),
          strip.text         = element_text(face = "bold", size = 9),
          legend.position    = "right",
          legend.key.size    = unit(0.4, "cm"),
          legend.text        = element_text(size = 7.5),
          legend.title       = element_text(size = 8),
          plot.title         = element_text(face = "bold", hjust = 0.5, size = 9),
          plot.margin        = margin(4, 4, 4, 4)
        )
      print(p)
      n_pages <- n_pages + 1
    }
    dev.off()
    cat("written:", basename(pdf_file), "-", n_pages, "pages\n")
  }
}
cat("done.\n")
