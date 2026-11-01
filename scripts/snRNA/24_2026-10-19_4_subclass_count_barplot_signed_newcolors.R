# 测试 2026-09-01_4_subclass_count_barplot_signed.pdf 的新配色版本
# 代码基于 2026-09-01_4_subclass_count_barplots.ipynb CELL 4（signed 上下对称版），
# 仅替换 subclass 配色，其余样式不变；只输出 PDF，不输出 PNG。
#
# 当前版本颜色记录（原版，用于回滚对比）：
#   DM  = #3C5488 (Developmental maintained)
#   DAL = #E64B35 (Developmental aging-lost)
#   AAM = #00A087 (Adult acquired and maintained)
#
# 本次测试新颜色（引自 /project/work/00use_scripts/00_colors.R 的 cols.subclass）：
#   DM = #66C2A5, DAL = #FC8D62, AAM = #8DA0CB
# 输出：Plot/2026-10-19_4_subclass_count_barplot_signed_newcolors.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

date_str <- "2026-10-19_4"
data_dir <- "/project/work/2026-09-01_DEG_summary/Data/"
plot_dir <- "/project/work/2026-09-01_DEG_summary/Plot/"
refined_file <- file.path(data_dir, "2026-09-01_1_five_stage_trajectory_refined.csv")

source("/project/work/00use_scripts/00_colors.R")

celltype_order <- c("Rod", "Cone", "AC", "BC", "HC", "MG")
subclass_order <- c(
  "Developmental maintained",
  "Developmental aging-lost",
  "Adult acquired and maintained"
)
subclass_labels <- c(
  "Developmental maintained" = "DM",
  "Developmental aging-lost" = "DAL",
  "Adult acquired and maintained" = "AAM"
)
cols.subclass.full <- c(
  "Developmental maintained" = unname(cols.subclass["DM"]),
  "Developmental aging-lost" = unname(cols.subclass["DAL"]),
  "Adult acquired and maintained" = unname(cols.subclass["AAM"])
)

# ---- 1. 统计 celltype × direction × subclass 基因数 -------------------------
plot_df <- read_csv(refined_file, show_col_types = FALSE) %>%
  filter(
    celltype %in% celltype_order,
    subclass_refined %in% subclass_order,
    direction_class %in% c("Macular direction", "Peripheral direction")
  ) %>%
  group_by(celltype, direction_class, subclass = subclass_refined) %>%
  summarise(n_gene = n_distinct(gene), .groups = "drop") %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    subclass = factor(subclass, levels = subclass_order)
  )

# ---- 2. signed 上下对称版（Macular 向上 / Peripheral 向下）-------------------
plot_df_signed <- plot_df %>%
  mutate(n_gene_signed = if_else(direction_class == "Macular direction", n_gene, -n_gene))

p_signed <- ggplot(
  plot_df_signed,
  aes(x = celltype, y = n_gene_signed, fill = subclass)
) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(aes(label = n_gene), position = position_stack(vjust = 0.5), size = 3) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  scale_fill_manual(
    values = cols.subclass.full,
    breaks = names(subclass_labels),
    labels = unname(subclass_labels)
  ) +
  labs(
    x = NULL, y = "Number of genes", fill = "Subclass",
    title = "Refined subclass gene counts",
    subtitle = "Up = Macular direction; down = Peripheral direction (new color test)"
  ) +
  theme_classic(base_size = 11)

pdf_file <- file.path(plot_dir, paste0(date_str, "_subclass_count_barplot_signed_newcolors.pdf"))
ggsave(pdf_file, p_signed, width = 8, height = 6, bg = "white")
cat("已输出：", pdf_file, "\n")
