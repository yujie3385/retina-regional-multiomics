# 测试 2026-09-01_3_trajectory_alluvial_3subclass.pdf 的新配色版本
# 代码基于 2026-09-01_3_trajectory_alluvial_3subclass.ipynb，仅替换 subclass 配色，其余样式不变
#
# 当前版本颜色记录（原版，用于回滚对比）：
#   cols.subclass.short (原):
#     DM  = #3C5488 (Developmental maintained)
#     DAL = #E64B35 (Developmental aging-lost)
#     AAM = #00A087 (Adult acquired and maintained)
#     Sub-threshold = #999999 (仅 with_subthreshold 版使用)
#   cols.state (不变): C=#6A3D9A, P=#E69F00, N=#BAB0AC, NA=#6B6B6B
#
# 本次测试新颜色：
#   DM  = #66C2A5
#   DAL = #FC8D62
#   AAM = #8DA0CB
# 输出：Plot/2026-10-19_3_trajectory_alluvial_3subclass_newcolors.pdf (+png)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggalluvial)
})

date_str <- "2026-10-19_3"
data_dir <- "/project/work/2026-09-01_DEG_summary/Data/"
plot_dir <- "/project/work/2026-09-01_DEG_summary/Plot/"
refined_file <- file.path(data_dir, "2026-09-01_1_five_stage_trajectory_refined.csv")

celltype_order <- c("Rod", "Cone", "AC", "BC", "HC", "MG")
subclass_order <- c(
  "Developmental maintained",
  "Developmental aging-lost",
  "Adult acquired and maintained",
  "Maintained-effect sub-threshold"
)
state_order <- c("C", "P", "N", "NA")

cols.state <- c("C" = "#6A3D9A", "P" = "#E69F00", "N" = "#BAB0AC", "NA" = "#6B6B6B")
cols.subclass.old <- c(
  "DM" = "#3C5488", "DAL" = "#E64B35", "AAM" = "#00A087", "Sub-threshold" = "#999999"
)

## ---- 新配色（测试版）----
cols.group <- c("DM" = "#66C2A5", "DAL" = "#FC8D62", "AAM" = "#8DA0CB")
cols.subclass.short <- c(cols.group, "Sub-threshold" = "#999999")

# ---- 1. 读取 refined 表并聚合路径 -------------------------------------------
trajectory_df <- read_csv(refined_file, show_col_types = FALSE) %>%
  filter(celltype %in% celltype_order, subclass_refined %in% subclass_order) %>%
  mutate(across(
    c(Early_status, Mid_status, Late_status, Young_status, Old_status),
    ~ {
      v <- as.character(.x)
      v[is.na(v)] <- "NA"
      v[v == "M"] <- "C"
      factor(v, levels = state_order)
    }
  ))

alluvial_count <- trajectory_df %>%
  count(
    Early_status, Mid_status, Late_status, Young_status, Old_status,
    subclass_refined,
    name = "n_genes"
  ) %>%
  mutate(
    subclass_refined = factor(subclass_refined, levels = subclass_order),
    subclass_short = factor(
      recode(as.character(subclass_refined),
        "Developmental maintained" = "DM",
        "Developmental aging-lost" = "DAL",
        "Adult acquired and maintained" = "AAM",
        "Maintained-effect sub-threshold" = "Sub-threshold"
      ),
      levels = c("DM", "DAL", "AAM", "Sub-threshold")
    )
  )
cat("celltype × gene 记录：", sum(alluvial_count$n_genes), "; 路径数：", nrow(alluvial_count), "\n")

# ---- 2. 绘图（仅 3subclass 版）-----------------------------------------------
plot_count <- filter(alluvial_count, subclass_short != "Sub-threshold")

p <- ggplot(
  plot_count,
  aes(
    axis1 = Early_status, axis2 = Mid_status, axis3 = Late_status,
    axis4 = Young_status, axis5 = Old_status, axis6 = subclass_short,
    y = n_genes
  )
) +
  geom_alluvium(aes(fill = subclass_short), width = 0.13, alpha = 0.72) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.13, color = "white", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6) +
  scale_fill_manual(
    values = c(cols.state, cols.subclass.short),
    breaks = names(cols.subclass.short),
    name = "Subclass"
  ) +
  scale_x_discrete(
    limits = c("Early", "Mid", "Late", "Young", "Old", "Subclass"),
    expand = c(0.04, 0.04)
  ) +
  labs(
    y = "Number of celltype-gene records",
    title = "Five-stage regional trajectory of the three focus subclasses (refined)",
    subtitle = "DM / DAL / AAM only; stage axes show regional status (C/P/N/NA) — new color test"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold")
  )

pdf_file <- file.path(plot_dir, paste0(date_str, "_trajectory_alluvial_3subclass_newcolors.pdf"))
ggsave(pdf_file, p, width = 9, height = 5, bg = "white")
ggsave(sub("\\.pdf$", ".png", pdf_file), p, width = 9, height = 5, dpi = 150, bg = "white")
cat("已输出：", pdf_file, "\n")
