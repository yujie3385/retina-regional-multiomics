# 最终整理脚本 10 · Fig4 Step1 区域 chromVAR 差异热图 (v3 + v4) (2026-10-27)
# 来源: 00use_scripts/24_2026-10-18_fig4_step1_regional_chromVAR_heatmap.ipynb (仅 v3/v4 cell + 前置数据构建)
# 产出: 2026-09-16_ATAC_start/Plot/2026-10-18_fig4_step1_chromVAR_heatmap_full_v3.pdf (633 motif × 12 列, candidate 符号)
#       2026-09-16_ATAC_start/Plot/2026-10-18_fig4_step1_chromVAR_heatmap_full_v4.pdf (333 motif × 12 列, 去两期方向冲突 TF)
# 整理说明: 去除 v1/v2 测试绘图; 保留 v2 的行聚类基准计算 (v3 行序校验依赖);
#           Delta chromVAR = Central - Peripheral, 对称色标 ±1.5 截断
# v3: 6 ct × 2 期全 633 motif; candidate 符号 (dev = paired t FDR<0.05 & |delta|>=0.10;
#     adult = 2-donor 方向一致 & |mean_delta|>=0.10); 行聚类 ward.D2
# v4: v3 基础上去掉两期同 celltype 均显著但方向相反的 TF (333 行), 过滤后重聚类
suppressPackageStartupMessages({
  library(data.table); library(ComplexHeatmap); library(circlize); library(grid); library(jsonlite)
})
options(stringsAsFactors = FALSE, width = 150)

DATE_PREFIX <- "2026-10-18"
ROOT <- "/project/work"; OUT <- file.path(ROOT, "2026-09-16_ATAC_start")
DIR_TAB <- file.path(OUT, "Table"); DIR_PLOT <- file.path(OUT, "Plot")
dir.create(DIR_PLOT, recursive = TRUE, showWarnings = FALSE)
DEV_F <- file.path(DIR_TAB, "2026-09-16-2_dev_regional_chromVAR.tsv")
AD_F  <- file.path(DIR_TAB, "2026-09-16-2_adult_regional_chromVAR.tsv")
FM_F  <- file.path(DIR_TAB, "2026-09-16-1_motif_TF_family_map.tsv")
stopifnot(file.exists(DEV_F), file.exists(AD_F), file.exists(FM_F))

## 配色 (color.ipynb 现有色板)
source_color_ipynb <- function(path = file.path(ROOT, "color.ipynb")) {
  j <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  codes <- vapply(j$cells, function(cc)
    if (identical(cc$cell_type, "code")) paste(unlist(cc$source), collapse = "") else "", character(1))
  eval(parse(text = paste(codes, collapse = "\n")), envir = globalenv())
}
invisible(tryCatch(source_color_ipynb(), error = function(e) message("color.ipynb: ", conditionMessage(e))))
stopifnot(exists("cols.celltype"))

Z_LIM <- 1.5     # 对称色标上限; |Delta|>1.5 截断显示最深色
cat("参数就绪 | 色标 [-", Z_LIM, ", ", Z_LIM, "]\n", sep = "")

## ---- 1. 读入两期 chromVAR + family map; 构建长表 ----
dev <- fread(DEV_F)
ad  <- fread(AD_F)
fm  <- fread(FM_F, select = c("motif_id", "motif_name", "family_for_collapse"))
cat("dev:", nrow(dev), "行 |", uniqueN(dev$motif_id), "motif |", uniqueN(dev$celltype), "ct\n")
cat("adult:", nrow(ad), "行 |", uniqueN(ad$motif_id), "motif |", uniqueN(ad$celltype), "ct (无 RPC)\n")

## 长表: period, celltype, motif_id, delta, is_candidate
dv <- dev[, .(period = "Development", celltype, motif_id, motif_name,
              delta, is_cand = candidate)]
da <- ad[,  .(period = "Adult", celltype, motif_id, motif_name,
              delta = mean_delta,
              is_cand = direction_consistent & abs(mean_delta) >= 0.10)]
lng <- rbindlist(list(dv, da))
cat("合并长表:", nrow(lng), "行 | dev cand:", sum(lng$period=="Development" & lng$is_cand),
    "| adult cand:", sum(lng$period=="Adult" & lng$is_cand), "\n")
cat("delta 范围:", paste(round(range(lng$delta), 3), collapse=" ~ "),
    "| |delta|>1.5 的格数:", sum(abs(lng$delta) > 1.5), "(截断显示)\n")

## motif 行序 (全 633, 按 motif_id 排序)
motifs <- fm[order(motif_id)]$motif_id

## ---- 2. v2 行序基准 (ward.D2 on 633x12; 供 v3 行序校验; 不输出 v2 图) ----
CT6 <- c("Rod","Cone","AC","BC","HC","MG")
cols_v2 <- rbind(data.table(celltype = CT6, period = "Development"),
                 data.table(celltype = CT6, period = "Adult"))
col_key_v2 <- paste(cols_v2$celltype, cols_v2$period, sep = "||")
lng2 <- copy(lng)[, ck := paste(celltype, period, sep = "||")]
mm2 <- matrix(NA_real_, nrow = length(motifs), ncol = length(col_key_v2),
              dimnames = list(motifs, col_key_v2))
for (i in seq_len(nrow(lng2))) {
  mi <- lng2$motif_id[i]; ck <- lng2$ck[i]
  if (mi %in% motifs && ck %in% col_key_v2) mm2[mi, ck] <- lng2$delta[i]
}
stopifnot(sum(is.na(mm2)) == 0)
mm2_plot <- mm2
mm2_plot[mm2_plot >  Z_LIM] <-  Z_LIM
mm2_plot[mm2_plot < -Z_LIM] <- -Z_LIM
hc <- hclust(dist(mm2_plot), method = "ward.D2")
row_ord2 <- motifs[hc$order]
cat("行序基准计算完成 (ward.D2, 633 行)\n")

## ---- 3. v3: 恢复 candidate 符号; 6 ct × 2 期全 633 motif; 仅 PDF ----
## 布局: 去 RPC/RGC; 6 ct × 2 期 = 12 列; Period 分组; 行聚类 ward.D2 (行序与 v2 基准一致)
## 符号: dev = 黑盘叠白心小圆环, adult = 黑色实心小点 (cex 0.09/0.05, 相邻行不粘连)
cols_v3 <- rbind(data.table(celltype = CT6, period = "Development"),
                 data.table(celltype = CT6, period = "Adult"))
col_key_v3 <- paste(cols_v3$celltype, cols_v3$period, sep = "||")
pd_v3 <- cols_v3$period
ct_v3 <- cols_v3$celltype

## 矩阵: 633 × 12 (delta) + candidate 逻辑矩阵 (自包含, 由 lng 重建; lng$is_cand 已含两期口径)
lng3 <- copy(lng)[celltype %in% CT6][, ck := paste(celltype, period, sep = "||")]
mm3 <- matrix(NA_real_, nrow = length(motifs), ncol = length(col_key_v3),
              dimnames = list(motifs, col_key_v3))
cand3 <- matrix(FALSE, nrow = length(motifs), ncol = length(col_key_v3),
                dimnames = list(motifs, col_key_v3))
for (i in seq_len(nrow(lng3))) {
  mi <- lng3$motif_id[i]; ck <- lng3$ck[i]
  if (mi %in% motifs && ck %in% col_key_v3) {
    mm3[mi, ck]   <- lng3$delta[i]
    cand3[mi, ck] <- isTRUE(lng3$is_cand[i])
  }
}
cat("v3 矩阵:", paste(dim(mm3), collapse="x"), "| NA 格:", sum(is.na(mm3)),
    "| candidate 格: dev", sum(cand3[, pd_v3=="Development"]),
    "| adult", sum(cand3[, pd_v3=="Adult"]), "\n")
stopifnot(sum(is.na(mm3)) == 0)

## 色标 + 截断: 对称 ±1.5
z_fun3 <- colorRamp2(seq(-Z_LIM, Z_LIM, length.out = 256),
                     colorRampPalette(c("#2166AC", "white", "#B2182B"))(256))
mm3_plot <- mm3
mm3_plot[mm3_plot >  Z_LIM] <-  Z_LIM
mm3_plot[mm3_plot < -Z_LIM] <- -Z_LIM

## 行聚类: 与基准同算法同输入 → 行序校验
hc3 <- hclust(dist(mm3_plot), method = "ward.D2")
row_ord3 <- motifs[hc3$order]
mm3_plot <- mm3_plot[row_ord3, ]
cand3    <- cand3[row_ord3, ]
stopifnot(identical(rownames(mm3_plot), row_ord2))
cat("行聚类完成 (ward.D2, 行序 == v2 基准)\n")

## 列注释: Period 优先分组 + celltype
ta3 <- HeatmapAnnotation(
  `Period`    = factor(pd_v3, levels = c("Development","Adult")),
  `Cell type` = factor(ct_v3, levels = CT6),
  col = list(`Period` = c(Development = "#8DA0CB", Adult = "#FC8D62"),
             `Cell type` = cols.celltype),
  annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(0.3, "cm"))

## 格内符号 (cex 0.09, 约 0.6 行高, 相邻行不粘连)
CELL_CEX <- 0.09
CORE_CEX <- 0.05
cell_fun3 <- function(j, i, x, y, width, height, fill) {
  if (is.na(mm3_plot[i, j])) return(invisible())
  if (cand3[i, j]) {
    grid::grid.points(x, y, pch = 16, gp = gpar(col = "black", cex = CELL_CEX))
    if (pd_v3[j] == "Development")
      grid::grid.points(x, y, pch = 16, gp = gpar(col = "white", cex = CORE_CEX))
  }
}

h3 <- Heatmap(mm3_plot, name = "delta", col = z_fun3, na_col = "grey90",
              top_annotation = ta3,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_order = seq_len(nrow(mm3_plot)),
              column_order = seq_len(ncol(mm3_plot)),
              column_split = factor(pd_v3, levels = c("Development","Adult")),
              show_row_names = FALSE, show_column_names = FALSE,
              row_title = NULL,
              column_gap = unit(1.5, "mm"),
              cell_fun = cell_fun3,
              heatmap_legend_param = list(title = "Delta chromVAR\n(Central -\nPeripheral)",
                                          at = c(-1.5, -0.75, 0, 0.75, 1.5),
                                          labels = c("<= -1.5", "-0.75", "0", "0.75", ">= 1.5")),
              use_raster = TRUE, raster_quality = 5)

lgd3 <- list(
  Legend(title = "Delta chromVAR\n(Central -\nPeripheral)",
         col_fun = z_fun3, at = c(-1.5, -0.75, 0, 0.75, 1.5),
         labels = c("<= -1.5", "-0.75", "0", "0.75", ">= 1.5"),
         legend_height = unit(3, "cm")),
  Legend(title = "Candidate motif",
         labels = c("Dev (FDR<0.05, |D|>=0.10)",
                    "Adult (2-donor consistent, |D|>=0.10)"),
         type = "points", pch = c(21, 16),
         legend_gp = gpar(col = "black", fill = c("white", NA),
                          lwd = c(0.8, NA), cex = 0.6)))

pdf_f3 <- file.path(DIR_PLOT, sprintf("%s_fig4_step1_chromVAR_heatmap_full_v3.pdf", DATE_PREFIX))
cairo_pdf(pdf_f3, width = 6, height = 9)
grid.newpage(); draw(h3, newpage = FALSE, annotation_legend_list = lgd3, merge_legend = TRUE)
dev.off()
cat("saved:", basename(pdf_f3), "|", file.size(pdf_f3), "bytes\n")
cat("[v3 QC] 633x12 | candidate 格: dev", sum(cand3[, pd_v3=="Development"]),
    "| adult", sum(cand3[, pd_v3=="Adult"]), "\n")

## ---- 4. v4: v3 基础上去掉两期同 celltype 均显著但方向相反的 TF; 仅 PDF ----
## 过滤: 同一 celltype 内 dev candidate 且 adult candidate 且 delta 反号 → 该 motif 整行剔除
CT6 <- c("Rod","Cone","AC","BC","HC","MG")

## 冲突检测: 从 lng 取同 celltype 两期 candidate 且 delta 反号的 motif
dv4 <- lng[period == "Development" & celltype %in% CT6,
           .(celltype, motif_id, dev_delta = delta, dev_cand = is_cand)]
da4 <- lng[period == "Adult" & celltype %in% CT6,
           .(celltype, motif_id, ad_delta = delta, ad_cand = is_cand)]
pair4 <- merge(dv4, da4, by = c("celltype", "motif_id"))
pair4[, dev_c := !is.na(dev_cand) & dev_cand]
pair4[, ad_c  := !is.na(ad_cand) & ad_cand]
conflict4 <- pair4[dev_c == TRUE & ad_c == TRUE & sign(dev_delta) != sign(ad_delta)]
conflict_motifs4 <- unique(conflict4$motif_id)
motifs_v4 <- setdiff(motifs, conflict_motifs4)   # 保持原 motif_id 行序
cat("两期均 candidate 的 (ct,motif):", nrow(pair4[dev_c == TRUE & ad_c == TRUE]),
    "| 其中方向冲突:", nrow(conflict4),
    "| 剔除 motif:", length(conflict_motifs4), "| 保留:", length(motifs_v4), "行\n")
print(conflict4[, .(n_conflict = .N), by = celltype][order(celltype)])

## 矩阵: 333 × 12 (delta) + candidate 逻辑矩阵 (同 v3 构造, 行 = motifs_v4)
cols_v4 <- rbind(data.table(celltype = CT6, period = "Development"),
                 data.table(celltype = CT6, period = "Adult"))
col_key_v4 <- paste(cols_v4$celltype, cols_v4$period, sep = "||")
pd_v4 <- cols_v4$period
ct_v4 <- cols_v4$celltype
lng4 <- copy(lng)[celltype %in% CT6][, ck := paste(celltype, period, sep = "||")]
mm4 <- matrix(NA_real_, nrow = length(motifs_v4), ncol = length(col_key_v4),
              dimnames = list(motifs_v4, col_key_v4))
cand4 <- matrix(FALSE, nrow = length(motifs_v4), ncol = length(col_key_v4),
                dimnames = list(motifs_v4, col_key_v4))
for (i in seq_len(nrow(lng4))) {
  mi <- lng4$motif_id[i]; ck <- lng4$ck[i]
  if (mi %in% motifs_v4 && ck %in% col_key_v4) {
    mm4[mi, ck]   <- lng4$delta[i]
    cand4[mi, ck] <- isTRUE(lng4$is_cand[i])
  }
}
cat("v4 矩阵:", paste(dim(mm4), collapse="x"), "| NA 格:", sum(is.na(mm4)),
    "| candidate 格: dev", sum(cand4[, pd_v4=="Development"]),
    "| adult", sum(cand4[, pd_v4=="Adult"]), "\n")
stopifnot(sum(is.na(mm4)) == 0, !any(conflict_motifs4 %in% rownames(mm4)))

## 色标 + 截断 (同 v3): 对称 ±1.5
z_fun4 <- colorRamp2(seq(-Z_LIM, Z_LIM, length.out = 256),
                     colorRampPalette(c("#2166AC", "white", "#B2182B"))(256))
mm4_plot <- mm4
mm4_plot[mm4_plot >  Z_LIM] <-  Z_LIM
mm4_plot[mm4_plot < -Z_LIM] <- -Z_LIM

## 行聚类: 过滤后矩阵重新 ward.D2 (333 行)
hc4 <- hclust(dist(mm4_plot), method = "ward.D2")
row_ord4 <- motifs_v4[hc4$order]
mm4_plot <- mm4_plot[row_ord4, ]
cand4    <- cand4[row_ord4, ]
cat("行聚类完成 (ward.D2,", nrow(mm4_plot), "行, 过滤后重聚类)\n")

## 列注释: 同 v3 (Period 优先分组 + celltype)
ta4 <- HeatmapAnnotation(
  `Period`    = factor(pd_v4, levels = c("Development","Adult")),
  `Cell type` = factor(ct_v4, levels = CT6),
  col = list(`Period` = c(Development = "#8DA0CB", Adult = "#FC8D62"),
             `Cell type` = cols.celltype),
  annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(0.3, "cm"))

## 格内符号 (与 v3 同一物理尺寸 cex 0.09, 便于 v3/v4 直接对比)
CELL_CEX4 <- 0.09
CORE_CEX4 <- 0.05
cell_fun4 <- function(j, i, x, y, width, height, fill) {
  if (is.na(mm4_plot[i, j])) return(invisible())
  if (cand4[i, j]) {
    grid::grid.points(x, y, pch = 16, gp = gpar(col = "black", cex = CELL_CEX4))
    if (pd_v4[j] == "Development")
      grid::grid.points(x, y, pch = 16, gp = gpar(col = "white", cex = CORE_CEX4))
  }
}

h4 <- Heatmap(mm4_plot, name = "delta", col = z_fun4, na_col = "grey90",
              top_annotation = ta4,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_order = seq_len(nrow(mm4_plot)),
              column_order = seq_len(ncol(mm4_plot)),
              column_split = factor(pd_v4, levels = c("Development","Adult")),
              show_row_names = FALSE, show_column_names = FALSE,
              row_title = NULL,
              column_gap = unit(1.5, "mm"),
              cell_fun = cell_fun4,
              heatmap_legend_param = list(title = "Delta chromVAR\n(Central -\nPeripheral)",
                                          at = c(-1.5, -0.75, 0, 0.75, 1.5),
                                          labels = c("<= -1.5", "-0.75", "0", "0.75", ">= 1.5")),
              use_raster = TRUE, raster_quality = 5)

lgd4 <- list(
  Legend(title = "Delta chromVAR\n(Central -\nPeripheral)",
         col_fun = z_fun4, at = c(-1.5, -0.75, 0, 0.75, 1.5),
         labels = c("<= -1.5", "-0.75", "0", "0.75", ">= 1.5"),
         legend_height = unit(3, "cm")),
  Legend(title = "Candidate motif",
         labels = c("Dev (FDR<0.05, |D|>=0.10)",
                    "Adult (2-donor consistent, |D|>=0.10)"),
         type = "points", pch = c(21, 16),
         legend_gp = gpar(col = "black", fill = c("white", NA),
                          lwd = c(0.8, NA), cex = 0.6)))

pdf_f4 <- file.path(DIR_PLOT, sprintf("%s_fig4_step1_chromVAR_heatmap_full_v4.pdf", DATE_PREFIX))
cairo_pdf(pdf_f4, width = 6, height = 9)
grid.newpage(); draw(h4, newpage = FALSE, annotation_legend_list = lgd4, merge_legend = TRUE)
dev.off()
cat("saved:", basename(pdf_f4), "|", file.size(pdf_f4), "bytes\n")
cat("[v4 QC]", nrow(mm4_plot), "x12 | 冲突 (ct,motif) 组合", nrow(conflict4),
    "| candidate 格: dev", sum(cand4[, pd_v4=="Development"]),
    "| adult", sum(cand4[, pd_v4=="Adult"]), "\n")

cat("done: 2 chromVAR heatmap PDFs\n")
