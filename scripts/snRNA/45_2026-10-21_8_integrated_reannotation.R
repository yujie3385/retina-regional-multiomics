#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-21_8  整合对象（retina7）重新注释 + UMAP（用户指定配色, 仅PDF）
# 阶段: markers  —— FindAllMarkers(res 0.4) + 经典marker平均表达 + 聚类构成
#       annotate —— 应用 cluster->celltype 映射, 保存重注释对象, 出图
# 输出:
#   Data/2026-10-21_8_tmp_markers_res0.4.csv
#   Data/2026-10-21_8_canonical_marker_avgexp.csv
#   Data/2026-10-21_8_cluster_composition.csv
#   Data/2026-10-21_8_adult_integrated_retina7_reannotated.rds   (annotate)
#   Data/2026-10-21_8_reannotation_crosstab.csv                  (annotate)
#   Plot/2026-10-21_8_adult_integrated_retina7_UMAP_celltype.pdf (annotate, 单独一张)
#   Plot/2026-10-21_8_adult_integrated_retina7_UMAP_overview.pdf (annotate, 四联)
#   Plot/2026-10-21_8_adult_integrated_retina7_DotPlot_markers.pdf (annotate)
# =====================================================================
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "markers"

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
plot_dir <- file.path(proj_dir, "Plot")
in_rds   <- file.path(data_dir, "2026-10-21_7_adult_integrated_retina7.rds")
anno_res <- "0.4"

log_msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n", sep = " ")

canonical_markers <- c(
  "RHO","NRL","NR2E3","RCVRN","GNAT1","PDE6B",            # Rod
  "ARR3","OPN1SW","GNAT2","RXRG","CNGA3","GUCA1A",        # Cone
  "VSX2","ISL1","CABP5","GRM6","GRM8","PRKCA","LHX3","SCGN", # BC
  "LHX1","PVALB","CALB1","OTX2",                          # HC
  "TFAP2A","TFAP2B","GAD1","GAD2","SLC32A1","SLC6A9","TH","EBF3", # AC
  "RBPMS","POU4F2","ISL2","SNCG","NEFM",                  # RGC
  "RLBP1","GLUL","SOX9","CLDN1","KCNJ10","SLC1A3","AQP4", # MG
  "GFAP","ALDH1L1",                                       # Astro
  "PTPRC","AIF1","CX3CR1","C1QA","C1QB","CSF1R",          # Micro
  "RPE65","BEST1","PMEL","MLANA"                          # RPE 检查
)

## ================= stage: markers =================
if (stage %in% c("markers", "all")) {
  log_msg("读取整合对象 ...")
  obj <- readRDS(in_rds)
  Idents(obj) <- paste0("RNA_snn_res.", anno_res)
  cat("cells:", ncol(obj), "| clusters:", length(levels(Idents(obj))), "\n")

  ## 聚类构成: cluster x 原celltype x dataset
  md <- obj[[]]
  comp <- as.data.frame.matrix(table(Idents(obj), md$celltype_major))
  comp_ds <- as.data.frame.matrix(table(Idents(obj), md$dataset))
  comp_out <- cbind(cluster = rownames(comp), n = comp_ds$GSE155288 + comp_ds$GSE277326,
                    frac_g277 = round(comp_ds$GSE277326 / (comp_ds$GSE155288 + comp_ds$GSE277326), 3),
                    as.data.frame(comp), as.data.frame(comp_ds))
  write.csv(comp_out, file.path(data_dir, "2026-10-21_8_cluster_composition.csv"), row.names = FALSE)
  cat("\n==== cluster x 原celltype 构成 ====\n"); print(comp_out, row.names = FALSE)

  ## FindAllMarkers
  log_msg("FindAllMarkers (res", anno_res, ") ...")
  options(future.globals.maxSize = 24 * 1024^3)
  future::plan("multisession", workers = 12)
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  future::plan("sequential")
  write.csv(markers, file.path(data_dir, "2026-10-21_8_tmp_markers_res0.4.csv"), row.names = FALSE)

  top10 <- markers %>% group_by(.data$cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
  cat("\n==== top10 markers per cluster ====\n")
  print(as.data.frame(top10 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2)), row.names = FALSE)

  ## 经典marker平均表达
  mk <- intersect(canonical_markers, rownames(obj))
  log_msg("canonical markers found:", length(mk), "/", length(canonical_markers))
  avg <- AverageExpression(obj, features = mk, layer = "data")$RNA
  avg_df <- as.data.frame(avg); avg_df$gene <- rownames(avg_df)
  avg_df <- avg_df[, c("gene", setdiff(colnames(avg_df), "gene"))]
  write.csv(avg_df, file.path(data_dir, "2026-10-21_8_canonical_marker_avgexp.csv"), row.names = FALSE)
  log_msg("markers stage done.")
}

## ================= stage: annotate =================
if (stage %in% c("annotate", "all")) {
  ## ---- cluster -> celltype 映射（依据 markers + 原注释构成, 人工核对）----
  cluster_map <- c(
    "0" = "PLACEHOLDER"
  )
  source("/project/work/2026-10-20_other_dataset/scripts/2026-10-21_8_cluster_map.R"
         , local = TRUE, echo = FALSE)

  log_msg("读取整合对象 ...")
  obj <- readRDS(in_rds)
  Idents(obj) <- paste0("RNA_snn_res.", anno_res)
  cl <- as.character(Idents(obj))
  stopifnot(all(cl %in% names(cluster_map)))
  obj$celltype_reanno <- unname(cluster_map[cl])
  ct_levels <- intersect(c("Rod","Cone","BC","HC","AC","RGC","MG","RPC","Micro","Astro"),
                         unique(obj$celltype_reanno))
  obj$celltype_reanno <- factor(obj$celltype_reanno, levels = ct_levels)

  ## 与原注释交叉表
  xt <- table(obj$celltype_reanno, obj$celltype_major)
  write.csv(as.data.frame(xt), file.path(data_dir, "2026-10-21_8_reannotation_crosstab.csv"))
  cat("==== 重注释 x 原注释 ====\n"); print(xt)
  cat("\n重注释构成:\n"); print(table(obj$celltype_reanno))

  out_rds <- file.path(data_dir, "2026-10-21_8_adult_integrated_retina7_reannotated.rds")
  stopifnot(!file.exists(out_rds))
  saveRDS(obj, out_rds)
  log_msg("重注释对象已保存:", basename(out_rds))

  ## ---- 绘图（用户指定配色, 仅PDF）----
  cols.celltype <- c(
    Rod   = "#A6CEE3", Cone  = "#1F78B4", AC    = "#B2DF8A", BC    = "#33A02C",
    HC    = "#FB9A99", MG    = "#E31A1C", RGC   = "#FDBF6F", RPC   = "#FF7F00",
    Micro = "#66C2A5", Astro = "#FC8D62"
  )
  types_present <- levels(factor(obj$celltype_reanno))
  stopifnot(all(types_present %in% names(cols.celltype)))
  cols_use <- cols.celltype[types_present]
  cols_ds <- c(GSE277326 = "#4E79A7", GSE155288 = "#E15759")
  cols_ti <- c(central = "#6A3D9A", peripheral = "#E69F00")
  cols_age <- c("21y" = "#3288BD", "27y" = "#66C2A5", "43y" = "#A6D96A",
                "78y" = "#FDAE61", "90y" = "#D73027")

  ## 单独一张 celltype UMAP
  p_ct <- DimPlot(obj, reduction = "umap", group.by = "celltype_reanno",
                  cols = cols_use, raster = TRUE, raster.dpi = c(300, 300), label = TRUE,
                  label.size = 5, pt.size = 0.1) +
    ggtitle("Integrated adult retina (neurons + MG) - re-annotated") + theme_bw()
  pdf(file.path(plot_dir, "2026-10-21_8_adult_integrated_retina7_UMAP_celltype.pdf"),
      width = 10, height = 8, useDingbats = FALSE)
  print(p_ct)
  dev.off()

  ## 四联 overview（celltype 用重注释 + 指定配色）
  p1 <- DimPlot(obj, reduction = "umap", group.by = "celltype_reanno",
                cols = cols_use, raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle("celltype (re-annotated)") + theme_bw()
  p2 <- DimPlot(obj, reduction = "umap", group.by = "dataset", cols = cols_ds,
                raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("dataset") + theme_bw()
  p3 <- DimPlot(obj, reduction = "umap", group.by = "donor_age", cols = cols_age,
                raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("donor age") + theme_bw()
  p4 <- DimPlot(obj, reduction = "umap", group.by = "tissue", cols = cols_ti,
                raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("tissue") + theme_bw()
  pdf(file.path(plot_dir, "2026-10-21_8_adult_integrated_retina7_UMAP_overview.pdf"),
      width = 15, height = 12, useDingbats = FALSE)
  print((p1 | p2) / (p3 | p4))
  dev.off()

  ## DotPlot 经典marker（按重注释）
  mk <- intersect(canonical_markers, rownames(obj))
  p_dot <- DotPlot(obj, features = mk, group.by = "celltype_reanno") +
    coord_flip() + theme_bw() +
    theme(axis.text.y = element_text(size = 7)) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggtitle("canonical markers by re-annotated celltype")
  pdf(file.path(plot_dir, "2026-10-21_8_adult_integrated_retina7_DotPlot_markers.pdf"),
      width = 10, height = 12, useDingbats = FALSE)
  print(p_dot)
  dev.off()
  log_msg("plots done (PDF only).")
}
