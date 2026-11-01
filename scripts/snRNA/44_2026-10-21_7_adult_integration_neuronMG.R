#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-21_7  成年数据跨数据集整合 v2（仅视网膜神经元 + Müller glia）
# 在 2026-10-21_6 基础上按用户要求收紧细胞范围:
#   保留: Rod, Cone, BC, HC, AC, RGC, MG（7 类）
#   去除: Endothelium, Pericytes, RPE_contamination, Other_unresolved,
#         Doublet_BC_MG, Astrocyte, Microglia
# 其余方法与 _6 完全一致（harmony by dataset; 不覆盖任何既有文件）
# 输出:
#   Data/2026-10-21_7_adult_integrated_retina7.rds
#   Data/2026-10-21_7_integration_celltype_dataset_counts.csv
#   Data/2026-10-21_7_integration_cluster_dataset_frac.csv
#   Data/2026-10-21_7_integration_parameters_sessionInfo.txt
#   Plot/2026-10-21_7_adult_integrated_retina7_UMAP_overview.pdf/.png
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
plot_dir <- file.path(proj_dir, "Plot")
out_rds  <- file.path(data_dir, "2026-10-21_7_adult_integrated_retina7.rds")
stopifnot(!file.exists(out_rds))   # 不覆盖保护

KEEP_TYPES <- c("Rod", "Cone", "BC", "HC", "AC", "RGC", "MG")
n_hvg <- 3000
npcs  <- 50
harmony_dims <- 1:30
cluster_res  <- c(0.2, 0.4, 0.6, 0.8)
default_res  <- "0.4"
n_workers    <- 12

log_msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n", sep = " ")

## ---- 1. 读入并按细胞类型收紧 ----
log_msg("读取 GSE277326 ...")
o1 <- readRDS(file.path(data_dir, "2026-10-20_1_GSE277326_annotated.rds"))
log_msg("读取 GSE155288 ...")
o2 <- readRDS(file.path(data_dir, "2026-10-20_2_GSE155288_merged.rds"))

m1 <- o1[[]]; keep1 <- as.character(m1$celltype_major) %in% KEEP_TYPES
m2 <- o2[[]]; keep2 <- as.character(m2$celltype_major) %in% KEEP_TYPES
cat("GSE277326 去除的细胞类型及数量:\n"); print(table(as.character(m1$celltype_major[!keep1])))
cat("GSE155288 去除的细胞类型及数量:\n"); print(table(as.character(m2$celltype_major[!keep2])))
o1 <- subset(o1, cells = rownames(m1)[keep1])
o2 <- subset(o2, cells = rownames(m2)[keep2])
cat("保留后 GSE277326:", ncol(o1), "| GSE155288:", ncol(o2), "\n")

## ---- 2. 共同基因 ----
cg <- intersect(rownames(o1), rownames(o2))
cat("共同基因:", length(cg), "\n")
c1 <- GetAssayData(o1, assay = "RNA", layer = "counts")[cg, ]
c2 <- GetAssayData(o2, assay = "RNA", layer = "counts")[cg, ]

## ---- 3. 统一meta + 细胞名前缀, 重建干净对象 ----
m1 <- o1[[]]; m2 <- o2[[]]
meta1 <- data.frame(
  row.names = colnames(o1),
  dataset = "GSE277326", sample_id = as.character(m1$sample_id),
  region_raw = as.character(m1$region), tissue = as.character(m1$tissue),
  donor_age = paste0(m1$age_years, "y"), age_group = as.character(m1$age_group),
  celltype_major = as.character(m1$celltype_major),
  annotation_source = "self_annotation",
  nCount_RNA = m1$nCount_RNA, nFeature_RNA = m1$nFeature_RNA, percent.mt = m1$percent.mt,
  stringsAsFactors = FALSE)
meta2 <- data.frame(
  row.names = colnames(o2),
  dataset = "GSE155288", sample_id = as.character(m2$sample_id),
  region_raw = as.character(m2$region), tissue = as.character(m2$tissue),
  donor_age = paste0(m2$age_years, "y"), age_group = as.character(m2$age_group),
  celltype_major = as.character(m2$celltype_major),
  annotation_source = "official",
  nCount_RNA = m2$nCount_RNA, nFeature_RNA = m2$nFeature_RNA, percent.mt = m2$percent.mt,
  stringsAsFactors = FALSE)
stopifnot(all(meta1$celltype_major %in% KEEP_TYPES),
          all(meta2$celltype_major %in% KEEP_TYPES),
          !any(duplicated(c(rownames(meta1), rownames(meta2)))))

n1 <- CreateSeuratObject(counts = c1, meta.data = meta1)
n2 <- CreateSeuratObject(counts = c2, meta.data = meta2)
rm(o1, o2, c1, c2, m1, m2); invisible(gc())

log_msg("合并 ...")
merged <- merge(n1, y = n2)
rm(n1, n2); invisible(gc())
merged <- JoinLayers(merged)
merged$donor_age <- factor(merged$donor_age, levels = c("21y","27y","43y","78y","90y"))
merged$celltype_major <- factor(merged$celltype_major, levels = KEEP_TYPES)
cat("合并后:", ncol(merged), "cells x", nrow(merged), "genes\n")
print(table(merged$dataset, merged$celltype_major))

## ---- 4. 标准整合流程（与 _6 一致）----
options(future.globals.maxSize = 24 * 1024^3)
future::plan("multisession", workers = n_workers)
log_msg("NormalizeData ...")
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
log_msg("FindVariableFeatures ...")
merged <- FindVariableFeatures(merged, nfeatures = n_hvg, verbose = FALSE)
log_msg("ScaleData (regress percent.mt) ...")
merged <- ScaleData(merged, vars.to.regress = "percent.mt", verbose = TRUE)
log_msg("RunPCA ...")
merged <- RunPCA(merged, npcs = npcs, verbose = FALSE)

future::plan("sequential")
log_msg("RunHarmony (batch = dataset) ...")
merged <- RunHarmony(merged, group.by.vars = "dataset", reduction.use = "pca",
                     dims.use = seq_len(npcs), plot.convergence = FALSE)
log_msg("RunUMAP ...")
merged <- RunUMAP(merged, reduction = "harmony", dims = harmony_dims, verbose = FALSE)
log_msg("FindNeighbors ...")
merged <- FindNeighbors(merged, reduction = "harmony", dims = harmony_dims, verbose = FALSE)
log_msg("FindClusters ...")
merged <- FindClusters(merged, resolution = cluster_res, verbose = FALSE)
Idents(merged) <- paste0("RNA_snn_res.", default_res)

## ---- 5. 保存 ----
saveRDS(merged, out_rds)
log_msg("整合对象已保存:", basename(out_rds))

## ---- 6. QC ----
ct_tab <- as.data.frame(table(merged$celltype_major, merged$dataset))
colnames(ct_tab)[1:2] <- c("celltype_major", "dataset")
write.csv(ct_tab, file.path(data_dir, "2026-10-21_7_integration_celltype_dataset_counts.csv"),
          row.names = FALSE)
cat("\n==== celltype x dataset ====\n")
print(table(merged$celltype_major, merged$dataset))

mix_ct <- merged[[]] %>% group_by(celltype_major) %>%
  summarise(n = n(), frac_g277 = round(mean(dataset == "GSE277326"), 3), .groups = "drop")
cat("\n==== 每细胞类型 GSE277326 占比（总体期望",
    round(mean(merged$dataset == "GSE277326"), 3), "）====\n")
print(as.data.frame(mix_ct), row.names = FALSE)

mix_cl <- merged[[]] %>% group_by(cluster = RNA_snn_res.0.4) %>%
  summarise(n = n(), frac_g277 = round(mean(dataset == "GSE277326"), 3), .groups = "drop")
write.csv(mix_cl, file.path(data_dir, "2026-10-21_7_integration_cluster_dataset_frac.csv"),
          row.names = FALSE)
cat("\n==== 每聚类(res0.4) GSE277326 占比 ====\n")
print(as.data.frame(mix_cl), row.names = FALSE)

## 同类型细胞跨数据集聚类混合指数
cat("\n==== 同类型细胞混合指数（1=同聚类完全共存）====\n")
cl_frac <- tapply(merged$dataset == "GSE277326", merged$RNA_snn_res.0.4, mean)
md <- merged[[]]
md$cl_frac <- cl_frac[as.character(md$RNA_snn_res.0.4)]
for (ct in KEEP_TYPES) {
  s <- md[md$celltype_major == ct, ]
  a <- mean(s$cl_frac[s$dataset == "GSE277326"]); b <- mean(s$cl_frac[s$dataset == "GSE155288"])
  cat(sprintf("%-4s: 混合指数 = %.3f\n", ct, 1 - abs(a - b)))
}

## ---- 7. UMAP 总览图 ----
cols_ct <- c(Rod = "#1F77B4", Cone = "#FF7F0E", BC = "#2CA02C", AC = "#D62728",
             HC = "#9467BD", MG = "#8C564B", RGC = "#E377C2")
cols_ds <- c(GSE277326 = "#4E79A7", GSE155288 = "#E15759")
cols_ti <- c(central = "#6A3D9A", peripheral = "#E69F00")
cols_age <- c("21y" = "#3288BD", "27y" = "#66C2A5", "43y" = "#A6D96A",
              "78y" = "#FDAE61", "90y" = "#D73027")

p1 <- DimPlot(merged, reduction = "umap", group.by = "celltype_major",
              cols = cols_ct[levels(merged$celltype_major)], raster = TRUE,
              raster.dpi = c(300, 300)) + ggtitle("celltype (neurons + MG)") + theme_bw()
p2 <- DimPlot(merged, reduction = "umap", group.by = "dataset", cols = cols_ds,
              raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("dataset") + theme_bw()
p3 <- DimPlot(merged, reduction = "umap", group.by = "donor_age", cols = cols_age,
              raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("donor age") + theme_bw()
p4 <- DimPlot(merged, reduction = "umap", group.by = "tissue", cols = cols_ti,
              raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("tissue") + theme_bw()
p_all <- (p1 | p2) / (p3 | p4)

pdf_file <- file.path(plot_dir, "2026-10-21_7_adult_integrated_retina7_UMAP_overview.pdf")
pdf(pdf_file, width = 15, height = 12, useDingbats = FALSE)
print(p_all)
dev.off()
ggsave(file.path(plot_dir, "2026-10-21_7_adult_integrated_retina7_UMAP_overview.png"),
       p_all, width = 15, height = 12, dpi = 130)
log_msg("UMAP 总览图已输出")

## ---- 8. 参数与 sessionInfo ----
writeLines(c(
  "pipeline: 2026-10-21_7_adult_integration_neuronMG.R (v2 of 2026-10-21_6)",
  "cell scope: neurons + Muller glia only (Rod, Cone, BC, HC, AC, RGC, MG)",
  "removed: Endothelium, Pericytes, RPE_contamination, Other_unresolved, Doublet_BC_MG, Astrocyte, Microglia",
  paste0("common_genes: ", length(cg)),
  paste0("n_hvg: ", n_hvg, "; npcs: ", npcs, "; harmony_dims: 1:30; cluster_res: ",
         paste(cluster_res, collapse = ","), "; default_res: ", default_res),
  "integration: harmony by dataset (2 levels)",
  "new files only; no existing object overwritten",
  "", capture.output(sessionInfo())
), file.path(data_dir, "2026-10-21_7_integration_parameters_sessionInfo.txt"))
log_msg("all done.")
