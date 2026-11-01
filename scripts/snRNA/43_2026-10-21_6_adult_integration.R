#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-21_6  成年数据跨数据集整合（GSE277326 snRNA 21–43y + GSE155288 scRNA 78–90y）
# - 平台: snRNA(multiome GEX) vs scRNA —— 以 dataset 为批次变量做 harmony 校正
#   （保留供体/样本水平生物学变异, 不过度校正, 便于后续检查21y异常）
# - 基因: 两对象共同基因（~19.6k）
# - 注释: 各自保留（GSE277326=自行注释, GSE155288=官方）; Astrocytes→Astrocyte 统一
# - 细胞名加数据集前缀（g277_/g155_）避免barcode碰撞
# - 输出全部为新文件, 不覆盖任何既有对象
# 输出:
#   Data/2026-10-21_6_adult_integrated.rds
#   Data/2026-10-21_6_integration_celltype_dataset_counts.csv
#   Data/2026-10-21_6_integration_cluster_dataset_frac.csv
#   Data/2026-10-21_6_integration_parameters_sessionInfo.txt
#   Plot/2026-10-21_6_adult_integrated_UMAP_overview.pdf/.png
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
out_rds  <- file.path(data_dir, "2026-10-21_6_adult_integrated.rds")
stopifnot(!file.exists(out_rds))   # 不覆盖保护

n_hvg <- 3000
npcs  <- 50
harmony_dims <- 1:30
cluster_res  <- c(0.2, 0.4, 0.6, 0.8)
default_res  <- "0.4"
n_workers    <- 12

log_msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n", sep = " ")

## ---- 1. 读入两个对象 ----
log_msg("读取 GSE277326 (snRNA, 21-43y)...")
o1 <- readRDS(file.path(data_dir, "2026-10-20_1_GSE277326_annotated.rds"))
log_msg("读取 GSE155288 (scRNA, 78-90y)...")
o2 <- readRDS(file.path(data_dir, "2026-10-20_2_GSE155288_merged.rds"))
cat("GSE277326:", ncol(o1), "cells x", nrow(o1), "genes\n")
cat("GSE155288:", ncol(o2), "cells x", nrow(o2), "genes\n")

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
  celltype_major = sub("Astrocytes$", "Astrocyte", as.character(m1$celltype_major)),
  annotation_source = "self_annotation",
  nCount_RNA = m1$nCount_RNA, nFeature_RNA = m1$nFeature_RNA, percent.mt = m1$percent.mt,
  stringsAsFactors = FALSE)
meta2 <- data.frame(
  row.names = colnames(o2),
  dataset = "GSE155288", sample_id = as.character(m2$sample_id),
  region_raw = as.character(m2$region), tissue = as.character(m2$tissue),
  donor_age = paste0(m2$age_years, "y"), age_group = as.character(m2$age_group),
  celltype_major = sub("Astrocytes$", "Astrocyte", as.character(m2$celltype_major)),
  annotation_source = "official",
  nCount_RNA = m2$nCount_RNA, nFeature_RNA = m2$nFeature_RNA, percent.mt = m2$percent.mt,
  stringsAsFactors = FALSE)
stopifnot(!any(duplicated(c(rownames(meta1), rownames(meta2)))))

n1 <- CreateSeuratObject(counts = c1, meta.data = meta1)
n2 <- CreateSeuratObject(counts = c2, meta.data = meta2)
rm(o1, o2, c1, c2, m1, m2); invisible(gc())

log_msg("合并...")
merged <- merge(n1, y = n2)
rm(n1, n2); invisible(gc())
merged <- JoinLayers(merged)
merged$donor_age <- factor(merged$donor_age, levels = c("21y","27y","43y","78y","90y"))
cat("合并后:", ncol(merged), "cells x", nrow(merged), "genes\n")
print(table(merged$dataset, merged$celltype_major))

## QC 概览（各数据集 nFeature / mt 分布）
cat("\n各数据集 QC 概览:\n")
print(merged[[]] %>% group_by(dataset) %>%
  summarise(n_cells = n(), med_feature = median(nFeature_RNA),
            q05_feature = quantile(nFeature_RNA, .05), q95_feature = quantile(nFeature_RNA, .95),
            med_mt = round(median(percent.mt), 2), q95_mt = round(quantile(percent.mt, .95), 2),
            .groups = "drop") %>% as.data.frame(), row.names = FALSE)

## ---- 4. 标准整合流程 ----
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

## ---- 5. 保存整合对象 ----
saveRDS(merged, out_rds)
log_msg("整合对象已保存:", basename(out_rds))

## ---- 6. 整合质量 QC ----
## 6a. celltype x dataset 计数
ct_tab <- as.data.frame(table(merged$celltype_major, merged$dataset))
colnames(ct_tab)[1:2] <- c("celltype_major", "dataset")
write.csv(ct_tab, file.path(data_dir, "2026-10-21_6_integration_celltype_dataset_counts.csv"),
          row.names = FALSE)
cat("\n==== celltype x dataset ====\n")
print(table(merged$celltype_major, merged$dataset))

## 6b. 每细胞类型的数据集混合度（g277占比; 总体期望 44,712/(44,712+92,385) = 32.6%）
mix_ct <- merged[[]] %>% group_by(celltype_major) %>%
  summarise(n = n(), frac_g277 = round(mean(dataset == "GSE277326"), 3), .groups = "drop")
cat("\n==== 每细胞类型 GSE277326 占比（总体期望 0.326; 接近则混合良好）====\n")
print(as.data.frame(mix_ct), row.names = FALSE)

## 6c. 每个harmony聚类的数据集占比（混合度）
mix_cl <- merged[[]] %>% group_by(cluster = RNA_snn_res.0.4) %>%
  summarise(n = n(), frac_g277 = round(mean(dataset == "GSE277326"), 3), .groups = "drop")
write.csv(mix_cl, file.path(data_dir, "2026-10-21_6_integration_cluster_dataset_frac.csv"),
          row.names = FALSE)
cat("\n==== 每聚类(res0.4) GSE277326 占比 ====\n")
print(as.data.frame(mix_cl), row.names = FALSE)

## ---- 7. UMAP 总览图 ----
cols_ct <- c(Rod = "#1F77B4", Cone = "#FF7F0E", BC = "#2CA02C", AC = "#D62728",
             HC = "#9467BD", MG = "#8C564B", RGC = "#E377C2", Astrocyte = "#17BECF",
             Microglia = "#7F7F7F", Pericytes = "#8C6D31", Endothelium = "#BC80BD",
             Doublet_BC_MG = "#FFD700", RPE_contamination = "#BCBD22", Other_unresolved = "#000000")
cols_ds <- c(GSE277326 = "#4E79A7", GSE155288 = "#E15759")
cols_ti <- c(central = "#6A3D9A", peripheral = "#E69F00")
cols_age <- c("21y" = "#3288BD", "27y" = "#66C2A5", "43y" = "#A6D96A",
              "78y" = "#FDAE61", "90y" = "#D73027")

p1 <- DimPlot(merged, reduction = "umap", group.by = "celltype_major",
              cols = cols_ct[levels(factor(merged$celltype_major))], raster = TRUE,
              raster.dpi = c(300, 300)) + ggtitle("celltype (harmonized)") + theme_bw()
p2 <- DimPlot(merged, reduction = "umap", group.by = "dataset", cols = cols_ds,
              raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("dataset") + theme_bw()
p3 <- DimPlot(merged, reduction = "umap", group.by = "donor_age", cols = cols_age,
              raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("donor age") + theme_bw()
p4 <- DimPlot(merged, reduction = "umap", group.by = "tissue", cols = cols_ti,
              raster = TRUE, raster.dpi = c(300, 300)) + ggtitle("tissue") + theme_bw()
p_all <- (p1 | p2) / (p3 | p4)

pdf_file <- file.path(plot_dir, "2026-10-21_6_adult_integrated_UMAP_overview.pdf")
pdf(pdf_file, width = 15, height = 12, useDingbats = FALSE)
print(p_all)
dev.off()
ggsave(file.path(plot_dir, "2026-10-21_6_adult_integrated_UMAP_overview.png"),
       p_all, width = 15, height = 12, dpi = 130)
log_msg("UMAP 总览图已输出")

## ---- 8. 参数与 sessionInfo ----
writeLines(c(
  "pipeline: 2026-10-21_6_adult_integration.R",
  "inputs: Data/2026-10-20_1_GSE277326_annotated.rds (snRNA, self-annotated) | Data/2026-10-20_2_GSE155288_merged.rds (scRNA, official)",
  paste0("common_genes: ", length(cg)),
  paste0("n_hvg: ", n_hvg, "; npcs: ", npcs, "; harmony_dims: 1:30; cluster_res: ",
         paste(cluster_res, collapse = ","), "; default_res: ", default_res),
  "integration: harmony by dataset (2 levels) — removes snRNA/scRNA platform shift, preserves donor/sample biology",
  "cell names prefixed g277_/g155_; Astrocytes unified to Astrocyte",
  "new files only; no existing object overwritten",
  "", capture.output(sessionInfo())
), file.path(data_dir, "2026-10-21_6_integration_parameters_sessionInfo.txt"))
log_msg("all done.")
