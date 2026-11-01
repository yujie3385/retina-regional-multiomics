#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-20_1  GSE277326 snRNA 注释 pipeline（读入 + QC + 聚类 + marker）
# 数据: 10x multiome GEX（Macula/Retina × 21y/27y/43y 供体; RPE 已排除）
# 平台: Seurat 5.2.1 / harmony 1.2.3
# 用法: Rscript 2026-10-20_1_GSE277326_annotation_pipeline.R <read|cluster|markers|all>
# 输出:
#   Data/2026-10-20_1_GSE277326_QC_before_filter.csv   QC汇总(过滤前)
#   Data/2026-10-20_1_GSE277326_QC_after_filter.csv    QC汇总(过滤后)
#   Data/2026-10-20_1_tmp_merged_qc.rds                过滤后合并对象(检查点)
#   Data/2026-10-20_1_tmp_clustered.rds                聚类后对象(检查点)
#   Data/2026-10-20_1_GSE277326_markers_res0.8.csv     全部marker
#   Data/2026-10-20_1_GSE277326_canonical_marker_avgexp.csv  经典marker平均表达
#   Data/2026-10-20_1_GSE277326_cluster_sample_counts.csv    cluster×sample计数
# =====================================================================
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 1) args[1] else "all"

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
rna_dir  <- file.path(proj_dir, "00data/2026-10-20_adult_new/GSE277326/RNA_matrices")
stopifnot(dir.exists(rna_dir))

## ---------------- 参数 ----------------
qc_min_feature <- 300
qc_max_feature <- 8000
qc_max_mt      <- 5
n_hvg          <- 3000
npcs           <- 50
harmony_dims   <- 1:30
cluster_res    <- seq(0.2, 1.2, by = 0.2)
default_res    <- "0.8"
n_workers      <- 12

samples <- data.frame(
  gsm = c("GSM8519431", "GSM8519432", "GSM8519437", "GSM8519433", "GSM8519434", "GSM8519438"),
  prefix = c("QY_1840_1_2_QY_1839_1_2", "QY_1952_1_2_QY_1951_1_2", "QY_1993_1_2_QY_1992_1_2",
             "QY_1954_1_2_QY_1953_1_2", "QY_1956_1_2_QY_1955_1_2", "QY_1995_1_2_QY_1994_1_2"),
  region = c("Macula", "Macula", "Macula", "Retina", "Retina", "Retina"),
  donor  = c("donor1", "donor2", "donor3", "donor1", "donor2", "donor3"),
  age_years = c(21, 27, 43, 21, 27, 43),
  stringsAsFactors = FALSE
)
samples$sample_id <- paste0(samples$region, "_", samples$age_years, "y")
samples$age_group <- ifelse(samples$age_years <= 39, "Young_20_39y", "40y_plus")
samples$tissue    <- ifelse(samples$region == "Macula", "central", "peripheral")

log_msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n", sep = " ")

get_counts <- function(obj) {
  tryCatch(
    LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

## ================= stage: read + QC =================
run_read <- function() {
  so_list <- list()
  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]
    base <- file.path(rna_dir, paste0(s$gsm, "_", s$prefix))
    log_msg("reading", s$sample_id, "(", i, "/", nrow(samples), ") ...")
    feats_con <- gzfile(paste0(base, ".features.tsv.gz"))
    feats <- read.delim(feats_con, header = FALSE,
                        col.names = c("ensembl", "symbol", "type", "chr", "start", "end"),
                        stringsAsFactors = FALSE)
    tryCatch(close(feats_con), error = function(e) NULL)
    keep <- feats$type == "Gene Expression"
    mtx_con <- gzfile(paste0(base, ".matrix.mtx.gz"))
    mtx <- readMM(mtx_con)
    tryCatch(close(mtx_con), error = function(e) NULL)
    bc_con <- gzfile(paste0(base, ".barcodes.tsv.gz"))
    bc <- readLines(bc_con)
    tryCatch(close(bc_con), error = function(e) NULL)
    stopifnot(nrow(mtx) == nrow(feats), ncol(mtx) == length(bc))
    m <- as(mtx[keep, ], "CsparseMatrix")
    rownames(m) <- make.unique(feats$symbol[keep])
    colnames(m) <- paste0(s$sample_id, "_", bc)
    n_cells <- ncol(m)
    so <- CreateSeuratObject(
      counts = m,
      meta.data = data.frame(
        sample_id = rep(s$sample_id, n_cells), gsm = rep(s$gsm, n_cells),
        region = rep(s$region, n_cells), donor = rep(s$donor, n_cells),
        age_years = rep(s$age_years, n_cells), age_group = rep(s$age_group, n_cells),
        tissue = rep(s$tissue, n_cells),
        row.names = colnames(m)
      )
    )
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    so_list[[s$sample_id]] <- so
    rm(mtx, m); gc()
  }

  ## ---- QC 汇总（过滤前）----
  qc_before <- bind_rows(lapply(so_list, function(so) {
    md <- so[[]]
    data.frame(
      sample_id = md$sample_id[1], n_cells = nrow(md),
      med_feature = round(median(md$nFeature_RNA)),
      q05_feature = round(quantile(md$nFeature_RNA, 0.05)),
      q95_feature = round(quantile(md$nFeature_RNA, 0.95)),
      max_feature = max(md$nFeature_RNA),
      med_count = round(median(md$nCount_RNA)),
      med_mt = round(median(md$percent.mt), 2),
      q95_mt = round(quantile(md$percent.mt, 0.95), 2)
    )
  }))
  cat("==== QC before filter ====\n"); print(qc_before, row.names = FALSE)
  write.csv(qc_before, file.path(data_dir, "2026-10-20_1_GSE277326_QC_before_filter.csv"), row.names = FALSE)

  ## ---- 细胞过滤 ----
  so_list <- lapply(so_list, function(so) {
    md <- so[[]]
    keep <- md$nFeature_RNA >= qc_min_feature & md$nFeature_RNA <= qc_max_feature & md$percent.mt <= qc_max_mt
    log_msg(so$sample_id[1], ":", nrow(md), "->", sum(keep), "cells")
    subset(so, cells = rownames(md)[keep])
  })

  qc_after <- bind_rows(lapply(so_list, function(so) {
    md <- so[[]]
    data.frame(sample_id = md$sample_id[1], n_cells = nrow(md),
               med_feature = round(median(md$nFeature_RNA)), med_mt = round(median(md$percent.mt), 2))
  }))
  cat("==== QC after filter ====\n"); print(qc_after, row.names = FALSE)
  write.csv(qc_after, file.path(data_dir, "2026-10-20_1_GSE277326_QC_after_filter.csv"), row.names = FALSE)

  merged <- merge(so_list[[1]], y = so_list[-1])
  rm(so_list); gc()
  merged <- JoinLayers(merged)

  ## ---- 基因过滤: >=3 细胞检出 ----
  cm <- get_counts(merged)
  genes_keep <- rownames(cm)[Matrix::rowSums(cm > 0) >= 3]
  merged <- subset(merged, features = genes_keep)
  log_msg("merged:", ncol(merged), "cells x", nrow(merged), "genes")
  rm(cm); gc()

  saveRDS(merged, file.path(data_dir, "2026-10-20_1_tmp_merged_qc.rds"))
  log_msg("read stage done.")
  merged
}

## ================= stage: cluster =================
run_cluster <- function(merged) {
  options(future.globals.maxSize = 20 * 1024^3)
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
  log_msg("RunHarmony (by sample_id) ...")
  merged <- RunHarmony(merged, group.by.vars = "sample_id", reduction.use = "pca",
                       dims.use = seq_len(npcs), plot.convergence = FALSE)
  log_msg("RunUMAP ...")
  merged <- RunUMAP(merged, reduction = "harmony", dims = harmony_dims, verbose = FALSE)
  log_msg("FindNeighbors ...")
  merged <- FindNeighbors(merged, reduction = "harmony", dims = harmony_dims, verbose = FALSE)
  log_msg("FindClusters (res", paste(cluster_res, collapse = ","), ") ...")
  merged <- FindClusters(merged, resolution = cluster_res, verbose = FALSE)
  Idents(merged) <- paste0("RNA_snn_res.", default_res)

  saveRDS(merged, file.path(data_dir, "2026-10-20_1_tmp_clustered.rds"))
  log_msg("cluster stage done.")
  merged
}

## ================= stage: markers =================
canonical_markers <- c(
  # Rod
  "RHO", "NRL", "NR2E3", "RCVRN", "GNAT1", "PDE6B",
  # Cone
  "ARR3", "OPN1SW", "GNAT2", "RXRG", "CNGA3", "GUCA1A",
  # BC
  "VSX2", "ISL1", "CABP5", "GRM6", "GRM8", "PRKCA", "LHX3", "SCGN",
  # HC
  "LHX1", "PVALB", "CALB1", "OTX2",
  # AC
  "TFAP2A", "TFAP2B", "GAD1", "GAD2", "SLC32A1", "SLC6A9", "TH", "EBF3",
  # RGC
  "RBPMS", "POU4F2", "ISL2", "SNCG", "NEFM",
  # MG
  "RLBP1", "GLUL", "SOX9", "CLDN1", "KCNJ10", "SLC1A3", "AQP4",
  # Astro
  "GFAP", "ALDH1L1",
  # Oligo / Schwann
  "SOX10", "MPZ", "PLP1", "S100B",
  # Microglia / immune
  "PTPRC", "AIF1", "CX3CR1", "C1QA", "C1QB", "CSF1R", "P2RY12", "CD68",
  # Pericyte / VSMC
  "PDGFRB", "RGS5", "NOTCH3", "ACTA2",
  # Endothelial
  "PECAM1", "VWF", "KDR", "CLDN5", "FLT1",
  # Fibroblast
  "COL1A1", "COL1A2", "DCN", "LUM",
  # T / NK
  "TRAC", "CD3E", "CD3D", "NKG7", "GNLY",
  # Mast
  "TPSAB1", "CPA3", "KIT",
  # RPE / melanocyte 污染检查
  "RPE65", "BEST1", "PMEL", "MLANA", "TYRP1", "DCT"
)

run_markers <- function() {
  merged <- readRDS(file.path(data_dir, "2026-10-20_1_tmp_clustered.rds"))
  Idents(merged) <- paste0("RNA_snn_res.", default_res)
  log_msg("FindAllMarkers (res", default_res, ") ...")
  options(future.globals.maxSize = 20 * 1024^3)
  future::plan("multisession", workers = n_workers)
  markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  future::plan("sequential")
  write.csv(markers, file.path(data_dir, "2026-10-20_1_GSE277326_markers_res0.8.csv"), row.names = FALSE)
  log_msg("markers written:", nrow(markers), "rows")

  ## 每个 cluster top10 marker 打到日志
  top10 <- markers %>%
    group_by(.data$cluster) %>%
    slice_max(avg_log2FC, n = 10) %>%
    ungroup()
  cat("\n==== top10 markers per cluster ====\n")
  print(as.data.frame(top10 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2)), row.names = FALSE)

  ## 经典 marker 平均表达
  mk <- intersect(canonical_markers, rownames(merged))
  log_msg("canonical markers found:", length(mk), "/", length(canonical_markers))
  avg <- AverageExpression(merged, features = mk, layer = "data")$RNA
  avg_df <- as.data.frame(avg)
  avg_df$gene <- rownames(avg_df)
  avg_df <- avg_df[, c("gene", setdiff(colnames(avg_df), "gene"))]
  write.csv(avg_df, file.path(data_dir, "2026-10-20_1_GSE277326_canonical_marker_avgexp.csv"), row.names = FALSE)

  ## cluster × sample 计数
  tab <- as.data.frame.matrix(table(Idents(merged), merged$sample_id))
  tab$cluster <- rownames(tab)
  write.csv(tab, file.path(data_dir, "2026-10-20_1_GSE277326_cluster_sample_counts.csv"), row.names = FALSE)
  cat("\n==== cluster x sample counts ====\n"); print(tab, row.names = FALSE)

  ## 参数与 sessionInfo
  writeLines(
    c(
      "pipeline: 2026-10-20_1_GSE277326_annotation_pipeline.R",
      paste0("samples: ", paste(samples$sample_id, collapse = ", ")),
      paste0("qc_min_feature: ", qc_min_feature, "; qc_max_feature: ", qc_max_feature, "; qc_max_mt: ", qc_max_mt),
      paste0("n_hvg: ", n_hvg, "; npcs: ", npcs, "; harmony_dims: 1:30"),
      paste0("cluster_res: ", paste(cluster_res, collapse = ","), "; default_res: ", default_res),
      paste0("integration: harmony by sample_id"),
      "",
      capture.output(sessionInfo())
    ),
    file.path(data_dir, "2026-10-20_1_annotation_parameters_sessionInfo.txt")
  )
  log_msg("markers stage done.")
}

## ================= main =================
t_start <- Sys.time()
if (stage %in% c("read", "all")) {
  merged <- run_read()
  if (stage == "all") merged <- run_cluster(merged)
} else if (stage == "cluster") {
  merged <- readRDS(file.path(data_dir, "2026-10-20_1_tmp_merged_qc.rds"))
  merged <- run_cluster(merged)
} else if (stage == "markers") {
  run_markers()
}
log_msg("total elapsed (min):", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1))
