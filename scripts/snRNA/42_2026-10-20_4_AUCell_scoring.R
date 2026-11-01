#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-20_4  三个数据集 AUCell 打分（方法与 2026-09-01_2 完全一致）
#   - 基因集: 31 signature（celltype × DM/DAL/AAM × Macular/Peripheral）
#   - 每数据集独立 ranking; aucMaxRank = ceiling(0.05*nGenes); ties.method="random"; set.seed(1)
#   - 基因匹配: symbol 直接匹配 + ENSG 形式经 Ensembl 列兜底
# 输出:
#   Data/2026-10-20_4_AUCell_GSE277326_meta_data.rds   含 AUCell_* 列
#   Data/2026-10-20_4_AUCell_GSE155288_meta_data.rds
#   Data/2026-10-20_4_AUCell_GSE142526_meta_data.rds
#   Data/2026-10-20_4_AUCell_signature_QC.csv          每signature每数据集检出QC
#   Data/2026-10-20_4_AUCell_parameters_sessionInfo.txt
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(DelayedArray)
  library(matrixStats)
  library(SummarizedExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
})
set.seed(1)

date_str <- "2026-10-20_4"
proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
sig_file <- "/project/work/2026-09-01_DEG_summary/Data/2026-09-01_2_AUCell_signature_definition.csv"

auc_max_rank_pct <- 0.05
n_cores <- 4L

log_msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n", sep = " ")

## ---- 服务器兼容版 ranking（与 2026-09-01_2 一致，勿改）----
.AUCell_buildRankings_new <- function(expr_mat) {
  expr_delayed <- DelayedArray::DelayedArray(-expr_mat)
  rankings <- do.call(cbind, DelayedArray::blockApply(
    expr_delayed,
    FUN = function(x, ...) matrixStats::colRanks(x, ..., useNames = TRUE),
    ties.method = "random", preserveShape = TRUE,
    grid = DelayedArray::colAutoGrid(expr_delayed)
  ))
  rownames(rankings) <- rownames(expr_mat)
  colnames(rankings) <- colnames(expr_mat)
  names(dimnames(rankings)) <- c("genes", "cells")
  new("aucellResults",
      SummarizedExperiment::SummarizedExperiment(assays = list(ranking = rankings)),
      nGenesDetected = numeric(0))
}

get_counts <- function(obj) {
  tryCatch(
    LayerData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

## ---- signature 定义 ----
signature_key <- read.csv(sig_file, stringsAsFactors = FALSE) %>%
  mutate(signature_genes = strsplit(signature_gene_names, ";", fixed = TRUE))
cat("signatures:", nrow(signature_key), "\n")

## =====================================================================
## 数据集 1: GSE277326（snRNA, 注释对象）
## =====================================================================
prep_GSE277326 <- function() {
  obj <- readRDS(file.path(data_dir, "2026-10-20_1_GSE277326_annotated.rds"))
  expr_mat <- get_counts(obj)
  meta_df <- obj[[]]
  rm(obj); gc()
  stopifnot(identical(colnames(expr_mat), rownames(meta_df)))

  ## ENSG 兜底映射：6 个样本 features 应一致，取并集验证
  rna_dir <- file.path(proj_dir, "00data/2026-10-20_adult_new/GSE277326/RNA_matrices")
  feat_files <- list.files(rna_dir, pattern = "\\.features\\.tsv\\.gz$")
  feats <- lapply(feat_files, function(f) {
    con <- gzfile(file.path(rna_dir, f))
    d <- read.delim(con, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("ensembl", "symbol", "type", "chr", "start", "end"))
    tryCatch(close(con), error = function(e) NULL)
    d <- d[d$type == "Gene Expression", c("ensembl", "symbol")]
  })
  ref <- feats[[1]]
  same <- all(sapply(feats[-1], function(x) identical(x, ref)))
  cat("GSE277326 features identical across samples:", same, "\n")
  if (!same) {
    ## 取并集构建映射
    ref <- unique(bind_rows(feats))
  }
  ensg_map <- data.frame(ensembl = ref$ensembl, rowname = make.unique(ref$symbol),
                         stringsAsFactors = FALSE)
  list(expr_mat = expr_mat, meta_df = meta_df, ensg_map = ensg_map,
       meta_std = meta_df %>%
         tibble::rownames_to_column("cell_id") %>%
         transmute(
           cell_id, dataset = "GSE277326",
           sample_id, region, tissue, donor,
           age_group = age_group, celltype = celltype_major
         ))
}

## =====================================================================
## 数据集 2: GSE155288（scRNA, 官方注释）
## =====================================================================
prep_GSE155288 <- function() {
  obj <- readRDS(file.path(data_dir, "2026-10-20_2_GSE155288_merged.rds"))
  expr_mat <- get_counts(obj)
  meta_df <- obj[[]]
  rm(obj); gc()
  stopifnot(identical(colnames(expr_mat), rownames(meta_df)))

  raw_dir <- file.path(proj_dir, "00data/2026-10-20_adult_new/GSE155288/RAW")
  gene_files <- sort(list.files(raw_dir, pattern = "_genes\\.tsv\\.gz$"))
  g1 <- read.delim(gzfile(file.path(raw_dir, gene_files[1])), header = FALSE,
                   col.names = c("ensembl", "symbol"), stringsAsFactors = FALSE)
  ensg_map <- data.frame(ensembl = g1$ensembl, rowname = make.unique(g1$symbol),
                         stringsAsFactors = FALSE)
  list(expr_mat = expr_mat, meta_df = meta_df, ensg_map = ensg_map,
       meta_std = meta_df %>%
         transmute(
           cell_id, dataset = "GSE155288",
           sample_id, region, tissue, donor,
           age_group = age_group, celltype = celltype_major
         ))
}

## =====================================================================
## 数据集 3: GSE142526（scRNA, F6 注释; 直接读矩阵）
## =====================================================================
prep_GSE142526 <- function() {
  dev_dir <- file.path(proj_dir, "00data/2026-10-20_development_new/GSE142526")
  meta_std <- read.csv(file.path(data_dir, "2026-10-20_3_GSE142526_metadata.csv"),
                       stringsAsFactors = FALSE) %>%
    transmute(cell_id, dataset = "GSE142526", sample_id, region, tissue,
              donor = "fetal_donor", age_group = "FD125",
              celltype = celltype_major)

  mats <- list(); ensg_maps <- list()
  for (smp in c("D125C", "D125P")) {
    genes <- read.delim(file.path(dev_dir, smp, "GRCh38", "genes.tsv"), header = FALSE,
                        col.names = c("ensembl", "symbol"), stringsAsFactors = FALSE)
    m <- readMM(file.path(dev_dir, smp, "GRCh38", "matrix.mtx"))
    bc <- readLines(file.path(dev_dir, smp, "GRCh38", "barcodes.tsv"))
    stopifnot(nrow(m) == nrow(genes), ncol(m) == length(bc))
    m <- as(m, "CsparseMatrix")
    rownames(m) <- make.unique(genes$symbol)
    colnames(m) <- paste0(smp, "_", sub("-1$", "", bc))
    mats[[smp]] <- m
    ensg_maps[[smp]] <- data.frame(ensembl = genes$ensembl,
                                   rowname = make.unique(genes$symbol),
                                   stringsAsFactors = FALSE)
  }
  ## 两样本基因注释相同（同一参考），直接拼接
  stopifnot(identical(ensg_maps$D125C, ensg_maps$D125P))
  expr_mat <- cbind(mats$D125C, mats$D125P)
  meta_std <- meta_std[match(colnames(expr_mat), meta_std$cell_id), ]
  stopifnot(!any(is.na(meta_std$celltype)))
  list(expr_mat = expr_mat, meta_df = meta_std, ensg_map = ensg_maps$D125C,
       meta_std = meta_std)
}

## =====================================================================
## 基因匹配 + 打分
## =====================================================================
match_signature_genes <- function(sig_genes, rownames_vec, ensg_map) {
  vapply(sig_genes, function(g) {
    if (g %in% rownames_vec) return(g)
    hit <- ensg_map$rowname[ensg_map$ensembl == g & ensg_map$rowname %in% rownames_vec]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }, character(1), USE.NAMES = FALSE)
}

score_dataset <- function(prep, dataset_name) {
  expr_mat <- prep$expr_mat
  meta_std <- prep$meta_std
  ## 统一列结构（GSE142526 的 meta_df 即 meta_std）
  rownames_vec <- rownames(expr_mat)

  matched <- lapply(signature_key$signature_genes, match_signature_genes,
                    rownames_vec = rownames_vec, ensg_map = prep$ensg_map)
  names(matched) <- signature_key$signature_id
  qc <- signature_key %>%
    select(signature_id, celltype, subclass, direction_class, n_signature_genes) %>%
    mutate(
      dataset = dataset_name,
      n_detected = vapply(matched, function(x) sum(!is.na(x)), integer(1)),
      detected_gene_fraction = n_detected / n_signature_genes,
      missing_genes = vapply(seq_len(nrow(signature_key)), function(i) {
        miss <- signature_key$signature_genes[[i]][is.na(matched[[i]])]
        if (length(miss) == 0) "" else paste(miss, collapse = ";")
      }, character(1))
    )
  gene_sets <- lapply(matched, function(x) unique(na.omit(x)))
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  cat("\n====", dataset_name, "signature QC ====\n")
  print(qc %>% select(signature_id, n_signature_genes, n_detected, detected_gene_fraction) %>%
          as.data.frame(), row.names = FALSE)

  log_msg(dataset_name, "building rankings (", nrow(expr_mat), "genes x", ncol(expr_mat), "cells) ...")
  rankings <- .AUCell_buildRankings_new(expr_mat)
  log_msg(dataset_name, "calculating AUC (aucMaxRank =", ceiling(auc_max_rank_pct * nrow(expr_mat)), ") ...")
  auc <- AUCell::AUCell_calcAUC(
    gene_sets, rankings,
    aucMaxRank = ceiling(auc_max_rank_pct * nrow(expr_mat)), nCores = n_cores
  )
  auc_mat <- t(AUCell::getAUC(auc))
  colnames(auc_mat) <- paste0("AUCell_", colnames(auc_mat))
  cells <- colnames(expr_mat)
  meta_std <- meta_std[match(cells, meta_std$cell_id), ]
  stopifnot(identical(meta_std$cell_id, cells))
  out_meta <- cbind(meta_std, auc_mat[cells, , drop = FALSE])
  saveRDS(out_meta, file.path(data_dir, paste0(date_str, "_AUCell_", dataset_name, "_meta_data.rds")))
  rm(expr_mat, rankings, auc, auc_mat); gc()
  log_msg(dataset_name, "done. cells:", nrow(out_meta))
  qc
}

## ================= main =================
t0 <- Sys.time()
qc_all <- list()

log_msg("preparing GSE277326 ...")
p <- prep_GSE277326()
qc_all$GSE277326 <- score_dataset(p, "GSE277326"); rm(p); gc()

log_msg("preparing GSE155288 ...")
p <- prep_GSE155288()
qc_all$GSE155288 <- score_dataset(p, "GSE155288"); rm(p); gc()

log_msg("preparing GSE142526 ...")
p <- prep_GSE142526()
qc_all$GSE142526 <- score_dataset(p, "GSE142526"); rm(p); gc()

qc_df <- bind_rows(qc_all)
write.csv(qc_df, file.path(data_dir, paste0(date_str, "_AUCell_signature_QC.csv")), row.names = FALSE)

## 汇总: 哪些 signature 在哪些数据集检出不足
low <- qc_df %>% filter(detected_gene_fraction < 0.5)
cat("\n==== 检出比例 < 0.5 的 signature ====\n")
if (nrow(low) == 0) cat("无\n") else print(as.data.frame(low), row.names = FALSE)

writeLines(
  c(
    paste0("date_str: ", date_str),
    paste0("sig_file: ", sig_file),
    paste0("auc_max_rank_pct: ", auc_max_rank_pct),
    "ranking: .AUCell_buildRankings_new, ties.method=random, set.seed(1)",
    "gene matching: symbol direct + ENSG via ensembl column",
    "independent ranking per dataset; scores not comparable across datasets",
    "",
    capture.output(sessionInfo())
  ),
  file.path(data_dir, paste0(date_str, "_AUCell_parameters_sessionInfo.txt"))
)
log_msg("all done. elapsed (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
