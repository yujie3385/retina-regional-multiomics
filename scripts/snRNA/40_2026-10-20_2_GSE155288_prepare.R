#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-20_2  GSE155288 老年组准备（78y/90y, scRNA）
# 读入16样本RAW counts -> join官方注释 -> 保存merged对象 + SingleR参考集
# 输出:
#   Data/2026-10-20_2_GSE155288_merged.rds          注释后合并对象(仅注释细胞)
#   Data/2026-10-20_2_GSE155288_metadata.csv        细胞级metadata
#   Data/2026-10-20_2_GSE155288_SingleR_ref.rds     SingleR参考(每类≤500细胞)
#   Data/2026-10-20_2_GSE155288_join_QC.csv          join汇总
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
raw_dir  <- file.path(proj_dir, "00data/2026-10-20_adult_new/GSE155288/RAW")
ident_file <- file.path(proj_dir, "00data/2026-10-20_adult_new/GSE155288/GSE155288_ident.csv.gz")

log_msg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n", sep = " ")

## ---- 样本信息解析 ----
bc_files <- sort(list.files(raw_dir, pattern = "_barcodes\\.tsv\\.gz$"))
sample_ids <- sub("^GSM[0-9]+_", "", sub("_barcodes\\.tsv\\.gz$", "", bc_files))
sample_info <- data.frame(
  sample_id = sample_ids,
  region = ifelse(substr(sample_ids, 1, 1) == "M", "Macula", "Periphery"),
  donor  = ifelse(as.integer(substr(sample_ids, 2, 2)) == 1, "donor78y", "donor90y"),
  age_years = ifelse(as.integer(substr(sample_ids, 2, 2)) == 1, 78, 90),
  stringsAsFactors = FALSE
)
sample_info$age_group <- "Old_78_90y"
sample_info$tissue <- ifelse(sample_info$region == "Macula", "central", "peripheral")
stopifnot(nrow(sample_info) == 16)

## ---- 官方注释 ----
ident <- read.csv(gzfile(ident_file), stringsAsFactors = FALSE)
log_msg("ident rows:", nrow(ident), "; cell types:", length(unique(ident$cell_type)))
tryCatch(close(gzfile(ident_file)), error = function(e) NULL)

## ---- 读入16样本 ----
so_list <- list()
for (i in seq_len(nrow(sample_info))) {
  s <- sample_info[i, ]
  base <- file.path(raw_dir, sub("_barcodes\\.tsv\\.gz$", "", bc_files[i]))
  genes_con <- gzfile(paste0(base, "_genes.tsv.gz"))
  genes <- read.delim(genes_con, header = FALSE, col.names = c("ensembl", "symbol"),
                      stringsAsFactors = FALSE)
  tryCatch(close(genes_con), error = function(e) NULL)
  mtx_con <- gzfile(paste0(base, "_matrix.mtx.gz"))
  mtx <- readMM(mtx_con)
  tryCatch(close(mtx_con), error = function(e) NULL)
  bc_con <- gzfile(paste0(base, "_barcodes.tsv.gz"))
  bc <- readLines(bc_con)
  tryCatch(close(bc_con), error = function(e) NULL)
  stopifnot(nrow(mtx) == nrow(genes), ncol(mtx) == length(bc))

  m <- as(mtx, "CsparseMatrix")
  rownames(m) <- make.unique(genes$symbol)
  cell_names <- paste0(s$sample_id, "_", sub("-1$", "", bc))
  colnames(m) <- cell_names
  n_cells <- ncol(m)
  so <- CreateSeuratObject(
    counts = m,
    meta.data = data.frame(
      sample_id = rep(s$sample_id, n_cells), region = rep(s$region, n_cells),
      donor = rep(s$donor, n_cells), age_years = rep(s$age_years, n_cells),
      age_group = rep(s$age_group, n_cells), tissue = rep(s$tissue, n_cells),
      row.names = cell_names
    )
  )
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so_list[[s$sample_id]] <- so
  log_msg(s$sample_id, ":", n_cells, "cells read")
  rm(mtx, m); gc()
}

## ---- 合并 + join 注释 ----
merged <- merge(so_list[[1]], y = so_list[-1])
rm(so_list); gc()
merged <- JoinLayers(merged)
log_msg("merged:", ncol(merged), "cells x", nrow(merged), "genes")

md <- merged[[]]
md$cell_id <- rownames(md)
md$cell_type_official <- ident$cell_type[match(md$cell_id, ident$cell_id)]
join_rate <- sum(!is.na(md$cell_type_official)) / nrow(md)
log_msg("annotation join rate:", round(join_rate, 4))

major_map <- c(
  Rod = "Rod", Cone = "Cone", Bipolar = "BC", Amacrine = "AC",
  Horizontal = "HC", Muller = "MG",
  Microglia = "Microglia", Astrocytes = "Astrocytes", Pericytes = "Pericytes",
  Ganglion = "RGC", Endothelium = "Endothelium"
)
md$celltype_major <- unname(major_map[md$cell_type_official])

## ---- 仅保留注释细胞 ----
keep_cells <- rownames(md)[!is.na(md$celltype_major)]
merged <- subset(merged, cells = keep_cells)
merged <- AddMetaData(merged, md[keep_cells, c("cell_id", "cell_type_official", "celltype_major")])
log_msg("after annotation filter:", ncol(merged), "cells")

## ---- join QC 汇总 ----
qc <- merged[[]] %>%
  group_by(sample_id, region, donor, age_years) %>%
  summarise(n_cells = n(), n_types = n_distinct(cell_type_official), .groups = "drop")
type_tab <- as.data.frame(table(merged$cell_type_official, merged$sample_id))
write.csv(qc, file.path(data_dir, "2026-10-20_2_GSE155288_join_QC.csv"), row.names = FALSE)
cat("==== join QC by sample ====\n"); print(as.data.frame(qc), row.names = FALSE)
cat("==== cell type totals ====\n"); print(table(merged$cell_type_official))

## ---- 保存合并对象 ----
saveRDS(merged, file.path(data_dir, "2026-10-20_2_GSE155288_merged.rds"))
write.csv(merged[[]], file.path(data_dir, "2026-10-20_2_GSE155288_metadata.csv"))
log_msg("merged object saved.")

## ---- SingleR 参考集（每类 ≤500 细胞, logcounts）----
merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 1e4, verbose = FALSE)
ref_md <- merged[[]]
set.seed(1)
ref_cells <- unlist(lapply(split(rownames(ref_md), ref_md$cell_type_official), function(x) {
  if (length(x) > 500) sample(x, 500) else x
}))
log_msg("SingleR ref cells:", length(ref_cells))
logcounts <- tryCatch(
  LayerData(merged, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(merged, assay = "RNA", slot = "data")
)
ref_obj <- list(
  logcounts = logcounts[, ref_cells],
  labels = ref_md$cell_type_official[ref_cells],
  cell_ids = ref_cells
)
saveRDS(ref_obj, file.path(data_dir, "2026-10-20_2_GSE155288_SingleR_ref.rds"))
log_msg("SingleR reference saved. done.")
