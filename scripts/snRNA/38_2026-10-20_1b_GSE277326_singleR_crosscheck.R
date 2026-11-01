#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-20_1b  GSE277326 cluster × SingleR 交叉验证
# 参考: GSE155288 官方注释（同组织, 11 类, 每类≤500细胞）
# 输出: Data/2026-10-20_1_GSE277326_singleR_cluster.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(dplyr)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")

obj <- readRDS(file.path(data_dir, "2026-10-20_1_tmp_clustered.rds"))
Idents(obj) <- "RNA_snn_res.0.8"
cluster_ids <- as.character(Idents(obj))

query <- tryCatch(
  LayerData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
)
ref_obj <- readRDS(file.path(data_dir, "2026-10-20_2_GSE155288_SingleR_ref.rds"))

cat("[", format(Sys.time(), "%H:%M:%S"), "] running SingleR (per cluster) ...\n")
res <- SingleR(
  test = query, ref = ref_obj$logcounts, labels = ref_obj$labels,
  clusters = cluster_ids
)
out <- as.data.frame(res$by.cluster[, c("labels", "delta.next")])
out$cluster <- rownames(out)
out$n_cells <- as.integer(table(cluster_ids)[out$cluster])
out <- out[order(as.integer(out$cluster)), c("cluster", "n_cells", "labels", "delta.next")]
write.csv(out, file.path(data_dir, "2026-10-20_1_GSE277326_singleR_cluster.csv"), row.names = FALSE)
cat("==== SingleR per-cluster labels ====\n")
print(out, row.names = FALSE)
cat("done.\n")
