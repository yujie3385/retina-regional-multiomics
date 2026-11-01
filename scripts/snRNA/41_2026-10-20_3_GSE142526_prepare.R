#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-20_3  GSE142526 胚胎期准备（D125C/D125P, FD125, scRNA）
# 验证 F6 注释与矩阵 barcode 的 join，生成细胞级 metadata
# 输出:
#   Data/2026-10-20_3_GSE142526_metadata.csv      细胞级metadata
#   Data/2026-10-20_3_GSE142526_join_QC.csv       join QC汇总
# =====================================================================
suppressPackageStartupMessages({ library(dplyr) })
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
dev_dir  <- file.path(proj_dir, "00data/2026-10-20_development_new/GSE142526")

## ---- F6 注释（D125C/D125P, 100% 覆盖）----
f6 <- read.csv(file.path(dev_dir, "metadata/GSE142526_metadata/F6/cca_fetalvsorg_125CP_205_metadata.csv"),
               stringsAsFactors = FALSE)
f6 <- f6[f6$orig.ident %in% c("D125C", "D125P"), ]
f6$bc16 <- sub("_\\d+$", "", f6[[1]])
stopifnot(!any(duplicated(paste(f6$orig.ident, f6$bc16))))
cat("F6 rows (D125C+D125P):", nrow(f6), "\n")

## ---- 矩阵 barcodes ----
meta_list <- list()
join_qc <- list()
for (smp in c("D125C", "D125P")) {
  bc <- readLines(file.path(dev_dir, smp, "GRCh38", "barcodes.tsv"))
  bc16 <- sub("-1$", "", bc)
  f6s <- f6[f6$orig.ident == smp, ]
  hit <- match(bc16, f6s$bc16)
  cat(smp, ": matrix barcodes =", length(bc16),
      "; matched in F6 =", sum(!is.na(hit)),
      "; F6 not in matrix =", sum(!f6s$bc16 %in% bc16), "\n")
  join_qc[[smp]] <- data.frame(
    sample = smp, n_matrix_barcodes = length(bc16),
    n_matched = sum(!is.na(hit)),
    n_f6_cells = nrow(f6s),
    n_f6_unmatched = sum(!f6s$bc16 %in% bc16)
  )
  cells <- paste0(smp, "_", bc16)
  f6_idx <- f6s[hit, ]
  meta_list[[smp]] <- data.frame(
    cell_id = cells,
    sample_id = smp,
    region = ifelse(smp == "D125C", "Central", "Peripheral"),
    tissue = ifelse(smp == "D125C", "central", "peripheral"),
    donor_age = "FD125",
    type = ifelse(is.na(hit), NA, f6_idx$type),
    stringsAsFactors = FALSE
  )
}
md <- bind_rows(meta_list)

## ---- type -> major 映射（PR=rod/cone未拆分, iMG=未成熟Müller）----
major_map <- c(PR = "PR", BC = "BC", AC = "AC", AC1 = "AC", AC2 = "AC",
               HC = "HC", iMG = "MG", Prog = "Prog", `Prog/Glia` = "Prog/Glia",
               `T1/T3` = "T1/T3", T2 = "T2", RGC = "RGC", miG = "miG", Midbrain = "Midbrain")
md$celltype_major <- unname(major_map[md$type])

cat("\n==== type x sample ====\n")
print(table(md$type, md$sample_id, useNA = "ifany"))
cat("\n==== celltype_major x sample ====\n")
print(table(md$celltype_major, md$sample_id, useNA = "ifany"))

write.csv(md, file.path(data_dir, "2026-10-20_3_GSE142526_metadata.csv"), row.names = FALSE)
write.csv(bind_rows(join_qc), file.path(data_dir, "2026-10-20_3_GSE142526_join_QC.csv"), row.names = FALSE)
cat("\nmetadata saved:", nrow(md), "cells\n")
