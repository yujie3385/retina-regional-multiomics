# 最终整理脚本 08 · motif↔TF↔family 映射表构建 (2026-10-29)
# 来源: 00use_scripts/26_2026-09-16-1_Phase0_environment_QC.ipynb (仅映射表构建主链;
#       cell-count QC 主表/配对 donor 功效/坐标核查/QC tile 图等诊断段已去除)
# 产出: 2026-09-16_ATAC_start/Table/2026-09-16-1_motif_TF_family_map.tsv
#       (下游: 07 chromVAR 统计挂注释 / 08 热图 / 11 候选 TF family 汇总)
# 输入: candidate_CRE 项目的 motif inventory 与 DM motif family 注释表
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(RColorBrewer); library(jsonlite); library(IRdisplay)
})
options(stringsAsFactors = FALSE, width = 140)

# ---- 1. Parameters ----
DATE_PREFIX <- "2026-09-16-1"
ROOT     <- "/project/work"
CRE_DATA <- file.path(ROOT, "2026-08-08_candidate_CRE/Data")
OUT      <- file.path(ROOT, "2026-09-16_ATAC_start")
DIR_DATA <- file.path(OUT, "Data"); DIR_TAB <- file.path(OUT, "Table")
DIR_PLOT <- file.path(OUT, "Plot"); DIR_LOG <- file.path(OUT, "Log")
for (d in c(DIR_DATA, DIR_TAB, DIR_PLOT, DIR_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOWPOWER_NCELL <- 50L   # celltype×donor×region 细胞数 < 此值 → 标 low-power，不下阳性结论

f_dev_ss <- file.path(CRE_DATA, "2026-08-12_4_develop_chromVAR_sample_summary.csv")
f_ad_ss  <- file.path(CRE_DATA, "2026-08-12_4_adult_chromVAR_sample_summary.csv")
f_motinv <- file.path(CRE_DATA, "2026-08-12_4_chromVAR_motif_inventory.csv")
f_fam    <- file.path(CRE_DATA, "2026-08-15_3_DM_motif_family_annotation.csv")
f_devDAR <- file.path(ROOT, "2026-07-06_development_DAR/Data/2026-07-06_develop_3_regional_dars_core_filtered_merged_ChIPseeker_annotated.csv")
f_adDAR  <- file.path(ROOT, "2026-08-03_adult_DAR/Data/2026-08-03_adult_shared_fovea_peripheral_DAR.csv")
f_peakM  <- file.path(CRE_DATA, "2026-08-12_DEG_DAR_ArchR_link_peak_level.csv")

check_required <- function(dt, cols, nm) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(sprintf("[%s] 缺列: %s", nm, paste(miss, collapse=", ")))
  invisible(TRUE)
}
celltype_order <- c("Rod","Cone","AC","BC","HC","MG","RGC","RPC")
cat("Params set. R", R.version.string, "\n")
motinv <- fread(f_motinv)
fam    <- fread(f_fam)
check_required(motinv, c("motif_id","motif_name","TF_annotation"), "motif_inventory")
check_required(fam, c("motif_id","JASPAR_family","family_for_collapse"), "family_annotation")

motmap <- merge(motinv, unique(fam[, .(motif_id, JASPAR_family, JASPAR_class, family_for_collapse)]),
                by = "motif_id", all.x = TRUE)
motmap[, is_multi_TF := grepl("::", TF_annotation)]
motmap[, family := fifelse(!is.na(family_for_collapse) & family_for_collapse != "", family_for_collapse, NA_character_)]

n_mot <- nrow(motmap)
n_fam_cov <- motmap[!is.na(family), .N]
n_multiTF <- motmap[is_multi_TF == TRUE, .N]
fam_sizes <- motmap[!is.na(family), .N, by = family][order(-N)]
cat(sprintf("motifs: %d | family 覆盖: %d (%.1f%%) | 多TF(::) motif: %d | 不同 family 数: %d\n",
            n_mot, n_fam_cov, 100*n_fam_cov/n_mot, n_multiTF, nrow(fam_sizes)))
cat("Top 8 冗余 family (motif 数):\n"); print(head(fam_sizes, 8))

fwrite(motmap, file.path(DIR_TAB, sprintf("%s_motif_TF_family_map.tsv", DATE_PREFIX)), sep = "\t")
cat("→ Table/", sprintf("%s_motif_TF_family_map.tsv", DATE_PREFIX), "\n", sep="")
