# 最终整理脚本 09 · 区域差异 chromVAR 统计 (每 celltype 内 central vs peripheral) (2026-10-29)
# 来源: 00use_scripts/27_2026-09-16-2_Phase1_regional_chromVAR.ipynb (仅统计主链;
#       volcano/dumbbell/PoC 速览/候选数量柱状图等诊断绘图已去除)
# 产出: 2026-09-16_ATAC_start/Table/2026-09-16-2_dev_regional_chromVAR.tsv
#       2026-09-16_ATAC_start/Table/2026-09-16-2_adult_regional_chromVAR.tsv
# 口径: region 侧细胞数 ≥ 50 (加权 pseudobulk); 胚胎配对 donor ≥ 6 (paired t + wilcox, BH FDR);
#       成年 = 21y/67y 共享 ΔchromVAR (无 P/FDR, 两 donor 方向一致); candidate = FDR<0.05 且 |Δ|≥0.10
# 输入: candidate_CRE 的 per-sample chromVAR 汇总 (2026-08-12_4_*) + 06 的 family map
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(jsonlite)
  library(RColorBrewer); library(IRdisplay); library(ggrepel)
})
options(stringsAsFactors = FALSE, width = 140)

# ---- 1. Parameters ----
DATE_PREFIX <- "2026-09-16-2"
ROOT     <- "/project/work"
CRE_DATA <- file.path(ROOT, "2026-08-08_candidate_CRE/Data")
OUT      <- file.path(ROOT, "2026-09-16_ATAC_start")
DIR_TAB  <- file.path(OUT, "Table"); DIR_PLOT <- file.path(OUT, "Plot"); DIR_LOG <- file.path(OUT, "Log")
P0_TAB   <- file.path(DIR_TAB, "2026-09-16-1_motif_TF_family_map.tsv")

LOWPOWER_NCELL <- 50L    # region 侧细胞数下限（与 Phase 0 一致）
MIN_PAIRS_DEV  <- 6L     # 发育配对检验最少配对 donor 数
SIG_FDR        <- 0.05
SIG_DELTA      <- 0.10   # 最小 |mean Δ chromVAR| 效应量（叫 candidate 用）
POC            <- c("MG","Cone")
ADULT_21Y <- "BCM_22_0183"; ADULT_67Y <- "BCM_22_0389"
celltype_order <- c("Rod","Cone","AC","BC","HC","MG","RGC","RPC")

f_dev_ss <- file.path(CRE_DATA, "2026-08-12_4_develop_chromVAR_sample_summary.csv")
f_ad_ss  <- file.path(CRE_DATA, "2026-08-12_4_adult_chromVAR_sample_summary.csv")

# colors from color.ipynb
source_color_ipynb <- function(path = file.path(ROOT, "color.ipynb")) {
  j <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  codes <- vapply(j$cells, function(cc)
    if (identical(cc$cell_type,"code")) paste(unlist(cc$source), collapse="") else "", character(1))
  eval(parse(text = paste(codes, collapse="\n")), envir = globalenv())
}
invisible(tryCatch(source_color_ipynb(), error=function(e) message("color.ipynb: ", conditionMessage(e))))
if (!exists("cols.direction")) cols.direction <- c("Macular"="#D55E00","Peripheral"="#0072B2")
DIRCOL <- c("central-high"=unname(cols.direction["Macular"]),
            "peripheral-high"=unname(cols.direction["Peripheral"]), "ns"="grey80")

save_plot <- function(p, name, w=8, h=6) {
  suppressMessages(ggsave(file.path(DIR_PLOT, sprintf("%s_%s.pdf", DATE_PREFIX, name)), p, width=w, height=h))
  suppressMessages(ggsave(file.path(DIR_PLOT, sprintf("%s_%s.png", DATE_PREFIX, name)), p, width=w, height=h, dpi=300))
  cat("saved:", sprintf("%s_%s", DATE_PREFIX, name), "(pdf+png)\n"); invisible(p)
}
motmap <- fread(P0_TAB)[, .(motif_id, motif_name, TF_annotation, family, is_multi_TF)]
cat("Params ready.\n")
load_ss <- function(f) {
  dt <- fread(f, select = c("dataset_period","celltype","donor_id","tissue","tissue_region",
                            "formal_regional_pair","sample","motif_id","motif_name","mean_chromVAR","n_cells"))
  dt[!is.na(celltype) & celltype != "" & celltype %in% celltype_order]
}
dev <- load_ss(f_dev_ss); ad <- load_ss(f_ad_ss)

# 重复样本诊断
dup_report <- function(dt, nm) {
  d <- dt[, .N, by = .(celltype, donor_id, tissue_region, motif_id)][N > 1]
  cat(sprintf("[%s] 每 (celltype×donor×region×motif) >1 样本的组合数: %d（将被加权合并）\n", nm, nrow(d)))
}
dup_report(dev, "dev"); dup_report(ad, "adult")

# region 级细胞数（跨样本求和）+ 加权 pseudobulk
region_ncell <- function(dt) {
  smeta <- unique(dt[, .(celltype, donor_id, tissue_region, sample, n_cells)])
  smeta[, .(n_cells_region = sum(n_cells)), by = .(celltype, donor_id, tissue_region)]
}
agg_pb <- function(dt) dt[, .(mean_chromVAR = weighted.mean(mean_chromVAR, n_cells)),
                          by = .(celltype, donor_id, tissue_region, motif_id, motif_name)]

elig_pairs <- function(dt) {
  rn <- region_ncell(dt)[n_cells_region >= LOWPOWER_NCELL]
  rn[, .(nreg = uniqueN(tissue_region)), by = .(celltype, donor_id)][nreg == 2L, .(celltype, donor_id)]
}
dev_pair <- elig_pairs(dev); ad_pair <- elig_pairs(ad)
pb_dev <- merge(agg_pb(dev), dev_pair, by = c("celltype","donor_id"))
pb_ad  <- merge(agg_pb(ad),  ad_pair,  by = c("celltype","donor_id"))[donor_id %in% c(ADULT_21Y, ADULT_67Y)]
# 唯一性校验
stopifnot(nrow(pb_dev[, .N, by=.(celltype,donor_id,tissue_region,motif_id)][N>1]) == 0)
stopifnot(nrow(pb_ad[,  .N, by=.(celltype,donor_id,tissue_region,motif_id)][N>1]) == 0)
cat("发育合格配对(celltype×donor):", nrow(dev_pair), " | 成年:", nrow(ad_pair),
    "(限21y/67y后:", nrow(unique(pb_ad[,.(celltype,donor_id)])), ") | 唯一性校验通过 ✓\n")
wdev <- dcast(pb_dev, celltype + motif_id + motif_name + donor_id ~ tissue_region, value.var = "mean_chromVAR")
wdev <- wdev[!is.na(central) & !is.na(peripheral)]
dev_res <- wdev[, {
  n <- .N
  if (n >= MIN_PAIRS_DEV && sd(central - peripheral) > 0) {
    tp <- tryCatch(t.test(central, peripheral, paired = TRUE)$p.value, error = function(e) NA_real_)
    wp <- tryCatch(suppressWarnings(wilcox.test(central, peripheral, paired = TRUE)$p.value), error = function(e) NA_real_)
  } else { tp <- NA_real_; wp <- NA_real_ }
  .(n_pairs = n, central_mean = mean(central), peripheral_mean = mean(peripheral),
    delta = mean(central - peripheral), t_p = tp, wilcox_p = wp)
}, by = .(celltype, motif_id, motif_name)]
dev_res[, t_fdr := p.adjust(t_p, "BH"), by = celltype]
dev_res[, wilcox_fdr := p.adjust(wilcox_p, "BH"), by = celltype]
dev_res[, direction := fifelse(delta > 0, "central-high", "peripheral-high")]
dev_res[, candidate := !is.na(t_fdr) & t_fdr < SIG_FDR & abs(delta) >= SIG_DELTA]
dev_res <- merge(dev_res, motmap[, .(motif_id, TF_annotation, family, is_multi_TF)], by = "motif_id", all.x = TRUE)
setcolorder(dev_res, c("celltype","motif_id","motif_name","TF_annotation","family",
                       "n_pairs","central_mean","peripheral_mean","delta",
                       "t_p","t_fdr","wilcox_p","wilcox_fdr","direction","candidate","is_multi_TF"))
setorder(dev_res, celltype, t_fdr, -delta)
fwrite(dev_res, file.path(DIR_TAB, sprintf("%s_dev_regional_chromVAR.tsv", DATE_PREFIX)), sep = "\t")

cat("发育 candidate motif（FDR<", SIG_FDR, "& |Δ|>=", SIG_DELTA, "）计数:\n")
print(dcast(dev_res[candidate == TRUE], celltype ~ direction, fun.aggregate = length, value.var = "motif_id"))
wad <- dcast(pb_ad, celltype + motif_id + motif_name + donor_id ~ tissue_region, value.var = "mean_chromVAR")
wad <- wad[!is.na(central) & !is.na(peripheral)]
wad[, delta := central - peripheral]
ad_res <- dcast(wad, celltype + motif_id + motif_name ~ donor_id, value.var = "delta")
setnames(ad_res, ADULT_21Y, "delta_21y", skip_absent = TRUE)
setnames(ad_res, ADULT_67Y, "delta_67y", skip_absent = TRUE)
for (cc in c("delta_21y","delta_67y")) if (!cc %in% names(ad_res)) ad_res[[cc]] <- NA_real_
ad_res[, n_donor := as.integer(!is.na(delta_21y)) + as.integer(!is.na(delta_67y))]
ad_res[, mean_delta := rowMeans(cbind(delta_21y, delta_67y), na.rm = TRUE)]
ad_res[, direction_consistent := !is.na(delta_21y) & !is.na(delta_67y) & sign(delta_21y) == sign(delta_67y)]
ad_res[, direction := fifelse(mean_delta > 0, "central-high", "peripheral-high")]
ad_res <- merge(ad_res, motmap[, .(motif_id, TF_annotation, family, is_multi_TF)], by = "motif_id", all.x = TRUE)
setcolorder(ad_res, c("celltype","motif_id","motif_name","TF_annotation","family",
                      "n_donor","delta_21y","delta_67y","mean_delta","direction_consistent","direction","is_multi_TF"))
setorder(ad_res, celltype, -mean_delta)
fwrite(ad_res, file.path(DIR_TAB, sprintf("%s_adult_regional_chromVAR.tsv", DATE_PREFIX)), sep = "\t")
cat("成年表已存（无 P/FDR 列）。方向一致(top |Δ|)占比（n_donor==2）:\n")
print(ad_res[n_donor == 2, .(n_motif = .N, n_consistent = sum(direction_consistent),
                             pct = round(100*mean(direction_consistent),1)), by = celltype][order(celltype)])
