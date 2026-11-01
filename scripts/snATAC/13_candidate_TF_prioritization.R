# 最终整理脚本 13 · 候选 TF 筛选与优先级 (Version B 主表) (2026-10-29)
# 来源: 00use_scripts/30_2026-10-09_Phase4_rev2_candidate_TF_prioritization.ipynb (全部主链 cell)
# 产出: 2026-09-16_ATAC_start/Table/2026-10-09_rev2_candidate_TF_by_DEG_group.tsv (Version B; 12 的直接输入)
#       + Version A sensitivity 表 + TF×DEG group 形态/family 汇总表
# 输入: 10 的 direction_group_overlap + 07 的 chromVAR (annotation) + 06 的 family map
#       + DESeq2 DEG 全表 (RNA 侧 macula/Central vs peripheral) + DAL/AAM 功能基因表 (annotation)
## ---- 1. 输入: Phase 3 direction-group overlap (真实 DEG direction 层, 已含 deg_group_n/coverage) ----
## direction-group 表即任务所需结构 (same/opposite 已映射回 DEG_direction, 且带 deg_group_n 与
## deg_group_coverage_fraction) — 但其 version A/B 都在; 本版以 B 为主表, A 为 sensitivity.
ds_files <- c(dev   = file.path(DIR_TAB, sprintf("%s_dev_TF_target_DM_DAL_AAM_direction_group_overlap_test.tsv", P3_PREFIX)),
              adult = file.path(DIR_TAB, sprintf("%s_adult_TF_target_DM_DAL_AAM_direction_group_overlap_test.tsv", P3_PREFIX)))
stopifnot(all(file.exists(ds_files)))
dg <- rbindlist(lapply(names(ds_files), function(p) fread(ds_files[[p]])))
dg <- dg[!is.na(deg_group_n)]                    # DEG group 不存在的行 (空组) 不参与
cat("direction-group overlap rows (dev+adult, group exists):", nrow(dg), "\n")
cat("version 分布:\n"); print(dg[, .N, by = .(period, version)])
## family_for_collapse 从 Phase 0 map 挂上
motmap <- fread(file.path(DIR_TAB, "2026-09-16-1_motif_TF_family_map.tsv"), select = c("motif_id","family_for_collapse"))
dg <- merge(dg, motmap, by = "motif_id", all.x = TRUE, sort = FALSE)
## ---- 2. chromVAR 背景 (为什么该 TF 被定义为 regional candidate; 只作 annotation) ----
cv_dev <- fread(file.path(DIR_TAB, sprintf("%s_dev_regional_chromVAR.tsv", P1_PREFIX)))
cv_ad  <- fread(file.path(DIR_TAB, sprintf("%s_adult_regional_chromVAR.tsv", P1_PREFIX)))
cv <- rbindlist(list(
  cv_dev[, .(period="dev",   celltype, motif_id, chromVAR_effect=delta,    chromVAR_FDR=t_fdr,
             chromVAR_flag=candidate, chromVAR_n=n_pairs)],
  cv_ad [, .(period="adult", celltype, motif_id, chromVAR_effect=mean_delta, chromVAR_FDR=NA_real_,
             chromVAR_flag=direction_consistent & abs(mean_delta)>=0.10, chromVAR_n=n_donor,
             adult_dir_consistent=direction_consistent)]
), fill = TRUE)
cv[, chromVAR_direction := fifelse(chromVAR_effect > 0, "central", "peripheral")]
cv[period=="dev", adult_dir_consistent := NA]
cat("chromVAR background rows:", nrow(cv), "\n")
## ---- 3. TF RNA 背景 (biological plausibility; 只 annotation, 不筛除, 不要求方向一致) ----
read_deg <- function(files, period) {
  rbindlist(lapply(files, function(f) {
    d <- suppressWarnings(fread(f))
    if (!all(c("celltype","stage","gene","baseMean","log2FoldChange","padj") %in% names(d))) return(NULL)
    d[, .(period = period, celltype, stage, gene, baseMean, log2FoldChange, padj)]
  }), fill = TRUE)
}
dev_deg_files <- list.files(DEV_DEG_DIR, pattern = "^2026-07-06_2_.*_macula_vs_peripheral\\.csv$", full.names = TRUE)
ad_deg_files  <- list.files(AD_DEG_DIR,  pattern = "^2026-07-28_adult_.*_Central_vs_peripheral\\.csv$", full.names = TRUE)
deg <- rbindlist(list(read_deg(dev_deg_files, "dev"), read_deg(ad_deg_files, "adult")), fill = TRUE)
## per (period, celltype, TF-gene): pooled baseMean (max across stages) + mature-stage regional log2FC/padj
expr <- deg[, .(rna_baseMean = max(baseMean, na.rm = TRUE)), by = .(period, celltype, gene)]
mature <- deg[(period == "dev" & stage == "Late") | (period == "adult" & stage == "Old"),
              .(period, celltype, gene, rna_log2FC = log2FoldChange, rna_padj = padj)]
rna <- merge(expr, mature, by = c("period","celltype","gene"), all.x = TRUE)
rna[, rna_expressed := rna_baseMean >= RNA_EXPR_MIN]
rna[, rna_direction := fifelse(is.na(rna_log2FC), NA_character_,
                          fifelse(rna_log2FC > 0, "central", "peripheral"))]
cat("TF RNA background rows:", nrow(rna), "| dev/adult files:", length(dev_deg_files), "/", length(ad_deg_files), "\n")
## ---- 4. 主表: candidate TF × DEG group (Version B 为主; A 附列) ----
mk_master <- function(ver){
  ## 2026-10-09 修复: 形参原名 version 与 dt 列名冲突, version == version 恒 TRUE → A/B 混入;
  ## 改用 ver 形参, i 中 version 列 == ver (data.table i 中 version 指列, ver 指函数形参)
  d <- dg[version == ver & overlap_n > 0]            # 只保留有交集的行 (0 行无候选意义, 不进主表)
  ## chromVAR 背景
  d <- merge(d, cv[, .(period, celltype, motif_id, chromVAR_direction, chromVAR_effect,
                       chromVAR_FDR, chromVAR_flag, chromVAR_n, adult_dir_consistent)],
             by = c("period","celltype","motif_id"), all.x = TRUE)
  ## TF RNA 背景 (motif_name 作为 TF 基因名映射 — 沿用旧 Phase 4 约定: motif_name 即 TF symbol)
  rsub <- rna[, .(period, celltype, gene, rna_baseMean, rna_expressed, rna_log2FC, rna_direction, rna_padj)]
  setnames(rsub, "gene", "motif_name")   # RNA 基因名对齐 motif_name (TF symbol)
  d <- merge(d, rsub, by = c("period","celltype","motif_name"), all.x = TRUE)
  setcolorder(d, c("period","version","celltype","subclass","DEG_direction","motif_id","motif_name",
                   "family_for_collapse","TF_direction","relation_to_TF_direction",
                   "deg_group_n","overlap_n","deg_group_coverage_fraction","overlap_genes",
                   "chromVAR_direction","chromVAR_effect","chromVAR_FDR","chromVAR_flag","chromVAR_n",
                   "adult_dir_consistent","rna_baseMean","rna_expressed","rna_log2FC","rna_direction","rna_padj"))
  d[order(period, celltype, subclass, DEG_direction, -deg_group_coverage_fraction, -overlap_n)]
}
master_B <- mk_master("B")
master_A <- mk_master("A")
## 强制 QC: A/B 严格隔离
stopifnot(
  all(master_B$version == "B"),
  all(master_A$version == "A")
)
cat("master_B versions:", paste(unique(master_B$version), collapse=","), "\n")
cat("master_A versions:", paste(unique(master_A$version), collapse=","), "\n")
cat("master rows: Version B", nrow(master_B), "| Version A", nrow(master_A), "\n")
## 主结果 = Version B; Version A 单独存 sensitivity (不与 B 混入同一主表)
fwrite(master_B, file.path(DIR_TAB, sprintf("%s_rev2_candidate_TF_by_DEG_group.tsv", DATE_PREFIX)), sep = "\t")
fwrite(master_A, file.path(DIR_TAB, sprintf("%s_rev2_candidate_TF_by_DEG_group_versionA.tsv", DATE_PREFIX)), sep = "\t")
cat("wrote:", sprintf("%s_rev2_candidate_TF_by_DEG_group.tsv", DATE_PREFIX), "(Version B main) +",
    sprintf("%s_rev2_candidate_TF_by_DEG_group_versionA.tsv", DATE_PREFIX), "(Version A sensitivity)\n")
cat("\n=== 分布概览 (Version B): overlap_n / coverage 分位数 by period ===\n")
print(master_B[, .(n_rows=.N, n_TF=uniqueN(motif_id),
                   med_overlap=median(overlap_n), q90_overlap=quantile(overlap_n, .9), max_overlap=max(overlap_n),
                   med_cov=round(median(deg_group_coverage_fraction),3),
                   q90_cov=round(quantile(deg_group_coverage_fraction,.9),3)), by=period])
## ---- 5. functional target hits (对 overlap_genes 与 celltype×subclass 对应 functional set 取交集) ----
dal_func <- unique(fread(DAL_FUNC_F, select = c("celltype","gene"))[, .(celltype, gene, subclass = "DAL")])
aam_func <- unique(fread(AAM_FUNC_F, select = c("celltype","gene"))[, .(celltype, gene, subclass = "AAM")])
func_genes <- unique(rbind(dal_func, aam_func))
func_by_ct_sub <- func_genes[, .(fgenes = list(unique(gene))), by = .(celltype, subclass)]
cat("curated functional-module genes: DAL", nrow(dal_func), " AAM", nrow(aam_func), " (DM 无 curated set → NA)\n")

add_func <- function(m){
  mm <- merge(m, func_by_ct_sub, by = c("celltype","subclass"), all.x = TRUE)
  mm[, `:=`(functional_hit_n = NA_integer_, functional_hit_genes = NA_character_)]
  has_set <- !is.na(mm$fgenes) & sapply(mm$fgenes, function(x) !is.null(x) && length(x) > 0)
  mm[has_set, functional_hit_n := {
    og <- strsplit(overlap_genes, ";")[[1]]
    hit <- intersect(og, fgenes[[1]])
    length(hit)
  }, by = seq_len(sum(has_set))]
  ## 逐行 safe 版本 (上面的 by=seq_len 写法在子集上易错, 改为显式循环)
  mm[, functional_hit_genes := NA_character_]
  for (i in which(has_set)) {
    og <- strsplit(mm$overlap_genes[i], ";")[[1]]
    hit <- intersect(og, mm$fgenes[[i]])
    mm$functional_hit_n[i] <- length(hit)
    mm$functional_hit_genes[i] <- paste(sort(hit), collapse = ";")
  }
  mm[, fgenes := NULL]
  mm
}
fh_B <- add_func(master_B)
fwrite(fh_B[, .(period, version, celltype, subclass, DEG_direction, motif_id, motif_name,
                overlap_genes, functional_hit_n, functional_hit_genes)],
       file.path(DIR_TAB, sprintf("%s_rev2_candidate_TF_functional_hits.tsv", DATE_PREFIX)), sep = "\t")
cat("wrote:", sprintf("%s_rev2_candidate_TF_functional_hits.tsv", DATE_PREFIX),
    "| DAL/AAM 命中行:", sum(!is.na(fh_B$functional_hit_n) & fh_B$functional_hit_n>0), "| DM 全 NA:", all(is.na(fh_B[subclass=="DM"]$functional_hit_n)), "\n")
## ---- 6. shared pattern 第一层: 同 celltype×DEG_direction 跨 subclass ----
## 对 TF × celltype × DEG_direction: 覆盖了几个 subclass (overlap_n>0 定义 hit)
sh1 <- master_B[overlap_n > 0, .(subclass_hit_list = paste(sort(unique(subclass)), collapse="+"),
                                  n_subclass_hit = uniqueN(subclass)),
                by = .(period, celltype, DEG_direction, motif_id, motif_name)]
sh1[, subclass_pattern := fifelse(n_subclass_hit == 3L, "DM+DAL+AAM",
                           fifelse(n_subclass_hit == 2L, subclass_hit_list,
                             paste0(subclass_hit_list, " only")))]
fwrite(sh1[order(period, celltype, DEG_direction, -n_subclass_hit)],
       file.path(DIR_TAB, sprintf("%s_rev2_candidate_TF_shared_subclass.tsv", DATE_PREFIX)), sep = "\t")
cat("wrote:", sprintf("%s_rev2_candidate_TF_shared_subclass.tsv", DATE_PREFIX), "\n")
cat("\n=== subclass pattern 分布 (Version B, by period) ===\n")
print(dcast(sh1[, .N, by = .(period, subclass_pattern)], period~subclass_pattern, value.var = "N", fill = 0L))
## ---- 7. shared pattern 第二层: 跨 celltype (同 DEG_direction); 第三层: dev/adult 保留 period ----
sh2 <- master_B[overlap_n > 0, .(
    n_celltype_hit = uniqueN(celltype),
    celltype_hit_list = paste(sort(unique(celltype)), collapse="+"),
    celltype_subclass_hit_list = paste(sort(unique(paste0(celltype, ":", subclass))), collapse="+")),
  by = .(period, DEG_direction, motif_id, motif_name)]
fwrite(sh2[order(period, DEG_direction, -n_celltype_hit)],
       file.path(DIR_TAB, sprintf("%s_rev2_candidate_TF_shared_celltype.tsv", DATE_PREFIX)), sep = "\t")
cat("wrote:", sprintf("%s_rev2_candidate_TF_shared_celltype.tsv", DATE_PREFIX), "\n")
cat("\n=== 跨 celltype 分布 (Version B): n_celltype_hit 分布 by period ===\n")
print(sh2[, .N, by = .(period, n_celltype_hit)][order(period, n_celltype_hit)])
## dev/adult presence comparison (descriptive only, 不定义 final 分类)
pres <- dcast(master_B[overlap_n>0, .N, by=.(motif_id, period)], motif_id~period, value.var="N", fill=0L)
cat("\n[dev/adult presence] both:", sum(pres$dev>0 & pres$adult>0),
    "| dev-only:", sum(pres$dev>0 & pres$adult==0),
    "| adult-only:", sum(pres$dev==0 & pres$adult>0), "(descriptive, 不作最终结论)\n")
## ---- 8. family-level summary (annotation; 不指定 family 唯一代表) ----
fam <- master_B[overlap_n > 0, .(n_motif = uniqueN(motif_id), motif_list = paste(sort(unique(motif_name)), collapse=","),
                                 max_overlap = max(overlap_n), max_cov = round(max(deg_group_coverage_fraction),3)),
                by = .(period, celltype, subclass, DEG_direction, family_for_collapse)][order(-max_overlap)]
fwrite(fam, file.path(DIR_TAB, sprintf("%s_rev2_candidate_TF_family_summary.tsv", DATE_PREFIX)), sep = "\t")
cat("wrote:", sprintf("%s_rev2_candidate_TF_family_summary.tsv", DATE_PREFIX), "|\n")
cat("\n=== top families (Version B, by max_overlap) ===\n")
print(head(fam, 12))
