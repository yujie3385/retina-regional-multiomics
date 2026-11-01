# 最终整理脚本 14 · TF-DEG supporting peaks 与 DAR 的 overlap 及 support 汇总 (2026-10-29)
# 来源: 00use_scripts/31_2026-10-09_Phase5_1_TF_DEG_supporting_peak_DAR_overlap.ipynb (主链 + QC cell;
#       分布图 (观察用, 无 cutoff) 已去除)
# 产出: Table/2026-10-09_Phase5_TF_DEG_supporting_peak_DAR_overlap_peaklevel.tsv (peak 级精确交集)
#       Table/2026-10-09_Phase5_TF_DEG_DAR_support_summary.tsv (TF-target DEG 级; 13 的直接输入)
#       + TF×DEG group 描述性汇总 + unmatched 关系 QC 表
# 输入: 11 的 Version B 主表 + 09 的 targets/links/motif presence
#       + candidate_CRE 的 DEG-DAR pairing 主表 (含各期 DAR 状态) + ArchR strict links
## ---- 1. 读入 Phase 4 主表 (Version B) 并展开到 TF × target DEG ----
## 分析单位从 "TF × DEG group" 细化为 "TF × target DEG": overlap_genes 按 ";" 拆开
p4 <- fread(P4_MASTER_F)
stopifnot(all(p4$version == "B"))     # QC: 主表必须纯 Version B
p4g <- p4[, .(period, celltype, subclass, DEG_direction, motif_id, motif_name, TF_direction,
              deg_group_n, overlap_n, target_gene = strsplit(overlap_genes, ";")[[1]]),
          by = seq_len(nrow(p4))]
p4g[, seq_len := NULL]
cat("QC1: Phase 4 展开的 TF × target DEG 关系数:", nrow(p4g), "\n")
print(p4g[, .(n_relations = .N, n_TF = uniqueN(motif_id), n_genes = uniqueN(target_gene)), by = period])
print(p4g[, .(n_relations = .N), by = .(period, celltype)][order(period, celltype)])
## ---- 2. 回查 Phase 2: 为每个 (period, celltype, motif, target_gene) 提取 supporting peaks ----
## Phase 2 targets 表是 gene 级 (无 peak 列); peak 级关系 = P2G links × motif-in-peak × region-accessibility
## 过滤 (与 Phase 2 assemble 完全同口径): 含该 motif ∧ P2G 指向该 gene ∧ acc_pct>=0.05 ∧ acc_ndet>=阈值
extract_supporting_peaks <- function(period){
  tf   <- fread(P2_TARGETS_F[[period]], select = c("celltype","motif_id","motif_name","direction","gene","versionB"))
  tf   <- unique(tf[versionB == TRUE, .(celltype, motif_id, motif_name, direction, gene)])   # Version B 高置信关系
  links <- fread(P2_LINKS_F[[period]], select = c("celltype","peak","gene","score","FDR","seqnames","start","end"))
  mip  <- readRDS(MIP_F[[period]])
  pct  <- readRDS(PCT_F[[period]]); ndet <- readRDS(NDET_F[[period]])
  mipcols <- colnames(mip); peaks <- rownames(mip)
  stopifnot(identical(rownames(pct), peaks))
  dmin <- DONOR_DET[[period]]

  outl <- list()
  cts <- intersect(unique(tf$celltype), unique(links$celltype))
  for(ct in cts){
    tfc <- tf[celltype == ct]
    if(!nrow(tfc)) next
    linkc <- links[celltype == ct, .(peak, gene, score, FDR)]
    grp <- function(reg) paste(ct, reg, sep = "__")
    if(!(grp("central") %in% colnames(pct) && grp("peripheral") %in% colnames(pct))) next
    accdt <- data.table(peak = peaks,
                        central_pct = pct[, grp("central")],  central_ndet = ndet[, grp("central")],
                        periph_pct  = pct[, grp("peripheral")], periph_ndet  = ndet[, grp("peripheral")])
    ## motif-in-peak 限 P2G peaks × 该 ct 的候选 motif
    lk_peaks <- unique(linkc$peak); pidx <- match(lk_peaks, peaks)
    midx <- match(unique(tfc$motif_id), mipcols)
    sub <- as.matrix(mip[pidx, midx, drop = FALSE])
    w <- which(sub, arr.ind = TRUE)
    if(!nrow(w)) next
    dtmp <- data.table(peak = lk_peaks[w[, 1]], motif_id = unique(tfc$motif_id)[w[, 2]])
    dtmp <- merge(dtmp, tfc, by = "motif_id", allow.cartesian = TRUE)
    dtmp <- merge(dtmp, accdt, by = "peak", allow.cartesian = TRUE)
    ## region 特异 accessibility 过滤 (与 Phase 2 完全一致)
    dtmp[, acc_pct  := fifelse(direction == "central-high", central_pct,  periph_pct)]
    dtmp[, acc_ndet := fifelse(direction == "central-high", central_ndet, periph_ndet)]
    dtmp <- dtmp[acc_pct >= ACCESS_MIN & acc_ndet >= dmin]
    if(!nrow(dtmp)) next
    ## 连 P2G gene (peak 多对多) — 只保留指向候选 (motif, gene) 组合的行
    tl <- merge(dtmp, linkc, by = "peak", allow.cartesian = TRUE)
    tl <- tl[gene.x == gene.y]                      # motif 挂的 target gene == P2G 指向 gene
    tl[, `:=`(target_gene = gene.x, gene.x = NULL, gene.y = NULL, period = period)]
    outl[[ct]] <- tl[, .(period, celltype, motif_id, motif_name, direction, target_gene,
                         supporting_peak = peak, P2G_cor = score, P2G_FDR = FDR)]
  }
  rbindlist(outl)
}
sup <- rbindlist(lapply(c("dev", "adult"), extract_supporting_peaks))
cat("QC3 (预): 提取的 supporting peak 行数:", nrow(sup), "\n")
## ---- 3. QC2: Phase 4 的 TF-DEG 关系能多少回查到 Phase 2 supporting peaks ----
## 原则上 Phase 4 来源于 Phase 2 (同 Version B), 应接近全部匹配
key_p4 <- unique(p4g[, .(period, celltype, motif_id, target_gene)])
key_sup <- unique(sup[, .(period, celltype, motif_id, target_gene)])
matched <- key_p4[key_sup, on = names(key_p4), nomatch = 0]
unmatched <- key_p4[!key_sup, on = names(key_p4)]
cat("QC2: Phase 4 关系数:", nrow(key_p4),
    "| 成功回查 supporting peaks:", nrow(matched),
    "| 未匹配:", nrow(unmatched), "\n")
if(nrow(unmatched)){
  cat("\n[检查] 未匹配的 TF-DEG 关系 (前 20):\n")
  print(head(unmatched, 20))
  ## 可能原因: (a) TF 同时挂多 gene 时 pmap 的 gene 与 motif 组合在 Phase 2 基因级表为 TRUE
  ##           但 peak 级重建时 P2G link 不含该 (motif, gene) 组合; (b) motif 在该 ct 的 link peaks 无 motif 命中.
  fwrite(unmatched, file.path(DIR_TAB, sprintf("%s_Phase5_unmatched_TF_DEG_relations.tsv", DATE_PREFIX)), sep = "\t")
  cat("已写出 unmatched 清单供核查:", sprintf("%s_Phase5_unmatched_TF_DEG_relations.tsv", DATE_PREFIX), "\n")
}
## ---- 4. 读入正式 DEG-DAR pairing 主表, 整理两期完整 DAR 信息 (2026-10-09 rev2) ----
## dev: 完整保留 Early/Mid/Late 逐 stage DAR (status/direction/effect/padj), 不压缩为单一方向;
##      is_DAR = 任一 stage 为正式 DAR; 方向冲突 (DAR_direction_conflict) 保留不删.
## adult: 正式 DAR = Adult_shared (core DAR/direction/log2FC/padj); 21y/67y 仅 donor-level annotation.
## P2G distance: 从 ArchR 原始 strict links 按 (celltype,peak,gene) 关联 (正式字段, 不自行计算).
pair <- fread(PAIRING_F, select = c("celltype","gene","peak",
  "Early_DAR_status","Mid_DAR_status","Late_DAR_status",
  "Early_log2FoldChange","Mid_log2FoldChange","Late_log2FoldChange",
  "Early_padj","Mid_padj","Late_padj",
  "Adult_shared_core_DAR","Adult_shared_direction","Adult_shared_log2FC","Adult_shared_padj",
  "Adult_21y_effect","Adult_21y_effect_direction","Adult_67y_effect","Adult_67y_effect_direction",
  "Adult_donor_match_n","adult_grade"))
pair <- unique(pair, by = c("celltype","gene","peak"))
cat("pairing 主表:", nrow(pair), "行 (celltype × gene × peak)\n")

## --- dev: 逐 stage 展开为正式 DAR 行 (status ∈ {M,P}) ---
st_cols <- c("Early","Mid","Late")
dev_rows <- list()
for(st in st_cols){
  s <- paste0(st, "_DAR_status")
  d <- pair[get(s) %in% c("M","P"), .(celltype, gene, peak,
      direction = fifelse(get(s) == "M", "central", "peripheral"),
      effect = get(paste0(st, "_log2FoldChange")),
      padj   = get(paste0(st, "_padj")), stage = st)]
  if(nrow(d)) dev_rows[[st]] <- d
}
dev_dar <- rbindlist(dev_rows)
## dev 汇总: 每 (ct,gene,peak) 完整 stage 级信息
dir_map <- c(central = "central", peripheral = "peripheral")
## 2026-10-09 rev3: stage 按生物学时间序排列 (Early → Mid → Late), 不用字母序
stage_order <- c("Early", "Mid", "Late")
dev_dar[, stage := factor(stage, levels = stage_order)]
setkey(dev_dar, celltype, gene, peak, stage)   # 按 key 顺序即时间序
dev_sum <- dev_dar[, .(
  DAR_stage_list        = paste(unique(as.character(stage)), collapse = ";"),
  DAR_direction_list    = paste(as.character(direction), collapse = ";"),
  DAR_stage_direction_list = paste(paste0(as.character(stage), ":", as.character(direction)), collapse = ";"),
  DAR_effect_list       = paste(round(effect, 4), collapse = ";"),
  DAR_padj_list         = paste(signif(padj, 3), collapse = ";"),
  DAR_direction_conflict= length(unique(direction)) > 1L,
  ## representative: |effect| 最大的 stage — direction/effect/padj 同源该 stage
  representative_DAR_stage = as.character(stage[which.max(abs(effect))]),
  DAR_direction = direction[which.max(abs(effect))],
  DAR_effect    = effect[which.max(abs(effect))],
  DAR_padj      = padj[which.max(abs(effect))]),
  by = .(celltype, gene, peak)]
dev_sum[, period := "dev"]

## --- adult: 正式 DAR = Adult_shared; 21y/67y 为 donor-level annotation (不定义 21y/67y DAR) ---
ad_sum <- pair[Adult_shared_core_DAR == TRUE, .(period = "adult", celltype, gene, peak,
  DAR_stage_list = "adult", DAR_direction = fifelse(Adult_shared_direction == "M", "central", "peripheral"),
  DAR_direction_list = fifelse(Adult_shared_direction == "M", "central", "peripheral"),
  DAR_stage_direction_list = paste0("adult:", fifelse(Adult_shared_direction == "M", "central", "peripheral")),
  DAR_effect_list = as.character(round(Adult_shared_log2FC, 4)),
  DAR_padj_list = as.character(signif(Adult_shared_padj, 3)),
  DAR_direction_conflict = FALSE,
  representative_DAR_stage = "adult",
  Adult_shared_core_DAR = TRUE, Adult_shared_direction,
  Adult_shared_log2FC, Adult_shared_padj,
  Adult_21y_effect, Adult_21y_effect_direction,
  Adult_67y_effect, Adult_67y_effect_direction,
  Adult_donor_match_n, adult_grade)]

## --- 逐 stage 宽列 (dev): Early/Mid/Late_is_DAR + direction/effect/padj (从 pair 主表直接取) ---
stage_wide <- pair[, .(celltype, gene, peak,
  Early_is_DAR = Early_DAR_status %in% c("M","P"), Early_DAR_direction = fifelse(Early_DAR_status %in% c("M","P"), fifelse(Early_DAR_status=="M","central","peripheral"), NA_character_),
    Early_DAR_effect = Early_log2FoldChange, Early_DAR_padj = Early_padj,
  Mid_is_DAR   = Mid_DAR_status   %in% c("M","P"), Mid_DAR_direction   = fifelse(Mid_DAR_status   %in% c("M","P"), fifelse(Mid_DAR_status=="M","central","peripheral"),   NA_character_),
    Mid_DAR_effect = Mid_log2FoldChange, Mid_DAR_padj = Mid_padj,
  Late_is_DAR  = Late_DAR_status  %in% c("M","P"), Late_DAR_direction  = fifelse(Late_DAR_status  %in% c("M","P"), fifelse(Late_DAR_status=="M","central","peripheral"),  NA_character_),
    Late_DAR_effect = Late_log2FoldChange, Late_DAR_padj = Late_padj,
  Adult_shared_core_DAR, Adult_shared_direction, Adult_shared_log2FC, Adult_shared_padj,
  Adult_21y_effect, Adult_21y_effect_direction, Adult_67y_effect, Adult_67y_effect_direction,
  Adult_donor_match_n, adult_grade)]

## --- pair_dar 合并 (dev_sum + ad_sum): peak-level 判定行 ---
pair_dar <- rbindlist(list(dev_sum, ad_sum), use.names = TRUE, fill = TRUE)
cat("判为 DAR 的 (ct,gene,peak): dev", nrow(dev_sum), "/ adult", nrow(ad_sum),
    "| dev 方向冲突 peak:", sum(dev_sum$DAR_direction_conflict), "\n")

## --- P2G distance: ArchR 原始 strict links (正式字段) 按 (celltype,peak,gene) 关联 ---
p2g_dist <- rbindlist(lapply(names(ARCHR_LINKS_F), function(p)
  fread(ARCHR_LINKS_F[[p]], select = c("celltype","peak","gene","distance"))[, period := p]))
setnames(p2g_dist, "gene", "target_gene")
cat("ArchR links distance 行:", nrow(p2g_dist), "\n")
## ---- 5. QC7: peak 坐标格式统一性检查 ----
## 所有来源 (Phase 2 links / pairing 表) 的 peak 均应为 "chrN-start-end" 格式; 检查三处是否一致
pat <- "^chr[0-9XYM]+-[0-9]+-[0-9]+$"
cat("QC7 峰坐标格式检查 (chr-seq-start-end):\n")
cat("  supporting peaks 符合格式:", mean(grepl(pat, sup$supporting_peak)),
    "| 非标准例:", if(any(!grepl(pat, sup$supporting_peak))) paste(head(sup$supporting_peak[!grepl(pat, sup$supporting_peak)], 3), collapse=" ") else "无", "\n")
cat("  pairing peaks  符合格式:", mean(grepl(pat, pair$peak)), "\n")
cat("  坐标集合交集示例 (supporting ∩ pairing peak, 前匹配方式=精确字符串):",
    length(intersect(unique(sup$supporting_peak), unique(pair$peak))), "个共有 peak\n")
## ---- 6. Peak-level 精确交集: supporting peak ∩ 该 target DEG 对应的 DAR (2026-10-09 rev2) ----
## sup 全部 supporting peaks 为底座; 左连接 pair_dar (period×ct×gene×peak 精确键) 与 stage_wide (完整逐期信息).
setnames(sup, "supporting_peak", "peak")     # 临时统一列名; 稍后恢复 supporting_peak
setnames(pair_dar, "gene", "target_gene")    # 统一列名后四键左连
pk <- merge(sup, pair_dar, by = c("period","celltype","peak","target_gene"), all.x = TRUE)
## 未匹配 → DAR_* 为 NA → is_DAR=FALSE (非 DAR supporting peaks 完整保留, 不按基因过滤丢弃)
pk[, supporting_peak := peak][, peak := NULL]
pk[, TF_direction := direction][, direction := NULL]
pk[, is_DAR := !is.na(DAR_direction)]

## 挂逐 stage 宽列 (dev 的 Early/Mid/Late + adult 的 shared/donor annotation; 按 ct+gene+peak)
stage_wide_key <- copy(stage_wide[, .(celltype, gene, peak,
  Early_is_DAR, Early_DAR_direction, Early_DAR_effect, Early_DAR_padj,
  Mid_is_DAR, Mid_DAR_direction, Mid_DAR_effect, Mid_DAR_padj,
  Late_is_DAR, Late_DAR_direction, Late_DAR_effect, Late_DAR_padj)])   # Adult_* 已由 ad_sum 提供, 避免 merge 冲突
setnames(stage_wide_key, c("gene","peak"), c("target_gene","supporting_peak"))
pk <- merge(pk, unique(stage_wide_key, by = c("celltype","target_gene","supporting_peak")),
            by = c("celltype","target_gene","supporting_peak"), all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)

## 挂 P2G distance (period × ct × peak × target_gene 精确键, 来自 ArchR 正式链接表)
setnames(p2g_dist, "peak", "supporting_peak")
pk <- merge(pk, p2g_dist, by = c("period","celltype","supporting_peak","target_gene"), all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
pk[, P2G_distance := distance][, distance := NULL]
if(all(is.na(pk$P2G_distance))) cat("[note] P2G distance unavailable in source table\n")

## 期别字段隔离 (2026-10-09 rev2): dev 期字段 (Early/Mid/Late_*) 仅 dev 行有效; adult 行置 NA (反之亦然)
## — 避免 adult 行出现 dev 期交叉信息造成误读
dev_only_cols  <- c("Early_is_DAR","Early_DAR_direction","Early_DAR_effect","Early_DAR_padj",
                    "Mid_is_DAR","Mid_DAR_direction","Mid_DAR_effect","Mid_DAR_padj",
                    "Late_is_DAR","Late_DAR_direction","Late_DAR_effect","Late_DAR_padj",
                    "DAR_direction_conflict","representative_DAR_stage",
                    "DAR_stage_vs_DEG_direction_list","DAR_stage_vs_TF_direction_list")
## 注: DAR_stage_list/DAR_direction_list/DAR_stage_direction_list/DAR_effect_list/DAR_padj_list 为两期通用
## (dev = stage 列表; adult = "adult" 单条), 不做期别隔离
adult_only_cols <- c("Adult_shared_core_DAR","Adult_shared_direction","Adult_shared_log2FC","Adult_shared_padj",
                     "Adult_21y_effect","Adult_21y_effect_direction","Adult_67y_effect","Adult_67y_effect_direction",
                     "Adult_donor_match_n","adult_grade","Adult_shared_direction_list")
for(cc in dev_only_cols) pk[period == "adult", (cc) := NA]
for(cc in adult_only_cols) if(cc %in% colnames(pk)) pk[period == "dev", (cc) := NA]
## 注意: 上方 pair_dar 的 ad_sum 已带 Adult_* 列, stage_wide merge 不再提供 Adult_*; 二者来源统一

## 挂 subclass / DEG_direction (来自 Phase 4 展开表) — 先挂方向信息再做 annotation
p4key <- unique(p4g[, .(period, celltype, subclass, DEG_direction, motif_id, target_gene)])
pk <- merge(pk, p4key, by = c("period","celltype","motif_id","target_gene"), all.x = TRUE, allow.cartesian = TRUE)
pk <- pk[!is.na(subclass)]     # 只保留 Phase 4 关系内的行

## 方向 annotation (不做筛选): DAR_vs_DEG / DAR_vs_TF + dev 逐 stage 版
tf_dir_tok <- fifelse(pk$TF_direction == "central-high", "central", "peripheral")
pk[, DAR_vs_DEG_direction := fifelse(is_DAR, fifelse(DAR_direction == DEG_direction, "same", "opposite"), NA_character_)]
pk[, DAR_vs_TF_direction  := fifelse(is_DAR, fifelse(DAR_direction == tf_dir_tok, "same", "opposite"), NA_character_)]
## dev 逐 stage vs DEG/TF (annotation): "Early:same;Mid:opposite"
mk_stage_list <- function(sdir_col, ref){
  out <- character(nrow(pk))
  for(st in c("Early","Mid","Late")){
    sd <- pk[[paste0(st, "_DAR_direction")]]
    hit <- !is.na(sd)
    if(!any(hit)) next
    cmp <- fifelse(sd[hit] == ref[hit], "same", "opposite")
    out[hit] <- ifelse(out[hit] == "", paste0(st, ":", cmp), paste0(out[hit], ";", paste0(st, ":", cmp)))
  }
  out[out == ""] <- NA_character_
  out
}
pk[, DAR_stage_vs_DEG_direction_list := mk_stage_list("dir", DEG_direction)]
pk[, DAR_stage_vs_TF_direction_list  := mk_stage_list("dir", tf_dir_tok)]

setcolorder(pk, c("period","celltype","subclass","motif_id","motif_name","TF_direction",
                  "target_gene","DEG_direction","supporting_peak","P2G_cor","P2G_FDR","P2G_distance",
                  "is_DAR","DAR_direction","DAR_effect","DAR_padj",
                  "DAR_stage_list","DAR_direction_list","DAR_stage_direction_list",
                  "DAR_effect_list","DAR_padj_list","DAR_direction_conflict","representative_DAR_stage",
                  "Early_is_DAR","Early_DAR_direction","Early_DAR_effect","Early_DAR_padj",
                  "Mid_is_DAR","Mid_DAR_direction","Mid_DAR_effect","Mid_DAR_padj",
                  "Late_is_DAR","Late_DAR_direction","Late_DAR_effect","Late_DAR_padj",
                  "Adult_shared_core_DAR","Adult_shared_direction","Adult_shared_log2FC","Adult_shared_padj",
                  "Adult_21y_effect","Adult_21y_effect_direction","Adult_67y_effect","Adult_67y_effect_direction",
                  "Adult_donor_match_n","adult_grade",
                  "DAR_vs_DEG_direction","DAR_vs_TF_direction",
                  "DAR_stage_vs_DEG_direction_list","DAR_stage_vs_TF_direction_list"))
## 规范化: DAR 相关列表/方向字段中空字符串统一为 NA (merge 过程 NA→"" 的清理)
for(cc in c("DAR_stage_list","DAR_direction_list","DAR_stage_direction_list","DAR_effect_list","DAR_padj_list",
            "DAR_direction","DAR_stage_vs_DEG_direction_list","DAR_stage_vs_TF_direction_list",
            "DAR_vs_DEG_direction","DAR_vs_TF_direction","representative_DAR_stage",
            "Early_DAR_direction","Mid_DAR_direction","Late_DAR_direction","Adult_shared_direction",
            "Adult_21y_effect_direction","Adult_67y_effect_direction"))
  if(cc %in% colnames(pk)) pk[get(cc) == "", (cc) := NA_character_]

fwrite(pk, file.path(DIR_TAB, sprintf("%s_Phase5_TF_DEG_supporting_peak_DAR_overlap_peaklevel.tsv", DATE_PREFIX)), sep = "\t")
cat("输出 A: peak-level 表", nrow(pk), "行 (含逐 stage DAR + adult donor annotation + P2G_distance)\n")
## ---- 7. QC4-6 + rev2 新 QC (peak-level 核查) ----
cat("QC4: supporting peaks 总数:", nrow(pk), "| 其中 DAR:", sum(pk$is_DAR),
    "| 比例:", round(mean(pk$is_DAR), 3), "\n")
print(pk[, .(n_peaks = .N, n_DAR = sum(is_DAR), frac = round(mean(is_DAR), 3)), by = period])
print(pk[, .(n_peaks = .N, n_DAR = sum(is_DAR), frac = round(mean(is_DAR), 3)), by = .(period, celltype)][order(period, celltype)])
cat("\nP2G_distance 覆盖:", sum(!is.na(pk$P2G_distance)), "/", nrow(pk), "\n")

## dev 新 QC
devpk <- pk[period == "dev"]
cat("\n[dev QC-1] Early/Mid/Late 中有 DAR 的 peak 数:", sum(devpk$is_DAR), "\n")
for(st in c("Early","Mid","Late")) cat("  ", st, "_is_DAR TRUE:", sum(devpk[[paste0(st,"_is_DAR")]] %in% TRUE, na.rm=TRUE), "\n")
## 2026-10-09 rev3 修复: 逐 stage 计数 (不用总体 is_DAR 混入)
n_stage_dar <-
    as.integer(devpk$Early_is_DAR %in% TRUE) +
    as.integer(devpk$Mid_is_DAR %in% TRUE) +
    as.integer(devpk$Late_is_DAR %in% TRUE)
cat("[dev QC-2] 同一 peak 在 >=2 个 stage 为 DAR:", sum(n_stage_dar >= 2), "\n")
cat("[dev QC-2b] n_stage_dar 分布 (supporting peaks 在 1/2/3 个 stage 为 DAR 的数量):\n")
print(table(n_stage_dar))
cat("[dev QC-3] DAR_direction_conflict == TRUE 的 DAR peak:", sum(devpk$is_DAR & devpk$DAR_direction_conflict %in% TRUE), "\n")
## QC-4/5: representative 同源验证 — 用 stage_wide 独立重算对比
chk <- devpk[is_DAR == TRUE]
chk_ok4 <- all(chk$representative_DAR_stage %in% c("Early","Mid","Late"))
cat("[dev QC-4] representative_DAR_stage ∈ 真实 stage:", chk_ok4, "\n")
## 每行验证: representative 的 effect == 该 stage 的 effect
rep_eff_ok <- vapply(seq_len(nrow(chk)), function(i){
  st <- chk$representative_DAR_stage[i]
  unname(chk[[paste0(st, "_DAR_effect")]][i] == chk$DAR_effect[i])
}, logical(1))
rep_dir_ok <- vapply(seq_len(nrow(chk)), function(i){
  st <- chk$representative_DAR_stage[i]
  unname(chk[[paste0(st, "_DAR_direction")]][i] == chk$DAR_direction[i])
}, logical(1))
cat("[dev QC-5] representative direction/effect 与该 stage 逐行一致:", all(rep_dir_ok), "/", all(rep_eff_ok), "\n")

## adult 新 QC
adpk <- pk[period == "adult"]
cat("\n[adult QC-1] is_DAR 与 Adult_shared_core_DAR 完全一致:",
    all(adpk$is_DAR == (adpk$Adult_shared_core_DAR %in% TRUE)), "\n")
cat("[adult QC-2] DAR_direction 与 Adult_shared_direction 完全一致:",
    all(adpk[is_DAR == TRUE]$DAR_direction == fifelse(adpk[is_DAR == TRUE]$Adult_shared_direction == "M", "central", "peripheral")), "\n")
cat("[adult QC-3] 21y/67y 仅 annotation (不参与 is_DAR 定义):",
    "已按 Adult_shared_core_DAR 定义 is_DAR — 21y/67y 列独立保留\n")

cat("\n方向 annotation 分布 (无筛选):\n")
print(pk[is_DAR == TRUE, .N, by = .(DAR_vs_DEG_direction, DAR_vs_TF_direction)])
## ---- 8. TF-target DEG 层面汇总 ----
tfdeg <- pk[, .(
  n_supporting_peaks       = .N,
  n_DAR_supporting_peaks   = sum(is_DAR),
  DAR_support_fraction     = round(sum(is_DAR) / .N, 4),
  has_DAR_support          = sum(is_DAR) > 0,
  supporting_peak_list     = paste(sort(unique(supporting_peak)), collapse = ";"),
  DAR_supporting_peak_list = paste(sort(unique(supporting_peak[is_DAR])), collapse = ";"),
  DAR_direction_list       = paste(unique(DAR_direction[is_DAR]), collapse = ";"),
  ## rev2 新增: dev 逐 stage 方向轨迹 + 冲突计数; adult shared 方向列表 (不筛选)
  DAR_stage_direction_list = paste(unique(DAR_stage_direction_list[is_DAR]), collapse = ";"),
  n_DAR_direction_conflict_peaks = sum(DAR_direction_conflict %in% TRUE),
  Adult_shared_direction_list = paste(unique(Adult_shared_direction[!is.na(Adult_shared_direction)]), collapse = ";"),
  TF_direction             = TF_direction[1],
  DEG_direction            = DEG_direction[1]),
  by = .(period, celltype, subclass, motif_id, motif_name, target_gene)]
setcolorder(tfdeg, c("period","celltype","subclass","motif_id","motif_name","TF_direction",
                     "target_gene","DEG_direction","n_supporting_peaks","n_DAR_supporting_peaks",
                     "DAR_support_fraction","has_DAR_support",
                     "supporting_peak_list","DAR_supporting_peak_list","DAR_direction_list",
                     "DAR_stage_direction_list","n_DAR_direction_conflict_peaks","Adult_shared_direction_list"))
fwrite(tfdeg, file.path(DIR_TAB, sprintf("%s_Phase5_TF_DEG_DAR_support_summary.tsv", DATE_PREFIX)), sep = "\t")
cat("输出 B: TF-DEG summary", nrow(tfdeg), "行 →",
    sprintf("%s_Phase5_TF_DEG_DAR_support_summary.tsv", DATE_PREFIX), "\n")
cat("\nhas_DAR_support 分布:\n"); print(tfdeg[, .N, by = .(period, has_DAR_support)])
cat("\nQC5: DAR_support_fraction ∈ [0,1]:", all(tfdeg$DAR_support_fraction >= 0 & tfdeg$DAR_support_fraction <= 1), "\n")
cat("QC6: n_DAR <= n_support 全部成立:", all(tfdeg$n_DAR_supporting_peaks <= tfdeg$n_supporting_peaks), "\n")
## ---- 9. TF × DEG group 层面描述性汇总 (无 cutoff) ----
deggrp <- tfdeg[, .(
  n_target_DEGs               = .N,
  n_DAR_supported_target_DEGs = sum(has_DAR_support),
  DAR_supported_target_fraction = round(sum(has_DAR_support) / .N, 4),
  DAR_supported_target_genes  = paste(sort(unique(target_gene[has_DAR_support])), collapse = ";")),
  by = .(period, celltype, subclass, DEG_direction, motif_id, motif_name)]
fwrite(deggrp, file.path(DIR_TAB, sprintf("%s_Phase5_TF_DEGgroup_DAR_support_summary.tsv", DATE_PREFIX)), sep = "\t")
cat("输出 C: DEG-group summary", nrow(deggrp), "行 →",
    sprintf("%s_Phase5_TF_DEGgroup_DAR_support_summary.tsv", DATE_PREFIX), "\n")
cat("\n=== TF × DEG group 的 DAR-supported target fraction 分布 (descriptive) ===\n")
print(deggrp[, .(n_groups = .N, med_frac = round(median(DAR_supported_target_fraction), 3),
                 q90_frac = round(quantile(DAR_supported_target_fraction, .9), 3)), by = period])
cat("\n=== PoC (MG/Cone) — Version B, 按 fraction 降序 top ===\n")
print(head(deggrp[celltype %in% POC][order(-DAR_supported_target_fraction, -n_target_DEGs),
      .(period, celltype, subclass, DEG_direction, TF = motif_name,
        n_target_DEGs, n_DAR_supported_target_DEGs, DAR_supported_target_fraction)], 15))
## ---- 11. 收尾: QC 总览 + 解释边界 ----
cat("\n==== 2026-10-09 Phase 5 第一步 QC 总览 ====\n")
cat("QC1  TF-DEG 关系数 (Phase 4 展开):", nrow(p4g), "\n")
cat("QC2  成功回查 Phase 2 supporting peaks:", nrow(matched), "/", nrow(key_p4),
    "| 未匹配:", nrow(unmatched), "\n")
cat("QC3  总 supporting peaks:", nrow(pk), "\n")
cat("QC4  DAR peaks:", sum(pk$is_DAR), "| 比例:", round(mean(pk$is_DAR), 3),
    "(dev:", sum(pk[period=="dev"]$is_DAR), "/", nrow(pk[period=="dev"]),
    "; adult:", sum(pk[period=="adult"]$is_DAR), "/", nrow(pk[period=="adult"]), ")\n")
cat("QC5  fraction ∈ [0,1]:", all(tfdeg$DAR_support_fraction >= 0 & tfdeg$DAR_support_fraction <= 1), "\n")
cat("QC6  n_DAR <= n_support:", all(tfdeg$n_DAR_supporting_peaks <= tfdeg$n_supporting_peaks), "\n")
cat("QC7  坐标格式统一 (chr-start-end):", mean(grepl("^chr[0-9XYM]+-[0-9]+-[0-9]+$", pk$supporting_peak)) == 1, "\n")
cat("QC8  P2G_distance 覆盖:", sum(!is.na(pk$P2G_distance)), "/", nrow(pk), " (ArchR 正式字段)\n")
cat("QC9  dev representative 同源:", all(rep_dir_ok), "/", all(rep_eff_ok),
    "| 方向冲突 peak:", sum(devpk$is_DAR & devpk$DAR_direction_conflict %in% TRUE), "\n")
cat("QC10 adult is_DAR==Adult_shared_core_DAR:", all(adpk$is_DAR == (adpk$Adult_shared_core_DAR %in% TRUE)), "\n")
cat("\nrev2 补充: dev 完整保留 Early/Mid/Late 逐 stage DAR (含 direction 冲突标记) + representative 同源四字段;\n")
cat("adult 正式 DAR = Adult_shared (21y/67y 仅 donor-level annotation); P2G_distance 来自 ArchR 正式链接表.\n")
cat("\n解释边界: 本步只增加 DAR-supported CRE evidence 层 — '部分 supporting peak 同时表现为 regional DAR';\n")
cat("不表述 TF 已证明调控该 DEG; 不据 TF/DAR/DEG direction 判断 activation/repression; 未做筛选 (无 cutoff).\n")
cat("下一步 (未执行): direction pattern → representative CRE → locus track → footprint → virtual perturbation.\n")
writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, sprintf("%s_Phase5_1_sessionInfo.txt", DATE_PREFIX)))
cat("\nPhase 5 第一步完成.\n")
