# 最终整理脚本 12 · TF targets 与 DM/DAL/AAM DEG 的方向分层 overlap (2026-10-29)
# 来源: 00use_scripts/29_2026-09-16-4_Phase3_TF_target_DEG_enrichment.ipynb
#       (仅 direct overlap + direction-stratified/group overlap 两段 — 其输出是候选 TF 链 (11) 的直接输入;
#        超几何富集主分析/PoC bubble/诊断/TFE3 小结/方向无关 [TEST] 段未进入最终链, 已去除)
# 产出: Table/2026-09-16-4_{dev,adult}_TF_target_DM_DAL_AAM_direct_overlap_test.tsv
#       Table/2026-09-16-4_{dev,adult}_TF_target_DM_DAL_AAM_direction_stratified_overlap_test.tsv
#       Table/2026-09-16-4_{dev,adult}_TF_target_DM_DAL_AAM_direction_group_overlap_test.tsv (11 的输入)
# 口径: 纯描述 overlap (无检验); TF regional direction 与 DEG direction 映射 same/opposite; 完整骨架
# 输入: 09 的 peak2gene targets + 五期轨迹 refined DEG 表 (DM70/DAL126/AAM169)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales) })
options(stringsAsFactors=FALSE, width=160)
DATE_PREFIX <- "2026-09-16-4"; UP_PREFIX <- "2026-09-16-3"
ROOT <- "/project/work"; OUT <- file.path(ROOT,"2026-09-16_ATAC_start")
DIR_TAB<-file.path(OUT,"Table"); DIR_DATA<-file.path(OUT,"Data"); DIR_PLOT<-file.path(OUT,"Plot"); DIR_LOG<-file.path(OUT,"Log")
DEG_DIR <- "/project/resources/2026-09-01_DEG_summary/Data"
DEG_TRAJ <- file.path(DEG_DIR,"2026-09-01_1_five_stage_trajectory_refined.csv")

# ---- Parameters ----
MIN_T <- 5L      # target set 在 universe 内 >=5 gene 才检验
MIN_S <- 3L      # subclass set 在 universe 内 >=3 gene 才检验
SIG_FDR <- 0.10
POC <- c("MG","Cone")
SUBORD <- c("DM","DAL","AAM")

source_color_ipynb <- function(path=file.path(ROOT,"color.ipynb")){
  j<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  codes<-vapply(j$cells,function(cc) if(identical(cc$cell_type,"code")) paste(unlist(cc$source),collapse="") else "",character(1))
  eval(parse(text=paste(codes,collapse="\n")),envir=globalenv())
}
invisible(tryCatch(suppressWarnings(source_color_ipynb()),error=function(e) message("color.ipynb: ",conditionMessage(e))))
save_plot<-function(p,name,w=8,h=6){
  suppressMessages(ggsave(file.path(DIR_PLOT,sprintf("%s_%s.pdf",DATE_PREFIX,name)),p,width=w,height=h))
  suppressMessages(ggsave(file.path(DIR_PLOT,sprintf("%s_%s.png",DATE_PREFIX,name)),p,width=w,height=h,dpi=300))
  cat("saved:",sprintf("%s_%s",DATE_PREFIX,name),"\n"); invisible(p)}
cat("Params: MIN_T=",MIN_T," MIN_S=",MIN_S," SIG_FDR=",SIG_FDR,"\n")
## ---- 9. [TEST 2026-10-08 rev2] Direct overlap: TF targets vs DM/DAL/AAM (无检验, 纯描述; 完整骨架) ----
# 分析目的: 查看三组 DEG (DM/DAL/AAM) 的交集候选 — 直接比较 Phase 2 TF candidate target genes 与
#           DM70/DAL126/AAM169 的重叠, 用于筛选三组 DEG 潜在的调控 TF.
# 方法定义 (与富集检验严格区分):
#   - 不区分 DEG central-up/peripheral-up 方向, 不做方向匹配 (无 matched/opposite);
#   - 不计算 Fisher / 超几何 / OR / P / FDR — 只作描述性 overlap analysis;
#   - DM/DAL/AAM 使用完整 canonical refined sets (70/126/169);
#   - Version A/B 分开; target 基因直取 Phase 2 输出原值 (不交 universe);
#   - TF regional direction (central-high/peripheral-high) 仅作为 annotation 保留, 不参与交集筛选;
#   - 完整组合骨架 period×version×celltype×TF×subclass(该 ct 实际存在的 DM/DAL/AAM), 无交集填
#     overlap_n=0 / overlap_genes="" / fraction=0 (分母缺失才 NA) — 便于完整比较同一 TF 对三组的覆盖模式;
#   - 排名三指标并列 (overlap_n / target_overlap_fraction / subclass_coverage_fraction).
map9 <- c("Developmental maintained"="DM","Developmental aging-lost"="DAL","Adult acquired and maintained"="AAM")
tr9 <- fread(DEG_TRAJ, select=c("celltype","gene","subclass_refined"))
tr9[, subclass := map9[subclass_refined]]
deg9 <- unique(tr9[!is.na(subclass), .(celltype, subclass, gene)])
cat("canonical refined sets (celltype,gene 行数, 应 DM70/DAL126/AAM169):\n")
print(deg9[, .N, by=subclass][order(match(subclass, SUBORD))])
stopifnot(deg9[subclass=="DM", .N]==70L, deg9[subclass=="DAL", .N]==126L, deg9[subclass=="AAM", .N]==169L)
## 每 celltype 实际存在的 subclass
ct_sub <- deg9[, .(subclasses=paste(sort(unique(subclass)), collapse=",")), by=celltype]

overlap_period <- function(period){
  outl <- list()
  for(version in c("A","B")){
    vcol <- if(version=="A") "versionA" else "versionB"
    tf <- fread(file.path(DIR_TAB, sprintf("%s_%s_TF_peak2gene_targets.tsv", UP_PREFIX, period)),
                select=c("celltype","motif_id","motif_name","direction","gene","versionA","versionB"))
    tv <- unique(tf[get(vcol)==TRUE, .(celltype,motif_id,motif_name,direction,gene)])
    ## TF 级 target 数与 direction (annotation; 同 motif 在同 ct 的 direction 唯一, Phase 1 定义)
    Tn <- tv[, .(target_gene_n=uniqueN(gene)), by=.(celltype,motif_id,motif_name,direction)]
    ## 完整骨架: TF × 该 celltype 实际存在的 subclass
    base <- merge(Tn, ct_sub, by="celltype", allow.cartesian=TRUE)
    base[, subclass := NULL]
    base <- base[, .(celltype,motif_id,motif_name,direction,target_gene_n,
                     subclass=unlist(strsplit(subclasses, ","))), by=seq_len(nrow(base))]
    base <- unique(base[, .(celltype,motif_id,motif_name,direction,target_gene_n,subclass)])
    ## subclass set 大小
    Ssz <- deg9[, .(subclass_gene_n=uniqueN(gene)), by=.(celltype,subclass)]
    base <- merge(base, Ssz, by=c("celltype","subclass"), all.x=TRUE)
    ## 实际 overlap left join
    ov <- merge(tv, deg9, by=c("celltype","gene"), allow.cartesian=TRUE)
    ovg <- ov[, .(overlap_n=uniqueN(gene), overlap_genes=paste(sort(unique(gene)),collapse=";")),
              by=.(celltype,motif_id,motif_name,subclass)]
    res <- merge(base, ovg, by=c("celltype","motif_id","motif_name","subclass"), all.x=TRUE)
    res[is.na(overlap_n), overlap_n:=0L]
    res[is.na(overlap_genes), overlap_genes:=""]
    res[, target_overlap_fraction  := fifelse(is.na(target_gene_n), NA_real_,
                                        round(overlap_n/pmax(target_gene_n,1L),4))]
    res[, subclass_coverage_fraction := fifelse(is.na(subclass_gene_n), NA_real_,
                                        round(overlap_n/pmax(subclass_gene_n,1L),4))]
    res[, `:=`(period=period, version=version)]
    setcolorder(res, c("period","version","celltype","motif_id","motif_name","direction","subclass",
                       "target_gene_n","subclass_gene_n","overlap_n",
                       "target_overlap_fraction","subclass_coverage_fraction","overlap_genes"))
    outl[[version]] <- res
  }
  rbindlist(outl)
}
res_ol <- rbindlist(lapply(c("dev","adult"), overlap_period))
for(p in c("dev","adult")){
  fwrite(res_ol[period==p], file.path(DIR_TAB, sprintf("%s_%s_TF_target_DM_DAL_AAM_direct_overlap_test.tsv", DATE_PREFIX, p)), sep="\t")
}
cat("\n==== [TEST] direct overlap (完整骨架) 已写出 (_test): dev", nrow(res_ol[period=="dev"]),
    "行 / adult", nrow(res_ol[period=="adult"]), "行 ====\n")

## ---- QC ----
cat("\n=== [TEST][QC] 每 celltype×version 的 TF 数 / 组合数 / overlap=0 组合数 / 完整性 ===\n")
qc <- res_ol[, .(n_TF=uniqueN(motif_id),
                 n_rows=.N,
                 n_subclass=uniqueN(subclass),
                 n_zero=sum(overlap_n==0)), by=.(period,version,celltype)][order(period,version,celltype)]
print(qc)
## 完整性: 每 (period,version,celltype,TF) 的 subclass 组合数 == 该 ct 的 subclass 数?
chk <- merge(res_ol[, .(n_comb=uniqueN(subclass)), by=.(period,version,celltype,motif_id)],
             ct_sub[, .(n_expect=as.integer(length(strsplit(subclasses,",")[[1]]))), by=celltype],
             by="celltype")
cat("\n[QC] 缺失组合数 (n_comb < n_expect):", sum(chk$n_comb < chk$n_expect), " (应 0)\n")
cat("[QC] 超额组合数 (n_comb > n_expect):", sum(chk$n_comb > chk$n_expect), " (应 0)\n")
stopifnot(sum(chk$n_comb != chk$n_expect) == 0L)
cat("[QC] direction 列取值:", paste(sort(unique(res_ol$direction)), collapse=","), "(annotation only)\n")
cat("[QC] fraction 无 NA (分母均存在):", sum(is.na(res_ol$target_overlap_fraction)), "+", sum(is.na(res_ol$subclass_coverage_fraction)), "NA 格\n")

## summary: 每 celltype 的 DM/DAL/AAM DEG 数
cat("\n=== [TEST] DM/DAL/AAM DEG 数 by celltype (canonical) ===\n")
print(dcast(deg9[, .(n_genes=uniqueN(gene)), by=.(celltype,subclass)], celltype~subclass, value.var="n_genes", fill=0L))

## PoC top (三指标并列)
for(ct in POC) for(p in c("dev","adult")) for(v in c("A","B")){
  d <- res_ol[period==p & version==v & celltype==ct & overlap_n>0]
  if(!nrow(d)) next
  cat(sprintf("\n=== [TEST][PoC] %s %s Version %s — top overlap TF (by subclass_coverage_fraction) ===\n", p, ct, v))
  print(head(d[order(-subclass_coverage_fraction, -overlap_n),
      .(TF=motif_name, direction, subclass, target_gene_n, subclass_gene_n, overlap_n,
        target_overlap_fraction, subclass_coverage_fraction, overlap_genes=substr(overlap_genes,1,50))], 8))
}
cat("\n[TEST direct overlap rev2 done] 完整骨架 + direction annotation + 0 填充; 无统计检验; 原富集结果未动。\n")
## ---- 10. [TEST 2026-10-08] Direction-stratified direct overlap: TF targets vs DM/DAL/AAM by DEG direction ----
# 分析目的: 在 direction-agnostic direct overlap (cell 9) 基础上, 按 DEG regional direction 拆分, 判断
#           同一 celltype×TF×subclass 的 P2G targets 主要覆盖与 TF regional direction 同向还是反向的 DEG,
#           用于筛选调控模式相对单一、后续更易解释的 candidate TF.
# 定义: TF central-high → same=central-up DEG, opposite=peripheral-up; TF peripheral-high → 反之.
#       仅描述方向一致/相反, 不解释为 activation/repression.
# 方法学: 完整骨架 (period×version×celltype×TF×subclass, 同 cell 9); 无交集填 0; total=0 时 fraction=NA;
#         dominant_mode: same / opposite / balanced(相等且>0) / no_overlap; 不做任何统计检验.
# 复用: deg9 (cell 9 已建 canonical DM/DAL/AAM + direction), ct_sub; Phase 2 正式 targets (含 direction).
## 自包含修复 (2026-10-08): cell 9 的 tr9 未含 direction_class, 本 cell 自建 4 列版 (不依赖 cell 9 的 tr9),
## 保证 notebook 全量顺序执行不报错; 口径与已产出 tsv 完全一致 (独立脚本同 select).
tr9d <- fread(DEG_TRAJ, select=c("celltype","gene","subclass_refined","direction_class"))
tr9d[, subclass := map9[subclass_refined]]
dirmap9 <- c("Macular direction"="central","Peripheral direction"="peripheral")
deg9d <- unique(tr9d[!is.na(subclass), .(celltype, subclass, gene, deg_dir=dirmap9[direction_class])])
ct_sub9 <- deg9d[, .(subclasses=paste(sort(unique(subclass)), collapse=",")), by=celltype]


strat_overlap <- function(period){
  outl <- list()
  for(version in c("A","B")){
    vcol <- if(version=="A") "versionA" else "versionB"
    tf <- fread(file.path(DIR_TAB, sprintf("%s_%s_TF_peak2gene_targets.tsv", UP_PREFIX, period)),
                select=c("celltype","motif_id","motif_name","direction","gene","versionA","versionB"))
    tv <- unique(tf[get(vcol)==TRUE, .(celltype,motif_id,motif_name,direction,gene)])
    ## TF-level direction (annotation; 同 motif 同 ct 唯一)
    TFdir <- tv[, .(direction=direction[1L]), by=.(celltype,motif_id,motif_name)]
    ## 完整骨架 (同 cell 9): TF × 该 ct 实际存在的 subclass
    base <- merge(TFdir, ct_sub9, by="celltype", allow.cartesian=TRUE)
    base <- base[, .(celltype,motif_id,motif_name,direction,
                     subclass=unlist(strsplit(subclasses, ","))), by=seq_len(nrow(base))]
    base <- unique(base)
    ## overlap: target × deg (celltype,gene) → 按 TF 方向归类 same/opposite
    ov <- merge(tv, deg9d, by=c("celltype","gene"), allow.cartesian=TRUE)
    ## direction 来自 tv (TF 方向); deg_dir 来自 DEG. 显式取列避免 merge 后缀歧义
    tf_dir_col <- if("direction.x" %in% colnames(ov)) "direction.x" else "direction"
    ov[, rel := fifelse((get(tf_dir_col)=="central-high" & deg_dir=="central") |
                        (get(tf_dir_col)=="peripheral-high" & deg_dir=="peripheral"), "same", "opposite")]
    ovg <- ov[, .(n=uniqueN(gene), genes=paste(sort(unique(gene)),collapse=";")),
              by=.(celltype,motif_id,motif_name,subclass,rel)]
    same <- ovg[rel=="same", .(celltype,motif_id,motif_name,subclass, same_overlap_n=n, same_overlap_genes=genes)]
    oppo <- ovg[rel=="opposite", .(celltype,motif_id,motif_name,subclass, opposite_overlap_n=n, opposite_overlap_genes=genes)]
    res <- merge(base, same, by=c("celltype","motif_id","motif_name","subclass"), all.x=TRUE)
    res <- merge(res, oppo, by=c("celltype","motif_id","motif_name","subclass"), all.x=TRUE)
    res[is.na(same_overlap_n), same_overlap_n:=0L][is.na(same_overlap_genes), same_overlap_genes:=""]
    res[is.na(opposite_overlap_n), opposite_overlap_n:=0L][is.na(opposite_overlap_genes), opposite_overlap_genes:=""]
    res[, total_overlap_n := same_overlap_n + opposite_overlap_n]
    res[, same_fraction := fifelse(total_overlap_n==0L, NA_real_, round(same_overlap_n/total_overlap_n,4))]
    res[, opposite_fraction := fifelse(total_overlap_n==0L, NA_real_, round(opposite_overlap_n/total_overlap_n,4))]
    res[, dominant_mode := fifelse(total_overlap_n==0L, "no_overlap",
                             fifelse(same_overlap_n > opposite_overlap_n, "same",
                              fifelse(opposite_overlap_n > same_overlap_n, "opposite", "balanced")))]
    res[, dominance_fraction := fifelse(total_overlap_n==0L, NA_real_,
                                  round(pmax(same_overlap_n, opposite_overlap_n)/total_overlap_n, 4))]
    res[, `:=`(period=period, version=version)]
    setnames(res, "direction", "TF_direction")
    setcolorder(res, c("period","version","celltype","motif_id","motif_name","TF_direction","subclass",
                       "same_overlap_n","opposite_overlap_n","total_overlap_n",
                       "same_fraction","opposite_fraction","dominant_mode","dominance_fraction",
                       "same_overlap_genes","opposite_overlap_genes"))
    outl[[version]] <- res
  }
  rbindlist(outl)
}
res_ds <- rbindlist(lapply(c("dev","adult"), strat_overlap))
for(p in c("dev","adult")){
  fwrite(res_ds[period==p], file.path(DIR_TAB, sprintf("%s_%s_TF_target_DM_DAL_AAM_direction_stratified_overlap_test.tsv", DATE_PREFIX, p)), sep="\t")
}
cat("\n==== [TEST] direction-stratified overlap 已写出 (_test): dev", nrow(res_ds[period=="dev"]),
    "行 / adult", nrow(res_ds[period=="adult"]), "行 ====\n")

## ---- QC ----
cat("\n=== [TEST][QC] dominant_mode 分布 (period×version) ===\n")
print(res_ds[, .N, by=.(period,version,dominant_mode)][order(period,version,dominant_mode)])
chk <- merge(res_ds[, .(n_comb=uniqueN(subclass)), by=.(period,version,celltype,motif_id)],
             ct_sub9[, .(n_expect=as.integer(length(strsplit(subclasses,",")[[1]]))), by=celltype], by="celltype")
cat("[QC] 缺失/超额组合:", sum(chk$n_comb < chk$n_expect), "/", sum(chk$n_comb > chk$n_expect), " (应 0/0)\n")
stopifnot(sum(chk$n_comb != chk$n_expect) == 0L)
cat("[QC] total=0 行的 fraction/dominance 全 NA:", all(is.na(res_ds[total_overlap_n==0L]$same_fraction)),
    "| no_overlap 行数:", res_ds[dominant_mode=="no_overlap", .N], "\n")
cat("[QC] 分层加和一致性 (sum vs cell9 overlap_n): 待对比\n")
## 与 direction-agnostic (cell 9) 数值互检 (容错: tsv 不存在则跳过)
ag_f <- file.path(DIR_TAB, sprintf("%s_dev_TF_target_DM_DAL_AAM_direct_overlap_test.tsv", DATE_PREFIX))
if(file.exists(ag_f)){
  res_ol9 <- fread(ag_f)
  m <- merge(res_ds[period=="dev" & version=="A", .(celltype,motif_id,subclass, t_ds=total_overlap_n)],
             res_ol9[version=="A", .(celltype,motif_id,subclass, t_ag=overlap_n)], by=c("celltype","motif_id","subclass"))
  cat("[QC] dev-A 分层总和 == agnostic overlap_n:", all(m$t_ds == m$t_ag), "\n")
} else cat("[QC] agnostic tsv 不存在, 跳过互检\n")

## ---- Summary 1: 每 celltype×subclass 的 dominant 计数 ----
cat("\n=== [TEST] dominant_mode 计数 by celltype×subclass (仅 Version B) ===\n")
print(dcast(res_ds[version=="B" & dominant_mode!="no_overlap", .N, by=.(period,celltype,subclass,dominant_mode)],
            period+celltype+subclass~dominant_mode, value.var="N", fill=0L))

## ---- Summary 2: Version B 重点 — total 大且 dominance 高 ----
cat("\n=== [TEST] Version B — high-priority TF (total>=5 & dominance>=0.8) by period×celltype ===\n")
print(head(res_ds[version=="B" & total_overlap_n>=5 & dominance_fraction>=0.8][order(-dominance_fraction,-total_overlap_n),
    .(period,celltype,TF=motif_name,TF_direction,subclass,same_overlap_n,opposite_overlap_n,
      total_overlap_n,dominant_mode,dominance_fraction)], 20))

## ---- Summary 3: PoC (MG/Cone) 单独汇总 (Version B + Version A) ----
for(ct in POC) for(p in c("dev","adult")) for(v in c("A","B")){
  d <- res_ds[period==p & version==v & celltype==ct & total_overlap_n>0]
  if(!nrow(d)) next
  cat(sprintf("\n=== [TEST][PoC] %s %s Version %s — top by total×dominance ===\n", p, ct, v))
  print(head(d[order(-total_overlap_n, -dominance_fraction),
      .(TF=motif_name,TF_direction,subclass,same=same_overlap_n,opp=opposite_overlap_n,
        total=total_overlap_n,dominant_mode,dominance_fraction)], 8))
}

## ---- Summary 4: 明显 mixed pattern (same 与 opposite 都多, 后续可降优先级) ----
cat("\n=== [TEST] mixed-pattern TF (same>=5 & opposite>=5, Version B) — 建议降优先级 ===\n")
mx <- res_ds[version=="B" & same_overlap_n>=5 & opposite_overlap_n>=5]
print(head(mx[order(-(same_overlap_n+opposite_overlap_n)),
    .(period,celltype,TF=motif_name,TF_direction,subclass,same=same_overlap_n,opp=opposite_overlap_n,
      total=total_overlap_n,same_fraction,dominant_mode,dominance_fraction)], 20))
cat("\n[TEST direction-stratified done] 描述性分析: 无检验; same/opposite 仅指方向一致/相反, 不作 activation/repression 解释; 输出 *_direction_stratified_overlap_test.tsv; 原分析未动。\n")

## ---- 10b. [TEST 2026-10-08] Direction-group 展开: 结果整理到真实 DEG direction 层面 + 分母/coverage ----
# 将 same/opposite 按 TF 方向映射回真实 DEG_direction (TF central-high: same→central, opposite→peripheral;
# TF peripheral-high 反之), 每 period×version×celltype×TF×subclass×DEG_direction 输出:
#   overlap_n / overlap_genes / deg_group_n (该 ct×subclass×deg_dir 的 unique DEG 数) /
#   deg_group_coverage_fraction = overlap_n / deg_group_n / relation_to_TF_direction (same|opposite)
# 零检验; 无 overlap → 0/""/0; DEG group 不存在 → deg_group_n=NA, coverage=NA; Version A/B 分开.
## DEG 分组大小 (celltype × subclass × deg_dir)
deg_grp <- deg9d[, .(deg_group_n=uniqueN(gene)), by=.(celltype, subclass, deg_dir)]

## 长表展开: same 行 + opposite 行 → 真实 DEG_direction
ss <- res_ds[, .(period, version, celltype, motif_id, motif_name, TF_direction, subclass,
                 overlap_n=same_overlap_n, overlap_genes=same_overlap_genes, relation_to_TF_direction="same")]
oo <- res_ds[, .(period, version, celltype, motif_id, motif_name, TF_direction, subclass,
                 overlap_n=opposite_overlap_n, overlap_genes=opposite_overlap_genes, relation_to_TF_direction="opposite")]
res_dg <- rbindlist(list(ss, oo))
res_dg[, DEG_direction := fifelse(TF_direction=="central-high",
                           fifelse(relation_to_TF_direction=="same", "central", "peripheral"),
                           fifelse(relation_to_TF_direction=="same", "peripheral", "central"))]
## 挂分母 (left join; 该 ct×subclass×deg_dir 组不存在 → NA) — 先统一列名避免 merge 键歧义
setnames(deg_grp, "deg_dir", "DEG_direction")
res_dg <- merge(res_dg, deg_grp, by=c("celltype","subclass","DEG_direction"), all.x=TRUE)
res_dg[, deg_group_coverage_fraction := fifelse(is.na(deg_group_n), NA_real_,
                                          round(overlap_n/pmax(deg_group_n,1L), 4))]
setcolorder(res_dg, c("period","version","celltype","motif_id","motif_name","TF_direction","subclass",
                      "DEG_direction","relation_to_TF_direction","overlap_n","overlap_genes",
                      "deg_group_n","deg_group_coverage_fraction"))
for(p in c("dev","adult")){
  fwrite(res_dg[period==p], file.path(DIR_TAB, sprintf("%s_%s_TF_target_DM_DAL_AAM_direction_group_overlap_test.tsv", DATE_PREFIX, p)), sep="\t")
}
cat("\n==== [TEST] direction-group overlap 已写出 (_test): dev", nrow(res_dg[period=="dev"]),
    "行 / adult", nrow(res_dg[period=="adult"]), "行 ====\n")

## ---- QC ----
cat("\n=== [TEST][QC-1] deg_group_n by celltype×subclass×DEG_direction (canonical, 不分 period/version) ===\n")
print(dcast(deg_grp, celltype+subclass~DEG_direction, value.var="deg_group_n", fill=0L))
## QC-2: coverage ∈ [0,1]
cv <- res_dg$deg_group_coverage_fraction
cat("\n[QC-2] coverage 越界 (0 或 <0 或 >1) 的非 NA 格数:", sum(!is.na(cv) & (cv<0 | cv>1)), " (应 0)\n")
stopifnot(sum(!is.na(cv) & (cv<0 | cv>1)) == 0L)
## QC-3: 同一 TF×subclass 的 central+peripheral overlap == direction-agnostic overlap_n (QC-4: dev/adult × A/B 全查)
ag_f <- function(p) file.path(DIR_TAB, sprintf("%s_%s_TF_target_DM_DAL_AAM_direct_overlap_test.tsv", DATE_PREFIX, p))
all_ok <- TRUE
for(p in c("dev","adult")) for(v in c("A","B")){
  if(!file.exists(ag_f(p))){ cat(sprintf("[QC-3] %s-%s: agnostic tsv 缺失, 跳过\n", p, v)); next }
  ag <- fread(ag_f(p))[version==v]
  m <- merge(res_dg[period==p & version==v, .(s=sum(overlap_n)), by=.(celltype,motif_id,subclass)],
             ag[, .(celltype,motif_id,subclass, t_ag=overlap_n)], by=c("celltype","motif_id","subclass"))
  ok <- all(m$s == m$t_ag)
  all_ok <- all_ok && ok
  cat(sprintf("[QC-3/4] %s Version %s: central+peripheral == agnostic overlap_n: %s (%d 组合)\n", p, v, ok, nrow(m)))
}
cat("[QC-3] 全部一致:", all_ok, "\n")
cat("[QC] DEG group 不存在 (deg_group_n=NA) 行数:", sum(is.na(res_dg$deg_group_n)),
    "| 其 coverage 全 NA:", all(is.na(res_dg[is.na(res_dg$deg_group_n)]$deg_group_coverage_fraction)), "\n")
cat("\n[TEST direction-group done] 输出 *_TF_target_DM_DAL_AAM_direction_group_overlap_test.tsv; same_fraction/opposite_fraction/dominance 等既有列未动; 零统计检验。\n")
