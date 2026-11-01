# 最终整理脚本 11 · 候选 TF → 含 motif 且区域可及的 peaks → peak2gene targets (2026-10-29)
# 来源: 00use_scripts/28_2026-09-16-3_Phase2_TF_peak_gene.ipynb (仅组装主链;
#       PoC funnel/TFE3 小结/ndonor 敏感性分析(历史记录)等诊断段已去除)
# 产出: 2026-09-16_ATAC_start/Table/2026-09-16-3_{dev,adult}_TF_peak2gene_targets.tsv (Version A/B)
#       2026-09-16_ATAC_start/Table/2026-09-16-3_{dev,adult}_P2G_universe_by_celltype.tsv
# 口径: region %cells-accessible ≥ 0.05; per-donor 检出 dev ≥6 / adult ≥2; |ΔchromVAR| ≥ 0.10
# 输入: 07 的 chromVAR 统计 (候选 TF) + ArchR peak2gene strict links + motif presence 对象
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2); library(jsonlite); library(IRdisplay)
})
options(stringsAsFactors=FALSE, width=150)
DATE_PREFIX <- "2026-09-16-3"
ROOT <- "/project/work"; OUT <- file.path(ROOT,"2026-09-16_ATAC_start")
DIR_TAB<-file.path(OUT,"Table"); DIR_DATA<-file.path(OUT,"Data"); DIR_PLOT<-file.path(OUT,"Plot"); DIR_LOG<-file.path(OUT,"Log")

# ---- Parameters ----
ACCESS_MIN     <- 0.05          # region %cells-accessible 下限
DONOR_DET_MIN  <- c(dev=6L, adult=2L)   # per-donor detection 下限：发育要求多数 donor(≥6/12)；成年佐证(≥2, 2026-10-08 定稿；敏感性分析证实 >=1 vs >=2 结果高度稳健)
SIG_DELTA      <- 0.10
POC            <- c("MG","Cone")
celltype_order <- c("Rod","Cone","AC","BC","HC","MG","RGC","RPC")
DIRCOL <- c("central-high"="#D55E00","peripheral-high"="#0072B2")

source_color_ipynb <- function(path=file.path(ROOT,"color.ipynb")){
  j<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  codes<-vapply(j$cells,function(cc) if(identical(cc$cell_type,"code")) paste(unlist(cc$source),collapse="") else "",character(1))
  eval(parse(text=paste(codes,collapse="\n")),envir=globalenv())
}
invisible(tryCatch(source_color_ipynb(),error=function(e) message("color.ipynb: ",conditionMessage(e))))
if(exists("cols.direction")) DIRCOL <- c("central-high"=unname(cols.direction["Macular"]),"peripheral-high"=unname(cols.direction["Peripheral"]))
save_plot<-function(p,name,w=8,h=6){
  suppressMessages(ggsave(file.path(DIR_PLOT,sprintf("%s_%s.pdf",DATE_PREFIX,name)),p,width=w,height=h))
  suppressMessages(ggsave(file.path(DIR_PLOT,sprintf("%s_%s.png",DATE_PREFIX,name)),p,width=w,height=h,dpi=300))
  cat("saved:",sprintf("%s_%s",DATE_PREFIX,name),"\n"); invisible(p)}
cat("Params ready. ACCESS_MIN=",ACCESS_MIN," DONOR_DET_MIN dev/adult=",DONOR_DET_MIN["dev"],"/",DONOR_DET_MIN["adult"],"\n")
load_candidates <- function(period){
  f <- file.path(DIR_TAB, sprintf("2026-09-16-2_%s_regional_chromVAR.tsv", if(period=="dev")"dev" else "adult"))
  p1 <- fread(f)
  if(period=="dev"){
    cand <- p1[candidate==TRUE, .(celltype,motif_id,motif_name,TF_annotation,family,direction, effect=delta)]
  } else {
    cand <- p1[direction_consistent==TRUE & abs(mean_delta)>=SIG_DELTA,
               .(celltype,motif_id,motif_name,TF_annotation,family,direction, effect=mean_delta)]
  }
  cand[]
}
cat("dev candidates:", nrow(load_candidates("dev")), "\n")
assemble_period <- function(period){
  linksf <- file.path(DIR_TAB, sprintf("%s_%s_ArchR_P2G_links.tsv", DATE_PREFIX, period))
  mipf   <- file.path(DIR_DATA, sprintf("%s_%s_motif_in_peak_sparse.rds", DATE_PREFIX, period))
  pctf   <- file.path(DIR_DATA, sprintf("%s_%s_region_accessibility_pct.rds", DATE_PREFIX, period))
  ndetf  <- file.path(DIR_DATA, sprintf("%s_%s_region_accessibility_ndonor_detected.rds", DATE_PREFIX, period))
  if(!all(file.exists(linksf,mipf,pctf,ndetf))){ cat("[",period,"] 中间件不全，跳过\n"); return(NULL) }
  links <- fread(linksf); mip <- readRDS(mipf); pct <- readRDS(pctf); ndet <- readRDS(ndetf)
  mipcols <- colnames(mip); peaks <- rownames(mip)
  stopifnot(identical(rownames(pct), peaks))
  cand <- load_candidates(period)
  dmin <- DONOR_DET_MIN[[period]]

  tgt_list <- list(); map_list <- list(); funnel <- list()
  cts <- intersect(celltype_order, intersect(unique(cand$celltype), unique(links$celltype)))
  for(ct in cts){
    candc <- cand[celltype==ct]; if(!nrow(candc)) next
    linkc <- links[celltype==ct, .(peak,gene,score)]
    grp <- function(reg) paste(ct,reg,sep="__")
    have <- grp("central") %in% colnames(pct) & grp("peripheral") %in% colnames(pct)
    if(!have){ cat("[",period,ct,"] 缺 accessibility 组，跳过\n"); next }
    accdt <- data.table(peak=peaks,
                        central_pct=pct[,grp("central")],  central_ndet=ndet[,grp("central")],
                        periph_pct =pct[,grp("peripheral")],periph_ndet =ndet[,grp("peripheral")])
    # motif-in-peak 限制在 P2G peaks × 候选 motif
    lk_peaks <- unique(linkc$peak); pidx <- match(lk_peaks, peaks)
    midx <- match(candc$motif_id, mipcols)
    sub <- as.matrix(mip[pidx, midx, drop=FALSE])          # (P2G peaks) x (候选 motif) 逻辑
    w <- which(sub, arr.ind=TRUE)
    if(!nrow(w)){ next }
    dtmp <- data.table(peak=lk_peaks[w[,1]], motif_id=candc$motif_id[w[,2]])
    dtmp <- merge(dtmp, candc, by="motif_id", allow.cartesian=TRUE)
    dtmp <- merge(dtmp, accdt, by="peak", allow.cartesian=TRUE)
    # region 特异 accessibility 过滤（§0.4-5）
    dtmp[, acc_pct  := fifelse(direction=="central-high", central_pct,  periph_pct)]
    dtmp[, acc_ndet := fifelse(direction=="central-high", central_ndet, periph_ndet)]
    dtmp <- dtmp[acc_pct>=ACCESS_MIN & acc_ndet>=dmin]
    if(!nrow(dtmp)){ next }
    # 连 P2G gene（peak 多对多）
    tl <- merge(dtmp, linkc, by="peak", allow.cartesian=TRUE)
    tl[, period := period]   # 物化为列：by=.() 内标量不会被回收
    # peak-level map（仅 PoC 存盘，控体积）
    if(ct %in% POC) map_list[[ct]] <- tl[, .(period, celltype, motif_id, motif_name, TF_annotation, family,
                                             direction, peak, gene, corr=score, acc_pct=round(acc_pct,3), acc_ndet)]
    # gene 级聚合
    genes <- tl[, .(n_supporting_peaks=uniqueN(peak), max_corr=max(score), mean_corr=mean(score),
                    effect=effect[1], TF_annotation=TF_annotation[1], family=family[1]),
                by=.(period,celltype,motif_id,motif_name,direction,gene)]
    genes[, versionA:=TRUE][, versionB:=n_supporting_peaks>=2L]
    tgt_list[[ct]] <- genes
    funnel[[ct]] <- data.table(period,celltype=ct,
      n_candidate_motif=nrow(candc),
      n_motif_peak_in_P2G=uniqueN(dtmp$peak),
      n_targets_vA=uniqueN(genes[versionA==TRUE]$gene),
      n_targets_vB=uniqueN(genes[versionB==TRUE]$gene),
      n_motif_target_pairs=nrow(genes))
  }
  targets <- rbindlist(tgt_list); pmap <- rbindlist(map_list)
  fun <- rbindlist(funnel)
  # P2G universe（Phase 3 用）：该期每 celltype 的 linked genes
  universe <- unique(links[, .(celltype, gene)])[, period:=period]
  list(targets=targets, pmap=pmap, funnel=fun, universe=universe)
}
res <- list()
for(period in c("dev","adult")){
  cat("\n==== assemble", period, "====\n"); r <- assemble_period(period)
  if(is.null(r)) next
  res[[period]] <- r
  fwrite(r$targets,  file.path(DIR_TAB, sprintf("%s_%s_TF_peak2gene_targets.tsv", DATE_PREFIX, period)), sep="\t")
  fwrite(r$universe, file.path(DIR_TAB, sprintf("%s_%s_P2G_universe_by_celltype.tsv", DATE_PREFIX, period)), sep="\t")
  if(nrow(r$pmap)) fwrite(r$pmap, file.path(DIR_TAB, sprintf("%s_%s_TF_motif_peak_map_PoC.tsv", DATE_PREFIX, period)), sep="\t")
  cat("targets rows:", nrow(r$targets), "| PoC peak-map rows:", nrow(r$pmap), "\n")
  cat("funnel:\n"); print(r$funnel)
}
