# 最终整理脚本 20 · Fig4 DEGlink snATAC DAR 可及性热图 V7 (linkDAR_peakDEG, 列按 region 集中) (2026-10-27)
# 来源: 00use_scripts/20_2026-10-20_Fig4_DEGlink_heatmap_snATAC_DAR_merged_v6.R (仅 V7 段 + 前置行集/可及性构建)
# 产出: /project/work/2026-09-04_Fig4_heatmap/Plot/2026-10-20_Fig4_DEGlink_heatmap_CRE_access_linkDAR_peakDEG_merged_v7.pdf
# 行集: snRNA peakDEG 图中全部 DEG (121 行/116 基因) 的同期配对 DAR∩link CRE, 249 行 (central 90 上 + peripheral 159 下);
#       行序与 DEG 顺序一致: dir(central 上) → group → DEG 行序(order) → 同基因 peak 坐标
# 列: 5 celltype × 10 期-区 (50 列全保留); V7 列序 = central 25 列集中在前, peripheral 25 列在后
# 数值: ATAC TF-IDF data 层假体积 → 胚胎 30 列块 / 成年 20 列块各自行 z-score, 色阶 [-2,2] solarExtra
# linkDAR 口径 (与 v5 notebook cell 7 完全一致):
#   emb_ok = 胚胎任一 stage DAR(任意方向) ∧ develop_ArchR_link; adu_ok = 成年 core DAR ∧ adult_ArchR_link;
#   filter(emb_ok | adu_ok)
# 整理说明: 去除 V6 图 (linkDAR_MGonly / linkDAR_peakDEG v6) 绘制段, 仅保留最终收录的 peakDEG V7;
#           可及性新缓存优先, 缺失 peak 需加载 12G+31G snATAC 对象增量补抽 (~25 min, 首次运行后走缓存)
suppressMessages({ library(Seurat); library(Signac); library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(grid); library(data.table); library(ArchR); library(RColorBrewer) })

proj_dir <- "/project/work/2026-09-04_Fig4_heatmap/"
data_dir <- paste0(proj_dir, "Data/")
plot_dir <- paste0(proj_dir, "Plot/")
prefix_v5 <- "2026-09-28_Fig4_DEGlink_heatmap"
master_csv <- "/project/work/2026-08-08_candidate_CRE/Data/2026-08-12_DEG_DAR_ArchR_link_peak_level.csv"
atac_dev_rds <- "/project/work/2026-08-03_develop_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds"
atac_adu_rds <- "/project/work/2026-08-03_adult_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds"
CACHE_PEAK <- paste0(data_dir, "2026-10-20_Fig4_DEGlink_heatmap_ACC_linkDAR_peakDEG_matrix_v6.rds")

## ---- 列定义与配色 (hex 与 00_colors.R 逐项一致) ----
cts      <- c("Rod","Cone","AC","BC","MG")
emb_sr   <- c("Early_central","Early_peripheral","Mid_central","Mid_peripheral","Late_central","Late_peripheral")
adu_sr   <- c("Young_central","Young_peripheral","Old_central","Old_peripheral")
all_sr   <- c(emb_sr, adu_sr)
all_cols <- as.vector(t(outer(all_sr, cts, paste, sep="__")))
sr_vec <- setNames(rep(all_sr, each=5), all_cols)
ct_vec <- setNames(rep(cts, 10), all_cols)
rg_vec <- setNames(rep(c("central","peripheral")[match(sub(".*_","",all_sr), c("central","peripheral"))], each=5), all_cols)
st_vec <- setNames(rep(sub("_.*","",all_sr), each=5), all_cols)

cols.group    <- c(DM="#66C2A5", DAL="#FC8D62", AAM="#8DA0CB")        # = 00_colors.R cols.subclass
cols.celltype <- c(Rod="#A6CEE3", Cone="#1F78B4", AC="#B2DF8A", BC="#33A02C", HC="#FB9A99",
                   MG="#E31A1C", RGC="#FDBF6F", RPC="#FF7F00", Micro="#66C2A5", Astro="#FC8D62")  # = 00_colors.R cols.celltype
cols.stage    <- c(Early="#59A14F", Mid="#EDC948", Late="#B07AA1", Young="#4E79A7", Old="#E15759") # = 00_colors.R cols.stage
cols.region   <- c(central="#6A3D9A", peripheral="#E69F00")            # = 00_colors.R cols.region
z_fun <- colorRamp2(seq(-2, 2, length.out=256), paletteContinuous("solarExtra"))
group_order <- c("DM","DAL","AAM")
emb_cols_v3 <- intersect(all_cols[1:30], all_cols); adu_cols_v3 <- intersect(all_cols[31:50], all_cols)
zscore_block <- function(m){
  z <- m
  for(blk in list(emb_cols_v3, adu_cols_v3)){
    z[,blk] <- t(apply(m[,blk,drop=FALSE], 1, function(x){ s<-sd(x); if(is.na(s)||s==0) return(rep(0,length(x))); (x-mean(x))/s }))
  }
  z
}

## ---- 1. linkDAR_peakDEG 行集 (peakDEG 图全部 DEG 的同期配对 DAR∩link) ----
ro <- do.call(rbind, lapply(c("central","peripheral"), function(dr)
  transform(read.csv(paste0(data_dir, prefix_v5, "_gene_row_order_peak_v5_", dr, ".csv"), stringsAsFactors=FALSE), dir=dr)))
names(ro)[names(ro)=="celltype"] <- "celltype_sort"; ro$order <- ro$go; ro$go <- NULL
master <- as.data.frame(fread(master_csv, select=c(
  "celltype","gene","peak","develop_ArchR_link","adult_ArchR_link",
  "Early_DAR_status","Mid_DAR_status","Late_DAR_status","Adult_shared_core_DAR")))
ld_tab <- master %>% filter(gene %in% unique(ro$gene), celltype %in% cts) %>%
  inner_join(ro %>% select(dir, group, celltype_sort, gene, order), by=c("gene","celltype"="celltype_sort"))
ld_tab <- ld_tab %>% mutate(
  is_dev_dar = (!is.na(Early_DAR_status)) | (!is.na(Mid_DAR_status)) | (!is.na(Late_DAR_status)),
  is_adu_dar = !is.na(Adult_shared_core_DAR) & Adult_shared_core_DAR,
  emb_ok = is_dev_dar & (!is.na(develop_ArchR_link) & develop_ArchR_link),
  adu_ok = is_adu_dar & (!is.na(adult_ArchR_link) & adult_ArchR_link)) %>%
  filter(emb_ok | adu_ok)
pc <- t(sapply(strsplit(ld_tab$peak, "-"), function(v) c(v[1], as.integer(v[2]))))
ld_tab$chr <- pc[,1]; ld_tab$pstart <- as.integer(pc[,2])
ld_tab <- as.data.frame(ld_tab %>% arrange(match(dir, c("central","peripheral")), match(group, group_order), order, chr, pstart) %>% distinct())
cat("linkDAR (peakDEG 全基因) 行:", nrow(ld_tab), "| 唯一 peak:", length(unique(ld_tab$peak)), "\n")

## ---- 2. 可及性: 新缓存优先; 否则 v5 缓存并集 + 缺失 peak 增量补抽 (补抽需加载 12G+31G snATAC 对象, ~25 min) ----
if(file.exists(CACHE_PEAK)){
  cache_pk <- readRDS(CACHE_PEAK); ACC_pk <- cache_pk$mat; keep_cols_ref <- cache_pk$keep_cols
  cat("使用已有新缓存:", paste(dim(ACC_pk), collapse="x"), "\n")
} else {
pk_need <- unique(ld_tab$peak); have <- character(0); mats <- list(); keep_cols_ref <- NULL
for(cf in c("_ACC_linkDAR_matrix_v5.rds","_ACC_matrix_v5.rds","_ACC_allDAR_matrix_v5.rds","_ACC_DAR_matrix_v5.rds")){
  f <- paste0(data_dir, prefix_v5, cf)
  if(file.exists(f)){
    cc <- readRDS(f)
    if(is.null(keep_cols_ref)) keep_cols_ref <- cc$keep_cols
    hit <- setdiff(intersect(pk_need, rownames(cc$mat)), have)
    if(length(hit)){ mats[[length(mats)+1]] <- cc$mat[hit, keep_cols_ref, drop=FALSE]; have <- c(have, hit) }
  }
}
cat("缓存并集覆盖:", length(have), "/", length(pk_need), "| 需补抽:", length(setdiff(pk_need, have)), "\n")
pk_miss <- setdiff(pk_need, have)
pk_none <- character(0)
if(length(pk_miss)){
  slim_atac <- function(obj, peaks){
    keep <- intersect(names(obj@assays), "ATAC")
    obj@assays <- obj@assays[keep]
    DefaultAssay(obj) <- "ATAC"
    hit <- intersect(peaks, rownames(obj[["ATAC"]]))
    obj <- obj[hit, ]
    invisible(gc())
    obj
  }
  pb_acc_full <- function(obj, peaks, stage_map, region_map){
    m <- obj[[]]
    mg <- data.frame(cell=rownames(m), celltype=as.character(m$celltype),
                     stage=stage_map(m), region=region_map(m), stringsAsFactors=FALSE) %>%
          filter(celltype %in% cts, !is.na(stage), !is.na(region))
    dat <- GetAssayData(obj, assay="ATAC", layer="data")
    hit <- intersect(unique(peaks), rownames(dat))
    res <- list(); nc <- list()
    for(ct in cts){ for(sr in all_sr){
      st <- sub("_.*","",sr); rg <- sub(".*_","",sr)
      cells <- mg$cell[mg$celltype==ct & mg$stage==st & mg$region==rg]
      if(length(cells)==0){ nc[[paste(ct,sr)]] <- data.frame(celltype=ct, column=sr, n_cells=0); next }
      res[[paste(ct,sr)]] <- data.frame(feature=hit, column=sr, celltype=ct,
        value=as.numeric(Matrix::rowMeans(dat[hit,cells,drop=FALSE])), stringsAsFactors=FALSE)
      nc[[paste(ct,sr)]] <- data.frame(celltype=ct, column=sr, n_cells=length(cells))
    }}
    list(long=do.call(rbind,res), ncell=do.call(rbind,nc), hit=hit)
  }
  cat("补抽: 载入成年对象 (12G)...\n")
  obj <- readRDS(atac_adu_rds); obj <- slim_atac(obj, pk_miss)
  a <- pb_acc_full(obj, pk_miss,
    function(m){ da<-gsub(" ","",m$donor_age); dplyr::case_when(da %in% c("21years","25years")~"Young",
                                 da %in% c("67years","70years","75years","86years")~"Old") },
    function(m) ifelse(m$tissue %in% c("fovea","macula"),"central","peripheral"))
  rm(obj); invisible(gc())
  cat("成年: 命中", length(a$hit), "/", length(pk_miss), "\n")
  ## 注意: 胚胎对象负责填 Early/Mid/Late 列 — 即使成年已命中全部 peak 也必须加载胚胎对象抽取 (v5 同模式)
  cat("补抽: 载入胚胎对象 (31G)...\n")
  obj <- readRDS(atac_dev_rds); obj <- slim_atac(obj, pk_miss)
  e <- pb_acc_full(obj, pk_miss,
    function(m) dplyr::case_when(m$donor_age %in% c("70days","79days","87days","91days")~"Early",
                                 m$donor_age %in% c("100days","103days","116days")~"Mid",
                                 m$donor_age %in% c("137days","141days","142days","162days","165days")~"Late"),
    function(m) ifelse(m$tissue=="macula","central","peripheral"))
  rm(obj); invisible(gc())
  cat("胚胎: 命中", length(e$hit), "\n")
  pk_found <- union(a$hit, e$hit)
  pk_none <- setdiff(pk_miss, pk_found)
  if(length(pk_none)) cat("WARNING: 以下 peak 不在任何 snATAC 对象中, 对应行剔除:", paste(pk_none, collapse=", "), "\n")
  if(length(pk_found)){
    d <- bind_rows(a$long, e$long); d$column <- paste(d$column, d$celltype, sep="__")
    mm <- as.data.frame(d %>% select(feature, column, value) %>% pivot_wider(names_from=column, values_from=value))
    rownames(mm) <- mm$feature; mm$feature <- NULL
    miss_cols <- setdiff(keep_cols_ref, colnames(mm)); for(mc in miss_cols) mm[[mc]] <- NA_real_
    mats[[length(mats)+1]] <- as.matrix(mm[intersect(pk_found, rownames(mm)), keep_cols_ref, drop=FALSE])
  }
}
ACC_pk <- do.call(rbind, mats)
na_cnt <- sum(is.na(ACC_pk))
if(na_cnt > 0){
  na_pk <- unique(rownames(ACC_pk)[apply(is.na(ACC_pk), 1, any)])
  cat("WARNING:", na_cnt, "NA 值, 涉及 peak:", paste(na_pk, collapse=", "), "(单侧对象缺失, 图中灰色)\n")
}
saveRDS(list(mat=ACC_pk, keep_cols=keep_cols_ref), CACHE_PEAK)
cat("新缓存:", CACHE_PEAK, "|", paste(dim(ACC_pk), collapse="x"), "\n")
}

## ---- 3. 行表存档 (剔除无数据的 peak 行) ----
ld_pk <- ld_tab %>% filter(peak %in% rownames(ACC_pk))
write.csv(ld_pk %>% select(dir,group,celltype,gene,order,peak,chr,pstart,emb_ok,adu_ok),
          paste0(data_dir, "2026-10-20_Fig4_DEGlink_heatmap_CRE_row_meta_linkDAR_peakDEG_v6.csv"), row.names=FALSE)

## ---- 4. V7: 仅调整列序 — central 列集中在前, peripheral 在后; 行名列名不显示 ----
## 区内 stage 顺序 (Early→Mid→Late→Young→Old) 与 celltype 顺序 (Rod→Cone→AC→BC→MG) 同 v6;
## 行序/left_anno/top_anno/色阶/行高全部同 v6
all_sr_v7 <- c("Early_central","Mid_central","Late_central","Young_central","Old_central",
               "Early_peripheral","Mid_peripheral","Late_peripheral","Young_peripheral","Old_peripheral")
keep_cols_v7 <- intersect(as.vector(t(outer(all_sr_v7, cts, paste, sep="__"))), keep_cols_ref)
z7 <- zscore_block(ACC_pk[ld_pk$peak, keep_cols_v7, drop=FALSE])
rownames(z7) <- make.unique(paste0(ld_pk$gene, " | ", ld_pk$peak))
colnames(z7) <- unname(ct_vec[keep_cols_v7])
ta7 <- HeatmapAnnotation(`stage`=factor(unname(st_vec[keep_cols_v7]), levels=c("Early","Mid","Late","Young","Old")),
                         `region`=factor(unname(rg_vec[keep_cols_v7]), levels=c("central","peripheral")),
                         `celltype`=factor(unname(ct_vec[keep_cols_v7]), levels=cts),
                         col=list(`stage`=cols.stage, `region`=cols.region, `celltype`=cols.celltype),
                         na_col="grey90", show_legend=c(TRUE,TRUE,TRUE),
                         annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(0.32,"cm"))
la7 <- rowAnnotation(region=factor(ld_pk$dir, levels=c("central","peripheral")),
                     group=factor(ld_pk$group, levels=group_order),
                     celltype=ld_pk$celltype,
                     col=list(region=cols.region, group=cols.group, celltype=cols.celltype),
                     annotation_name_gp=gpar(fontsize=8))
h7 <- Heatmap(z7, name="zblock", col=z_fun, na_col="grey90",
              column_title="Fig4 DEG-link | central(top) + peripheral(bottom) | peakDEG genes' same-period DAR-link CRE | snATAC accessibility (v7: columns grouped central|peripheral)",
              column_title_gp=gpar(fontsize=11), top_annotation=ta7, left_annotation=la7,
              row_split=factor(ld_pk$dir, levels=c("central","peripheral")),
              column_split=factor(unname(sr_vec[keep_cols_v7]), levels=all_sr_v7),
              cluster_columns=FALSE, cluster_rows=FALSE, row_order=seq_len(nrow(z7)),
              show_row_names=FALSE, show_column_names=FALSE,
              row_title_gp=gpar(fontsize=10), use_raster=TRUE, raster_quality=5,
              heatmap_legend_param=list(title="row z\n(block)"))
hh7 <- max(6, nrow(ld_pk)*0.05+3)
out7 <- paste0(plot_dir, "2026-10-20_Fig4_DEGlink_heatmap_CRE_access_linkDAR_peakDEG_merged_v7.pdf")
cairo_pdf(out7, width=13, height=hh7); grid.newpage(); draw(h7, newpage=FALSE); dev.off()
cat("done v7 linkDAR_peakDEG (列序 central→peripheral) | 行:", nrow(ld_pk),
    "| 高:", round(hh7,1), "in | pdf:", file.size(out7), "bytes\n")
