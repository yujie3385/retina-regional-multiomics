# 最终整理脚本 19 · Fig4 DEGlink snRNA 表达热图 V7 (胚胎+成年合并, 列按 region 集中) (2026-10-27)
# 来源: 00use_scripts/19_2026-10-20_Fig4_DEGlink_heatmap_snRNA_expr_fullDEG_merged_v6.R (仅 V7 段)
# 产出: /project/work/2026-09-04_Fig4_heatmap/Plot/2026-10-20_Fig4_DEGlink_heatmap_snRNA_expr_fullDEG_merged_v7.pdf
# 行集: 三类 DEG (DM/DAL/AAM) 全部行 (dir,group,celltype,gene 粒度), central 169 行在上 + peripheral 196 行在下 (V5 固化行序)
# 列: 5 celltype × 10 期-区 (48 列, 2 个零细胞列已剔除); V7 列序 = central 24 列集中在前, peripheral 24 列在后
# 数值: log-normalized 表达假体积 → 胚胎 30 列块 / 成年 20 列块各自行 z-score, 色阶 [-2,2] blueYellow
# 整理说明: 去除 V6 三图 (fullDEG/peakDEG/MGonly) 绘制段, 仅保留最终收录的 fullDEG V7
suppressMessages({ library(ComplexHeatmap); library(circlize); library(grid); library(ArchR); library(RColorBrewer) })

proj_dir <- "/project/work/2026-09-04_Fig4_heatmap/"
data_dir <- paste0(proj_dir, "Data/")
plot_dir <- paste0(proj_dir, "Plot/")
prefix_v5 <- "2026-09-28_Fig4_DEGlink_heatmap"

## ---- 1. v5 缓存 (表达矩阵 + 行序 CSV, 最简输入) ----
cache <- readRDS(paste0(data_dir, prefix_v5, "_RNA_matrix_full_v5.rds"))
RNA <- cache$mat; keep_cols <- cache$keep_cols

## ---- 2. 列定义与配色 (hex 与 00_colors.R 逐项一致) ----
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
z_fun <- colorRamp2(seq(-2, 2, length.out=256), paletteContinuous("blueYellow"))
group_order <- c("DM","DAL","AAM")
emb_cols_v3 <- intersect(all_cols[1:30], keep_cols); adu_cols_v3 <- intersect(all_cols[31:50], keep_cols)
zscore_block <- function(m, ecols, acols){
  z <- m
  for(blk in list(ecols, acols)){
    z[,blk] <- t(apply(m[,blk,drop=FALSE], 1, function(x){ s<-sd(x); if(is.na(s)||s==0) return(rep(0,length(x))); (x-mean(x))/s }))
  }
  z
}

## ---- 3. V7: fullDEG 仅调整列序 — central 列集中在前, peripheral 在后; 列名不显示 ----
## 区内 stage 顺序 (Early→Mid→Late→Young→Old) 与 celltype 顺序 (Rod→Cone→AC→BC→MG) 同 v6;
## 行序/left_anno/top_anno/blueYellow 色阶/行高同 v6
all_sr_v7 <- c("Early_central","Mid_central","Late_central","Young_central","Old_central",
               "Early_peripheral","Mid_peripheral","Late_peripheral","Young_peripheral","Old_peripheral")
keep_cols_v7 <- intersect(as.vector(t(outer(all_sr_v7, cts, paste, sep="__"))), keep_cols)
oc7 <- read.csv(paste0(data_dir, prefix_v5, "_gene_row_order_full_v5_central.csv"))
op7 <- read.csv(paste0(data_dir, prefix_v5, "_gene_row_order_full_v5_peripheral.csv"))
rows7 <- rbind(oc7, op7)
z7 <- zscore_block(RNA[rows7$gene, keep_cols_v7, drop=FALSE], emb_cols_v3, adu_cols_v3)
rownames(z7) <- make.unique(rows7$gene)
colnames(z7) <- unname(ct_vec[keep_cols_v7])
ta7 <- HeatmapAnnotation(`stage`=factor(unname(st_vec[keep_cols_v7]), levels=c("Early","Mid","Late","Young","Old")),
                         `region`=factor(unname(rg_vec[keep_cols_v7]), levels=c("central","peripheral")),
                         `celltype`=factor(unname(ct_vec[keep_cols_v7]), levels=cts),
                         col=list(`stage`=cols.stage, `region`=cols.region, `celltype`=cols.celltype),
                         na_col="grey90", show_legend=c(TRUE,TRUE,TRUE),
                         annotation_name_gp=gpar(fontsize=8), simple_anno_size=unit(0.32,"cm"))
la7 <- rowAnnotation(region=factor(rows7$dir, levels=c("central","peripheral")),
                     group=factor(rows7$group, levels=group_order),
                     celltype=rows7$celltype,
                     col=list(region=cols.region, group=cols.group, celltype=cols.celltype),
                     annotation_name_gp=gpar(fontsize=8))
h7 <- Heatmap(z7, name="zblock", col=z_fun, na_col="grey90",
              column_title="Fig4 DEG-link | central(top) + peripheral(bottom) | full DEG | Gene snRNA expression (v7: columns grouped central|peripheral)",
              column_title_gp=gpar(fontsize=11), top_annotation=ta7, left_annotation=la7,
              row_split=factor(rows7$dir, levels=c("central","peripheral")),
              column_split=factor(unname(sr_vec[keep_cols_v7]), levels=all_sr_v7),
              cluster_columns=FALSE, cluster_rows=FALSE, row_order=seq_len(nrow(z7)),
              show_row_names=FALSE, show_column_names=FALSE,
              row_title_gp=gpar(fontsize=10), use_raster=TRUE, raster_quality=5,
              heatmap_legend_param=list(title="row z\n(block)"))
hh7 <- max(6, nrow(rows7)*0.05+3)
out7 <- paste0(plot_dir, "2026-10-20_Fig4_DEGlink_heatmap_snRNA_expr_fullDEG_merged_v7.pdf")
cairo_pdf(out7, width=13, height=hh7); grid.newpage(); draw(h7, newpage=FALSE); dev.off()
cat("done v7 fullDEG (列序 central→peripheral, 无列名) | 行:", nrow(rows7),
    "| 高:", round(hh7,1), "in | pdf:", file.size(out7), "bytes\n")
