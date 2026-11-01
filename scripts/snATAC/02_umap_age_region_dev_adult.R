# 最终整理脚本 02 · snATAC UMAP (donor_age 着色 + central/peripheral region 着色) · 胚胎+成年 (2026-10-27)
# 来源: 00use_scripts/04_2026-05-29_胚胎全部样本summary.ipynb (umap_age + 末尾 umap_region 段)
#       00use_scripts/08_2026-05-29_成年summary.ipynb (同结构)
# 产出: Plot/2026-06-06_develop_umap_age.pdf        (8x8)
#       Plot/2026-06-06_develop_umap_region.pdf     (8x8)
#       Plot/2026-06-06_adult_umap_age.pdf          (10x8)
#       Plot/2026-06-06_adult_umap_region.pdf       (10x8)
# 注: 成年原 summary 读 ..._anno_sub_1.rds (该文件已更名为 ..._anno_sub.rds), 此处用现存文件名
suppressMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(RColorBrewer)
})

data_dir <- "/project/work/2026-05-08_ATAC_signac_all_sample/Data/"
plot_dir <- "/project/work/2026-05-08_ATAC_signac_all_sample/Plot/"

cols.palette1 <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Dark2"),brewer.pal(8, "Set2"),brewer.pal(8, "Accent"))
## age colors
cols.age <- c(cols.palette1[18:29],cols.palette1[6],cols.palette1[1:5],cols.palette1[8])
names(cols.age) <- c('70days','79days','87days','91days','100days','103days',
                     '116days','137days','141days','142days','162days','165days',
                     '21years','25years','67years','70years','75years','86years')

## region colors (fovea/macula 统一为 central)
cols.region <- c(
  central    = "#6A3D9A",
  peripheral = "#E69F00"
)

## ---- 成年 ----
date <- "2026-06-06_adult"
combined.atac <- readRDS(paste0(data_dir, "2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds"))

p_umap <- DimPlot(combined.atac, reduction = "umap.atac.harmony", group.by = "donor_age", label = FALSE , cols = cols.age, raster= TRUE)
pdf(paste0(plot_dir,date,"_umap_age",".pdf"), width = 10, height = 8)
  print(p_umap)
dev.off()

combined.atac$region <- as.character(combined.atac$tissue)
combined.atac$region[combined.atac$tissue %in% c('fovea','macula')] = 'central'

p_umap <- DimPlot(combined.atac, reduction = "umap.atac.harmony", group.by = "region", label = FALSE , cols = cols.region, raster=TRUE)
pdf(paste0(plot_dir,date,"_umap_region",".pdf"), width = 10, height = 8)
  print(p_umap)
dev.off()

rm(combined.atac); invisible(gc())

## ---- 胚胎 ----
date <- "2026-06-06_develop"
combined.atac <- readRDS(paste0(data_dir, "2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds"))

p_umap <- DimPlot(combined.atac, reduction = "umap.atac.harmony", group.by = "donor_age", label = FALSE , cols = cols.age, raster=TRUE)
pdf(paste0(plot_dir,date,"_umap_age",".pdf"), width = 8, height = 8)
  print(p_umap)
dev.off()

combined.atac$region <- as.character(combined.atac$tissue)
combined.atac$region[combined.atac$tissue %in% c('fovea','macula')] = 'central'

p_umap <- DimPlot(combined.atac, reduction = "umap.atac.harmony", group.by = "region", label = FALSE , cols = cols.region, raster=TRUE)
pdf(paste0(plot_dir,date,"_umap_region",".pdf"), width = 8, height = 8)
  print(p_umap)
dev.off()

rm(combined.atac); invisible(gc())
cat("done: 4 UMAP PDFs\n")
