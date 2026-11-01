#!/usr/bin/env Rscript
# =====================================================================
# 2026-10-20_1c  GSE277326 应用注释 + 保存对象 + QC 图
# cluster(res0.8) -> celltype_fine -> celltype_major
# 输出:
#   Data/2026-10-20_1_GSE277326_annotated.rds          注释后完整对象
#   Data/2026-10-20_1_GSE277326_celltype_counts.csv    celltype×sample计数
#   Plot/2026-10-20_1_GSE277326_UMAP_*.pdf/.png
#   Plot/2026-10-20_1_GSE277326_DotPlot_canonical.pdf/.png
#   Plot/2026-10-20_1_GSE277326_celltype_proportion.pdf/.png
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
set.seed(1)

proj_dir <- "/project/work/2026-10-20_other_dataset"
data_dir <- file.path(proj_dir, "Data")
plot_dir <- file.path(proj_dir, "Plot")

## ---- cluster -> fine -> major 映射（依据 FindAllMarkers + canonical avgexp）----
fine_map <- c(
  "0" = "Rod", "1" = "Rod", "3" = "Rod", "4" = "Rod", "23" = "Rod_MT_hi",
  "6" = "Cone", "19" = "Cone_S_enriched",
  "5" = "Rod_bipolar", "7" = "OFF_cone_BC", "10" = "ON_cone_BC_GRIK1",
  "14" = "BC_KIT_pos", "17" = "ON_cone_BC_NEUROG1", "18" = "ON_cone_BC_KIRREL3",
  "24" = "BC_LHX3_pos", "25" = "BC_CALB1_IRX6",
  "12" = "HC", "28" = "HC_CALB1_hi",
  "8" = "AC_GAD2_hi", "9" = "AC_glycinergic", "11" = "AC_GABAergic",
  "13" = "AC_AII_like", "21" = "AC_TFAP2A_hi", "22" = "AC_VGluT3",
  "26" = "AC_dopaminergic", "29" = "AC_CARTPT", "30" = "AC_starburst_CHAT",
  "2" = "MG", "20" = "RGC", "27" = "Astrocyte", "32" = "Microglia",
  "31" = "RPE_contamination", "15" = "Doublet_BC_MG",
  "16" = "Other_unresolved_43yR", "33" = "Other_unresolved_21yM"
)
fine_to_major <- c(
  Rod = "Rod", Rod_MT_hi = "Rod", Cone = "Cone", Cone_S_enriched = "Cone",
  Rod_bipolar = "BC", OFF_cone_BC = "BC", ON_cone_BC_GRIK1 = "BC",
  BC_KIT_pos = "BC", ON_cone_BC_NEUROG1 = "BC", ON_cone_BC_KIRREL3 = "BC",
  BC_LHX3_pos = "BC", BC_CALB1_IRX6 = "BC",
  HC = "HC", HC_CALB1_hi = "HC",
  AC_GAD2_hi = "AC", AC_glycinergic = "AC", AC_GABAergic = "AC",
  AC_AII_like = "AC", AC_TFAP2A_hi = "AC", AC_VGluT3 = "AC",
  AC_dopaminergic = "AC", AC_CARTPT = "AC", AC_starburst_CHAT = "AC",
  MG = "MG", RGC = "RGC", Astrocyte = "Astrocyte", Microglia = "Microglia",
  RPE_contamination = "RPE_contamination", Doublet_BC_MG = "Doublet_BC_MG",
  Other_unresolved_43yR = "Other_unresolved", Other_unresolved_21yM = "Other_unresolved"
)
major_of_fine <- function(f) unname(fine_to_major[f])

## ---- 读入并应用 ----
obj <- readRDS(file.path(data_dir, "2026-10-20_1_tmp_clustered.rds"))
Idents(obj) <- "RNA_snn_res.0.8"
cl <- as.character(Idents(obj))
stopifnot(all(cl %in% names(fine_map)))
obj$celltype_fine <- unname(fine_map[cl])
obj$celltype_major <- major_of_fine(obj$celltype_fine)

major_order <- c("Rod", "Cone", "BC", "HC", "AC", "MG", "RGC", "Astrocyte",
                 "Microglia", "Doublet_BC_MG", "RPE_contamination", "Other_unresolved")
obj$celltype_major <- factor(obj$celltype_major, levels = major_order)

cat("==== celltype_major x sample ====\n")
tab <- table(obj$celltype_major, obj$sample_id)
print(tab)
write.csv(as.data.frame(tab), file.path(data_dir, "2026-10-20_1_GSE277326_celltype_counts.csv"),
          row.names = TRUE)

saveRDS(obj, file.path(data_dir, "2026-10-20_1_GSE277326_annotated.rds"))
cat("annotated object saved.\n")

## ---- 配色 ----
cols_major <- c(
  Rod = "#1F77B4", Cone = "#FF7F0E", BC = "#2CA02C", AC = "#D62728",
  HC = "#9467BD", MG = "#8C564B", RGC = "#E377C2", Astrocyte = "#17BECF",
  Microglia = "#7F7F7F", Doublet_BC_MG = "#FFD700", RPE_contamination = "#BCBD22",
  Other_unresolved = "#000000"
)

## ---- UMAP 图 ----
p_major <- DimPlot(obj, reduction = "umap", group.by = "celltype_major",
                   cols = cols_major, raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("GSE277326 annotated (res 0.8)") + theme_bw()
p_cluster <- DimPlot(obj, reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE,
                     raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("clusters res 0.8") + theme_bw()
p_sample <- DimPlot(obj, reduction = "umap", group.by = "sample_id",
                    raster = TRUE, raster.dpi = c(300, 300)) + theme_bw()
p_region <- DimPlot(obj, reduction = "umap", group.by = "region",
                    cols = c(Macula = "#6A3D9A", Retina = "#E69F00"),
                    raster = TRUE, raster.dpi = c(300, 300)) + theme_bw()

pdf(file.path(plot_dir, "2026-10-20_1_GSE277326_UMAP_overview.pdf"), width = 14, height = 12, useDingbats = FALSE)
print((p_major | p_cluster) / (p_sample | p_region))
dev.off()
ggsave(file.path(plot_dir, "2026-10-20_1_GSE277326_UMAP_overview.png"),
       (p_major | p_cluster) / (p_sample | p_region), width = 14, height = 12, dpi = 130)

## ---- DotPlot 经典 marker ----
canonical <- c(
  "RHO","NRL","GNAT1","ARR3","GNAT2","OPN1SW","RXRG",
  "VSX2","CABP5","GRM6","GRM8","PRKCA","LHX3","PCP4L1","SCGN",
  "LHX1","PVALB","CALB1",
  "TFAP2A","TFAP2B","GAD1","GAD2","SLC6A9","SLC32A1","SLC17A8","CHAT",
  "RBPMS","POU4F2","NEFM","SNCG",
  "RLBP1","GLUL","SOX9","SLC1A3","GFAP","AQP4",
  "PTPRC","AIF1","C1QA",
  "PDGFRB","RGS5","PECAM1","VWF","COL1A1","DCN","LUM",
  "TPSAB1","TRAC","RPE65","BEST1","DCT","MLANA","KIT"
)
canonical <- canonical[canonical %in% rownames(obj)]
obj_dot <- obj
Idents(obj_dot) <- droplevels(factor(obj$celltype_major))
p_dot <- DotPlot(obj_dot, features = canonical, group.by = "celltype_major") +
  coord_flip() + theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle("canonical markers by celltype_major")
pdf(file.path(plot_dir, "2026-10-20_1_GSE277326_DotPlot_canonical.pdf"), width = 10, height = 12, useDingbats = FALSE)
print(p_dot)
dev.off()
ggsave(file.path(plot_dir, "2026-10-20_1_GSE277326_DotPlot_canonical.png"), p_dot, width = 10, height = 12, dpi = 130)

## ---- 细胞类型比例 ----
prop_df <- as.data.frame(table(obj$sample_id, obj$celltype_major)) %>%
  group_by(Var1) %>% mutate(pct = Freq / sum(Freq) * 100) %>% ungroup()
colnames(prop_df)[1:2] <- c("sample_id", "celltype_major")
p_prop <- ggplot(prop_df, aes(x = sample_id, y = pct, fill = celltype_major)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cols_major) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "cell fraction (%)", fill = "celltype") +
  ggtitle("celltype composition per sample")
pdf(file.path(plot_dir, "2026-10-20_1_GSE277326_celltype_proportion.pdf"), width = 9, height = 6, useDingbats = FALSE)
print(p_prop)
dev.off()
ggsave(file.path(plot_dir, "2026-10-20_1_GSE277326_celltype_proportion.png"), p_prop, width = 9, height = 6, dpi = 130)

cat("all plots written.\n")
