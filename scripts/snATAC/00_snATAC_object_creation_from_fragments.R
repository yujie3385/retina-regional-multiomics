# 最终整理脚本 00 · snATAC 对象构建完整流程 (fragment 文件 → 注释后 RDS + meta RDS) · 胚胎+成年 (2026-10-30)
# 来源: 00use_scripts/01/02/03/04 (胚胎链: create / SVD 离线版 / intergrate_anno / summary)
#       00use_scripts/05/06/07/08 (成年链: 9donors 每样本构建 / create / intergrate_anno / summary)
# 产出: Data/2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds (+ _meta.rds)  胚胎 282,291 细胞
#       Data/2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds   (+ _meta.rds)  成年 126,192 细胞
#       (下游全部最终脚本的起点对象; 中间检查点 _1~_5 保留原历史文件名)
# 流程: fragments → MACS2 峰集 → 每样本对象 → merge + 样本元数据 → QC(TSS>2 & nCount>500)
#       → TF-IDF/LSI → Harmony/UMAP/聚类 → GeneActivity → snRNA 标签迁移(CCA, 置信度≥0.3) → 过滤后重整合
#       → celltype 注释 (胚胎=迁移直取; 成年=两轮 cluster 人工核对) → anno_sub.rds + meta.rds
# 整理说明: ① 去掉 notebook 中的查看/绘图/summary 显示/注释掉的旧命令;
#           ② 胚胎链的每样本构建代码在原 notebook 中已注释 (产物来自早期执行的同款代码), 此处按 05 号活跃代码原样恢复;
#           ③ 成年链的 GeneActivity 段在原 notebook 中为注释状态 (产物来自早期运行), 此处按胚胎同款参数补全;
#           ④ 标签迁移的 RNA 参考对象: 胚胎 develop_all / 成年 6-donor 版 (各降采样每类 500 细胞后找锚点)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(future)
  library(magrittr)
  library(Matrix)
})
library(harmony)
library(EnsDb.Hsapiens.v86)
set.seed(1)

data_dir <- "/project/work/2026-05-08_ATAC_signac_all_sample/Data/"
input_dir <- "/project/work/00data/ChenLab_human_atac/Fragments/"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
macs2 <- Sys.which("macs2")
plan("multicore", workers = 32)
options(future.globals.maxSize = 256 * 1024^3)

meta.sample.age.info <- readxl::read_excel("/project/resources/2026-04-30_human_multiome/GSE265801_ATAC_with_fragments_0430.xlsx")
annotations <- readRDS("/project/resources/2026-04-30_human_multiome/EnsDb.Hsapiens.v86.annotations.rds")
blacklist_regions <- readRDS("/project/resources/2026-04-30_human_multiome/blacklist_regions.rds")

## ---- 0. 统一峰集 (MACS2 per-sample pseudobulk call peaks → reduce 合并 → 20-10,000 bp; 跨样本共享缓存) ----
peak_rds <- file.path(data_dir, "combined_peaks_from_macs.rds")
if (file.exists(peak_rds)) {
  message("读取已保存的 peak：", peak_rds)
  combined.peaks <- readRDS(peak_rds)
} else {
  peak_list <- vector("list", length(frag_paths))
  names(peak_list) <- names(frag_paths)
  for (nm in names(frag_paths)) {
    message("CallPeaks: ", nm)
    tmp.dir <- file.path(data_dir, "macs2_per_sample", nm)
    dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
    peak_list[[nm]] <- CallPeaks(object = frag_paths[[nm]], macs2.path = macs2,
                                 name = nm, effective.genome.size = 2.913e9)
    saveRDS(peak_list[[nm]], paste0(tmp.dir, "/", nm, "_peaks_from_macs.rds"))
  }
  gr_concat <- do.call(c, as.list(unname(peak_list)))
  combined.peaks <- GenomicRanges::reduce(gr_concat)
  pw <- width(combined.peaks)
  combined.peaks <- combined.peaks[pw < 10000L & pw > 20L]
  saveRDS(combined.peaks, peak_rds)
  message("已保存 combined.peaks -> ", peak_rds, " （n=", length(combined.peaks), "）")
}

# ======================================================================
# 第一部分: 胚胎链 (25 个样本, 12 名供体; 来源 01/02/03/04)
# ======================================================================
date <- "2026-06-01_develop_all"
target_samples <- c('GSM8295567','GSM8295571','GSM8295552','GSM8295556','GSM8295558','GSM8295562','GSM8295563',
                    'GSM8295588','GSM8295570','GSM8295586','GSM8295565','GSM8295574','GSM8295576','GSM8295596',
                    'GSM8295592','GSM8295597','GSM8295580','GSM8295594','GSM8295550','GSM8295559','GSM8295577',
                    'GSM8295582','GSM8295584','GSM8295553','GSM8295590')

## ---- 1. fragments 定位 + 每样本对象 + merge + 样本元数据 ----
all_fragments <- list.files(input_dir, pattern = "fragments\\.tsv\\.gz$", full.names = TRUE)
target_fragments <- sapply(target_samples, function(id) all_fragments[grepl(id, all_fragments)])
frag_paths <- unname(target_fragments)
frag_paths <- frag_paths[!grepl("\\.tbi\\.gz$", frag_paths, ignore.case = TRUE)]
frag_paths <- sort(frag_paths)
names(frag_paths) <- unname(sapply(frag_paths, function(x) strsplit(basename(x), "_")[[1]][1]))

objs <- vector("list", length(frag_paths)); names(objs) <- names(frag_paths)
for (nm in names(frag_paths)) {
  message("CreateFragmentObject + FeatureMatrix: ", nm)
  frags <- CreateFragmentObject(path = frag_paths[[nm]])
  cnts <- FeatureMatrix(fragments = frags, features = combined.peaks, verbose = TRUE)
  assay <- CreateChromatinAssay(counts = cnts, fragments = frags, annotation = NULL)
  objs[[nm]] <- CreateSeuratObject(assay, assay = "ATAC", project = nm)
}
combined <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = names(objs))
combined$sample <- gsub(rownames(CONTACT_EMAIL), pattern = "_\\S+", replacement = "")
combined$dataset <- combined$sample
meta.info <- subset(meta.sample.age.info, Sample_ID %in% names(frag_paths))
Idents(combined) <- combined$dataset
combined$age      <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$age)
combined$Sex      <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$Sex)
combined$ethnicity <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$ethnicity)
combined$tissue   <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$tissue)
combined$assay    <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$assay)
combined$donor_id <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$`donor id`)
saveRDS(combined, paste0(data_dir, date, "_1_merged.rds"))

## ---- 2. QC: 标准染色体过滤 + 基因注释 + fragments 重挂 + QC 指标 + 过滤 (TSS>2 & nCount>500) ----
combined.atac <- readRDS(paste0(data_dir, date, "_1_merged.rds"))
DefaultAssay(combined.atac) <- "ATAC"
peaks.keep <- seqnames(granges(combined.atac)) %in% standardChromosomes(granges(combined.atac))
combined.atac <- combined.atac[as.vector(peaks.keep), ]
Annotation(combined.atac) <- annotations
## merge 后 fragments 为空, 按 sample 前缀重挂 (validate 校验)
frag_paths_use <- frag_paths[target_samples]
frags <- lapply(names(frag_paths_use), function(sid) {
  obj_cells <- colnames(combined.atac)[as.character(combined.atac$sample) == sid]
  if (length(obj_cells) == 0) stop("No cells found in combined.atac for sample: ", sid)
  raw_barcodes <- sub(paste0("^", sid, "_"), "", obj_cells)
  cell_map <- raw_barcodes; names(cell_map) <- obj_cells
  CreateFragmentObject(path = frag_paths_use[[sid]], cells = cell_map, validate.fragments = TRUE)
})
names(frags) <- names(frag_paths_use)
Fragments(combined.atac[["ATAC"]]) <- frags
combined.atac <- NucleosomeSignal(object = combined.atac)
combined.atac <- TSSEnrichment(object = combined.atac)
combined.atac$blacklist_ratio <- FractionCountsInRegion(
  object = combined.atac, assay = 'ATAC', regions = blacklist_regions)
combined.atac <- subset(x = combined.atac, TSS.enrichment > 2 & nCount_ATAC > 500)
saveRDS(combined.atac, paste0(data_dir, date, "_2_merged_filtered.rds"))

## ---- 3. 降维整合: TF-IDF → TopFeatures → LSI → Harmony → UMAP → 聚类 (res 0.4) ----
combined.atac <- readRDS(paste0(data_dir, date, "_2_merged_filtered.rds"))
combined.atac <- RunTFIDF(combined.atac)
combined.atac <- FindTopFeatures(combined.atac)
combined.atac <- RunSVD(combined.atac)
saveRDS(combined.atac, paste0(data_dir, date, "_3_merged_RunTFIDF.rds"))

combined.atac <- readRDS(paste0(data_dir, date, "_3_merged_RunTFIDF.rds"))
combined.atac <- RunHarmony(object = combined.atac, group.by.vars = "sample", reduction = "lsi",
                            assay.use = "ATAC", dims.use = 2:30, theta = 2, lambda = 1,
                            max.iter.harmony = 30, project.dim = FALSE, reduction.save = "harmony.atac")
combined.atac <- RunUMAP(combined.atac, reduction = "harmony.atac", dims = 1:29,
                         reduction.name = "umap.atac.harmony", reduction.key = "atacHarmonyUMAP_")
combined.atac <- FindNeighbors(combined.atac, reduction = "harmony.atac", dims = 1:29,
                               graph.name = c("atac_harmony_nn", "atac_harmony_snn"))
combined.atac <- FindClusters(combined.atac, graph.name = "atac_harmony_snn",
                              resolution = 0.4, cluster.name = "harmony_ATAC_clusters_0.4")
saveRDS(combined.atac, paste0(data_dir, date, "_4_merged_harmony_integrated.rds"))

## ---- 4. GeneActivity → ACTIVITY assay (LogNormalize, median scale factor) ----
combined.atac <- readRDS(paste0(data_dir, date, "_4_merged_harmony_integrated.rds"))
gene.activities <- GeneActivity(combined.atac)
combined.atac[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
combined.atac <- NormalizeData(object = combined.atac, assay = 'ACTIVITY',
                               normalization.method = 'LogNormalize',
                               scale.factor = median(combined.atac$nCount_ACTIVITY))
saveRDS(combined.atac, paste0(data_dir, date, "_5_merged_harmony_integrated_add_ACTIVITY.rds"))

## ---- 5. snRNA 标签迁移 (CCA 锚点; 参考每类降采样 500; weight = harmony.atac, dims 1:29) ----
combined.atac_rna <- readRDS("/project/resources/Data/2026-06-05_develop_all_human_harmony_anno.rds")
Idents(combined.atac_rna) <- "celltype"
combined.atac_rna_sub <- subset(combined.atac_rna, downsample = 500)
combined.atac_rna <- combined.atac_rna_sub

combined.atac <- readRDS(paste0(data_dir, date, "_5_merged_harmony_integrated_add_ACTIVITY.rds"))
DefaultAssay(combined.atac) <- "ACTIVITY"
transfer.anchors <- FindTransferAnchors(reference = combined.atac_rna, query = combined.atac, reduction = 'cca')
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = combined.atac_rna$celltype,
                                 weight.reduction = combined.atac[['harmony.atac']], dims = 1:29)
predicted.labels$predicted.id[predicted.labels$prediction.score.max < 0.3] <- paste0(predicted.labels$predicted.id, "_", "Not_sure")
combined.atac <- AddMetaData(object = combined.atac, metadata = predicted.labels[, c(1, 10)])

## ---- 6. 置信度过滤 (prediction.score.max ≥ 0.3) + 过滤后重整合 ----
cells_keep <- WhichCells(combined.atac, expression = prediction.score.max >= 0.3)
combined.atac <- subset(combined.atac, cells = cells_keep)
combined.atac <- RunHarmony(object = combined.atac, group.by.vars = "sample", reduction = "lsi",
                            assay.use = "ATAC", dims.use = 2:30, theta = 2, lambda = 1,
                            max.iter.harmony = 30, project.dim = FALSE, reduction.save = "harmony.atac")
combined.atac <- RunUMAP(combined.atac, reduction = "harmony.atac", dims = 1:29,
                         reduction.name = "umap.atac.harmony", reduction.key = "atacHarmonyUMAP_")
combined.atac <- FindNeighbors(combined.atac, reduction = "harmony.atac", dims = 1:29,
                               graph.name = c("atac_harmony_nn", "atac_harmony_snn"))
combined.atac <- FindClusters(combined.atac, graph.name = "atac_harmony_snn",
                              resolution = 0.4, cluster.name = "harmony_ATAC_clusters_0.4")

## ---- 7. 注释 (胚胎: celltype = 迁移标签直取) + 最终对象 ----
combined.atac$tissue <- as.character(combined.atac$tissue)
combined.atac$tissue[combined.atac$tissue == 'fovea centralis'] = 'fovea'
combined.atac$tissue[combined.atac$tissue == 'macula lutea'] = 'macula'
combined.atac$tissue[combined.atac$tissue == 'peripheral region of retina'] = 'peripheral'
combined.atac$donor_age <- as.character(combined.atac$age)
combined.atac$donor_age[combined.atac$donor_id == 'Donor_1'] = '70days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_2'] = '79days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_3'] = '87days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_4'] = '91days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_5'] = '100days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_6'] = '103days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_7'] = '116days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_8'] = '137days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_9'] = '141days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_10'] = '142days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_11'] = '162days'
combined.atac$donor_age[combined.atac$donor_id == 'Donor_12'] = '165days'
combined.atac$celltype <- as.character(combined.atac$predicted.id)
combined.atac$celltype <- factor(combined.atac$celltype, level = c('Rod','Cone','AC','BC','HC','MG','RGC','RPC'))
combined.atac$tissue <- factor(combined.atac$tissue, level = c('fovea','macula','peripheral'))
combined.atac$tissue_age <- paste0(combined.atac$tissue, "_", combined.atac$donor_age)
combined.atac$tissue_age <- factor(combined.atac$tissue_age, level = c(
  'macula_70days','peripheral_70days','macula_79days','peripheral_79days',
  'macula_87days','peripheral_87days','macula_91days','peripheral_91days',
  'macula_100days','peripheral_100days','macula_103days','peripheral_103days',
  'macula_116days','peripheral_116days','macula_137days','peripheral_137days',
  'macula_141days','peripheral_141days','macula_142days','peripheral_142days',
  'macula_162days','peripheral_162days','macula_165days','peripheral_165days'))
saveRDS(combined.atac, file = paste0(data_dir, "2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds"))
saveRDS(CONTACT_EMAIL,
        file = paste0(data_dir, "2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds"))
rm(combined.atac, objs, combined); invisible(gc())

# ======================================================================
# 第二部分: 成年链 (18 个文库, 6 名供体; 来源 05/06/07/08)
# ======================================================================
date <- "2026-05-26_adult"
target_samples <- c('GSM8622737','GSM8622738','GSM8622739','GSM8622770','GSM8622772','GSM8622774','GSM8622759','GSM8622760','GSM8622730',
                    'GSM8622732','GSM8622734','GSM8622816','GSM8622817','GSM8622818','GSM8622819','GSM8622845','GSM8622846','GSM8622847')

## ---- 1. 每样本对象 + merge + 样本元数据 (同胚胎模式) ----
all_fragments <- list.files(input_dir, pattern = "fragments\\.tsv\\.gz$", full.names = TRUE)
target_fragments <- sapply(target_samples, function(id) all_fragments[grepl(id, all_fragments)])
frag_paths <- unname(target_fragments)
frag_paths <- frag_paths[!grepl("\\.tbi\\.gz$", frag_paths, ignore.case = TRUE)]
frag_paths <- sort(frag_paths)
names(frag_paths) <- unname(sapply(frag_paths, function(x) strsplit(basename(x), "_")[[1]][1]))

objs <- vector("list", length(frag_paths)); names(objs) <- names(frag_paths)
for (nm in names(frag_paths)) {
  message("CreateFragmentObject + FeatureMatrix: ", nm)
  frags <- CreateFragmentObject(path = frag_paths[[nm]])
  cnts <- FeatureMatrix(fragments = frags, features = combined.peaks, verbose = TRUE)
  assay <- CreateChromatinAssay(counts = cnts, fragments = frags, annotation = NULL)
  objs[[nm]] <- CreateSeuratObject(assay, assay = "ATAC", project = nm)
}
combined <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = names(objs))
combined$sample <- gsub(rownames(CONTACT_EMAIL), pattern = "_\\S+", replacement = "")
combined$dataset <- combined$sample
meta.info <- subset(meta.sample.age.info, Sample_ID %in% names(frag_paths))
Idents(combined) <- combined$dataset
combined$age      <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$age)
combined$Sex      <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$Sex)
combined$ethnicity <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$ethnicity)
combined$tissue   <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$tissue)
combined$assay    <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$assay)
combined$donor_id <- plyr::mapvalues(x = Idents(combined), from = meta.info$Sample_ID, to = meta.info$`donor id`)
saveRDS(combined, paste0(data_dir, date, "_1_merged.rds"))

## ---- 2. QC (同胚胎: 标准染色体 + 注释 + fragments 重挂 + QC 指标 + TSS>2 & nCount>500) ----
combined.atac <- readRDS(paste0(data_dir, date, "_1_merged.rds"))
DefaultAssay(combined.atac) <- "ATAC"
peaks.keep <- seqnames(granges(combined.atac)) %in% standardChromosomes(granges(combined.atac))
combined.atac <- combined.atac[as.vector(peaks.keep), ]
Annotation(combined.atac) <- annotations
frag_paths_use <- frag_paths[target_samples]
frags <- lapply(names(frag_paths_use), function(sid) {
  obj_cells <- colnames(combined.atac)[as.character(combined.atac$sample) == sid]
  if (length(obj_cells) == 0) stop("No cells found in combined.atac for sample: ", sid)
  raw_barcodes <- sub(paste0("^", sid, "_"), "", obj_cells)
  cell_map <- raw_barcodes; names(cell_map) <- obj_cells
  CreateFragmentObject(path = frag_paths_use[[sid]], cells = cell_map, validate.fragments = TRUE)
})
names(frags) <- names(frag_paths_use)
Fragments(combined.atac[["ATAC"]]) <- frags
combined.atac <- NucleosomeSignal(object = combined.atac)
combined.atac <- TSSEnrichment(object = combined.atac)
combined.atac$blacklist_ratio <- FractionCountsInRegion(
  object = combined.atac, assay = 'ATAC', regions = blacklist_regions)
combined.atac <- subset(x = combined.atac, TSS.enrichment > 2 & nCount_ATAC > 500)
saveRDS(combined.atac, paste0(data_dir, date, "_2_merged_filtered.rds"))

## ---- 3. 降维整合 (同胚胎参数) ----
combined.atac <- readRDS(paste0(data_dir, date, "_2_merged_filtered.rds"))
combined.atac <- RunTFIDF(combined.atac)
combined.atac <- FindTopFeatures(combined.atac)
combined.atac <- RunSVD(combined.atac)
saveRDS(combined.atac, paste0(data_dir, date, "_3_merged_RunTFIDF.rds"))

combined.atac <- readRDS(paste0(data_dir, date, "_3_merged_RunTFIDF.rds"))
combined.atac <- RunHarmony(object = combined.atac, group.by.vars = "sample", reduction = "lsi",
                            assay.use = "ATAC", dims.use = 2:30, theta = 2, lambda = 1,
                            max.iter.harmony = 30, project.dim = FALSE, reduction.save = "harmony.atac")
combined.atac <- RunUMAP(combined.atac, reduction = "harmony.atac", dims = 1:29,
                         reduction.name = "umap.atac.harmony", reduction.key = "atacHarmonyUMAP_")
combined.atac <- FindNeighbors(combined.atac, reduction = "harmony.atac", dims = 1:29,
                               graph.name = c("atac_harmony_nn", "atac_harmony_snn"))
combined.atac <- FindClusters(combined.atac, graph.name = "atac_harmony_snn",
                              resolution = 0.4, cluster.name = "harmony_ATAC_clusters_0.4")
saveRDS(combined.atac, paste0(data_dir, date, "_4_merged_harmony_integrated.rds"))

## ---- 4. GeneActivity → ACTIVITY assay (原 notebook 此段为注释状态, 产物来自早期运行; 参数同胚胎) ----
combined.atac <- readRDS(paste0(data_dir, date, "_4_merged_harmony_integrated.rds"))
gene.activities <- GeneActivity(combined.atac)
combined.atac[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
combined.atac <- NormalizeData(object = combined.atac, assay = 'ACTIVITY',
                               normalization.method = 'LogNormalize',
                               scale.factor = median(combined.atac$nCount_ACTIVITY))
saveRDS(combined.atac, paste0(data_dir, date, "_5_merged_harmony_integrated_add_ACTIVITY.rds"))

## ---- 5. snRNA 标签迁移 (参考: 成年 6-donor 版 RNA 对象, 每类降采样 500) ----
combined.atac_rna <- readRDS("/project/resources/Data/2026-06-05_adult_human_harmony_anno.rds")
Idents(combined.atac_rna) <- "celltype"
combined.atac_rna_sub <- subset(combined.atac_rna, downsample = 500)
combined.atac_rna <- combined.atac_rna_sub

combined.atac <- readRDS(paste0(data_dir, date, "_5_merged_harmony_integrated_add_ACTIVITY.rds"))
DefaultAssay(combined.atac) <- "ACTIVITY"
transfer.anchors <- FindTransferAnchors(reference = combined.atac_rna, query = combined.atac, reduction = 'cca')
predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = combined.atac_rna$celltype,
                                 weight.reduction = combined.atac[['harmony.atac']], dims = 1:29)
predicted.labels$predicted.id[predicted.labels$prediction.score.max < 0.3] <- paste0(predicted.labels$predicted.id, "_", "Not_sure")
combined.atac <- AddMetaData(object = combined.atac, metadata = predicted.labels[, c(1, 10)])

## ---- 6. 置信度过滤 (≥0.3) + 第一轮 cluster 注释 + 去未注释 cluster ----
cells_keep <- WhichCells(combined.atac, expression = prediction.score.max >= 0.3)
combined.atac <- subset(combined.atac, cells = cells_keep)
combined.atac$celltype <- paste0('UN_', combined.atac$harmony_ATAC_clusters_0.4)
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('2','19')] = 'MG'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('5')] = 'HC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('4','24')] = 'Cone'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('1','0')] = 'Rod'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('18','23','15','16','13','14','17','21')] = 'AC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('9','7','8','6','3','10','22')] = 'BC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('12')] = 'RGC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('20')] = 'Micro'
combined.atac <- subset(combined.atac, celltype != "UN_11")

## ---- 7. 过滤后重整合 + 第二轮 cluster 注释 ----
combined.atac <- RunHarmony(object = combined.atac, group.by.vars = "sample", reduction = "lsi",
                            assay.use = "ATAC", dims.use = 2:30, theta = 2, lambda = 1,
                            max.iter.harmony = 30, project.dim = FALSE, reduction.save = "harmony.atac")
combined.atac <- RunUMAP(combined.atac, reduction = "harmony.atac", dims = 1:29,
                         reduction.name = "umap.atac.harmony", reduction.key = "atacHarmonyUMAP_")
combined.atac <- FindNeighbors(combined.atac, reduction = "harmony.atac", dims = 1:29,
                               graph.name = c("atac_harmony_nn", "atac_harmony_snn"))
combined.atac <- FindClusters(combined.atac, graph.name = "atac_harmony_snn",
                              resolution = 0.4, cluster.name = "harmony_ATAC_clusters_0.4")
combined.atac$celltype <- paste0('UN_', combined.atac$harmony_ATAC_clusters_0.4)
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('2','18')] = 'MG'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('5')] = 'HC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('4','23')] = 'Cone'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('1','0')] = 'Rod'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('21','14','17','15','20','12','16','13')] = 'AC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('3','8','10','22','7','6','9')] = 'BC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('11')] = 'RGC'
combined.atac$celltype[combined.atac$harmony_ATAC_clusters_0.4 %in% c('19')] = 'Micro'

## ---- 8. 最终注释与对象保存 ----
combined.atac$tissue <- as.character(combined.atac$tissue)
combined.atac$tissue[combined.atac$tissue == 'fovea centralis'] = 'fovea'
combined.atac$tissue[combined.atac$tissue == 'macula lutea'] = 'macula'
combined.atac$tissue[combined.atac$tissue == 'peripheral region of retina'] = 'peripheral'
combined.atac$celltype <- factor(combined.atac$celltype, level = c('Rod','Cone','AC','BC','HC','MG','RGC','Micro'))
combined.atac$tissue <- factor(combined.atac$tissue, level = c('fovea','macula','peripheral'))
combined.atac$tissue_age <- paste0(combined.atac$tissue, "_", combined.atac$age)
combined.atac$tissue_age <- factor(combined.atac$tissue_age, level = c(
  'fovea_21years','macula_21years','peripheral_21years',
  'fovea_25years','macula_25years','peripheral_25years',
  'fovea_67years','peripheral_67years','macula_70years',
  'macula_75years','peripheral_75years','macula_86years','peripheral_86years'))
saveRDS(combined.atac, file = paste0(data_dir, "2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds"))
saveRDS(CONTACT_EMAIL,
        file = paste0(data_dir, "2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds"))
rm(combined.atac, objs, combined); invisible(gc())
cat("done: 2 anno_sub RDS + 2 meta RDS\n")
