# 最终整理脚本 18 · TF footprint 轨道图 V5 (Tn5 校正后 central vs peripheral 曲线, 分页合并 PDF) · 成年+胚胎 (2026-10-27)
# 来源: 00use_scripts/22_2026-10-14_Adult_footprint_plot_test.ipynb (仅 V5 cell)
#       00use_scripts/23_2026-10-16_Development_footprint_plot_test.ipynb (仅 V5 cell)
# 产出: 2026-10-11_ATAC_Footprint/Plot/2026-10-14_Adult_footprint_{KLF9,MAFK,NFE2L1,ONECUT1}_V5_age.pdf (各 2 页: 21y/67y)
#       2026-10-11_ATAC_Footprint/Plot/2026-10-16_Development_footprint_ONECUT1_V5_stage.pdf (3 页: Early/Mid/Late)
# 整理说明: 去除 V1-V4 测试版本; motif_target 收敛到最终收录的 motif (成年 4 个, 胚胎 1 个);
#           输入为轻量缓存 (Tn5 校正已含), 无重计算
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(GenomicRanges)
  library(jsonlite)
  library(data.table)
})

ROOT     <- "/project/work"
DIR_PLOT <- file.path(ROOT, "2026-10-11_ATAC_Footprint", "Plot")
DIR_DATA <- file.path(ROOT, "2026-10-11_ATAC_Footprint", "Data")

## 配色 (color.ipynb: cols.region = central 紫 / peripheral 橙; cols.celltype)
j <- fromJSON(file.path(ROOT, "color.ipynb"), simplifyVector = FALSE)
color_code <- vapply(j$cells, function(cc) {
  if (identical(cc$cell_type, "code")) paste(unlist(cc$source), collapse = "") else ""
}, character(1))
eval(parse(text = paste(color_code[nzchar(color_code)], collapse = "\n")), envir = globalenv())
stopifnot(all(c("central", "peripheral") %in% names(cols.region)))

## motif 名转文件安全名 (与计算脚本 SAFE_MOTIF 一致, 如 THRB(var.3) -> THRB_var.3_)
safe_name <- function(motif) gsub("[^A-Za-z0-9_.-]", "_", motif)
## 读取指定 motif 的轻量结果缓存 (只读, 无重计算; PERIOD 由各段设定)
load_mini <- function(motif) {
  f <- file.path(DIR_DATA, sprintf("2026-10-12_%s_footprint_%s.rds", PERIOD, safe_name(motif)))
  stopifnot(file.exists(f))   # 缓存不存在时需先用计算脚本跑该 motif
  readRDS(f)
}

## 统一 y 轴范围: 按 PlotFootprint 内部算法精确复现 (组内背景归一化均值 - expected)
sub_y_range <- function(obj_region, motif, group_by = "celltype_n") {
  pd  <- Signac::GetFootprintData(obj_region, features = motif, group.by = group_by)
  obs <- pd[pd$class == "Observed", ]
  ev  <- setNames(pd$norm.value[pd$class == "Expected"], pd$position[pd$class == "Expected"])
  range(obs$norm.value - ev[as.character(obs$position)], na.rm = TRUE)
}
## 仅替换主图 (footprint 面板) 的 y 轴，不影响下方 expected 面板
set_unified_ylim <- function(p, yr) {
  suppressWarnings({
    p$patches$plots[[1]]$patches$plots[[1]] <-
      p$patches$plots[[1]]$patches$plots[[1]] + ggplot2::ylim(yr)
  })
  p
}

# ============ 成年 V5: 目标细胞类型 × 年龄分页 (21y/67y 各一页), central vs peripheral ============
# 两页统一 y 轴范围 (跨年龄可比); 仅取 21y/67y, 其余年龄 (25/70/75/86 岁) 不纳入
DATE_PREFIX <- "2026-10-14"
PERIOD      <- "Adult"
## 最终收录 motif (均 MG; 完整 16 motif 登记见来源脚本 22)
motif_target <- c(ONECUT1 = "MG", KLF9 = "MG", MAFK = "MG", NFE2L1 = "MG")

## barcode -> 年龄 (2026-05-08 构建时的正式元数据, 与结果缓存 barcode 完全一致)
meta_age <- readRDS(file.path(ROOT, "2026-05-08_ATAC_signac_all_sample", "Data",
                               "2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds"))
age_vec  <- setNames(as.character(meta_age$age), rownames(meta_age))
ages_use <- c("21 years", "67 years")
stopifnot(all(ages_use %in% unique(age_vec)))

for (motif in names(motif_target)) {
  target_ct <- unname(motif_target[motif])
  mini <- load_mini(motif)
  motif_name <- mini@misc$footprint_params$motif_name
  m <- GetAssayData(mini[["ATAC"]], slot = "positionEnrichment")[[motif_name]]

  ## 目标细胞类型的细胞，按年龄取子集
  cells_all <- colnames(mini)[as.character(mini$celltype) == target_ct]
  age_of    <- age_vec[cells_all]

  objs_v5 <- list(); plots_v5 <- list()
  for (a in ages_use) {
    cells_keep <- cells_all[!is.na(age_of) & age_of == a]
    assay <- Signac::CreateChromatinAssay(
      counts = Matrix::Matrix(1, nrow = 2, ncol = length(cells_keep), sparse = TRUE,
                              dimnames = list(c("chr1-1-1", "chr1-3-3"), cells_keep)),
      ranges = GRanges("chr1", IRanges(start = c(1, 3), width = c(1, 1))),
      positionEnrichment = setNames(list(m[c(cells_keep, "expected", "motif"), ]), motif_name)
    )
    mini_a <- CreateSeuratObject(counts = assay, assay = "ATAC")
    rg  <- factor(as.character(mini$region[cells_keep]), levels = c("central", "peripheral"))
    tab <- table(droplevels(rg))
    mini_a$region_n <- factor(as.character(rg), levels = names(tab),
                              labels = sprintf("%s (n=%d)", names(tab), as.integer(tab)))
    stopifnot(nrow(GetAssayData(mini_a[["ATAC"]], slot = "positionEnrichment")[[motif_name]]) == ncol(mini_a) + 2)
    cat(sprintf("\n===== motif %s | %s | %s: central (n=%d) vs peripheral (n=%d) =====\n",
                motif, a, target_ct, tab[["central"]], tab[["peripheral"]]))

    p_a <- Signac::PlotFootprint(mini_a, features = motif, group.by = "region_n",
                                 normalization = "subtract", show.expected = TRUE)
    labs_a <- levels(droplevels(mini_a$region_n))
    cols_a <- setNames(cols.region[sub(" \\(n=.*$", "", labs_a)], labs_a)
    p_a <- p_a & scale_color_manual(values = cols_a)
    p_a <- p_a + plot_annotation(
      title = sprintf("%s motif footprint - Adult %s, %s central vs peripheral (Tn5-corrected)",
                      motif, target_ct, a),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
    objs_v5[[a]]  <- mini_a
    plots_v5[[a]] <- p_a
  }

  ## 两页统一 y 轴范围 (region 分组，跨两页取公共范围)
  y_v5 <- range(
    sub_y_range(objs_v5[["21 years"]], motif, group_by = "region_n"),
    sub_y_range(objs_v5[["67 years"]], motif, group_by = "region_n"))
  cat(sprintf("统一 y 轴范围: [%.3f, %.3f]\n", y_v5[1], y_v5[2]))
  for (a in ages_use) plots_v5[[a]] <- set_unified_ylim(plots_v5[[a]], y_v5)

  ## 两页合并输出一个 PDF (每页保持完整版式)
  f_v5 <- file.path(DIR_PLOT, sprintf("%s_%s_footprint_%s_V5_age.pdf", DATE_PREFIX, PERIOD, safe_name(motif)))
  pdf(f_v5, width = 8, height = 5.5)
  print(plots_v5[["21 years"]])
  print(plots_v5[["67 years"]])
  dev.off()
  cat("PDF 已导出 (2 页):\n ", f_v5, "\n")
  rm(mini, objs_v5, plots_v5); invisible(gc())
}

# ============ 胚胎 V5: 目标细胞类型 × stage 分页 (Early/Mid/Late 各一页), central vs peripheral ============
# 三页统一 y 轴范围 (跨 stage 可比); stage 由 donor_age 按已验证映射推导
DATE_PREFIX <- "2026-10-16"
PERIOD      <- "Development"
## 最终收录 motif (完整 10 motif 登记见来源脚本 23)
motif_target <- c(ONECUT1 = "MG")

## barcode -> stage (2026-05-08 构建时的正式元数据, 与结果缓存 barcode 完全一致)
meta_dev  <- readRDS(file.path(ROOT, "2026-05-08_ATAC_signac_all_sample", "Data",
                                "2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub_meta.rds"))
d_age     <- as.integer(sub("days.*", "", as.character(meta_dev$donor_age)))
stage_vec <- setNames(ifelse(d_age <= 91, "Early", ifelse(d_age <= 116, "Mid", "Late")), rownames(meta_dev))
stages_use <- c("Early", "Mid", "Late")
stopifnot(all(stages_use %in% unique(stage_vec)))
print(table(stage_vec))

for (motif in names(motif_target)) {
  target_ct <- unname(motif_target[motif])
  mini <- load_mini(motif)
  motif_name <- mini@misc$footprint_params$motif_name
  m <- GetAssayData(mini[["ATAC"]], slot = "positionEnrichment")[[motif_name]]

  ## 目标细胞类型的细胞，按 stage 取子集
  cells_all <- colnames(mini)[as.character(mini$celltype) == target_ct]
  st_of     <- stage_vec[cells_all]

  objs_v5 <- list(); plots_v5 <- list()
  for (s in stages_use) {
    cells_keep <- cells_all[!is.na(st_of) & st_of == s]
    assay <- Signac::CreateChromatinAssay(
      counts = Matrix::Matrix(1, nrow = 2, ncol = length(cells_keep), sparse = TRUE,
                              dimnames = list(c("chr1-1-1", "chr1-3-3"), cells_keep)),
      ranges = GRanges("chr1", IRanges(start = c(1, 3), width = c(1, 1))),
      positionEnrichment = setNames(list(m[c(cells_keep, "expected", "motif"), ]), motif_name)
    )
    mini_s <- CreateSeuratObject(counts = assay, assay = "ATAC")
    rg  <- factor(as.character(mini$region[cells_keep]), levels = c("central", "peripheral"))
    tab <- table(droplevels(rg))
    mini_s$region_n <- factor(as.character(rg), levels = names(tab),
                              labels = sprintf("%s (n=%d)", names(tab), as.integer(tab)))
    stopifnot(nrow(GetAssayData(mini_s[["ATAC"]], slot = "positionEnrichment")[[motif_name]]) == ncol(mini_s) + 2)
    cat(sprintf("\n===== motif %s | %s | %s: central (n=%d) vs peripheral (n=%d) =====\n",
                motif, s, target_ct, tab[["central"]], tab[["peripheral"]]))

    p_s <- Signac::PlotFootprint(mini_s, features = motif, group.by = "region_n",
                                 normalization = "subtract", show.expected = TRUE)
    labs_s <- levels(droplevels(mini_s$region_n))
    cols_s <- setNames(cols.region[sub(" \\(n=.*$", "", labs_s)], labs_s)
    p_s <- p_s & scale_color_manual(values = cols_s)
    p_s <- p_s + plot_annotation(
      title = sprintf("%s motif footprint - Development %s, %s central vs peripheral (Tn5-corrected)",
                      motif, target_ct, s),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
    objs_v5[[s]]  <- mini_s
    plots_v5[[s]] <- p_s
  }

  ## 三页统一 y 轴范围 (region 分组，跨三页取公共范围)
  y_v5 <- range(
    sub_y_range(objs_v5[["Early"]], motif, group_by = "region_n"),
    sub_y_range(objs_v5[["Mid"]],   motif, group_by = "region_n"),
    sub_y_range(objs_v5[["Late"]],  motif, group_by = "region_n"))
  cat(sprintf("统一 y 轴范围: [%.3f, %.3f]\n", y_v5[1], y_v5[2]))
  for (s in stages_use) plots_v5[[s]] <- set_unified_ylim(plots_v5[[s]], y_v5)

  ## 三页合并输出一个 PDF (每页保持完整版式)
  f_v5 <- file.path(DIR_PLOT, sprintf("%s_%s_footprint_%s_V5_stage.pdf", DATE_PREFIX, PERIOD, safe_name(motif)))
  pdf(f_v5, width = 8, height = 5.5)
  for (s in stages_use) print(plots_v5[[s]])
  dev.off()
  cat("PDF 已导出 (3 页):\n ", f_v5, "\n")
  rm(mini, objs_v5, plots_v5); invisible(gc())
}

cat("done: 5 footprint V5 PDFs\n")
