# 最终整理脚本 17 · 胚胎 TF footprint 计算 (Signac 官方 Footprint + Tn5 校正) (2026-10-29)
# 来源: 00use_scripts/34_2026-10-16_Footprint_Development.ipynb (计算主链 模块 1-3;
#       绘图模块 4 与 MAFK/ZNF460/LHX5 追加段 (非收录图) 已去除)
# 产出: 2026-10-11_ATAC_Footprint/Data/2026-10-12_Development_footprint_<MOTIF>.rds (16 的输入)
# 参数: 同 14 号 (窗口 ±250 bp / peaks 内位点 / bias 2 万细胞 seed=42 / 全细胞不降采样)
# 运行: 单 motif 参数式; 收录链条 (16) 需要的 motif = ONECUT1
#       (MOTIF 已由来源脚本的探索状态 NFE2L1 改为 ONECUT1)
# 输入: 胚胎 ArchR snATAC 对象 (含 motifs/chromVAR/P2G links, 只读, 31G)
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(jsonlite)
  library(data.table)
})

# ============ 模块 1a：参数（集中设置；换 TF 只改 MOTIF） ============
DATE_PREFIX <- "2026-10-12"
PERIOD      <- "Development"          # 本 notebook 时期
MOTIF       <- "ONECUT1"                # motif 名（JASPAR 基因名或 MA 编号）

ROOT   <- "/project/work"
F_OBJ  <- file.path(ROOT, "2026-08-03_develop_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds")   # 原始对象（只读，不修改不覆盖）
GENOME <- BSgenome.Hsapiens.UCSC.hg38   # 与 AddMotifs / chromVAR 一致（hg38）
UPSTREAM   <- 250                        # footprint 窗口
DOWNSTREAM <- 250
IN_PEAKS   <- TRUE                       # 仅统计 ATAC peaks 内的 motif 位点

BIAS_SEED  <- 42                         # bias 抽样（仅用于估计 bias，不影响全细胞计算）
BIAS_NCELL <- 20000

# 输出路径（缓存均持久化于 Data/，文件名含时期与 motif）
DIR_PLOT <- file.path(ROOT, "2026-10-11_ATAC_Footprint", "Plot")
DIR_DATA <- file.path(ROOT, "2026-10-11_ATAC_Footprint", "Data")
SAFE_MOTIF <- gsub("[^A-Za-z0-9_.-]", "_", MOTIF)
F_BIAS_CACHE <- file.path(DIR_DATA, sprintf("%s_%s_Tn5bias_%dcells_seed%d.rds",
                                            DATE_PREFIX, PERIOD, BIAS_NCELL, BIAS_SEED))
F_FP_CACHE <- file.path(DIR_DATA, sprintf("%s_%s_footprint_%s.rds",
                                          DATE_PREFIX, PERIOD, SAFE_MOTIF))
F_PLOT_PDF <- file.path(DIR_PLOT, sprintf("%s_%s_footprint_%s.pdf",
                                          DATE_PREFIX, PERIOD, SAFE_MOTIF))

for (d in c(DIR_PLOT, DIR_DATA)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
stopifnot(file.exists(F_OBJ))
t_start <- Sys.time()
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
log_msg("时期:", PERIOD, "| motif:", MOTIF)
cat("结果缓存:", F_FP_CACHE, "| 存在:", file.exists(F_FP_CACHE), "\n")
# ============ 模块 1b：配色与版本记录 ============
# 读取 color.ipynb 配色（cols.celltype / lighten_col）
j <- fromJSON(file.path(ROOT, "color.ipynb"), simplifyVector = FALSE)
color_code <- vapply(j$cells, function(cc) {
  if (identical(cc$cell_type, "code")) paste(unlist(cc$source), collapse = "") else ""
}, character(1))
eval(parse(text = paste(color_code[nzchar(color_code)], collapse = "\n")), envir = globalenv())

# 分组配色：celltype 本色 = central，lighten 0.45 = peripheral（沿用 Fig4 约定）
cols.fp_group <- c()
for (ct in names(cols.celltype)) {
  cols.fp_group[paste0(ct, "_central")]    <- unname(cols.celltype[ct])
  cols.fp_group[paste0(ct, "_peripheral")] <- lighten_col(cols.celltype[ct], 0.45)
}

# 软件版本（环境升级后用于检查缓存 RDS 兼容性）
info_versions <- data.table(
  item    = c("R", "Signac", "Seurat", "SeuratObject", "BSgenome.Hsapiens.UCSC.hg38"),
  version = c(as.character(getRversion()), as.character(packageVersion("Signac")),
              as.character(packageVersion("Seurat")), as.character(packageVersion("SeuratObject")),
              as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg38")))
)
print(info_versions)
fwrite(info_versions, file.path(DIR_DATA, sprintf("%s_%s_versions.csv", DATE_PREFIX, PERIOD)))
# ============ 模块 1c：检查 motif 结果缓存（命中则不加载原始大对象） ============
# 校验项：motif / 窗口 / in.peaks / 分组与配色完整 / positionEnrichment 数据完整
# 注意：只 readRDS 一次（有效即保留使用，无效则丢弃重新计算）
mini <- NULL
if (file.exists(F_FP_CACHE)) {
  cached <- tryCatch(readRDS(F_FP_CACHE), error = function(e) NULL)
  p  <- if (is.null(cached)) NULL else cached@misc$footprint_params
  pe <- tryCatch(
    if (is.null(p)) NULL else GetAssayData(cached[["ATAC"]], slot = "positionEnrichment")[[p$motif_name]],
    error = function(e) NULL
  )
  valid <- !is.null(pe) &&
    identical(p$motif_name, MOTIF) &&
    identical(p$upstream, UPSTREAM) && identical(p$downstream, DOWNSTREAM) &&
    identical(p$in_peaks, IN_PEAKS) &&
    !anyNA(cached$group) &&
    all(as.character(unique(cached$group)) %in% names(cols.fp_group))
  if (isTRUE(valid)) {
    mini <- cached
    log_msg("结果缓存有效，直接绘图（不加载原始大对象、不重跑计算）")
    cat(sprintf("缓存记录: motif %s | 细胞 %d | peak 内位点 %d | 计算耗时 %.1f min\n",
                p$motif_id, p$n_cells, p$n_sites_inpeaks, p$fp_mins))
  } else {
    log_msg("结果缓存无效或参数不匹配，将重新计算")
    mini <- NULL
  }
  rm(cached)
}
if (is.null(mini)) cat("进入路径: 完整计算\n") else cat("进入路径: 缓存绘图\n")
# ============ 模块 1d：[仅缓存未命中] 读取对象并整理注释 ============
if (is.null(mini)) {
  log_msg("读取原始对象（文件大，需数分钟）...")
  t_read <- Sys.time()
  obj <- readRDS(F_OBJ)
  log_msg(sprintf("读取完成: %.1f min", as.numeric(difftime(Sys.time(), t_read, units = "mins"))))

  DefaultAssay(obj) <- "ATAC"

  # ---- 输入检查：celltype 注释（缺颜色定义时报错补色，不删细胞） ----
  ct_missing <- setdiff(unique(as.character(CONTACT_EMAIL$celltype)), names(cols.celltype))
  if (length(ct_missing) > 0 || anyNA(CONTACT_EMAIL$celltype)) {
    stop("以下 celltype 在 color.ipynb 的 cols.celltype 中缺少颜色定义，请补充配色后重跑（不自动删细胞）: ",
         paste(ct_missing, collapse = ", "))
  }
  obj$celltype <- factor(as.character(CONTACT_EMAIL$celltype), levels = names(cols.celltype))

  # ---- region：central / peripheral ----
  # Development 无 region 列，从 tissue 列按明确规则映射；未知标签报错（不默认 peripheral）
  region_map <- function(x) {
    x0  <- tolower(gsub("\\s+", "", as.character(x)))
    out <- rep(NA_character_, length(x0))
    out[grepl("fovea|macula|central", x0)] <- "central"
    out[grepl("peripheral|periphery", x0)] <- "peripheral"
    out
  }
  obj$region <- region_map(CONTACT_EMAIL$tissue)
  print(table(tissue = as.character(CONTACT_EMAIL$tissue), region = obj$region, useNA = "ifany"))
  stopifnot(!anyNA(obj$region))

  # 分组：celltype × region
  obj$group <- droplevels(factor(paste(obj$celltype, obj$region, sep = "_")))
  stopifnot(!anyNA(obj$group))
  tab_group <- table(obj$group)
  print(tab_group)
  fwrite(data.table(group = names(tab_group), n_cells = as.integer(tab_group)),
         file.path(DIR_DATA, sprintf("%s_%s_footprint_%s_group_cell_counts.csv",
                                     DATE_PREFIX, PERIOD, SAFE_MOTIF)))
  cat("细胞数:", ncol(obj), "| peaks:", nrow(obj[["ATAC"]]), "
")
}
# ============ 模块 1e：[仅缓存未命中] 运行前预检查（耗时计算前确认全部输入可用） ============
if (is.null(mini)) {
  # ---- motif 存在性与 positions ----
  motif_pos <- Signac::GetMotifData(Signac::Motifs(obj[["ATAC"]]), slot = "positions")
  motif_tbl <- Signac::GetMotifData(Signac::Motifs(obj[["ATAC"]]), slot = "motif.names")
  stopifnot(!is.null(motif_pos))
  motif_id <- if (MOTIF %in% names(motif_pos)) MOTIF else
    if (MOTIF %in% motif_tbl) names(which(motif_tbl == MOTIF)) else
    stop("motif ", MOTIF, " 不在对象 JASPAR2020 注释中（可用基因名或 MA 编号）")
  regions_all <- motif_pos[[motif_id]]

  # in.peaks = TRUE 后的位点数（不足则停止，不进入耗时计算）
  n_motif_genome  <- length(regions_all)
  regions_inpeaks <- subsetByOverlaps(regions_all, granges(obj[["ATAC"]]))
  n_motif_inpeaks <- length(regions_inpeaks)
  cat(sprintf("motif: %s（%s）| 全基因组位点: %d | peak 内位点: %d（保留 %.1f%%）\n",
              MOTIF, motif_id, n_motif_genome, n_motif_inpeaks,
              100 * n_motif_inpeaks / n_motif_genome))
  stopifnot(n_motif_inpeaks > 0)

  # ---- fragment 文件与索引存在、可访问 ----
  frags      <- Signac::Fragments(obj[["ATAC"]])
  frag_paths <- vapply(frags, Signac::GetFragmentData, character(1), slot = "path")
  stopifnot(length(frags) > 0)                       # fragment 数量（括号修复）
  frag_ok <- file.exists(frag_paths) & file.exists(paste0(frag_paths, ".tbi"))
  print(data.table(fragment = basename(frag_paths), file_ok = frag_ok))
  stopifnot(all(frag_ok))                            # 任一缺失立即停止
  cat("fragment 文件数:", length(frags), "| 全部存在且带索引\n")

  # ---- barcode 与 fragment 的 cells 映射完整 ----
  n_mapped <- sum(vapply(frags, function(f) {
    cm <- Signac::GetFragmentData(f, slot = "cells")
    if (is.null(cm)) 0L else sum(names(cm) %in% colnames(obj))
  }, integer(1)))
  cat("fragment 映射细胞数:", n_mapped, "/ 对象细胞数:", ncol(obj), "\n")
  stopifnot(n_mapped == ncol(obj))
}
# ============ 模块 2：InsertionBias 缓存（时期级，换 motif 复用；仅估计时用 2 万抽样细胞） ============
# 注意：Signac 1.12.0 未导出 GetAssayData / SetAssayData（属 SeuratObject），
# 直接调用即可，Signac 提供 ChromatinAssay 方法分派（slot = "bias" 已实测可用）
if (is.null(mini)) {
  bias_check <- function(b) {          # bias 有效性：长度 4096 / 无 NA / 无 Inf / hexamer 名称有效
    ok_names <- !is.null(names(b)) && length(names(b)) == 4096 &&
      all(nchar(names(b)) == 6) && all(grepl("^[ACGT]+$", names(b)))
    length(b) == 4096 && !anyNA(b) && all(is.finite(b)) && all(b > 0) && ok_names
  }

  bias_from_cache <- FALSE
  bias_ref <- NULL
  if (file.exists(F_BIAS_CACHE)) {
    bias_env <- tryCatch(readRDS(F_BIAS_CACHE), error = function(e) NULL)
    if (!is.null(bias_env) &&
        identical(bias_env$period, PERIOD) &&
        identical(bias_env$genome, as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg38"))) &&
        identical(bias_env$ncell, BIAS_NCELL) && identical(bias_env$seed, BIAS_SEED) &&
        bias_check(bias_env$bias)) {
      bias_ref <- bias_env$bias
      bias_from_cache <- TRUE
      log_msg("bias 缓存命中且有效，直接读取（不重复 InsertionBias）")
    } else {
      log_msg("bias 缓存不存在 / 不匹配 / 无效，将重新估计")
      bias_env <- NULL
    }
  }

  if (!bias_from_cache) {
    # celltype × region 分层比例抽样（固定种子；元数据按 barcode 从 meta.data 提取）
    # 用 min() 截断，保证抽样总数不超过 BIAS_NCELL（ceiling 向上取整可能超出）
    set.seed(BIAS_SEED)
    strata <- droplevels(interaction(obj$celltype, obj$region, drop = TRUE))
    frac   <- BIAS_NCELL / sum(table(strata))
    alloc  <- pmin(ceiling(table(strata) * frac), table(strata))
    if (sum(alloc) > BIAS_NCELL) {           # 仍超出时随机去掉多余名额（只减实际抽中的名额，不影响其他分层）
      drop_n  <- sum(alloc) - BIAS_NCELL
      drop_tb <- table(sample(rep(names(alloc), pmin(alloc, drop_n)), drop_n))
      drop_count <- setNames(integer(length(alloc)), names(alloc))
      drop_count[names(drop_tb)] <- as.integer(drop_tb)
      alloc <- alloc - drop_count
    }
    idx <- unlist(mapply(function(g, n) if (n > 0) sample(which(strata == g), n),
                         names(alloc), alloc), use.names = FALSE)
    cells_bias <- colnames(obj)[idx]
    meta_bias <- CONTACT_EMAIL[cells_bias, , drop = FALSE]
    cat("bias 抽样:", length(cells_bias), "细胞（上限", BIAS_NCELL, "）；celltype × region 构成：\n")
    print(table(celltype = as.character(meta_bias$celltype),
                region = as.character(meta_bias$region)))

    log_msg("估计 Tn5 insertion bias（chr1）...")
    t_bias <- Sys.time()
    bias_sub <- subset(obj, cells = cells_bias)
    bias_sub[["ATAC"]] <- Signac::InsertionBias(
      object = bias_sub[["ATAC"]], genome = GENOME,
      region = paste0("chr1-1-", seqlengths(GENOME)[["chr1"]])
    )
    bias_ref <- GetAssayData(bias_sub[["ATAC"]], slot = "bias")
    rm(bias_sub); invisible(gc())
    stopifnot(bias_check(bias_ref))      # 估计结果同样校验后再保存

    # 原子保存 RDS（先 .tmp 再改名，避免中断产生损坏缓存）
    bias_env <- list(
      period = PERIOD, genome = as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg38")),
      ncell = BIAS_NCELL, seed = BIAS_SEED, n_cells_sampled = length(cells_bias),
      sampling = table(celltype = as.character(meta_bias$celltype),
                       region = as.character(meta_bias$region)),
      bias = bias_ref, elapsed_mins = round(as.numeric(
        difftime(Sys.time(), t_bias, units = "mins")), 1)
    )
    f_tmp <- paste0(F_BIAS_CACHE, ".tmp")
    saveRDS(bias_env, f_tmp)
    stopifnot(file.rename(f_tmp, F_BIAS_CACHE))
    log_msg("bias 已估计并保存缓存:", F_BIAS_CACHE)
  }

  # 写入 assay：Footprint() 检测到 bias 后不会再自动估计（官方逻辑）
  obj[["ATAC"]] <- SetAssayData(obj[["ATAC"]], slot = "bias", new.data = bias_ref)
  stopifnot(!is.null(GetAssayData(obj[["ATAC"]], slot = "bias")))
  cat("bias 来源:", ifelse(bias_from_cache, "缓存", "本次估计"), "\n")
}
# ============ 模块 3a：官方 Signac::Footprint() 计算（全部有效细胞，不降采样） ============
if (is.null(mini)) {
  invisible(gc())
  log_msg("Footprint 计算开始:", MOTIF, "|", n_motif_inpeaks, "个 peak 内位点 |",
          ncol(obj), "个细胞")
  t_fp <- Sys.time()
  obj <- Signac::Footprint(
    object     = obj,
    motif.name = MOTIF,
    genome     = GENOME,
    in.peaks   = IN_PEAKS,
    upstream   = UPSTREAM,
    downstream = DOWNSTREAM
  )
  fp_mins <- as.numeric(difftime(Sys.time(), t_fp, units = "mins"))
  log_msg(sprintf("Footprint 计算完成: %.1f min", fp_mins))
}
# ============ 模块 3b：验证后保存轻量结果缓存（时期 × motif 独立 RDS；不保存 12/31 GB 完整对象） ============
if (is.null(mini)) {
  fp_mat <- GetAssayData(obj[["ATAC"]], slot = "positionEnrichment")[[MOTIF]]

  # 结构检查：行 = 细胞 barcode + expected + motif；列 = motif 周围相对碱基位置
  # 细胞用 setequal 比较（避免顺序差异误报）；行尾两行必须是 expected / motif
  stopifnot(
    !is.null(fp_mat),
    nrow(fp_mat) == ncol(obj) + 2,
    identical(tail(rownames(fp_mat), 2), c("expected", "motif")),
    setequal(head(rownames(fp_mat), -2), colnames(obj)),
    !anyNA(fp_mat["expected", ])
  )
  cat("positionEnrichment:", paste(dim(fp_mat), collapse = " × "), "（细胞 + expected + motif）\n")

  # 轻量对象：分组元数据 + positionEnrichment + 参数快照（PlotFootprint 仅需这些）
  mini_assay <- Signac::CreateChromatinAssay(
    counts = Matrix::Matrix(1, nrow = 2, ncol = ncol(obj), sparse = TRUE,
                            dimnames = list(c("chr1-1-1", "chr1-3-3"), colnames(obj))),
    ranges = GRanges("chr1", IRanges(start = c(1, 3), width = c(1, 1))),
    positionEnrichment = setNames(list(fp_mat), MOTIF)
  )
  mini <- CreateSeuratObject(counts = mini_assay, assay = "ATAC")
  mini$group    <- obj$group
  mini$celltype <- obj$celltype
  mini$region   <- obj$region
  mini@misc$footprint_params <- list(
    motif_name = MOTIF, motif_id = motif_id,
    upstream = UPSTREAM, downstream = DOWNSTREAM, in_peaks = IN_PEAKS,
    n_cells = ncol(obj), n_sites_genome = n_motif_genome, n_sites_inpeaks = n_motif_inpeaks,
    fp_mins = round(fp_mins, 1),
    versions = info_versions, group_counts = tab_group
  )

  # 先验证轻量对象可正常 PlotFootprint()，验证通过才保存正式 RDS
  invisible(tryCatch(Signac::PlotFootprint(mini, features = MOTIF, group.by = "group"),
                     error = function(e) stop("轻量对象绘图验证失败，不保存缓存: ",
                                              conditionMessage(e))))
  cat("轻量对象绘图验证通过\n")

  # 原子保存（先 .tmp 再改名；不同 motif 文件名不同，互不覆盖；不自动删除历史缓存）
  f_tmp <- paste0(F_FP_CACHE, ".tmp")
  saveRDS(mini, f_tmp)
  stopifnot(file.rename(f_tmp, F_FP_CACHE))
  log_msg("轻量结果缓存已保存:", F_FP_CACHE)

  rm(fp_mat, mini_assay, obj); invisible(gc())   # 释放大对象
}
