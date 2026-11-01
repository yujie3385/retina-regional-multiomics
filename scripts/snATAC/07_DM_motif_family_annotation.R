# 最终整理脚本 07 · DM motif 的 JASPAR family 注释 (candidate_CRE 项目) (2026-10-29)
# 来源: 00use_scripts/36_codex_2026-08-15_3_plot_DM_motif_enrichment_chromVAR_integration.R
#       (仅保留 0-6 节至 family annotation 写出 — 该输出是 08 号 family map 的输入;
#        family_collapsed 折叠与 7-10 节 (chromVAR regional/adult central/integrated evidence/输出QC)
#        属 DM motif enrichment 自身分析, 未进入最终链条, 已去除)
# 产出: candidate_CRE/Data/2026-08-15_3_DM_motif_family_annotation.csv
#       (motif_id / motif_name / TF_name / JASPAR_family / JASPAR_class / family_for_collapse)
# 链路: pairing 主表 + DM foreground → motif object → DM-specific FindMotifs → JASPAR2020 family 映射
# ---- 0. 环境设置 -----------------------------------------------------------
Sys.setenv(SP_EVOLUTION_STATUS = "2")
set.seed(1)

required_pkgs <- c(
  "Signac", "SeuratObject", "Matrix", "dplyr", "tidyr", "readr",
  "tibble", "stringr", "TFBSTools", "JASPAR2020"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("缺少必要 R 包：", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(Signac)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(TFBSTools)
  library(JASPAR2020)
})

# ---- 1. 参数设置 ------------------------------------------------------------
date_str <- "2026-08-15_3"
data_dir <- "/project/work/2026-08-08_candidate_CRE/Data"
output_dir <- file.path(data_dir, "2026-08-15_3_DM_motif_chromVAR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pair_file <- file.path(data_dir, "2026-08-11_DEG_DAR_pairing",
                       "2026-08-11_DEG_DAR_pairing_peak_level.csv")
gene_summary_file <- file.path(data_dir, "2026-08-11_DEG_DAR_pairing",
                               "2026-08-11_DEG_DAR_pairing_gene_summary.csv")
develop_obj_file <- "/project/work/2026-08-03_develop_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds"
adult_obj_file <- "/project/work/2026-08-03_adult_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds"
develop_dar_dir <- "/project/work/2026-07-06_development_DAR/Data"
adult_dar_file <- "/project/work/2026-08-03_adult_DAR/Data/2026-08-03_adult_shared_fovea_peripheral_DAR.csv"
develop_chromvar_file <- file.path(data_dir, "2026-08-12_4_develop_chromVAR_sample_summary.csv")
adult_chromvar_file <- file.path(data_dir, "2026-08-12_4_adult_chromVAR_sample_summary.csv")

foreground_qc_out <- file.path(output_dir, paste0(date_str, "_DM_motif_enrichment_foreground_QC.csv"))
enrichment_full_out <- file.path(output_dir, paste0(date_str, "_DM_motif_enrichment_full.csv"))
family_collapsed_out <- file.path(output_dir, paste0(date_str, "_DM_motif_enrichment_family_collapsed.csv"))
mapping_out <- file.path(output_dir, paste0(date_str, "_DM_motif_family_annotation.csv"))
chromvar_regional_out <- file.path(output_dir, paste0(date_str, "_DM_motif_chromVAR_regional_summary.csv"))
adult_central_out <- file.path(output_dir, paste0(date_str, "_DM_motif_adult_central_chromVAR_summary.csv"))
integrated_out <- file.path(output_dir, paste0(date_str, "_DM_motif_integrated_evidence.csv"))
qc_out <- file.path(output_dir, paste0(date_str, "_DM_motif_QC.csv"))
session_out <- file.path(output_dir, paste0(date_str, "_DM_motif_parameters_sessionInfo.txt"))

celltype_order <- c("Rod", "Cone", "AC", "BC", "HC", "MG")
development_stages <- c("Early", "Mid", "Late")
direction_order <- c("M", "P")
min_foreground_peaks <- 10L
motif_assay <- "ATAC"

# ---- 2. 辅助函数 ------------------------------------------------------------
check_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " 不存在：", path, call. = FALSE)
}
check_columns <- function(x, required, label) {
  miss <- setdiff(required, colnames(x))
  if (length(miss)) stop(label, " 缺少必要列：", paste(miss, collapse = ", "), call. = FALSE)
}
norm_text <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(x) | y %in% c("", "NA", "NaN", "<NA>")] <- NA_character_
  y
}
as_num <- function(x) suppressWarnings(as.numeric(x))
is_mp <- function(x) !is.na(x) & x %in% direction_order
same_mp <- function(rna, dar) is_mp(rna) & is_mp(dar) & rna == dar
mean_or_na <- function(x) {
  x <- as_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}
median_or_na <- function(x) {
  x <- as_num(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}
parse_peak_width <- function(peaks) {
  mm <- stringr::str_match(as.character(peaks), "^(.+)-(\\d+)-(\\d+)$")
  ans <- as_num(mm[, 4]) - as_num(mm[, 3]) + 1
  ans[!is.finite(ans)] <- NA_real_
  ans
}
standardize_findmotifs <- function(x) {
  if (is.null(x) || !nrow(x)) return(tibble::tibble())
  y <- as.data.frame(x, stringsAsFactors = FALSE)
  y$motif_id <- if ("motif" %in% colnames(y)) as.character(y$motif) else rownames(y)
  y$motif_name <- if ("motif.name" %in% colnames(y)) as.character(y$motif.name) else NA_character_
  y$pvalue <- if ("pvalue" %in% colnames(y)) as_num(y$pvalue) else if ("p.value" %in% colnames(y)) as_num(y$p.value) else NA_real_
  y$padj <- if ("p.adjust" %in% colnames(y)) as_num(y$p.adjust) else if ("p_val_adj" %in% colnames(y)) as_num(y$p_val_adj) else stats::p.adjust(y$pvalue, "BH")
  y$fold_enrichment <- if ("fold.enrichment" %in% colnames(y)) as_num(y$fold.enrichment) else NA_real_
  y$foreground_motif_n <- if ("foreground.motif" %in% colnames(y)) as_num(y$foreground.motif) else NA_real_
  y$background_motif_n <- if ("background.motif" %in% colnames(y)) as_num(y$background.motif) else NA_real_
  tibble::as_tibble(y)
}
motif_id_of <- function(motif) {
  ans <- tryCatch(as.character(TFBSTools::ID(motif)), error = function(e) NA_character_)
  if (!length(ans)) NA_character_ else ans[[1]]
}
motif_name_of <- function(motif) {
  ans <- tryCatch(as.character(TFBSTools::name(motif)), error = function(e) NA_character_)
  if (!length(ans)) NA_character_ else ans[[1]]
}
extract_tag <- function(motif, key) {
  tg <- tryCatch(TFBSTools::tags(motif), error = function(e) NULL)
  if (is.null(tg) || is.null(tg[[key]])) return(NA_character_)
  ans <- as.character(tg[[key]])
  ans <- ans[!is.na(ans) & trimws(ans) != ""]
  if (!length(ans)) NA_character_ else paste(unique(ans), collapse = ";")
}
match_background <- function(foreground, eligible, feature_stats, context_id) {
  foreground <- unique(as.character(foreground))
  eligible <- unique(as.character(eligible))
  if (!length(foreground)) return(character())
  if (!all(foreground %in% eligible)) stop(context_id, "：foreground 不在 tested universe 中。", call. = FALSE)
  pool <- setdiff(eligible, foreground)
  if (length(pool) < length(foreground)) return(character())
  stats_tbl <- feature_stats[unique(c(foreground, pool)), , drop = FALSE]
  stats_tbl <- stats_tbl[stats::complete.cases(stats_tbl), , drop = FALSE]
  if (!all(foreground %in% rownames(stats_tbl))) return(character())
  hit <- Signac::MatchRegionStats(meta.feature = stats_tbl, target = foreground,
                                  features = pool, n = length(foreground))
  if (is.data.frame(hit) && "feature" %in% colnames(hit)) hit <- hit$feature
  unique(as.character(hit))[unique(as.character(hit)) %in% pool]
}
direction_status <- function(delta, direction) {
  delta <- delta[is.finite(delta)]
  if (!length(delta) || !direction %in% direction_order) return("Unavailable")
  same <- if (direction == "M") delta > 0 else delta < 0
  same_fraction <- mean(same)
  med <- stats::median(delta)
  if (same_fraction > 0.5 && ((direction == "M" && med > 0) || (direction == "P" && med < 0))) return("Supported")
  if (same_fraction > 0.5 && ((direction == "M" && med < 0) || (direction == "P" && med > 0))) return("Opposite")
  "Mixed"
}
get_motif_names <- function(object, motif_ids) {
  motif_obj <- Signac::Motifs(object, assay = motif_assay)
  out <- tibble::tibble(motif_id = as.character(motif_ids), motif_name = NA_character_)
  if ("motif.names" %in% slotNames(motif_obj)) {
    nm <- slot(motif_obj, "motif.names")
    if (!is.null(names(nm))) {
      out <- out |>
        dplyr::left_join(tibble::tibble(motif_id = names(nm), motif_name_new = as.character(nm)),
                         by = "motif_id") |>
        dplyr::mutate(motif_name = dplyr::coalesce(motif_name_new, motif_name)) |>
        dplyr::select(-dplyr::any_of("motif_name_new"))
    }
  }
  out |>
    dplyr::mutate(TF_name = motif_name)
}

# ---- 3. 读取 pairing 主表与 DM foreground -----------------------------------
for (f in c(pair_file, gene_summary_file, develop_obj_file, adult_obj_file,
            adult_dar_file, develop_chromvar_file, adult_chromvar_file)) {
  check_file(f, basename(f))
}
pair_master <- readr::read_csv(pair_file, show_col_types = FALSE, progress = FALSE)
gene_summary <- readr::read_csv(gene_summary_file, show_col_types = FALSE, progress = FALSE)
check_columns(pair_master, c(
  "celltype", "gene", "peak", "subclass", "direction_class",
  "Early_status", "Mid_status", "Late_status",
  "Early_DAR_status", "Mid_DAR_status", "Late_DAR_status",
  "embryonic_grade", "adult_grade", "Adult_RNA_direction", "Adult_shared_direction"
), "DEG-DAR peak-level 主表")
check_columns(gene_summary, c("celltype", "gene", "subclass"), "DEG-DAR gene summary")

pair_master <- pair_master |>
  dplyr::mutate(dplyr::across(dplyr::any_of(c(
    "celltype", "gene", "peak", "subclass", "direction_class",
    "Early_status", "Mid_status", "Late_status",
    "Early_DAR_status", "Mid_DAR_status", "Late_DAR_status",
    "embryonic_grade", "adult_grade", "Adult_RNA_direction", "Adult_shared_direction"
  )), norm_text))
pair_key <- c("celltype", "gene", "peak")
if (anyNA(pair_master[pair_key])) stop("DEG-DAR 主表关键字段包含 NA。", call. = FALSE)
if (anyDuplicated(pair_master[pair_key]) > 0L) stop("DEG-DAR 主表 celltype × gene × peak 不唯一。", call. = FALSE)

dm <- pair_master |>
  dplyr::filter(subclass == "Developmental maintained") |>
  dplyr::mutate(
    embryonic_gold_strong = embryonic_grade %in% c("Gold", "Strong"),
    adult_gold_strong = adult_grade %in% c("Gold", "Strong")
  )
if (!nrow(dm)) stop("主表中没有 Developmental maintained 行。", call. = FALSE)

development_source <- dplyr::bind_rows(lapply(development_stages, function(st) {
  rna_col <- paste0(st, "_status")
  dar_col <- paste0(st, "_DAR_status")
  dm |>
    dplyr::filter(
      embryonic_gold_strong,
      .data[[rna_col]] %in% direction_order,
      .data[[dar_col]] %in% direction_order,
      .data[[rna_col]] == .data[[dar_col]]
    ) |>
    dplyr::transmute(
      period = "Development", celltype, gene, peak, stage_period = st,
      direction = .data[[rna_col]], source_RNA_status = .data[[rna_col]],
      source_DAR_status = .data[[dar_col]], embryonic_grade, adult_grade
    )
}))
adult_source <- dm |>
  dplyr::filter(
    adult_gold_strong,
    Adult_RNA_direction %in% direction_order,
    Adult_shared_direction %in% direction_order,
    Adult_RNA_direction == Adult_shared_direction
  ) |>
  dplyr::transmute(
    period = "Adult", celltype, gene, peak, stage_period = "Adult_shared",
    direction = Adult_RNA_direction, source_RNA_status = Adult_RNA_direction,
    source_DAR_status = Adult_shared_direction, embryonic_grade, adult_grade
  )
foreground_source <- dplyr::bind_rows(development_source, adult_source) |>
  dplyr::distinct(period, celltype, direction, stage_period, gene, peak, .keep_all = TRUE)
if (nrow(foreground_source) && any(!foreground_source$source_RNA_status %in% direction_order)) {
  stop("foreground 中存在非 M/P direction，说明 N=N stage 被错误纳入。", call. = FALSE)
}
if (nrow(foreground_source) && any(foreground_source$source_RNA_status != foreground_source$source_DAR_status)) {
  stop("foreground 中存在 RNA/DAR direction 不一致的 pair。", call. = FALSE)
}

foreground <- foreground_source |>
  dplyr::group_by(period, celltype, direction, peak) |>
  dplyr::summarise(
    source_stages = paste(sort(unique(stage_period)), collapse = ";"),
    n_source_stages = dplyr::n_distinct(stage_period),
    source_genes = paste(sort(unique(gene)), collapse = ";"),
    n_source_genes = dplyr::n_distinct(gene),
    n_source_pairs = dplyr::n(),
    .groups = "drop"
  )
context_grid <- tidyr::expand_grid(
  period = c("Development", "Adult"),
  celltype = celltype_order,
  direction = direction_order
)
foreground_stats <- foreground_source |>
  dplyr::group_by(period, celltype, direction) |>
  dplyr::summarise(
    n_foreground_peaks = dplyr::n_distinct(peak),
    n_source_pairs = dplyr::n(),
    n_source_genes = dplyr::n_distinct(gene),
    .groups = "drop"
  )
foreground_qc <- context_grid |>
  dplyr::left_join(foreground_stats, by = c("period", "celltype", "direction")) |>
  dplyr::mutate(
    n_foreground_peaks = dplyr::coalesce(n_foreground_peaks, 0L),
    n_source_pairs = dplyr::coalesce(n_source_pairs, 0L),
    n_source_genes = dplyr::coalesce(n_source_genes, 0L),
    meets_minimum_foreground = n_foreground_peaks >= min_foreground_peaks,
    enrichment_status = dplyr::if_else(meets_minimum_foreground, "ready", "Insufficient_foreground")
  )

# ---- 4. motif object、peak features 和 tested universe -----------------------
develop_obj <- readRDS(develop_obj_file)
adult_obj <- readRDS(adult_obj_file)
if (!(motif_assay %in% SeuratObject::Assays(develop_obj)) ||
    !(motif_assay %in% SeuratObject::Assays(adult_obj))) {
  stop("Development/Adult 对象均必须包含 ATAC assay。", call. = FALSE)
}
develop_peaks <- rownames(develop_obj[[motif_assay]])
adult_peaks <- rownames(adult_obj[[motif_assay]])
if (!setequal(develop_peaks, adult_peaks)) stop("Development/Adult ATAC peak universe 不一致。", call. = FALSE)
unified_peaks <- sort(develop_peaks)

motif_matrix <- Signac::GetMotifData(develop_obj, assay = motif_assay, slot = "data")
adult_motif_matrix <- Signac::GetMotifData(adult_obj, assay = motif_assay, slot = "data")
if (!setequal(rownames(motif_matrix), rownames(adult_motif_matrix)) ||
    !setequal(colnames(motif_matrix), colnames(adult_motif_matrix))) {
  stop("Development/Adult motif presence matrix 的 peak/motif ID 集合不一致。", call. = FALSE)
}
motif_matrix <- motif_matrix[unified_peaks, sort(colnames(motif_matrix)), drop = FALSE]
adult_motif_matrix <- adult_motif_matrix[unified_peaks, sort(colnames(adult_motif_matrix)), drop = FALSE]
if (Matrix::nnzero(motif_matrix - adult_motif_matrix) != 0L) {
  stop("按 peak/motif ID 对齐后 motif presence matrix 不一致。", call. = FALSE)
}
meta_features <- develop_obj[[motif_assay]]@meta.features
if (!"GC.percent" %in% colnames(meta_features)) stop("ATAC assay 缺少 GC.percent。", call. = FALSE)
meta_features <- meta_features[unified_peaks, , drop = FALSE]
meta_features$GC.percent <- as_num(meta_features$GC.percent)
meta_features$peak_width <- parse_peak_width(rownames(meta_features))
if (any(!is.finite(meta_features$GC.percent)) || any(!is.finite(meta_features$peak_width))) {
  stop("GC.percent 或 peak width 存在非有限值。", call. = FALSE)
}
motif_map <- get_motif_names(develop_obj, colnames(motif_matrix))

develop_test_files <- list.files(
  develop_dar_dir,
  pattern = "^2026-07-06_develop_1_.+_(Early|Mid|Late)_macula_vs_peripheral\\.csv$",
  full.names = TRUE
)
if (!length(develop_test_files)) stop("未找到 Early/Mid/Late developmental DAR tested 文件。", call. = FALSE)
develop_tested <- dplyr::bind_rows(lapply(develop_test_files, function(path) {
  nm <- basename(path)
  mm <- stringr::str_match(nm, "^2026-07-06_develop_1_(.+)_(Early|Mid|Late)_macula_vs_peripheral\\.csv$")
  tab <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  check_columns(tab, c("peak", "baseMean", "log2FoldChange", "padj"), nm)
  tibble::tibble(
    celltype = mm[, 2], stage_period = mm[, 3], peak = norm_text(tab$peak),
    baseMean = as_num(tab$baseMean), log2FoldChange = as_num(tab$log2FoldChange),
    padj = as_num(tab$padj)
  ) |>
    dplyr::filter(!is.na(peak), is.finite(baseMean), is.finite(log2FoldChange),
                  is.finite(padj), peak %in% unified_peaks)
})) |>
  dplyr::group_by(celltype, stage_period, peak) |>
  dplyr::summarise(baseMean = mean(baseMean), .groups = "drop")

adult_raw <- readr::read_csv(adult_dar_file, show_col_types = FALSE, progress = FALSE)
check_columns(adult_raw, c("celltype", "peak", "baseMean", "log2FoldChange", "padj"), "Adult shared DAR")
adult_tested <- adult_raw |>
  dplyr::transmute(
    celltype = norm_text(celltype), stage_period = "Adult_shared", peak = norm_text(peak),
    baseMean = as_num(baseMean), log2FoldChange = as_num(log2FoldChange), padj = as_num(padj)
  ) |>
  dplyr::filter(!is.na(celltype), !is.na(peak), peak %in% unified_peaks,
                is.finite(baseMean), is.finite(log2FoldChange), is.finite(padj)) |>
  dplyr::group_by(celltype, stage_period, peak) |>
  dplyr::summarise(baseMean = mean(baseMean), .groups = "drop")

tested_by_period <- dplyr::bind_rows(
  develop_tested |> dplyr::mutate(period = "Development"),
  adult_tested |> dplyr::mutate(period = "Adult")
) |>
  dplyr::group_by(period, celltype, peak) |>
  dplyr::summarise(accessibility_proxy = mean(baseMean), .groups = "drop")
if (anyDuplicated(tested_by_period[c("period", "celltype", "peak")]) > 0L) {
  stop("tested background universe 键重复。", call. = FALSE)
}

# ---- 5. DM-specific FindMotifs ----------------------------------------------
enrichment_list <- list()
background_qc_list <- list()
for (i in seq_len(nrow(foreground_qc))) {
  period_i <- as.character(foreground_qc$period[i])
  celltype_i <- as.character(foreground_qc$celltype[i])
  direction_i <- as.character(foreground_qc$direction[i])
  context_id <- paste(period_i, celltype_i, direction_i, sep = "|")
  fg <- foreground |>
    dplyr::filter(period == period_i, celltype == celltype_i, direction == direction_i) |>
    dplyr::pull(peak) |>
    unique()
  eligible <- tested_by_period |>
    dplyr::filter(period == period_i, celltype == celltype_i) |>
    dplyr::pull(peak) |>
    unique()
  if (length(fg) && !all(fg %in% eligible)) stop(context_id, "：foreground 不在 tested universe 中。", call. = FALSE)
  stats_ctx <- meta_features[eligible, c("GC.percent", "peak_width"), drop = FALSE]
  proxy <- tested_by_period |>
    dplyr::filter(period == period_i, celltype == celltype_i) |>
    dplyr::select(peak, accessibility_proxy)
  stats_ctx$accessibility_proxy <- proxy$accessibility_proxy[match(rownames(stats_ctx), proxy$peak)]
  matched_bg <- if (length(fg) >= min_foreground_peaks) {
    match_background(fg, eligible, stats_ctx, context_id)
  } else character()
  background_qc_list[[i]] <- tibble::tibble(
    period = period_i, celltype = celltype_i, direction = direction_i,
    n_eligible_background_peaks = length(eligible),
    n_matched_background_peaks = length(matched_bg)
  )
  if (length(fg) < min_foreground_peaks || length(matched_bg) < min_foreground_peaks) next

  motifs <- tryCatch(
    Signac::FindMotifs(object = develop_obj, features = fg, background = matched_bg, assay = motif_assay),
    error = function(e) stop(context_id, " FindMotifs 失败：", conditionMessage(e), call. = FALSE)
  )
  one <- standardize_findmotifs(motifs)
  if (!nrow(one)) next
  fg_counts <- Matrix::colSums(motif_matrix[fg, , drop = FALSE])
  bg_counts <- Matrix::colSums(motif_matrix[matched_bg, , drop = FALSE])
  one <- one |>
    dplyr::mutate(
      period = period_i, celltype = celltype_i, direction = direction_i,
      enrichment_scope = "DM_Candidate_CRE",
      n_foreground_peaks = length(fg), n_background_peaks = length(matched_bg),
      foreground_motif_n = dplyr::coalesce(foreground_motif_n, as_num(fg_counts[match(motif_id, names(fg_counts))])),
      background_motif_n = dplyr::coalesce(background_motif_n, as_num(bg_counts[match(motif_id, names(bg_counts))]))
    ) |>
    dplyr::left_join(motif_map, by = "motif_id", suffix = c("", ".map")) |>
    dplyr::mutate(
      motif_name = dplyr::coalesce(motif_name, motif_name.map),
      TF_name = dplyr::coalesce(TF_name, motif_name),
      enrichment_positive = is.finite(padj) & padj < 0.05 &
        is.finite(fold_enrichment) & fold_enrichment > 1
    ) |>
    dplyr::select(-dplyr::any_of("motif_name.map"))
  enrichment_list[[i]] <- one
}
enrichment_full <- dplyr::bind_rows(enrichment_list)
if (!nrow(enrichment_full)) {
  enrichment_full <- tibble::tibble(
    period = character(), celltype = character(), direction = character(),
    enrichment_scope = character(), motif_id = character(), motif_name = character(),
    TF_name = character(), n_foreground_peaks = integer(), n_background_peaks = integer(),
    foreground_motif_n = numeric(), background_motif_n = numeric(),
    fold_enrichment = numeric(), pvalue = numeric(), padj = numeric(),
    enrichment_positive = logical()
  )
}
if (nrow(enrichment_full) && anyDuplicated(enrichment_full[c("period", "celltype", "direction", "motif_id")]) > 0L) {
  stop("motif enrichment context × motif_id 键重复。", call. = FALSE)
}
if (!nrow(enrichment_full)) warning("没有 context 返回 motif enrichment 结果。", call. = FALSE)

# ---- 6. JASPAR2020 family mapping 与 context-specific collapse ----------------
enrichment_ids <- sort(unique(enrichment_full$motif_id[!is.na(enrichment_full$motif_id)]))
jaspar_set <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606, collection = "CORE", all_versions = FALSE)
)
jaspar_map <- dplyr::bind_rows(lapply(jaspar_set, function(motif) {
  tibble::tibble(
    motif_id = motif_id_of(motif),
    jaspar_motif_name = motif_name_of(motif),
    jaspar_TF_name = motif_name_of(motif),
    JASPAR_family = extract_tag(motif, "family"),
    JASPAR_class = extract_tag(motif, "class")
  )
})) |>
  dplyr::filter(!is.na(motif_id)) |>
  dplyr::distinct(motif_id, .keep_all = TRUE)

mapping <- tibble::tibble(motif_id = enrichment_ids) |>
  dplyr::left_join(motif_map, by = "motif_id") |>
  dplyr::left_join(jaspar_map, by = "motif_id") |>
  dplyr::mutate(
    motif_name = dplyr::coalesce(jaspar_motif_name, motif_name),
    TF_name = dplyr::coalesce(jaspar_TF_name, TF_name, motif_name),
    family_annotation_status = dplyr::case_when(
      !is.na(jaspar_motif_name) & !is.na(JASPAR_family) & trimws(JASPAR_family) != "" ~ "Annotated",
      !is.na(jaspar_motif_name) ~ "Matched_no_family",
      TRUE ~ "Unmatched_ID"
    ),
    family_for_collapse = dplyr::case_when(
      family_annotation_status == "Annotated" ~ JASPAR_family,
      family_annotation_status == "Matched_no_family" ~ paste0("Singleton_", motif_id),
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::select(motif_id, motif_name, TF_name, JASPAR_family, JASPAR_class,
                family_annotation_status, family_for_collapse)
readr::write_csv(mapping, mapping_out, na = "")
if (anyDuplicated(mapping$motif_id) > 0L) stop("motif family mapping 不唯一。", call. = FALSE)
if (nrow(mapping) > 0L && any(mapping$family_annotation_status == "Unmatched_ID")) {
  stop("发现 Unmatched_ID motif；已写出 mapping 文件，请先核对 motif ID。", call. = FALSE)
}

cat("done: DM motif family annotation ->", mapping_out, "\n")
