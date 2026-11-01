# 最终整理脚本 06 · per-sample chromVAR 汇总计算 (胚胎+成年, candidate_CRE 项目) (2026-10-29)
# 来源: 00use_scripts/35_codex_2026-08-12_4_chromVAR_complete_summary.R (完整流水线执行脚本, 全部保留;
#       同名 ipynb 为查看/QC 版, 仅部分执行)
# 产出: candidate_CRE/Data/2026-08-12_4_chromVAR_motif_inventory.csv        (08 号 family map 的输入)
#       candidate_CRE/Data/2026-08-12_4_{develop,adult}_chromVAR_sample_summary.csv (09 号区域统计的输入)
#       + cell×motif RDS / formal regional summary / QC 表 / 参数记录
# 说明: 从两期 snATAC 对象的 chromVAR deviations 计算 per-sample (celltype×donor×region) 加权汇总;
#       本脚本是区域 chromVAR 分析链 (06→09→10) 的最上游
# ---- 0. 环境设置 -----------------------------------------------------------
Sys.setenv(SP_EVOLUTION_STATUS = "2")
set.seed(1)

required_pkgs <- c(
  "Signac", "SeuratObject", "Matrix", "dplyr", "tidyr", "readr", "tibble"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "缺少必要 R 包：", paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Signac)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ---- 1. 参数设置 ------------------------------------------------------------
date_str <- "2026-08-12_4"
output_dir <- "/project/work/2026-08-08_candidate_CRE/Data/"

develop_rds <- "/project/work/2026-08-03_develop_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds"
adult_rds <- "/project/work/2026-08-03_adult_ArchR/Data/2026-06-26_snATAC_with_motifs_chromVAR_ArchR_Peak2GeneLinks.rds"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

motif_inventory_out <- file.path(output_dir, paste0(date_str, "_chromVAR_motif_inventory.csv"))
develop_cell_motif_out <- file.path(output_dir, paste0(date_str, "_develop_chromVAR_cell_motif.rds"))
adult_cell_motif_out <- file.path(output_dir, paste0(date_str, "_adult_chromVAR_cell_motif.rds"))
develop_summary_out <- file.path(output_dir, paste0(date_str, "_develop_chromVAR_sample_summary.csv"))
adult_summary_out <- file.path(output_dir, paste0(date_str, "_adult_chromVAR_sample_summary.csv"))
formal_summary_out <- file.path(output_dir, paste0(date_str, "_chromVAR_formal_regional_summary.csv"))
qc_out <- file.path(output_dir, paste0(date_str, "_chromVAR_QC.csv"))
session_out <- file.path(output_dir, paste0(date_str, "_chromVAR_parameters_sessionInfo.txt"))

celltype_order <- c("Rod", "Cone", "AC", "BC", "HC", "MG", "RGC", "RPC")
develop_stage_order <- c("Early", "Mid", "Late")
adult_age_order <- c(
  "21years", "25years", "43years", "64years", "67years",
  "70years", "75years", "86years", "89years", ">89years"
)

develop_stage_map <- c(
  "70days" = "Early", "79days" = "Early", "87days" = "Early", "91days" = "Early",
  "100days" = "Mid", "103days" = "Mid", "116days" = "Mid",
  "137days" = "Late", "141days" = "Late", "142days" = "Late",
  "162days" = "Late", "165days" = "Late"
)

# ---- 2. 辅助函数 ------------------------------------------------------------
check_file <- function(path, label) {
  if (!file.exists(path)) stop(label, "不存在：", path, call. = FALSE)
}

check_columns <- function(x, required, label) {
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    stop(label, "缺少必要列：", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
}

norm_text <- function(x) {
  y <- trimws(as.character(x))
  y[y %in% c("", "NA", "NaN", "<NA>")] <- NA_character_
  y
}

first_metadata_column <- function(meta, candidates, label, required = FALSE) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) {
    if (required) {
      stop(label, "缺少以下任一列：", paste(candidates, collapse = ", "), call. = FALSE)
    }
    return(rep(NA_character_, nrow(meta)))
  }
  norm_text(meta[[hit[[1]]]])
}

get_assay_matrix <- function(object, assay, layer_name = "data") {
  mat <- tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = layer_name),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- SeuratObject::GetAssayData(object, assay = assay, slot = layer_name)
  }
  mat
}

get_motif_ids <- function(object, assay = "ATAC") {
  motif_matrix <- tryCatch(
    Signac::GetMotifData(object, assay = assay, slot = "data"),
    error = function(e) NULL
  )
  if (is.null(motif_matrix)) {
    stop("无法从 ", assay, " assay 获取 Motif data。", call. = FALSE)
  }
  colnames(motif_matrix)
}

get_motif_names <- function(object, motif_ids, assay = "ATAC") {
  out <- tibble::tibble(motif_id = as.character(motif_ids), motif_name = NA_character_)
  motif_obj <- tryCatch(Signac::Motifs(object, assay = assay), error = function(e) NULL)
  if (is.null(motif_obj)) return(out)

  slot_names <- methods::slotNames(motif_obj)
  if (!"motif.names" %in% slot_names) return(out)
  nm <- methods::slot(motif_obj, "motif.names")
  if (is.null(nm)) return(out)

  nm <- as.character(nm)
  if (!is.null(names(nm))) {
    name_tbl <- tibble::tibble(motif_id = names(nm), motif_name = unname(nm))
    return(out %>% dplyr::left_join(name_tbl, by = "motif_id", suffix = c("", ".new")) %>%
      dplyr::mutate(motif_name = dplyr::coalesce(motif_name, motif_name.new)) %>%
      dplyr::select(-dplyr::any_of("motif_name.new")))
  }
  if (length(nm) == nrow(out)) out$motif_name <- nm
  out
}

standardize_tissue <- function(x) {
  y <- tolower(gsub("\\s+", "", norm_text(x)))
  dplyr::case_when(
    grepl("fovea|macula|central|centralfovea", y) ~ "central",
    grepl("peripheral|periphery", y) ~ "peripheral",
    TRUE ~ NA_character_
  )
}

formal_pair_by_donor <- function(meta, dataset_period) {
  if (dataset_period == "Development") {
    central_pattern <- "macula|central"
    required_ages <- names(develop_stage_map)
  } else {
    central_pattern <- "fovea"
    required_ages <- c("21years", "67years")
  }

  donor_tissue <- meta %>%
    dplyr::transmute(
      donor_age,
      tissue_key = tolower(gsub("\\s+", "", norm_text(tissue)))
    ) %>%
    dplyr::filter(!is.na(donor_age), !is.na(tissue_key)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(donor_age) %>%
    dplyr::summarise(
      has_central = any(grepl(central_pattern, tissue_key)),
      has_peripheral = any(grepl("peripheral|periphery", tissue_key)),
      formal_donor_pair = has_central & has_peripheral,
      .groups = "drop"
    )

  meta %>%
    dplyr::left_join(donor_tissue, by = "donor_age") %>%
    dplyr::mutate(
      formal_regional_pair = !is.na(formal_donor_pair) &
        formal_donor_pair & donor_age %in% required_ages,
      regional_eligibility_reason = dplyr::case_when(
        dataset_period == "Adult" & donor_age %in% required_ages & formal_regional_pair ~ "required_adult_paired_donor",
        dataset_period == "Development" & donor_age %in% required_ages & formal_regional_pair ~ "paired_development_donor",
        dataset_period == "Adult" ~ "adult_not_formal_paired_donor",
        TRUE ~ "development_not_formal_paired_donor"
      )
    ) %>%
    dplyr::select(-formal_donor_pair)
}

make_cell_metadata <- function(object, dataset_period) {
  meta <- CONTACT_EMAIL
  meta$cell <- rownames(meta)
  check_columns(meta, c("celltype", "donor_age", "tissue"), paste0(dataset_period, " metadata"))
  meta <- meta %>%
    dplyr::transmute(
      cell,
      celltype = norm_text(celltype),
      sample = first_metadata_column(meta, c("sample", "dataset", "orig.ident"), "sample metadata"),
      donor_id = first_metadata_column(meta, c("donor_id", "donor", "donorID"), "donor metadata"),
      donor_age = gsub("\\s+", "", norm_text(donor_age)),
      tissue = norm_text(tissue),
      tissue_region = standardize_tissue(tissue),
      dataset_period = dataset_period
    ) %>%
    dplyr::mutate(
      stage_period = if (dataset_period == "Development") unname(develop_stage_map[donor_age]) else donor_age,
      celltype = as.character(celltype),
      stage_period = as.character(stage_period)
    )
  formal_pair_by_donor(meta, dataset_period)
}

matrix_qc <- function(mat, dataset_period, motif_ids, cell_ids) {
  feature_non_na_n <- Matrix::rowSums(!is.na(mat))
  feature_all_na <- rownames(mat)[feature_non_na_n == 0]
  values <- mat
  if (inherits(values, "sparseMatrix")) values <- values@x
  else values <- as.vector(values)
  tibble::tibble(
    dataset_period = dataset_period,
    n_cells = ncol(mat),
    n_motifs = nrow(mat),
    n_matrix_rows = nrow(mat),
    n_matrix_cols = ncol(mat),
    cell_ids_identical = identical(colnames(mat), cell_ids),
    n_motif_ids = length(motif_ids),
    n_na_values = sum(is.na(values)),
    n_nonfinite_values = sum(!is.finite(values)),
    n_all_na_features = length(feature_all_na),
    all_na_feature_ids = paste(feature_all_na, collapse = ";")
  )
}

sample_summary <- function(mat, cell_meta, motif_inventory, dataset_period) {
  if (!identical(colnames(mat), cell_meta$cell)) {
    stop(dataset_period, " chromVAR matrix 与 cell metadata 顺序不一致。", call. = FALSE)
  }

  group_cols <- c("dataset_period", "celltype", "stage_period", "donor_age", "tissue", "tissue_region", "sample", "donor_id", "formal_regional_pair")
  groups <- cell_meta %>% dplyr::select(dplyr::all_of(group_cols)) %>% dplyr::mutate(.cell_index = seq_len(n()))
  group_key <- do.call(paste, c(groups[group_cols], sep = "\r"))
  split_indices <- split(seq_len(nrow(groups)), group_key)

  rows <- lapply(split_indices, function(idx) {
    one_meta <- groups[idx, , drop = FALSE]
    one_mat <- mat[, idx, drop = FALSE]
    means <- Matrix::rowMeans(one_mat)
    medians <- apply(as.matrix(one_mat), 1, stats::median, na.rm = TRUE)
    sds <- apply(as.matrix(one_mat), 1, stats::sd, na.rm = TRUE)
    n_valid <- Matrix::rowSums(!is.na(one_mat))
    tibble::tibble(
      dataset_period = one_meta$dataset_period[[1]],
      celltype = one_meta$celltype[[1]],
      stage_period = one_meta$stage_period[[1]],
      donor_age = one_meta$donor_age[[1]],
      tissue = one_meta$tissue[[1]],
      tissue_region = one_meta$tissue_region[[1]],
      sample = one_meta$sample[[1]],
      donor_id = one_meta$donor_id[[1]],
      formal_regional_pair = one_meta$formal_regional_pair[[1]],
      motif_id = rownames(mat),
      mean_chromVAR = as.numeric(means),
      median_chromVAR = as.numeric(medians),
      sd_chromVAR = as.numeric(sds),
      n_cells = nrow(one_meta),
      n_valid_chromVAR = as.integer(n_valid)
    )
  })

  dplyr::bind_rows(rows) %>%
    dplyr::left_join(motif_inventory %>% dplyr::select(motif_id, motif_name), by = "motif_id")
}

# ---- 3. 读取对象并检查完整 chromVAR assay ----------------------------------
check_file(develop_rds, "development Signac 对象")
check_file(adult_rds, "adult Signac 对象")

develop_obj <- readRDS(develop_rds)
adult_obj <- readRDS(adult_rds)

validate_object <- function(object, dataset_period) {
  if (!("ATAC" %in% SeuratObject::Assays(object))) stop(dataset_period, " 对象缺少 ATAC assay。", call. = FALSE)
  if (!("chromvar" %in% SeuratObject::Assays(object))) stop(dataset_period, " 对象缺少现有 chromvar assay；本脚本不会重新计算 chromVAR。", call. = FALSE)
  motif_matrix <- tryCatch(Signac::GetMotifData(object, assay = "ATAC", slot = "data"), error = function(e) NULL)
  if (is.null(motif_matrix)) stop(dataset_period, " 对象缺少可读取的 ATAC Motif object。", call. = FALSE)
  chromvar_matrix <- get_assay_matrix(object, "chromvar", "data")
  if (is.null(chromvar_matrix)) stop(dataset_period, " 无法读取 chromvar data layer。", call. = FALSE)
  if (!identical(colnames(chromvar_matrix), colnames(object))) stop(dataset_period, " chromVAR cell IDs 与对象细胞 IDs 不一致。", call. = FALSE)
  if (!all(rownames(chromvar_matrix) %in% colnames(motif_matrix))) stop(dataset_period, " chromVAR motif IDs 无法对应 ATAC Motif object。", call. = FALSE)
  list(motif_matrix = motif_matrix, chromvar_matrix = chromvar_matrix)
}

develop_valid <- validate_object(develop_obj, "Development")
adult_valid <- validate_object(adult_obj, "Adult")

develop_motif_ids <- colnames(develop_valid$motif_matrix)
adult_motif_ids <- colnames(adult_valid$motif_matrix)
all_motif_ids <- union(develop_motif_ids, adult_motif_ids)

develop_motif_names <- get_motif_names(develop_obj, develop_motif_ids)
adult_motif_names <- get_motif_names(adult_obj, adult_motif_ids)
motif_inventory <- tibble::tibble(motif_id = all_motif_ids) %>%
  dplyr::left_join(develop_motif_names %>% dplyr::rename(motif_name_development = motif_name), by = "motif_id") %>%
  dplyr::left_join(adult_motif_names %>% dplyr::rename(motif_name_adult = motif_name), by = "motif_id") %>%
  dplyr::mutate(
    motif_name = dplyr::coalesce(motif_name_development, motif_name_adult),
    TF_annotation = motif_name,
    present_in_development = motif_id %in% develop_motif_ids,
    present_in_adult = motif_id %in% adult_motif_ids
  ) %>%
  dplyr::select(motif_id, motif_name, TF_annotation, present_in_development, present_in_adult)

readr::write_csv(motif_inventory, motif_inventory_out, na = "")

if (!all(rownames(develop_valid$chromvar_matrix) %in% all_motif_ids) ||
    !all(rownames(adult_valid$chromvar_matrix) %in% all_motif_ids)) {
  stop("chromVAR feature IDs 不在统一 motif inventory 中。", call. = FALSE)
}

# ---- 4. 整理 cell metadata 与完整 cell × motif RDS --------------------------
develop_meta <- make_cell_metadata(develop_obj, "Development")
adult_meta <- make_cell_metadata(adult_obj, "Adult")

if (!identical(colnames(develop_valid$chromvar_matrix), develop_meta$cell)) stop("Development cell metadata 顺序错误。", call. = FALSE)
if (!identical(colnames(adult_valid$chromvar_matrix), adult_meta$cell)) stop("Adult cell metadata 顺序错误。", call. = FALSE)

develop_cell_motif <- list(
  dataset_period = "Development",
  chromvar_matrix = develop_valid$chromvar_matrix,
  cell_metadata = develop_meta,
  motif_inventory = motif_inventory,
  note = "Complete existing chromVAR matrix; no motif pre-filter and no chromVAR re-computation."
)
adult_cell_motif <- list(
  dataset_period = "Adult",
  chromvar_matrix = adult_valid$chromvar_matrix,
  cell_metadata = adult_meta,
  motif_inventory = motif_inventory,
  note = "Complete existing chromVAR matrix; no motif pre-filter and no chromVAR re-computation."
)

saveRDS(develop_cell_motif, develop_cell_motif_out)
saveRDS(adult_cell_motif, adult_cell_motif_out)

# ---- 5. sample-level chromVAR summary --------------------------------------
develop_summary <- sample_summary(
  develop_valid$chromvar_matrix,
  develop_meta,
  motif_inventory,
  "Development"
)
adult_summary <- sample_summary(
  adult_valid$chromvar_matrix,
  adult_meta,
  motif_inventory,
  "Adult"
)

readr::write_csv(develop_summary, develop_summary_out, na = "")
readr::write_csv(adult_summary, adult_summary_out, na = "")

# ---- 6. formal regional comparison-ready summary ---------------------------
formal_summary <- dplyr::bind_rows(
  develop_summary %>% dplyr::filter(formal_regional_pair, tissue_region %in% c("central", "peripheral")),
  adult_summary %>% dplyr::filter(formal_regional_pair, donor_age %in% c("21years", "67years"), tissue_region %in% c("central", "peripheral"))
)
readr::write_csv(formal_summary, formal_summary_out, na = "")

# ---- 7. QC ------------------------------------------------------------------
develop_matrix_qc <- matrix_qc(
  develop_valid$chromvar_matrix,
  "Development",
  rownames(develop_valid$chromvar_matrix),
  colnames(develop_obj)
)
adult_matrix_qc <- matrix_qc(
  adult_valid$chromvar_matrix,
  "Adult",
  rownames(adult_valid$chromvar_matrix),
  colnames(adult_obj)
)

motif_qc <- tibble::tibble(
  dataset_period = "Combined",
  n_development_motifs = length(develop_motif_ids),
  n_adult_motifs = length(adult_motif_ids),
  n_shared_motifs = length(intersect(develop_motif_ids, adult_motif_ids)),
  n_development_only_motifs = length(setdiff(develop_motif_ids, adult_motif_ids)),
  n_adult_only_motifs = length(setdiff(adult_motif_ids, develop_motif_ids)),
  n_union_motifs = length(all_motif_ids)
)

cell_count_qc <- dplyr::bind_rows(develop_meta, adult_meta) %>%
  dplyr::count(dataset_period, celltype, stage_period, donor_age, tissue_region, formal_regional_pair, name = "n_cells")

formal_group_qc <- dplyr::bind_rows(develop_meta, adult_meta) %>%
  dplyr::filter(formal_regional_pair) %>%
  dplyr::count(dataset_period, donor_age, tissue_region, name = "n_cells") %>%
  dplyr::mutate(
    expected_formal_group = (dataset_period == "Development" & donor_age %in% names(develop_stage_map)) |
      (dataset_period == "Adult" & donor_age %in% c("21years", "67years"))
  )

required_adult_groups <- tibble::tibble(
  donor_age = c("21years", "21years", "67years", "67years"),
  tissue_region = c("central", "peripheral", "central", "peripheral")
) %>%
  dplyr::left_join(
    adult_meta %>% dplyr::count(donor_age, tissue_region, name = "n_cells"),
    by = c("donor_age", "tissue_region")
  ) %>%
  dplyr::mutate(
    n_cells = dplyr::coalesce(n_cells, 0L),
    present = n_cells > 0L
  )

if (!all(required_adult_groups$present)) {
  warning("Adult formal regional comparison 中至少一个 21y/67y 年龄-组织组合不存在；详见 QC。", call. = FALSE)
}

qc_table <- dplyr::bind_rows(
  develop_matrix_qc %>% dplyr::mutate(qc_type = "matrix"),
  adult_matrix_qc %>% dplyr::mutate(qc_type = "matrix")
) %>%
  dplyr::mutate(
    n_development_motifs = motif_qc$n_development_motifs[[1]],
    n_adult_motifs = motif_qc$n_adult_motifs[[1]],
    n_shared_motifs = motif_qc$n_shared_motifs[[1]],
    n_development_only_motifs = motif_qc$n_development_only_motifs[[1]],
    n_adult_only_motifs = motif_qc$n_adult_only_motifs[[1]],
    n_union_motifs = motif_qc$n_union_motifs[[1]]
  )

qc_table <- dplyr::bind_rows(
  qc_table,
  cell_count_qc %>% dplyr::mutate(qc_type = "cell_count") %>%
    dplyr::mutate(
      dataset_period = as.character(dataset_period),
      n_cells = as.numeric(n_cells)
    ),
  formal_group_qc %>% dplyr::mutate(qc_type = "formal_group"),
  required_adult_groups %>% dplyr::mutate(dataset_period = "Adult", qc_type = "required_adult_group")
)

readr::write_csv(qc_table, qc_out, na = "")

# ---- 8. 参数与运行记录 ------------------------------------------------------
parameter_lines <- c(
  paste0("date_str = ", date_str),
  paste0("development_rds = ", develop_rds),
  paste0("adult_rds = ", adult_rds),
  paste0("output_dir = ", output_dir),
  "Existing chromVAR deviations were read; chromVAR was not re-computed.",
  "All existing chromVAR motifs were retained; no motif enrichment or candidate CRE pre-filter was used.",
  "formal_regional_pair is donor-level: a donor_age is formal only when both central and peripheral tissues exist in that object.",
  "Development central tissue is macula; Adult central tissue is fovea.",
  "Adult formal regional groups are restricted to 21years and 67years; all other adult cells remain in complete summaries.",
  "No differential chromVAR, TF expression, ArchR Link filtering, or TF-CRE-gene analysis was performed.",
  capture.output(sessionInfo())
)
writeLines(parameter_lines, con = session_out, useBytes = TRUE)

message("完成完整 chromVAR 整理与 QC。输出目录：", output_dir)
