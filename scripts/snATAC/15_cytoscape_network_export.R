# 最终整理脚本 15 · Cytoscape 网络输入构建 (Development / Adult 两个 directed 网络) (2026-10-29)
# 来源: 00use_scripts/32_2026-10-10_plot_phase51_2_cytoscape.ipynb (全部 cell)
# 产出: 2026-09-16_ATAC_start/Table/2026-10-10_plot_phase51_2_cytoscape_{Development,Adult}_edges.tsv / _nodes.tsv / _network.graphml
# 原则: DAR 为 edge-level 证据标注 (保留 non-DAR edge, 不做 hard filter); node = gene symbol (不拆 TF&DEG);
#       同一 TF→DEG 跨 celltype 保留多 edge; subclass 冲突显式记录 (LRP2/PDZRN3)
# 输入: 12 的 TF_DEG_DAR_support_summary.tsv (Phase 5.1 Version B)
# 2026-10-10 · Phase 5.1 → Cytoscape 网络输入文件 (仅数据整理, 不重算任何 Phase 1-5.1 结果)
# =====================================================================
# 目的: 基于 Phase 5.1 正式 TF-target DEG-DAR 结果 (Version B), 整理 Cytoscape 可直接使用的
#       Development / Adult 两个独立 directed 网络输入 (edges.tsv + nodes.tsv + GraphML).
# 原则: DAR 是 edge-level evidence annotation (保留 non-DAR edge, 不做 hard filter);
#       node = gene symbol (不拆 TF&DEG); 同一 TF->DEG 跨 celltype 保留多 edge;
#       subclass 冲突显式记录 (LRP2/PDZRN3 已知跨 celltype 多 subclass — 记录不静默);
#       不筛 top TF / 不设计 Cytoscape layout.
suppressPackageStartupMessages({
  library(data.table); library(igraph); library(jsonlite)
})
options(stringsAsFactors = FALSE, width = 150)

DATE_PREFIX <- "2026-10-10"
ROOT <- "/project/work"; OUT <- file.path(ROOT, "2026-09-16_ATAC_start")
DIR_TAB <- file.path(OUT, "Table")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)

## 输入: Phase 5.1 正式结果 (summary 为 TF-DEG 级, 含 DAR 支持计数; 只用 Version B — 5.1 本身仅 B)
SUM_F <- file.path(DIR_TAB, "2026-10-09_Phase5_TF_DEG_DAR_support_summary.tsv")
PK_F  <- file.path(DIR_TAB, "2026-10-09_Phase5_TF_DEG_supporting_peak_DAR_overlap_peaklevel.tsv")   # 必要时回查
stopifnot(file.exists(SUM_F))

## celltype 配色 (color.ipynb 现有 cols.celltype; 不重新定义)
source_color_ipynb <- function(path = file.path(ROOT, "color.ipynb")) {
  j <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  codes <- vapply(j$cells, function(cc)
    if (identical(cc$cell_type, "code")) paste(unlist(cc$source), collapse = "") else "", character(1))
  eval(parse(text = paste(codes, collapse = "\n")), envir = globalenv())
}
invisible(tryCatch(source_color_ipynb(), error = function(e) message("color.ipynb: ", conditionMessage(e))))
stopifnot(exists("cols.celltype"))
cat("参数就绪\n")
## ---- 1. 读入 Phase 5.1 summary 并构建 edge table ----
s5 <- fread(SUM_F, na.strings = c("NA",""))
## relation_to_TF_direction: 沿用 Phase 3/4 定义 (TF_direction vs DEG_direction, same→match)
s5[, relation_to_TF_direction := fifelse(
  (TF_direction == "central-high" & DEG_direction == "central") |
  (TF_direction == "peripheral-high" & DEG_direction == "peripheral"), "match", "opposite")]
## 字段映射 (5.1 实际列名 → Cytoscape 输出列名; 数值不改):
## names(edge_cols) = 原始列名 (motif_name 等); values = 输出列名 (TF 等)
## [2026-10-10 修复] d[, ..edge_cols] 会用 values 选列 (原始表无 TF 列会报错) — 先按原始列名取, 再重命名
edge_cols <- c(
  motif_name = "TF", target_gene = "target_gene",
  period = "period", celltype = "celltype", subclass = "subclass",
  TF_direction = "TF_direction", DEG_direction = "DEG_direction",
  relation_to_TF_direction = "relation_to_TF_direction",
  n_supporting_peaks = "n_supporting_peaks",
  n_DAR_supporting_peaks = "n_DAR_supporting_peaks",
  DAR_support_fraction = "DAR_support_fraction",
  has_DAR_support = "has_DAR_support",
  supporting_peak_list = "supporting_peak_list",
  DAR_supporting_peak_list = "DAR_peak_list",
  DAR_direction_list = "DAR_direction_list")

build_edges <- function(per_lab){
  d <- s5[period == ifelse(per_lab == "Development", "dev", "adult")]
  input_cols <- names(edge_cols)                    # 原始列名 (motif_name 等)
  e <- d[, ..input_cols]                            # 按原始列名取子集
  setnames(e, old = names(edge_cols), new = unname(edge_cols))   # 重命名为输出列名
  e[, source := TF]
  e[, target := target_gene]
  e[, edge_id := paste(period, celltype, subclass, TF, target_gene, sep = "||")]
  e[, period := per_lab]                                        # 统一为 Development/Adult 标签
  ## has_DAR_support / fraction 复核 (不改数值, 只核验):
  stopifnot(all(e$has_DAR_support == (e$n_DAR_supporting_peaks >= 1)))
  ## celltype_color (color.ipynb cols.celltype)
  e[, celltype_color := unname(cols.celltype[celltype])]
  setcolorder(e, c("edge_id","source","target","TF","target_gene","period","celltype","subclass",
                   "TF_direction","DEG_direction","relation_to_TF_direction",
                   "n_supporting_peaks","n_DAR_supporting_peaks","DAR_support_fraction","has_DAR_support",
                   "supporting_peak_list","DAR_peak_list","DAR_direction_list","celltype_color"))
  e
}
edges_dev <- build_edges("Development")
edges_ad  <- build_edges("Adult")
cat("Development edges:", nrow(edges_dev), "| Adult edges:", nrow(edges_ad), "\n")
## ---- 2. Node table ----
## node_id = gene symbol (全局唯一); node_type = TF / DEG / TF&DEG;
## subclass: DEG node 的 subclass (冲突时显式记录; TF 纯节点 NA; TF&DEG 取其 DEG subclass)
build_nodes <- function(ed, per_lab){
  sources <- unique(ed$source); targets <- unique(ed$target)
  allg <- union(sources, targets)
  nd <- data.table(node_id = allg)
  nd[, node_type := fifelse(node_id %in% sources & node_id %in% targets, "TF&DEG",
                      fifelse(node_id %in% sources, "TF", "DEG"))]
  ## DEG subclass (含 TF&DEG 的 DEG 侧): 从 edge 表取该 target 的 subclass 集合
  deg_sub <- ed[, .(subclasses = paste(sort(unique(subclass)), collapse = ";")), by = .(target_gene)]
  nd <- merge(nd, deg_sub, by.x = "node_id", by.y = "target_gene", all.x = TRUE)
  setnames(nd, "subclasses", "subclass")
  ## 冲突检查: 同一 DEG 多 subclass → 打印并记录 (不静默)
  conflicts <- nd[grepl(";", subclass, fixed = TRUE)]
  if (nrow(conflicts)){
    cat(sprintf("[subclass 冲突 - %s] %d 个 DEG 跨多个 subclass (跨 celltype 归属):\n", per_lab, nrow(conflicts)))
    print(conflicts[, .(node_id, subclass)])
    ## 处理: subclass 列保留组合串 (如 "AAM;DM"), 并额外加 subclass_primary + 冲突标记
    ## 注意: subclass_primary 仅为字符串顺序第一个, 无生物学优先级 —
    ## 后续 Cytoscape 着色建议: subclass ∈ {DM,DAL,AAM} 正常染色; subclass_conflict==TRUE 用中性色单独处理
    nd[, subclass := as.character(subclass)]
    nd[, subclass_primary := sub(";.*", "", subclass)]
    nd[, subclass_conflict := grepl(";", subclass, fixed = TRUE)]
    nd[node_type == "TF", `:=`(subclass = NA_character_, subclass_primary = NA_character_, subclass_conflict = FALSE)]
  } else {
    nd[, subclass := as.character(subclass)]
    nd[node_type == "TF", subclass := NA_character_]
    nd[, subclass_primary := subclass]
    nd[, subclass_conflict := FALSE]
  }
  ## degree: directed network (该 period 网络)
  g <- graph_from_data_frame(d = ed[, .(source, target)], directed = TRUE)
  deg_all <- degree(g); deg_in <- degree(g, mode = "in"); deg_out <- degree(g, mode = "out")
  nd[, degree := deg_all[node_id]]
  nd[, in_degree := deg_in[node_id]]
  nd[, out_degree := deg_out[node_id]]
  nd[, label := node_id]
  setcolorder(nd, c("node_id","label","node_type","subclass","degree","in_degree","out_degree"))
  nd <- nd[order(node_id)]
  list(nodes = nd, conflicts = conflicts)
}
nodes_dev_l <- build_nodes(edges_dev, "Development")
nodes_ad_l  <- build_nodes(edges_ad, "Adult")
cat("Development nodes:", nrow(nodes_dev_l$nodes), "| Adult nodes:", nrow(nodes_ad_l$nodes), "\n")
## ---- 3. 写出 6 个核心文件 (edges.tsv / nodes.tsv / graphml × 2 period) ----
write_set <- function(ed, nd, per_lab){
  fe <- file.path(DIR_TAB, sprintf("%s_plot_phase51_2_cytoscape_%s_edges.tsv", DATE_PREFIX, per_lab))
  fn <- file.path(DIR_TAB, sprintf("%s_plot_phase51_2_cytoscape_%s_nodes.tsv", DATE_PREFIX, per_lab))
  fg <- file.path(DIR_TAB, sprintf("%s_plot_phase51_2_cytoscape_%s_network.graphml", DATE_PREFIX, per_lab))
  fwrite(ed, fe, sep = "\t")
  fwrite(nd, fn, sep = "\t")
  ## GraphML: igraph directed; edge 属性全部保留 (igraph 会将 data.frame 列转为 edge attr)
  ## [2026-10-10 修复] graph_from_data_frame 默认把【前两列】当 from/to — ed 前两列是 edge_id/source!
  ## 临时副本 setcolorder 把 source/target 排最前; edges.tsv 列顺序不变.
  edge_graph <- copy(ed)
  setcolorder(edge_graph, c("source", "target", setdiff(names(edge_graph), c("source", "target"))))
  g <- graph_from_data_frame(d = edge_graph, directed = TRUE, vertices = nd)
  write_graph(g, fg, format = "graphml")
  cat("written:", basename(fe), "|", basename(fn), "|", basename(fg), "\n")
  c(fe, fn, fg)
}
out_dev <- write_set(edges_dev, nodes_dev_l$nodes, "Development")
out_ad  <- write_set(edges_ad,  nodes_ad_l$nodes,  "Adult")
## ---- 4. QC ----
qc_period <- function(ed, nd, per_lab){
  cat(sprintf("\n===== [%s QC] =====\n", per_lab))
  cat("unique nodes:", uniqueN(nd$node_id), "\n")
  cat("unique TF nodes (纯 TF):", sum(nd$node_type == "TF"), "\n")
  cat("unique DEG nodes (纯 DEG):", sum(nd$node_type == "DEG"), "\n")
  cat("TF&DEG nodes:", sum(nd$node_type == "TF&DEG"), "\n")
  cat("total edges:", nrow(ed), "\n")
  cat("DAR-supported edges:", sum(ed$has_DAR_support), "\n")
  cat("non-DAR-supported edges:", sum(!ed$has_DAR_support), "\n")
}
qc_period(edges_dev, nodes_dev_l$nodes, "Development")
qc_period(edges_ad,  nodes_ad_l$nodes,  "Adult")

## 全局一致性检查
cat("\n===== [一致性 QC] =====\n")
cat("source/target 均在 node table (dev):",
    all(edges_dev$source %in% nodes_dev_l$nodes$node_id) && all(edges_dev$target %in% nodes_dev_l$nodes$node_id), "\n")
cat("source/target 均在 node table (adult):",
    all(edges_ad$source %in% nodes_ad_l$nodes$node_id) && all(edges_ad$target %in% nodes_ad_l$nodes$node_id), "\n")
cat("node_id 无重复 (dev):", !any(duplicated(nodes_dev_l$nodes$node_id)),
    "| (adult):", !any(duplicated(nodes_ad_l$nodes$node_id)), "\n")
cat("edge_id 无重复 (dev):", !any(duplicated(edges_dev$edge_id)),
    "| (adult):", !any(duplicated(edges_ad$edge_id)), "\n")
cat("n_DAR <= n_support (dev):", all(edges_dev$n_DAR_supporting_peaks <= edges_dev$n_supporting_peaks),
    "| (adult):", all(edges_ad$n_DAR_supporting_peaks <= edges_ad$n_supporting_peaks), "\n")
cat("fraction ∈ [0,1] (dev):", all(edges_dev$DAR_support_fraction >= 0 & edges_dev$DAR_support_fraction <= 1),
    "| (adult):", all(edges_ad$DAR_support_fraction >= 0 & edges_ad$DAR_support_fraction <= 1), "\n")
cat("has_DAR == (n_DAR>=1) (dev):", all(edges_dev$has_DAR_support == (edges_dev$n_DAR_supporting_peaks >= 1)),
    "| (adult):", all(edges_ad$has_DAR_support == (edges_ad$n_DAR_supporting_peaks >= 1)), "\n")
cat("period 不混合 (dev 只 Development):", all(edges_dev$period == "Development"),
    "| (adult 只 Adult):", all(edges_ad$period == "Adult"), "\n")
cat("无缺失 TF/target_gene/celltype/subclass (dev):",
    !any(is.na(edges_dev$TF) | is.na(edges_dev$target_gene) | is.na(edges_dev$celltype) | is.na(edges_dev$subclass)), "\n")
cat("无缺失 TF/target_gene/celltype/subclass (adult):",
    !any(is.na(edges_ad$TF) | is.na(edges_ad$target_gene) | is.na(edges_ad$celltype) | is.na(edges_ad$subclass)), "\n")
cat("\nsubclass 冲突记录: dev", nrow(nodes_dev_l$conflicts), "| adult", nrow(nodes_ad_l$conflicts),
    "(LRP2/PDZRN3 跨 celltype 多 subclass — 保留组合串 + subclass_conflict 标记, 不静默)\n")

cat("\n===== [输出文件] =====\n")
for (fp in c(out_dev, out_ad)) cat(fp, "| size:", file.size(fp), "bytes\n")
