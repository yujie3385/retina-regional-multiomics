# 最终整理脚本 21 · Fig4 Panel A 示例基因 track 图 V2 (覆盖度+DAR+Motif+links) · 胚胎+成年 (2026-10-27)
# 来源: 00use_scripts/25_2026-10-20_track_APOE_MG_link_motif_V1.ipynb (V2 绘图逻辑)
# 产出: /project/work/2026-10-20_track/Plot/2026-10-20_Fig4_A_{GENE}_{CT}_{dev,adult}_stage_tissue_link_motif_V2.pdf (8 张)
#       基因: ABCB11_MG / PDE10A_BC / TFRC_MG / TTC39B_MG × 胚胎(8x7.5)与成年(8x7)
# 轨道: 覆盖度 (stage×region 分组) + 基因 + DAR(central 紫/peripheral 橙) + Motif(堆叠 box) + ArchR links
# 整理说明: 去除 V1 输出与未收录基因; gene_motif 表保留完整登记但仅绘制最终收录的 4 个基因
suppressMessages({
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
})

date     <- "2026-10-20"
ver      <- "V2"
data_dir <- "/project/work/2026-09-03_CRE_new/Data/"
cand_dir <- "/project/work/2026-08-08_candidate_CRE/Data/"
plot_dir <- "/project/work/2026-10-20_track/Plot/"
frag_dir <- "/project/work/00data/ChenLab_human_atac/Fragments/"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

## 输入
support  <- read.csv(paste0(data_dir, "2026-09-04_Fig4_0_gene_support_master.csv"))
master   <- read.csv(paste0(cand_dir, "2026-08-12_DEG_DAR_ArchR_link_peak_level.csv"))
presence <- read.csv(paste0(cand_dir, "2026-08-12_3_candidate_CRE_motif_presence_long.csv"))

## 目标基因 × motif 表 (celltype + gene + motifs; 完整探索登记)
gene_motif <- data.frame(
  celltype = c("MG","MG","MG","MG","MG","MG","MG","MG","MG",
               "BC","BC","BC","BC","Cone","Cone","Cone","BC","AC"),
  gene     = c("APOE","ABCB11","TTC39B","TFRC","DIO2","PITPNC1","COTL1","ZBTB21","ADAMTS7",
               "NSG1","BHLHE22","PDE10A","PODXL","NOL4","ZIC1","PLXDC1","ZIC1","ZIC1"),
  motifs   = c("MEF2C","ONECUT1;MEF2C","ONECUT1","MAFK;NFE2L1","ONECUT1",
               "MAFK;NFE2L1","MAFK;NFE2L1","MAFK","NFE2L1",
               "EGR1;NEUROD1","NEUROD1;HAND2","KLF15;EGR1;EGR4","KLF9;KLF4;TCF3",
               "MEF2C;ASCL1(var.2);SP3;ARGFX","KLF2;E2F6;SP3","OTX2;KLF14;SP4",
               "KLF2;KLF9;SP3","KLF2;KLF9;SP3"),
  stringsAsFactors = FALSE)
subclass_group <- c("Developmental aging-lost" = "DAL",
                    "Developmental maintained" = "DM",
                    "Adult acquired and maintained" = "AAM")
gene_motif <- gene_motif %>%
  left_join(support %>% mutate(group = unname(subclass_group[subclass])) %>%
              select(celltype, gene, origin_class, group),
            by = c("celltype", "gene"))
stopifnot(!any(is.na(gene_motif$origin_class)),
          !anyDuplicated(gene_motif[c("celltype", "gene")]))
## 最终收录图集仅含 4 个基因 (ABCB11_MG / PDE10A_BC / TFRC_MG / TTC39B_MG); 其余为探索性示例, 不绘制
gene_motif <- gene_motif %>%
  filter((celltype == "MG" & gene %in% c("ABCB11","TTC39B","TFRC")) |
         (celltype == "BC" & gene == "PDE10A"))

## 配色 (与 2026-09-14 参考脚本一致)
cols.stage.full <- c('#59A14F','#EDC948','#B07AA1','#4E79A7','#E15759')
names(cols.stage.full) <- c('Early','Mid','Late','Young','Old')
cols.stage      <- cols.stage.full[c("Early","Mid","Late")]                       ## dev: Early/Mid/Late
cols.adultstage <- setNames(cols.stage.full[c("Young","Old")], c("young","old"))  ## adult: young/old
lighten_col <- function(col, f = 0.45) {
  m <- grDevices::col2rgb(col)
  grDevices::rgb(t(m + (255 - m) * f), maxColorValue = 255)
}
region_cols <- function(base) {
  v <- c(base, sapply(base, lighten_col))
  names(v) <- c(paste0(names(base), "_central"), paste0(names(base), "_peripheral"))
  v
}
cols.dar      <- c("central_open"="#6A3D9A", "peripheral_open"="#E69F00")
col.highlight <- "#FDE0A3"     ## highlight 背景色 (与参考图一致)
blank_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                 axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())

## ---- 辅助函数 ----
get_support_peaks <- function(ct, g) {
  r <- support %>% filter(celltype == ct, gene == g)
  unique(unlist(strsplit(c(r$embryonic_support_peaks, r$adult_support_peaks), ";")))
}
get_window <- function(obj, ct, g, pad = 15000) {
  gr <- StringToGRanges(get_support_peaks(ct, g))
  ann <- Annotation(obj); gb <- ann[!is.na(ann$gene_name) & ann$gene_name == g]
  if (length(gb) > 0) gr <- c(gr, GRanges(seqnames(gb)[1], IRanges(min(start(gb)), max(end(gb)))))
  w <- range(gr); start(w) <- start(w) - pad; end(w) <- end(w) + pad
  GRangesToString(w)
}
## DAR 轨: dev 取任一 stage 的方向一致 regional DAR; adult 取 shared core DAR (与参考一致)
get_dar_ranges <- function(ct, g, period) {
  m <- master %>% filter(celltype == ct, gene == !!g)
  if (period == "dev") {
    d <- m %>% filter(Early_DAR_status %in% c("M","P") | Mid_DAR_status %in% c("M","P") | Late_DAR_status %in% c("M","P")) %>%
      mutate(lab = ifelse(Early_DAR_status %in% "M" | Mid_DAR_status %in% "M" | Late_DAR_status %in% "M",
                          "central_open", "peripheral_open"))
  } else {
    d <- m %>% filter(Adult_shared_core_DAR) %>%
      mutate(lab = ifelse(Adult_shared_direction == "M", "central_open", "peripheral_open"))
  }
  if (nrow(d) == 0) return(NULL)
  gr <- StringToGRanges(d$peak); gr$DAR <- d$lab; gr
}
## links: 该期 ArchR peak2gene link (peak 中点 <-> TSS, score = 相关系数)
get_links <- function(obj, ct, g, period) {
  lcol <- ifelse(period == "dev", "develop_ArchR_link", "adult_ArchR_link")
  ccol <- ifelse(period == "dev", "develop_ArchR_correlation", "adult_ArchR_correlation")
  m <- master %>% filter(celltype == ct, gene == !!g, .data[[lcol]])
  if (nrow(m) == 0) return(NULL)
  ann <- Annotation(obj); gb <- ann[!is.na(ann$gene_name) & ann$gene_name == g]
  if (length(gb) == 0) return(NULL)
  tss <- ifelse(as.character(strand(gb))[1] == "+", min(start(gb)), max(end(gb)))
  pk <- StringToGRanges(m$peak); mid <- (start(pk) + end(pk)) %/% 2
  GRanges(as.character(seqnames(pk)), IRanges(pmin(mid, tss), pmax(mid, tss)),
          score = m[[ccol]], gene = g, peak = m$peak)
}
## motif peaks: presence 表 (peak 级), 取窗口内目标 TF 的 peak×TF 明细
get_motif_df <- function(region, tfs) {
  sub <- presence %>% filter(TF %in% tfs, motif_present == TRUE)
  if (nrow(sub) == 0) return(NULL)
  gr <- StringToGRanges(sub$peak)
  keep <- overlapsAny(gr, StringToGRanges(region))
  sub <- sub[keep, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  sub
}
get_motif_union <- function(motif.df) {
  if (is.null(motif.df)) return(GRanges())
  unique(StringToGRanges(motif.df$peak))
}
## highlight: motif peak 背景高亮
get_highlight <- function(peaks.gr) {
  if (length(peaks.gr) == 0) return(NULL)
  GRanges(seqnames(peaks.gr), ranges(peaks.gr), color = col.highlight)
}
## motif 轨 (样式参考 2026-08-15_6): 堆叠灰底 box + 斜体 motif 名标注
make_motif_track <- function(motif.df, region, title = "Motif") {
  win <- StringToGRanges(region)
  if (is.null(motif.df) || !nrow(motif.df)) {
    return(ggplot() +
             coord_cartesian(xlim = c(start(win), end(win)), ylim = c(0, 1)) +
             annotate("text", x = (start(win) + end(win)) / 2, y = 0.5,
                      label = "None", size = 3) +
             labs(y = title, x = NULL) + Signac::theme_browser())
  }
  df <- motif.df %>%
    group_by(peak) %>%
    summarise(motif_label = paste(sort(unique(TF)), collapse = " | "), .groups = "drop")
  gr <- StringToGRanges(df$peak)
  df$start <- start(gr); df$end <- end(gr); df$midpoint <- (df$start + df$end) / 2
  df <- df %>% arrange(start) %>% mutate(y = seq_len(n()))
  ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = y - 0.22, ymax = y + 0.22),
              fill = "grey65", colour = "grey30", linewidth = 0.25) +
    geom_text(aes(x = midpoint, y = y + 0.42, label = motif_label),
              size = 3, fontface = "italic", check_overlap = FALSE) +
    coord_cartesian(xlim = c(start(win), end(win)), clip = "off") +
    scale_y_continuous(breaks = NULL) +
    labs(y = title, x = NULL) +
    Signac::theme_browser() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(3, 8, 3, 5.5))
}

## ---- track 绘制: 覆盖度 + 基因 + DAR(紫/橙) + Motif(堆叠 box) + links ----
## 该期无 links 时省略 links 轨 (motif 轨置于最底并保留 x 轴)
plot_track_v2 <- function(obj, ct, g, tfs, group_col, track_cols, period, title) {
  region <- get_window(obj, ct, g)
  lk <- get_links(obj, ct, g, period)
  Links(obj) <- if (!is.null(lk)) lk else GRanges()
  ann <- Annotation(obj)
  Annotation(obj) <- ann[!is.na(ann$gene_name) & ann$gene_name == g]   ## 基因轨只留目标基因
  motif.df <- get_motif_df(region, tfs)
  motif.gr <- get_motif_union(motif.df)
  dar.gr   <- get_dar_ranges(ct, g, period)
  p <- CoveragePlot(obj, region = region, group.by = group_col,
                    ranges = dar.gr, ranges.group.by = "DAR", ranges.title = "DAR",
                    region.highlight = get_highlight(motif.gr),
                    peaks = FALSE, links = !is.null(lk)) &
    scale_fill_manual(values = track_cols)
  if (length(p$patches$plots) == 1L && inherits(p$patches$plots[[1]], "patchwork"))
    p <- p$patches$plots[[1]]
  n.sub <- length(p$patches$plots) + 1L
  if (n.sub < 2L) stop("CoveragePlot 轨道数异常: ", n.sub)
  tracks  <- vapply(seq_len(n.sub), function(i) paste(p[[i]]$labels$y, collapse = " "), "")
  i.dar   <- which(tracks == "DAR")
  i.links <- which(tracks == "Links")
  if (length(i.dar)) p[[i.dar]] <- p[[i.dar]] + scale_color_manual(values = cols.dar)
  motif.p <- make_motif_track(motif.df, region)
  mh <- if (is.null(motif.df) || !nrow(motif.df)) 1 else max(1, 0.8 * nrow(motif.df))
  if (length(i.links)) {
    links.p <- p[[i.links]]                    ## links 轨 (带 x 轴)
    p[[i.links]] <- motif.p + blank_x          ## motif 占原 links 位, 隐藏 x 轴
    p <- p / links.p + plot_layout(heights = c(10, 3, if (length(i.dar)) 1, mh, 3))
  } else {
    p[[n.sub]] <- p[[n.sub]] + blank_x         ## 原最后一轨 (DAR/基因) 隐藏 x 轴
    p <- p / motif.p + plot_layout(heights = c(10, 3, if (length(i.dar)) 1, mh))  ## motif 最底, 保留 x 轴
  }
  p + patchwork::plot_annotation(title = title,
                                 theme = theme(plot.title = element_text(size = 11)))
}

## ---- 胚胎 (dev) ----
dev <- readRDS("/project/work/2026-05-08_ATAC_signac_all_sample/Data/2026-06-06_develop_all_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds")
DefaultAssay(dev) <- "ATAC"
frags <- Fragments(dev[["ATAC"]]); paths <- sapply(frags, function(f) f@path)
if (!all(file.exists(paths))) {
  for (i in which(!file.exists(paths))) frags[[i]] <- UpdatePath(frags[[i]], file.path(frag_dir, basename(paths[i])))
  Fragments(dev[["ATAC"]]) <- NULL; Fragments(dev[["ATAC"]]) <- frags
}
cat("dev fragments ok:", all(file.exists(sapply(Fragments(dev[["ATAC"]]), function(f) f@path))), "\n")

## dev: stage x region 分组 (central = macula/fovea); 按细胞类型子集在绘图 cell 中进行
dev$stage <- dplyr::case_when(
  dev$donor_age %in% c("70days","79days","87days","91days") ~ "Early",
  dev$donor_age %in% c("100days","103days","116days") ~ "Mid",
  dev$donor_age %in% c("137days","141days","142days","162days","165days") ~ "Late")
dev$region <- ifelse(dev$tissue %in% c("macula","fovea"), "central", "peripheral")
lv <- as.vector(rbind(paste0(names(cols.stage),"_central"), paste0(names(cols.stage),"_peripheral")))
dev$stage_region <- factor(paste(dev$stage, dev$region, sep = "_"), levels = lv)
table(dev$celltype)

## dev 图: 按细胞类型循环绘制
cols_stage_track <- region_cols(cols.stage)
k <- 0
for (ct in unique(gene_motif$celltype)) {
  sub_tab <- gene_motif %>% filter(celltype == ct)
  obj_ct <- subset(dev, celltype == ct)
  if (ncol(obj_ct) == 0) { cat("!!", ct, "在 dev 对象中不存在, 跳过\n"); next }
  cat(sprintf("== dev %s: %d cells, %d genes ==\n", ct, ncol(obj_ct), nrow(sub_tab)))
  for (i in seq_len(nrow(sub_tab))) {
    g   <- sub_tab$gene[i]
    tfs <- unlist(strsplit(sub_tab$motifs[i], ";"))
    gm  <- sub_tab[i, ]
    rg    <- get_window(obj_ct, ct, g)
    dar.v <- get_dar_ranges(ct, g, "dev")
    lk.v  <- get_links(obj_ct, ct, g, "dev")
    mp.v  <- get_motif_df(rg, tfs)
    k <- k + 1
    cat(sprintf("[dev %d/%d] %-8s (%-4s) window=%s | DAR=%d (M=%d,P=%d) | links=%d | motif peaks=%d\n",
                k, nrow(gene_motif), g, ct, rg,
                ifelse(is.null(dar.v), 0, length(dar.v)),
                ifelse(is.null(dar.v), 0, sum(dar.v$DAR == "central_open")),
                ifelse(is.null(dar.v), 0, sum(dar.v$DAR == "peripheral_open")),
                ifelse(is.null(lk.v), 0, length(lk.v)),
                ifelse(is.null(mp.v), 0, nrow(mp.v))))
    p <- plot_track_v2(obj_ct, ct, g, tfs, "stage_region", cols_stage_track, "dev",
                       paste0(g, "  (", gm$origin_class, ", ", gm$group, ")  embryonic ATAC"))
    ggsave(paste0(plot_dir, date, "_Fig4_A_", g, "_", ct, "_dev_stage_tissue_link_motif_", ver, ".pdf"),
           p, width = 8, height = 7.5)
  }
  rm(obj_ct); invisible(gc())
}
rm(dev); invisible(gc())

## ---- 成年 (adult), 使用全部成年数据 ----
adult <- readRDS("/project/work/2026-05-08_ATAC_signac_all_sample/Data/2026-06-06_adult_merged_harmony_integrated_add_ACTIVITY_anno_sub.rds")
DefaultAssay(adult) <- "ATAC"
frags <- Fragments(adult[["ATAC"]]); paths <- sapply(frags, function(f) f@path)
if (!all(file.exists(paths))) {
  for (i in which(!file.exists(paths))) frags[[i]] <- UpdatePath(frags[[i]], file.path(frag_dir, basename(paths[i])))
  Fragments(adult[["ATAC"]]) <- NULL; Fragments(adult[["ATAC"]]) <- frags
}
cat("adult fragments ok:", all(file.exists(sapply(Fragments(adult[["ATAC"]]), function(f) f@path))), "\n")

## adult: young/old x region 分组 (young<40y, old>=40y); 按细胞类型子集在绘图 cell 中进行
adult$donor_age <- gsub(" ", "", adult$age)
age_num <- suppressWarnings(as.numeric(gsub("years", "", adult$donor_age)))
adult$adultstage <- ifelse(age_num < 40, "young", "old")
adult$region <- ifelse(adult$tissue %in% c("macula","fovea"), "central", "peripheral")
lv <- as.vector(rbind(paste0(names(cols.adultstage),"_central"), paste0(names(cols.adultstage),"_peripheral")))
adult$stage_region <- factor(paste(adult$adultstage, adult$region, sep = "_"), levels = lv)
table(adult$celltype)

## adult 图: 按细胞类型循环绘制
cols_adustage_track <- region_cols(cols.adultstage)
k <- 0
for (ct in unique(gene_motif$celltype)) {
  sub_tab <- gene_motif %>% filter(celltype == ct)
  obj_ct <- subset(adult, celltype == ct)
  if (ncol(obj_ct) == 0) { cat("!!", ct, "在 adult 对象中不存在, 跳过\n"); next }
  cat(sprintf("== adult %s: %d cells, %d genes ==\n", ct, ncol(obj_ct), nrow(sub_tab)))
  for (i in seq_len(nrow(sub_tab))) {
    g   <- sub_tab$gene[i]
    tfs <- unlist(strsplit(sub_tab$motifs[i], ";"))
    gm  <- sub_tab[i, ]
    rg    <- get_window(obj_ct, ct, g)
    dar.v <- get_dar_ranges(ct, g, "adult")
    lk.v  <- get_links(obj_ct, ct, g, "adult")
    mp.v  <- get_motif_df(rg, tfs)
    k <- k + 1
    cat(sprintf("[adult %d/%d] %-8s (%-4s) window=%s | DAR=%d (M=%d,P=%d) | links=%d | motif peaks=%d\n",
                k, nrow(gene_motif), g, ct, rg,
                ifelse(is.null(dar.v), 0, length(dar.v)),
                ifelse(is.null(dar.v), 0, sum(dar.v$DAR == "central_open")),
                ifelse(is.null(dar.v), 0, sum(dar.v$DAR == "peripheral_open")),
                ifelse(is.null(lk.v), 0, length(lk.v)),
                ifelse(is.null(mp.v), 0, nrow(mp.v))))
    p <- plot_track_v2(obj_ct, ct, g, tfs, "stage_region", cols_adustage_track, "adult",
                       paste0(g, "  (", gm$origin_class, ", ", gm$group, ")  adult ATAC"))
    ggsave(paste0(plot_dir, date, "_Fig4_A_", g, "_", ct, "_adult_stage_tissue_link_motif_", ver, ".pdf"),
           p, width = 8, height = 7)
  }
  rm(obj_ct); invisible(gc())
}
rm(adult); invisible(gc())
cat("done: 8 Fig4_A track V2 PDFs\n")
