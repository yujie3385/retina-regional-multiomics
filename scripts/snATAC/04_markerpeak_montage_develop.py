# 最终整理脚本 04 · 胚胎 marker peak 拼图 (Python 矢量拼合, 仅最终收录输出) (2026-10-27)
# 来源: 00use_scripts/12_build_atac_montage_develop_formal_0909.py
# 产出: marker_peak_track/output/pdf/2026-09-09_develop_markerpeak_formal_montage_editable_vector.pdf (AI 可编辑矢量)
#       + 各列素材来源清单 manifest CSV
# 整理说明: 去除 150dpi/600dpi/black_grid 等未收录变体与 readme 写出;
#           ROOT_DIR 由脚本自身位置改为原项目绝对路径 (移至 00final_scripts 后仍可运行)
from __future__ import annotations  # 启用现代类型注解，便于长期维护脚本。

import csv  # 导入 CSV 工具，用于读取列顺序表和保存素材来源清单。
import hashlib  # 导入哈希工具，用于记录每张输入 PDF 的内容指纹。
import re  # 导入正则表达式工具，用于解析文件名和 PDF 内的信号范围。
from dataclasses import dataclass  # 导入数据类，用于明确记录每一列的来源信息。
from pathlib import Path  # 导入路径工具，避免把本机绝对路径写死在脚本中。

import pymupdf as fitz  # 导入 PyMuPDF，用于读取单页 PDF 并按指定 dpi 渲染。
from PIL import Image, ImageDraw, ImageFont, ImageOps  # 导入 Pillow，用于裁剪、缩放、标注、拼接和导出。


ROOT_DIR      = Path("/project/work/2026-05-08_ATAC_signac_all_sample/Plot/marker_peak_track")  # 原项目根目录 (绝对路径, 保证脚本移位后仍可运行)。
INPUT_DIR     = ROOT_DIR / "columns_develop_formal_v2"  # 设置胚胎正式版单列 PDF 的素材目录。
CSV_PATH      = Path("/project/work/2026-05-08_ATAC_signac_all_sample/Data/2026-09-09_develop_markerpeak_formal.csv")  # 读取列顺序与 peak 坐标的来源表。
MANIFEST_PATH = ROOT_DIR / "2026-09-09_develop_formal_montage_manifest.csv"  # 保存各列素材来源与质控清单。
OUTPUT_STEM   = "2026-09-09_develop_markerpeak_formal_montage"  # 拼图输出文件的统一名称前缀。
OUTPUT_PDF_DIR = ROOT_DIR / "output/pdf"  # 适合 Adobe Illustrator 微调的矢量 PDF 输出目录。
FONT_PATH     = Path("/usr/share/fonts/dejavu/DejaVuSans.ttf")  # 本机可用字体（DejaVu Sans 替代 Arial）。

KNOWN_CELLTYPES = {'Rod','Cone','AC','BC','HC','MG','RGC','RPC','Micro','Astro'}  # 用于从 PDF 自动检测行顺序。

## ---- 实测裁剪坐标（来自 126x1440pt 的 CoveragePlot PDF, 72dpi 量测） ----
SOURCE_PAGE_WIDTH_PT  = 126   # 单列 PDF 页面宽度（pt）。
SOURCE_PAGE_HEIGHT_PT = 1440  # 单列 PDF 页面高度（pt）。
SRC_COVERAGE_Y0   = 4    # coverage 主体上边界（pt）。
SRC_COVERAGE_Y1   = 979  # coverage 主体下边界（pt），完整保留所有 celltype 行。
SRC_GENE_BODY_Y0  = 1216  # 蓝色 gene-body 轨道上边界（pt）。
SRC_GENE_BODY_Y1  = 1236  # 蓝色 gene-body 轨道下边界（pt）。
SRC_LABEL_X0      = 38   # celltype 标签区左边界（pt），排除源 PDF 自身的纵轴标题。
SRC_SIGNAL_X0     = 68   # signal 区左边界（pt）。
SRC_SIGNAL_X1     = 114  # signal 区右边界（pt）。

## ---- 矢量版专用: 每个网格格子(行×列)的正方形边长（pt），可调 ----
CELL_SIDE_PT = 46  # 矢量版格子边长(pt): 缩小整图使其更紧凑; 字体绝对大小保持可读, 不随格子等比缩小。


@dataclass(frozen=True)
class PanelSource:  # 定义单个拼图列的来源信息。
    order: int  # 保存最终顺序编号。
    celltype: str  # 保存当前列对应的 celltype。
    marker_rank: int  # 保存当前列在该 celltype 内的 marker 编号（1/2/3）。
    gene: str  # 保存当前列对应的 gene 标签。
    genomic_region: str  # 保存 peak 坐标（chr-start-end）。
    signal_max: int  # 保存当前 PDF 自身 normalized signal 的纵轴上限。
    pdf_path: Path  # 保存当前列的 PDF 素材路径。


def sha256_file(path: Path) -> str:  # 计算输入 PDF 的 SHA-256，便于确认素材内容是否变化。
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:  # 加载论文图常用的无衬线字体。
    for font_path in [str(FONT_PATH), "/opt/conda/fonts/DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(font_path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def safe_name(x: str) -> str:  # 将 celltype/gene/peak 转成合法文件名（与 R 端 safe() 函数一致）。
    return re.sub(r"[^A-Za-z0-9_.-]", "-", x)


def detect_row_order(pdf_path: Path) -> list[str]:  # 从第一张 PDF 自动检测 CoveragePlot 实际行顺序（top→bottom）。
    with fitz.open(pdf_path) as doc:
        page = doc[0]
        labels = []
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    text = span["text"].strip()
                    if text in KNOWN_CELLTYPES and span["bbox"][0] < SRC_SIGNAL_X0:
                        labels.append((span["bbox"][1], text))  # (y_top, celltype)
        labels.sort(key=lambda x: x[0])  # 按纵坐标从上到下排序。
        order = []
        for _, ct in labels:
            if ct not in order:  # 去重（同一 celltype 可能有多个文字碎片）。
                order.append(ct)
        return order


def discover_panel_sources(row_order: list[str]) -> list[PanelSource]:  # 按 CSV 列顺序发现并严格验证所有单列 PDF。
    sources: list[PanelSource] = []
    with CSV_PATH.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for order, row in enumerate(reader, start=1):
            ct    = row["celltype"]
            rank  = int(row["marker_rank"])
            gene  = row["gene"]
            peak  = row["peak"]
            pdf_path = INPUT_DIR / f"{safe_name(ct)}_{safe_name(gene)}_{safe_name(peak)}.pdf"
            if not pdf_path.exists():
                raise FileNotFoundError(f"Missing column PDF: {pdf_path}")
            with fitz.open(pdf_path) as doc:
                if doc.page_count != 1:
                    raise ValueError(f"Expected one-page PDF: {pdf_path}")
                page = doc[0]
                rect = page.rect
                if abs(rect.width - SOURCE_PAGE_WIDTH_PT) > 0.5 or abs(rect.height - SOURCE_PAGE_HEIGHT_PT) > 0.5:
                    raise ValueError(f"Unexpected page size {rect.width:.1f}x{rect.height:.1f} for {pdf_path.name}")
                page_text = " ".join(page.get_text().split())
            m = re.search(r"range\s*0\s*[−–-]\s*([0-9]+)", page_text)
            if m is None:
                raise ValueError(f"Cannot parse signal range from {pdf_path.name}")
            sources.append(PanelSource(order, ct, rank, gene, peak, int(m.group(1)), pdf_path))
    ## 检查所有行中的 celltype 都出现在检测到的行顺序里
    missing = {s.celltype for s in sources} - set(row_order)
    if missing:
        print(f"WARNING: celltypes in CSV but not detected in row order: {missing}")
    return sources


def check_marker_counts(sources: list[PanelSource]) -> None:  # Python montage 前核查每个 celltype 的列数，不足 3 列时打印警告。
    from collections import Counter
    counts = Counter(s.celltype for s in sources)
    print("\n=== MARKER COUNT CHECK (before montage) ===")
    all_ok = True
    for ct, n in sorted(counts.items()):
        flag = "" if n >= 3 else "  << WARNING: only {} markers".format(n)
        print(f"  {ct:<8} : {n} markers{flag}")
        if n < 3: all_ok = False
    if all_ok:
        print("All celltypes have >= 3 markers. Proceeding.\n")
    else:
        print("Some celltypes have fewer than 3 markers. See above. Proceeding anyway.\n")


def write_manifest(sources: list[PanelSource]) -> None:  # 写出可回溯的各列素材与质控清单。
    fieldnames = ["order","celltype","marker_rank","gene","genomic_region","signal_min","signal_max",
                  "page_width_pt","page_height_pt","source_pdf_relative","source_pdf_sha256"]
    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for s in sources:
            writer.writerow({"order":s.order,"celltype":s.celltype,"marker_rank":s.marker_rank,
                             "gene":s.gene,"genomic_region":s.genomic_region,"signal_min":0,
                             "signal_max":s.signal_max,"page_width_pt":SOURCE_PAGE_WIDTH_PT,
                             "page_height_pt":SOURCE_PAGE_HEIGHT_PT,
                             "source_pdf_relative":s.pdf_path.relative_to(ROOT_DIR).as_posix(),
                             "source_pdf_sha256":sha256_file(s.pdf_path)})


def render_pdf_to_image(pdf_path: Path, dpi: int) -> Image.Image:  # 将单页 PDF 按指定 dpi 渲染为 RGB 图片。
    with fitz.open(pdf_path) as doc:
        pix = doc.load_page(0).get_pixmap(matrix=fitz.Matrix(dpi/72, dpi/72), alpha=False)
    return Image.frombytes("RGB", [pix.width, pix.height], pix.samples)


def resize_panel(image: Image.Image, width: int, height: int) -> Image.Image:  # 标准化到统一单元格尺寸。
    return image.resize((width, height), resample=Image.Resampling.LANCZOS)


def paste_horizontal(parts: list[Image.Image]) -> Image.Image:  # 横向无间距拼接第一列的标签区与 signal 区。
    w = sum(p.width for p in parts); h = max(p.height for p in parts)
    canvas = Image.new("RGB", (w, h), "white"); x = 0
    for p in parts: canvas.paste(p, (x, 0)); x += p.width
    return canvas


def draw_rotated_gene_label(panel_width: int, gene: str, center_x: int, panel_height: int, font_size: int) -> Image.Image:  # 生成底部竖排 gene 标签面板。
    panel = Image.new("RGB", (panel_width, panel_height), "white")
    fs = font_size
    while fs > max(8, round(14 * panel_height / 96)):
        font_c = load_font(fs)
        bb = font_c.getbbox(gene)
        if (bb[2]-bb[0]+8) <= panel_height: break
        fs -= 1
    font = load_font(fs); bb = font.getbbox(gene)
    tw, th = bb[2]-bb[0], bb[3]-bb[1]
    ti = Image.new("RGBA", (tw+8, th+8), (255,255,255,0))
    ImageDraw.Draw(ti).text((4-bb[0], 4-bb[1]), gene, fill="black", font=font)
    rot = ti.rotate(90, expand=True)
    px = max(0, min(panel_width-rot.width, int(center_x-rot.width/2)))
    py = max(0, panel_height-rot.height-1)
    panel.paste(rot.convert("RGB"), (px, py), rot)
    return panel


def stack_vertical(parts: list[Image.Image]) -> Image.Image:  # 纵向拼接 coverage、gene 名和 gene-body 轨道。
    w = max(p.width for p in parts); h = sum(p.height for p in parts)
    canvas = Image.new("RGB", (w, h), "white"); y = 0
    for p in parts: canvas.paste(p, (0, y)); y += p.height
    return canvas


def paste_columns(columns: list[Image.Image]) -> Image.Image:  # 横向无间距拼接全部 gene 列。
    total_w = sum(c.width for c in columns); max_h = max(c.height for c in columns)
    canvas = Image.new("RGB", (total_w, max_h), "white"); x = 0
    for c in columns: canvas.paste(c, (x, 0)); x += c.width
    return canvas


def add_left_y_label_panel(montage: Image.Image, coverage_height: int, scale: float) -> Image.Image:  # 添加左侧 normalized signal 标题。
    pw = round(66*scale); panel = Image.new("RGB", (pw, montage.height), "white")
    font = load_font(max(8, round(16*scale))); text = "Normalized signal\n(column-specific scale)"
    tl = Image.new("RGBA", (round(280*scale), round(58*scale)), (255,255,255,0))
    td = ImageDraw.Draw(tl); bb = td.multiline_textbbox((0,0), text, font=font, spacing=round(3*scale), align="center")
    td.multiline_text((round(4*scale)-bb[0], round(4*scale)-bb[1]), text, fill="black", font=font, spacing=round(3*scale), align="center")
    rot = tl.rotate(90, expand=True); py = max(0, int((coverage_height-rot.height)/2))
    panel.paste(rot.convert("RGB"), (round(4*scale), py), rot)
    canvas = Image.new("RGB", (panel.width+montage.width, montage.height), "white")
    canvas.paste(panel, (0,0)); canvas.paste(montage, (panel.width,0))
    return canvas


def draw_black_grid(montage: Image.Image, column_widths: list[int], coverage_height: int, label_separator_x: int, n_rows: int, scale: float) -> Image.Image:  # 为主体及 celltype 标签边界添加黑色网格。
    gridded = montage.copy(); draw = ImageDraw.Draw(gridded)
    lw = max(2, round(1.8*scale)); row_h = coverage_height / n_rows
    for i in range(n_rows+1):
        y = round(i*row_h); draw.line((0, y, gridded.width, y), fill="black", width=lw)
    x = 0; draw.line((x, 0, x, gridded.height), fill="black", width=lw)
    draw.line((label_separator_x, 0, label_separator_x, gridded.height), fill="black", width=lw)
    for cw in column_widths:
        x += cw; draw.line((x, 0, x, gridded.height), fill="black", width=lw)
    draw.rectangle((0, 0, gridded.width-1, gridded.height-1), outline="black", width=lw)
    return gridded


def build_montage(sources: list[PanelSource], row_order: list[str], dpi: int, suffix: str) -> tuple[Path, Path]:  # 按指定 dpi 生成无网格与黑色网格两版拼图。
    scale = dpi / 72
    cov_y0 = round(SRC_COVERAGE_Y0*scale); cov_y1 = round(SRC_COVERAGE_Y1*scale)
    gb_y0  = round(SRC_GENE_BODY_Y0*scale); gb_y1  = round(SRC_GENE_BODY_Y1*scale)
    lbl_x0 = round(SRC_LABEL_X0*scale); sig_x0 = round(SRC_SIGNAL_X0*scale); sig_x1 = round(SRC_SIGNAL_X1*scale)
    cov_h  = cov_y1 - cov_y0; lbl_w = sig_x0 - lbl_x0; sig_w = sig_x1 - sig_x0; gb_h = gb_y1 - gb_y0
    gene_lbl_h = round(96*scale); gene_lbl_fs = max(8, round(22*scale))
    first_col_w = lbl_w + sig_w; columns = []; col_widths = []
    for idx, source in enumerate(sources):
        rendered = render_pdf_to_image(source.pdf_path, dpi=dpi)
        raw_sig  = rendered.crop((sig_x0, cov_y0, sig_x1, cov_y1))
        pk_sig   = resize_panel(raw_sig, sig_w, cov_h)
        raw_gb   = rendered.crop((sig_x0, gb_y0, sig_x1, gb_y1))
        gb_sig   = resize_panel(raw_gb, sig_w, gb_h)
        if idx == 0:
            raw_lbl  = rendered.crop((lbl_x0, cov_y0, sig_x0, cov_y1))
            lbl_pan  = resize_panel(raw_lbl, lbl_w, cov_h)
            pk_part  = paste_horizontal([lbl_pan, pk_sig])
            gb_part  = Image.new("RGB", (first_col_w, gb_h), "white")
            gb_part.paste(gb_sig, (lbl_w, 0))
            gene_cx  = lbl_w + sig_w//2
            raw_lbl.close(); lbl_pan.close()
        else:
            pk_part = pk_sig.copy(); gb_part = gb_sig.copy(); gene_cx = sig_w//2
        gene_lbl = draw_rotated_gene_label(pk_part.width, source.gene, gene_cx, gene_lbl_h, gene_lbl_fs)
        col = stack_vertical([pk_part, gene_lbl, gb_part])
        columns.append(col); col_widths.append(col.width)
        rendered.close(); raw_sig.close(); pk_sig.close(); raw_gb.close(); gb_sig.close()
        pk_part.close(); gene_lbl.close(); gb_part.close()
    core      = paste_columns(columns)
    core_grid = draw_black_grid(core, col_widths, cov_h, lbl_w, len(row_order), scale)
    montage   = add_left_y_label_panel(core, cov_h, scale)
    montage_g = add_left_y_label_panel(core_grid, cov_h, scale)
    border = round(18*scale)
    montage   = ImageOps.expand(montage,   border=border, fill="white")
    montage_g = ImageOps.expand(montage_g, border=border, fill="white")
    out   = ROOT_DIR / f"{OUTPUT_STEM}_{suffix}.pdf"
    out_g = ROOT_DIR / f"{OUTPUT_STEM}_black_grid_{suffix}.pdf"
    montage.save(out,   format="PDF", resolution=dpi)
    montage_g.save(out_g, format="PDF", resolution=dpi)
    for img in [montage, montage_g, core, core_grid] + columns: img.close()
    return out, out_g


def fit_pdf_font_size(font: fitz.Font, text: str, target: float, max_len: float) -> float:  # 为矢量 PDF 中的长文字计算可完整显示的字号。
    fs = target; fs_min = target * 14 / 22
    while fs > fs_min and font.text_length(text, fontsize=fs) > max_len: fs -= 0.25
    return max(fs_min, fs)


def build_editable_pdf(sources: list[PanelSource], row_order: list[str], with_grid: bool) -> Path:  # 直接使用原始 PDF 矢量对象构建适合 Adobe Illustrator 的单页拼图。
    OUTPUT_PDF_DIR.mkdir(parents=True, exist_ok=True)
    if not FONT_PATH.is_file(): raise FileNotFoundError(f"Font not found: {FONT_PATH}")
    lbl_w = CELL_SIDE_PT; sig_w = CELL_SIDE_PT  # 正方形: 标签列与每个信号列同宽, 均 = CELL_SIDE_PT
    cov_h = CELL_SIDE_PT * len(row_order)  # 正方形: coverage 高度 = 行数 × CELL_SIDE_PT, 使每行高 = 列宽
    gb_h = SRC_GENE_BODY_Y1 - SRC_GENE_BODY_Y0  # 底部蓝色 gene-body 轨道高度(注释条, 非网格格子)
    gene_lbl_h = 88; border = 14; title_w = 68; n_col = len(sources)  # 与缩小格子(46pt)匹配: 更紧凑, 字体绝对大小仍可读
    core_x0 = border + title_w; core_y0 = border
    cov_y1 = core_y0 + cov_h; gene_y1 = cov_y1 + gene_lbl_h; core_y1 = gene_y1 + gb_h
    page_w = core_x0 + lbl_w + n_col*sig_w + border
    page_h = core_y1 + border
    name = f"{OUTPUT_STEM}{'_black_grid' if with_grid else ''}_editable_vector.pdf"
    out_path = OUTPUT_PDF_DIR / name
    doc = fitz.open(); page = doc.new_page(width=page_w, height=page_h)
    fname = "MontageSans"
    page.insert_font(fontname=fname, fontfile=str(FONT_PATH))
    vec_font = fitz.Font(fontfile=str(FONT_PATH))
    sig_clip = fitz.Rect(SRC_SIGNAL_X0, SRC_COVERAGE_Y0, SRC_SIGNAL_X1, SRC_COVERAGE_Y1)
    gb_clip  = fitz.Rect(SRC_SIGNAL_X0, SRC_GENE_BODY_Y0, SRC_SIGNAL_X1, SRC_GENE_BODY_Y1)
    for idx, s in enumerate(sources):
        sx0 = core_x0 + lbl_w + idx*sig_w; sx1 = sx0 + sig_w
        with fitz.open(s.pdf_path) as src_doc:
            ## 去除原始文字，保留矢量图形
            tmp = fitz.open(); tmp.insert_pdf(src_doc)
            tp = tmp[0]
            tp.add_redact_annot(tp.rect, fill=False, cross_out=False)
            tp.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
            clean = fitz.open(stream=tmp.tobytes(garbage=4, clean=True, deflate=True), filetype="pdf")
            tmp.close()
            page.show_pdf_page(fitz.Rect(sx0, core_y0, sx1, cov_y1), clean, 0, clip=sig_clip, keep_proportion=False, overlay=True)
            page.show_pdf_page(fitz.Rect(sx0, gene_y1, sx1, core_y1), clean, 0, clip=gb_clip, keep_proportion=False, overlay=True)
            clean.close()
        gfs = fit_pdf_font_size(vec_font, s.gene, 19, gene_lbl_h-4)  # 基因名字号加大(和谐, 不显小)
        gl = vec_font.text_length(s.gene, fontsize=gfs)
        cx = (sx0+sx1)/2; cy = (cov_y1+gene_y1)/2
        page.insert_text((cx+gfs*0.34, cy+gl/2), s.gene, fontname=fname, fontsize=gfs, color=(0,0,0), rotate=90, overlay=True)
    ## celltype 行标签（重新绘制，用 ROW_CELLTYPE_ORDER）
    row_h = cov_h / len(row_order); ct_cx = core_x0 + lbl_w/2
    for i, ct in enumerate(row_order):
        ct_cy = core_y0 + (i+0.5)*row_h
        fs = fit_pdf_font_size(vec_font, ct, 17, lbl_w-4)  # celltype 行标签字号加大
        tl = vec_font.text_length(ct, fontsize=fs)
        page.insert_text((ct_cx-tl/2, ct_cy+fs*0.34), ct, fontname=fname, fontsize=fs, color=(0,0,0), overlay=True)
    ## 左侧纵轴标题
    tcy = core_y0 + cov_h/2
    for line, x_off in [("Normalized signal", border+22), ("(column-specific scale)", border+46)]:
        ll = vec_font.text_length(line, fontsize=15)
        page.insert_text((x_off+15*0.34, tcy+ll/2), line, fontname=fname, fontsize=15, color=(0,0,0), rotate=90, overlay=True)  # 纵轴标题字号加大
    if with_grid:
        lw = 0.8; total_w = lbl_w + n_col*sig_w
        for i in range(len(row_order)+1):
            y = core_y0 + i*row_h
            page.draw_line((core_x0, y), (core_x0+total_w, y), color=(0,0,0), width=lw, overlay=True)
        for x in [core_x0, core_x0+lbl_w] + [core_x0+lbl_w+j*sig_w for j in range(1, n_col+1)]:
            page.draw_line((x, core_y0), (x, core_y1), color=(0,0,0), width=lw, overlay=True)
        page.draw_rect(fitz.Rect(core_x0, core_y0, core_x0+total_w, core_y1), color=(0,0,0), width=lw, overlay=True)
    doc.subset_fonts()
    doc.save(out_path, garbage=4, deflate=True, clean=True)
    doc.close()
    return out_path


def main() -> None:  # 定义正式版拼图主流程。
    ## 1. 自动检测行顺序（从第一张 PDF）
    first_pdf = next(iter(sorted(INPUT_DIR.glob("*.pdf"))))
    row_order = detect_row_order(first_pdf)
    print(f"Detected row order ({len(row_order)} celltypes): {' > '.join(row_order)}")
    ## 2. 发现并验证所有列 PDF
    sources = discover_panel_sources(row_order)
    ## 3. montage 前核查每类 marker 数，不足 3 则打印警告
    check_marker_counts(sources)
    ## 4. 写 manifest
    write_manifest(sources)
    ## 5. 生成最终收录的矢量 PDF (无网格 editable_vector; 150dpi/600dpi/black_grid 变体为探索版, 已去除)
    vector_pdf = build_editable_pdf(sources, row_order, with_grid=False)
    print(f"validated_panels : {len(sources)}")
    print(f"signal_max_range : {min(s.signal_max for s in sources)}-{max(s.signal_max for s in sources)}")
    print(f"manifest         : {MANIFEST_PATH}")
    print(f"vector_pdf       : {vector_pdf}")
    print("columns:")
    for s in sources:
        print(f"  {s.order:03d}  {s.celltype:<7}  rank{s.marker_rank}  {s.gene:<10}  0-{s.signal_max}")


if __name__ == "__main__":
    main()
