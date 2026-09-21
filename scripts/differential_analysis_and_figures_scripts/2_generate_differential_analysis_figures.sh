#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# FIGURES: CUT&Tag vs ChIP-seq differential analysis
#
# Consumes the manifests, DESeq2 tables, feature tables and parameter file
# written by 01_differential_analysis.sh.
#
#   MAIN  violin + volcano, one column per epitope, coloured by ONE feature
#         (K562). One PDF per feature.
#
#   SUPP  violin + volcano + MA, one column per comparison, with all three
#         features stacked in a single PDF. Three files are produced:
#             SUPP_H3K4me3_violin_volcano_MA.pdf
#             SUPP_H3K27ac_violin_volcano_MA.pdf
#             SUPP_broadMarks_violin_volcano_MA.pdf   (H3K36me3 + H3K27me3)
#         Each holds nine panel rows: violin / volcano / MA coloured by GC
#         content, then the same three coloured by DNase-seq, then by ATAC-seq.
#         The intent is that a section can be lifted into an A4 Illustrator
#         canvas without hunting across separate PDFs for the matching panels,
#         so each feature block carries its own colour bar and legend and can
#         be cut out on its own.
#
# Orientation is fixed and identical in every panel: DESeq2 reports
# log2(CUT&Tag / ChIP-seq), so positive log2FC is CUT&Tag enrichment. The
# volcano puts CUT&Tag-enriched bins on the right and the violin puts the red
# CUT&Tag-enriched group on the right. Nothing is negated anywhere.
#
# Usage:
#   conda activate <env with R: ggplot2 dplyr readr tidyr patchwork ggrastr
#                               ggtext scales ggh4x>
#   bash 02_make_figures.sh
# =============================================================================

# --------------------------- USER CONFIG -------------------------------------
WORK_DIR="${WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output}"
DIFF_DIR="$WORK_DIR/differential"
RS_DIR="$WORK_DIR/Rscripts"
FIG_DIR="${FIG_DIR:-$WORK_DIR/figures_updatedSupp}"
mkdir -p "$FIG_DIR" "$RS_DIR"

MAIN_MANIFEST="$DIFF_DIR/manifest_main.tsv"
SUPP_MANIFEST="$DIFF_DIR/manifest_supp.tsv"
PARAMS="$DIFF_DIR/analysis_params.tsv"

# Features to colour by. In the supplementary figure the blocks are stacked
# TOP TO BOTTOM in this order; the main figure gets one PDF per entry.
COLOR_MODES=(GC DNASE ATAC)

# Global font multiplier applied to every text element in every panel.
FONT_SCALE="${FONT_SCALE:-1.5}"

# Which epitopes count as promoter (narrow) marks. Each of these gets its own
# supplementary file; everything else is pooled into "broadMarks", because
# H3K36me3 and H3K27me3 are read against each other and belong on one page
# whereas the two promoter marks are each wide enough to fill a sheet alone.
PROMOTER_MARKS="H3K4me3,H3K27ac"

for f in "$MAIN_MANIFEST" "$SUPP_MANIFEST" "$PARAMS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing input: $f" >&2; echo "Run 01_differential_analysis.sh first." >&2; exit 1; }
done
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found" >&2; exit 1; }

# =============================================================================
# SUPPLEMENTARY FIGURE SPLITTING
# =============================================================================
# Writes one manifest per epitope group, preserving the header and the original
# row order so panel order inside each file is unchanged. Prints "key<TAB>path"
# for each part, in the order the groups first appear. Because the comparison
# list is ordered by epitope, each group's rows are contiguous and the nested
# facet strips stay intact.
split_manifest() {   # in_manifest, out_dir
    local IN="$1" OUTDIR="$2"
    mkdir -p "$OUTDIR"
    awk -F'\t' -v d="$OUTDIR" -v prom="$PROMOTER_MARKS" '
      BEGIN { np = split(prom, pa, ","); for (i = 1; i <= np; i++) isprom[pa[i]] = 1 }
      NR == 1 { hdr = $0; next }
      {
        key = ($1 in isprom) ? $1 : "broadMarks"
        if (!(key in file)) { file[key] = d "/manifest_supp_" key ".tsv"; order[++n] = key; print hdr > file[key] }
        print > file[key]
      }
      END { for (i = 1; i <= n; i++) print order[i] "\t" file[order[i]] }
    ' "$IN"
}

# =============================================================================
# SHARED R CODE
# =============================================================================
cat > "$RS_DIR/fig_common.R" <<'RCOMMON'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(ggrastr); library(ggtext); library(scales)
})
options(warn = -1)

# --- Thresholds come from the analysis, never redeclared here ----------------
read_params <- function(path) {
  p <- suppressMessages(read_tsv(path, col_types = cols(.default = col_character())))
  setNames(as.list(p$value), p$key)
}
PARAMS      <- read_params(Sys.getenv("PARAMS_FILE"))
FC_THRESH   <- as.numeric(PARAMS$fc_thresh)
FDR_THRESH  <- as.numeric(PARAMS$fdr_thresh)
FONT_SCALE  <- as.numeric(Sys.getenv("FONT_SCALE", "1.5"))
fs <- function(x) x * FONT_SCALE          # one knob scales every font size

CHIP_COL    <- "#3936ff"
CNT_COL     <- "#ff3b3b"
NONDIFF_COL <- "grey70"
MEDIAN_COL  <- "#1a9d3b"

SIG_TEXT_SIZE    <- fs(3.6)   # adjusted p-values are long strings
# The green medians and the red/blue n= counts are the two numbers read off the
# figure most often, so they carry an extra 1.5x on top of the global scale.
MEDIAN_TEXT_SIZE <- fs(3.7 * 1.5)
COUNT_TEXT_SIZE  <- fs(4.0 * 1.5)

# Fixed display scales, identical in every PDF.
GC_COLOR_LIMITS <- c(30, 75)  # absolute % GC
Z_LIMITS        <- c(-3, 3)   # SD units within cell line x region set

# --- Orientation ------------------------------------------------------------
# log2FC is log2(CUT&Tag / ChIP-seq) as reported by DESeq2. Level 1 (left,
# blue) is ChIP-seq enriched and level 3 (right, red) is CUT&Tag enriched,
# matching the negative and positive halves of the volcano x-axis.
STATUS_LEVELS <- c("ChIP-seq Enriched", "Non-differential", "CUT&Tag Enriched")
STATUS_COLORS <- c("ChIP-seq Enriched" = CHIP_COL,
                   "Non-differential"  = NONDIFF_COL,
                   "CUT&Tag Enriched"  = CNT_COL)

REGION_SET_LABEL <- c(
  promoter_bins = paste0("active promoter, ", PARAMS$prom_bin_size, " bp bins"),
  genebody_bins = paste0("active gene body, ", as.numeric(PARAMS$broad_bin_size) / 1000, " kb bins"),
  polycomb_bins = paste0("Polycomb repressed, ", as.numeric(PARAMS$broad_bin_size) / 1000, " kb bins"),
  promoter_full = "whole active promoter"
)
region_set_label <- function(rs) {
  out <- unname(REGION_SET_LABEL[as.character(rs)])
  ifelse(is.na(out), as.character(rs), out)
}

# --- Feature readers --------------------------------------------------------
# deepTools writes a quoted comment header; the value column is the mean bigWig
# signal over the bin.
read_feature <- function(path, value_name) {
  if (!file.exists(path)) return(data.frame())
  df <- suppressMessages(read_tsv(path, comment = "#",
                                  col_names = c("Chr", "Start", "End", value_name),
                                  show_col_types = FALSE))
  df$Start <- as.integer(df$Start); df$End <- as.integer(df$End)
  df[[value_name]] <- suppressWarnings(as.numeric(df[[value_name]]))
  df
}

# GC unit detection is done ONCE per file from the file maximum, never per row,
# so a genuinely AT-poor bin can never be mistaken for a 0-1 fraction.
prep_gc <- function(path) {
  df <- read_feature(path, "GC_Raw")
  if (nrow(df) == 0) return(df)
  mx <- suppressWarnings(max(df$GC_Raw, na.rm = TRUE))
  df %>% mutate(GC_Percent = GC_Raw * if (is.finite(mx) && mx <= 1) 100 else 1)
}

prep_signal <- function(path, out_name) {
  df <- read_feature(path, "raw")
  if (nrow(df) == 0) return(df)
  df[[out_name]] <- log2(pmax(df$raw, 0) + 1)
  df[, c("Chr", "Start", "End", out_name)]
}

# --- Accessibility Z-score references ---------------------------------------
# A reference is built per cell line x region set, from that region set's own
# feature table. Bins of different size are therefore never pooled: promoter
# 500 bp bins are standardised against promoter bins only, gene-body 3 kb bins
# against gene-body bins only, and Polycomb bins against Polycomb bins only.
# The table defines the region set, so every comparison over the same bins uses
# an identical reference and repeated comparisons carry no extra weight.
build_reference <- function(manifest, path_col, value_name) {
  keys <- manifest %>% distinct(cell_line, region_set, .keep_all = TRUE)
  out <- list()
  for (i in seq_len(nrow(keys))) {
    d <- prep_signal(keys[[path_col]][i], value_name)
    if (nrow(d) == 0) next
    d <- d[is.finite(d[[value_name]]), ]
    d <- dplyr::distinct(d, Chr, Start, End, .keep_all = TRUE)
    if (nrow(d) < 3) next
    out[[length(out) + 1]] <- data.frame(
      cell_line = keys$cell_line[i], region_set = keys$region_set[i],
      n_bins = nrow(d), ref_mean = mean(d[[value_name]]), ref_sd = sd(d[[value_name]]))
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

# --- Loading ----------------------------------------------------------------
load_comparison <- function(row, col_order) {
  if (!file.exists(row$diff_path)) return(data.frame())
  d <- suppressMessages(read_tsv(row$diff_path, show_col_types = FALSE))
  d$Start <- as.integer(d$Start); d$End <- as.integer(d$End)
  gc <- prep_gc(row$gc_path)
  if (nrow(gc) == 0) return(data.frame())

  d <- d %>%
    inner_join(gc[, c("Chr", "Start", "End", "GC_Percent")], by = c("Chr", "Start", "End")) %>%
    filter(!is.na(GC_Percent), !is.na(log2FC), !is.na(padj)) %>%
    mutate(
      Status = case_when(
        log2FC >  FC_THRESH & padj < FDR_THRESH ~ "CUT&Tag Enriched",
        log2FC < -FC_THRESH & padj < FDR_THRESH ~ "ChIP-seq Enriched",
        TRUE ~ "Non-differential"),
      nLog10_FDR    = pmin(ifelse(padj == 0, 300, -log10(padj)), 300),
      BaseMean_Log2 = log2(baseMean + 1),
      Epitope = row$epitope, Cell_Line = row$cell_line, Region_Set = row$region_set,
      Region_Label = region_set_label(row$region_set),
      # Unique per comparison. Keying on cnt_group alone would collide as soon
      # as one CUT&Tag sample is compared against two ChIP-seq sources, and the
      # two comparisons would be silently pooled into a single panel.
      ColKey = row$comparison,
      ColOrder = col_order, Cnt_Full = row$cnt_label, Chip_Full = row$chip_label)

  for (spec in list(list("atac_path", "ATAC_log2"), list("dnase_path", "DNase_log2"))) {
    s <- prep_signal(row[[spec[[1]]]], spec[[2]])
    if (nrow(s) > 0) d <- left_join(d, s, by = c("Chr", "Start", "End")) else d[[spec[[2]]]] <- NA_real_
  }

  d %>% select(BinID, Chr, Start, End, GC_Percent, ATAC_log2, DNase_log2,
               log2FC, padj, nLog10_FDR, BaseMean_Log2, Status, Epitope, Cell_Line,
               Region_Set, Region_Label, ColKey, ColOrder, Cnt_Full, Chip_Full)
}

attach_z <- function(df, ref, value_col, z_col) {
  if (nrow(ref) == 0) { df[[z_col]] <- NA_real_; return(df) }
  df <- left_join(df, ref[, c("cell_line", "region_set", "ref_mean", "ref_sd")],
                  by = c("Cell_Line" = "cell_line", "Region_Set" = "region_set"))
  df[[z_col]] <- ifelse(is.na(df$ref_sd) | df$ref_sd == 0, NA_real_,
                        (df[[value_col]] - df$ref_mean) / df$ref_sd)
  df$ref_mean <- NULL; df$ref_sd <- NULL
  df
}

load_manifest <- function(path) {
  m <- suppressMessages(read_tsv(path, show_col_types = FALSE))
  out <- vector("list", nrow(m))
  for (i in seq_len(nrow(m))) out[[i]] <- load_comparison(m[i, ], i)
  big <- bind_rows(out)
  if (nrow(big) == 0) return(big)
  atac_ref  <- build_reference(m, "atac_path",  "ATAC_log2")
  dnase_ref <- build_reference(m, "dnase_path", "DNase_log2")
  big <- attach_z(big, atac_ref,  "ATAC_log2",  "ATAC_Z")
  big <- attach_z(big, dnase_ref, "DNase_log2", "DNase_Z")
  attr(big, "atac_ref") <- atac_ref; attr(big, "dnase_ref") <- dnase_ref
  big
}

# --- Metric specification ---------------------------------------------------
metric_spec <- function(mode) {
  switch(toupper(mode),
    ATAC  = list(var = "ATAC_Z", y_label = "ATAC-seq (Z-score)",
                 scale = scale_colour_viridis_c(name = "ATAC-seq\n(Z-score)", option = "viridis",
                                                limits = Z_LIMITS, oob = scales::squish)),
    DNASE = list(var = "DNase_Z", y_label = "DNase-seq (Z-score)",
                 scale = scale_colour_viridis_c(name = "DNase-seq\n(Z-score)", option = "mako",
                                                limits = Z_LIMITS, oob = scales::squish)),
    list(var = "GC_Percent", y_label = "GC content (%)",
         scale = scale_colour_viridis_c(name = "GC %", option = "magma",
                                        limits = GC_COLOR_LIMITS, oob = scales::squish)))
}

# --- Statistics -------------------------------------------------------------
# The adjusted p-value itself is printed, not a star code. wilcox.test uses the
# normal approximation at these sample sizes and pnorm() underflows to exactly
# zero for very large |z|, so an exact zero is reported as a bound rather than
# as "0" or "Inf".
format_padj <- function(p) {
  if (is.na(p)) return("n/a")
  if (p <= 0)     return(paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 0)))
  if (p >= 0.001) return(paste0("p = ", formatC(p, format = "f", digits = 3)))
  paste0("p = ", formatC(p, format = "e", digits = 1))
}

make_annotations <- function(big_df, median_y, b_y1, b_y2, metric_col) {
  stat_df <- big_df %>%
    group_by(ColKey, Status, .drop = FALSE) %>%
    summarise(n = n(), median_val = median(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(y_pos = median_y, med_label = ifelse(n > 0, sprintf("%.2f", median_val), ""))

  b_df <- data.frame()
  for (k in unique(as.character(big_df$ColKey))) {
    sub <- big_df %>% filter(as.character(ColKey) == k)
    # Groups are pulled left-to-right in the same order as STATUS_LEVELS.
    g1 <- sub[[metric_col]][sub$Status == "ChIP-seq Enriched"]
    g2 <- sub[[metric_col]][sub$Status == "Non-differential"]
    g3 <- sub[[metric_col]][sub$Status == "CUT&Tag Enriched"]
    # Benjamini-Hochberg across the three contrasts WITHIN this panel.
    padj <- p.adjust(c(
      tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA),
      tryCatch(wilcox.test(g2, g3)$p.value, error = function(e) NA),
      tryCatch(wilcox.test(g1, g3)$p.value, error = function(e) NA)), method = "BH")
    b_df <- rbind(b_df,
      data.frame(x = 1, xend = 2, y = b_y1, label = format_padj(padj[1]), ColKey = k),
      data.frame(x = 2, xend = 3, y = b_y1, label = format_padj(padj[2]), ColKey = k),
      data.frame(x = 1, xend = 3, y = b_y2, label = format_padj(padj[3]), ColKey = k))
  }
  km <- big_df %>% distinct(ColKey, Epitope, Cell_Line)
  b_df$ColKey <- factor(b_df$ColKey, levels = levels(big_df$ColKey))
  list(stat_df = left_join(stat_df, km, by = "ColKey"),
       b_df    = left_join(b_df,    km, by = "ColKey"))
}

# --- Panels -----------------------------------------------------------------
sub_theme <- function() theme_minimal(base_size = fs(18)) +
  theme(strip.text = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.title = element_text(size = fs(22), face = "bold"),
        axis.text  = element_text(size = fs(18), colour = "black"),
        panel.spacing = unit(0.9, "cm"))

# Strip label used by the MAIN figure: ChIP on top and CUT&Tag below, matching
# left-to-right on the x-axis, with the region set that defines the Z-score
# reference underneath.
label_check <- function(cnt, chip, region) paste0(
  "<span style='font-size:", fs(8),   "pt;color:grey35'>ChIP: ", chip, "</span><br>",
  "<span style='font-size:", fs(8),   "pt;color:grey35'>CnT: ",  cnt,  "</span><br>",
  "<span style='font-size:", fs(7.5), "pt;color:grey55'>", region, "</span>")

# --- Short source name ------------------------------------------------------
# Reduces a full sample label to "Author Year" for the supplementary strip
# header:
#   chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2   -> Bernstein 2017
#   CnT_H3K4me3_K562_KayaOkur_2020_rep1_2020_rep2     -> KayaOkur 2020
#   CnT_H3K27ac_K562_Abbasova_2025_ab177r1            -> Abbasova 2025
#
# Matched on a four-digit year and the token immediately before it, rather than
# on a fixed field position: group labels append replicate suffixes, so
# counting underscores from either end would land on a different token
# depending on how many replicates a group happens to have.
#
# EVERY distinct year in the label is collected, not just the first. Several
# K562 ChIP-seq groups pool replicates submitted in different years, e.g.
#   chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3  ->  Bernstein 2011/2022
# and reporting only "Bernstein 2011" there would name one submission for a
# group that contains two. The author is taken from the first match, since the
# replicate tokens ("rep2_2022") also satisfy the pattern and would otherwise
# be mistaken for author names.
#
# Anything with no year at all is returned unchanged, so a renamed sample
# degrades to its full label rather than to an empty header.
short_source <- function(label) {
  vapply(as.character(label), function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    m <- regmatches(x, gregexpr("[A-Za-z][A-Za-z0-9.-]*_(19|20)[0-9]{2}", x))[[1]]
    if (length(m) == 0) return(x)
    parts   <- strsplit(m, "_", fixed = TRUE)
    authors <- unique(vapply(parts, function(p) p[1], character(1)))
    years   <- unique(vapply(parts, function(p) p[2], character(1)))
    paste0(authors[1], " ", paste(years, collapse = "/"))
  }, character(1), USE.NAMES = FALSE)
}

# Strip label used by the SUPPLEMENTARY figures. A large coloured "Author Year"
# pair sits above the small grey lines, so a panel can be identified at a
# glance while the exact sample that produced it is still printed underneath
# and remains checkable.
#
# ChIP-seq is on the left in blue and CUT&Tag on the right in red: the same two
# colours as the violin groups and the volcano n= counts, in the same order as
# the violin x-axis, where the blue ChIP-seq-enriched group sits left and the
# red CUT&Tag-enriched group sits right.
label_check_supp <- function(cnt, chip, region) paste0(
  "<span style='font-size:", fs(13),  "pt;color:", CHIP_COL, "'>**", short_source(chip), "**</span>",
  "<span style='font-size:", fs(10),  "pt;color:grey45'> vs </span>",
  "<span style='font-size:", fs(13),  "pt;color:", CNT_COL,  "'>**", short_source(cnt),  "**</span><br>",
  "<span style='font-size:", fs(7.5), "pt;color:grey35'>ChIP: ", chip, "</span><br>",
  "<span style='font-size:", fs(7.5), "pt;color:grey35'>CnT: ",  cnt,  "</span><br>",
  "<span style='font-size:", fs(7),   "pt;color:grey55'>", region, "</span>")

side_legend_plot <- function() {
  d1 <- data.frame(x = 0, y = 0.80, label = paste0(
    "<span style='font-size:", fs(16), "pt'>**Enriched in**</span><br><br>",
    "<span style='font-size:", fs(18), "pt;color:", CHIP_COL, "'>**ChIP-seq**</span><br>",
    "<span style='font-size:", fs(14), "pt'>vs</span><br>",
    "<span style='font-size:", fs(18), "pt;color:", CNT_COL, "'>**CUT&Tag**</span>"))
  d2 <- data.frame(x = 0.16, y = 0.34, label = paste0(
    "<span style='font-size:", fs(16), "pt;color:", MEDIAN_COL, "'>**median values**</span>"))
  ggplot() + theme_void() + xlim(0, 1) + ylim(0, 1) +
    geom_richtext(data = d1, aes(x, y, label = label), hjust = 0, vjust = 1, fill = NA, label.color = NA) +
    annotate("point", x = 0.06, y = 0.34, colour = MEDIAN_COL, size = fs(4)) +
    geom_richtext(data = d2, aes(x, y, label = label), hjust = 0, vjust = 0.5, fill = NA, label.color = NA)
}

build_violin <- function(big_df, facet_layer, strip_element, metric_col, y_label) {
  big_df <- big_df %>% filter(!is.na(.data[[metric_col]]))
  if (nrow(big_df) == 0) stop("no non-missing values for metric: ", metric_col)
  max_val <- max(big_df[[metric_col]], na.rm = TRUE)
  min_val <- min(big_df[[metric_col]], na.rm = TRUE)
  rng <- ifelse(max_val - min_val == 0, 1, max_val - min_val)
  # Headroom above the data for the median labels and the two bracket rows.
  median_y <- max_val + rng * 0.08
  b_y1 <- max_val + rng * 0.22; b_y2 <- max_val + rng * 0.38; y_top <- max_val + rng * 0.52
  ann <- make_annotations(big_df, median_y, b_y1, b_y2, metric_col)

  ggplot(big_df, aes(x = Status, y = .data[[metric_col]], fill = Status)) +
    geom_violin(alpha = 0.7, colour = "black", width = 0.85) +
    geom_boxplot(width = 0.18, colour = "black", outlier.shape = NA) +
    geom_segment(data = ann$b_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
    geom_segment(data = ann$b_df, aes(x = x, xend = x, y = y - rng * 0.03, yend = y), inherit.aes = FALSE) +
    geom_segment(data = ann$b_df, aes(x = xend, xend = xend, y = y - rng * 0.03, yend = y), inherit.aes = FALSE) +
    geom_text(data = ann$b_df, aes(x = (x + xend) / 2, y = y + rng * 0.01, label = label),
              inherit.aes = FALSE, size = SIG_TEXT_SIZE, fontface = "bold", vjust = 0) +
    geom_text(data = dplyr::filter(ann$stat_df, n > 0),
              aes(x = Status, y = y_pos, label = med_label), inherit.aes = FALSE,
              size = MEDIAN_TEXT_SIZE, fontface = "bold.italic", colour = MEDIAN_COL) +
    scale_fill_manual(values = STATUS_COLORS, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE) + facet_layer +
    coord_cartesian(ylim = c(min_val - rng * 0.05, y_top)) +
    labs(y = y_label, x = NULL) +
    theme_minimal(base_size = fs(18)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = fs(20), colour = "black"),
          axis.title.y = element_text(size = fs(22), face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
          panel.spacing = unit(0.9, "cm"),
          strip.background = element_blank(), strip.text = strip_element)
}

build_volcano <- function(big_df, facet_layer, color_var, color_scale) {
  km <- big_df %>% distinct(ColKey, Epitope, Cell_Line)
  ann <- big_df %>% group_by(ColKey) %>%
    summarise(n_chip = sum(Status == "ChIP-seq Enriched", na.rm = TRUE),
              n_cnt  = sum(Status == "CUT&Tag Enriched",  na.rm = TRUE), .groups = "drop") %>%
    mutate(lbl_chip = paste0("n=", n_chip), lbl_cnt = paste0("n=", n_cnt)) %>%
    left_join(km, by = "ColKey")
  ggplot(big_df, aes(x = log2FC, y = nLog10_FDR, colour = .data[[color_var]])) +
    geom_point_rast(alpha = 0.5, size = 0.5, raster.dpi = 300) + color_scale +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(FDR_THRESH), linetype = "dashed", colour = "grey50") +
    geom_text(data = ann, aes(x =  Inf, y = Inf, label = lbl_cnt),  hjust = 1.1,  vjust = 1.5,
              size = COUNT_TEXT_SIZE, fontface = "bold", colour = CNT_COL,  inherit.aes = FALSE) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lbl_chip), hjust = -0.1, vjust = 1.5,
              size = COUNT_TEXT_SIZE, fontface = "bold", colour = CHIP_COL, inherit.aes = FALSE) +
    facet_layer + labs(x = "log2(CUT&Tag / ChIP-seq)", y = "-log10(FDR)") + sub_theme()
}

build_ma <- function(big_df, facet_layer, color_var, color_scale) {
  ggplot(big_df, aes(x = BaseMean_Log2, y = log2FC, colour = .data[[color_var]])) +
    geom_point_rast(alpha = 0.5, size = 0.5, raster.dpi = 300) + color_scale +
    geom_hline(yintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "black") +
    facet_layer + labs(x = "Mean normalised signal (log2)", y = "log2(CUT&Tag / ChIP-seq)") + sub_theme()
}

# --- Caption ----------------------------------------------------------------
caption_for <- function(mode) {
  quant <- paste0(
    "Quantification: fragments per bin counted from filtered, deduplicated BAMs (deepTools multiBamSummary). ",
    "Paired-end libraries are extended to their true fragment span and counted on first mates only, so each fragment contributes once; ",
    "single-end libraries are extended to ", PARAMS$se_fragment_length, " bp. ",
    "Differential enrichment by DESeq2, contrast CUT&Tag vs ChIP-seq; positive log2FC is CUT&Tag enrichment throughout. ",
    "Thresholds: FDR < ", FDR_THRESH, ", |log2FC| > ", FC_THRESH, ".")

  bins <- paste0(
    "Bins: active promoter = TSS \u00b1 ", PARAMS$prom_up, " bp of active genes tiled in ", PARAMS$prom_bin_size, " bp bins (H3K4me3, H3K27ac); ",
    "active gene body = ChromHMM Tx/TxWk bodies in ", as.numeric(PARAMS$broad_bin_size)/1000, " kb bins, TSS \u00b1 ", PARAMS$tss_exclude_window, " bp excluded (H3K36me3); ",
    "Polycomb repressed = ChromHMM ReprPC/ReprPCWk domains in ", as.numeric(PARAMS$broad_bin_size)/1000, " kb bins (H3K27me3).")

  gc_txt <- paste0(
    "GC content is the mean gc5Base value per bin, plotted on a fixed absolute scale (% GC, colour range ",
    GC_COLOR_LIMITS[1], "-", GC_COLOR_LIMITS[2], "%) with no standardisation, so GC panels are directly comparable across every figure.")

  acc_txt <- paste0(
    "Accessibility is the mean bigWig signal per bin, summarised as log2(signal + 1) and standardised as a Z-score ",
    "within each cell line \u00d7 region set: promoter bins are standardised against promoter bins only, gene-body bins against gene-body bins only ",
    "and Polycomb bins against Polycomb bins only, so bins of different size and chromatin class are never pooled. ",
    "Each reference mean and SD is computed once over the deduplicated bins of that region set, independently of which comparisons appear in a figure. ",
    "Z = 0 is the mean accessibility of that cell line and region set; colour range \u00b1", Z_LIMITS[2], " SD, values beyond the range are clamped. ",
    "ATAC-seq is the merged biological-replicate CPM track; DNase-seq is the matched ENCODE track and is Tn5-independent, ",
    "so it does not share the transposase chemistry that CUT&Tag and ATAC-seq have in common.")

  stats_txt <- paste0(
    "Violins: pairwise Wilcoxon rank-sum tests, Benjamini-Hochberg corrected across the three contrasts within each panel; ",
    "brackets show the adjusted p-value, and p-values below double precision are reported as a bound. ",
    "Green values are group medians of the plotted metric. ",
    "Strip headers give the ChIP-seq source in blue and the CUT&Tag source in red, with the exact sample identifiers printed beneath.")

  feature_txt <- if (toupper(mode) == "GC") paste(gc_txt, acc_txt, sep = "\n") else paste(acc_txt, gc_txt, sep = "\n")
  paste(quant, bins, feature_txt, stats_txt, sep = "\n")
}

# Written next to the figures so the exact reference used is recoverable.
write_reference_table <- function(big_df, out_tsv) {
  a <- attr(big_df, "atac_ref"); d <- attr(big_df, "dnase_ref")
  if (!is.null(a) && nrow(a) > 0) a$feature <- "ATAC-seq"
  if (!is.null(d) && nrow(d) > 0) d$feature <- "DNase-seq"
  ref <- bind_rows(a, d)
  if (nrow(ref) > 0) write_tsv(ref[, c("feature", "cell_line", "region_set", "n_bins", "ref_mean", "ref_sd")], out_tsv)
}
RCOMMON

# --- MAIN figure: violin + volcano, one column per epitope -------------------
cat > "$RS_DIR/fig_main_body.R" <<'RMAIN'
args <- commandArgs(trailingOnly = TRUE)
manifest <- args[1]; out_pdf <- args[2]; color_mode <- args[3]

big_df <- load_manifest(manifest)
if (nrow(big_df) == 0) { message("empty manifest; skipping ", out_pdf); q(save = "no") }

ms <- metric_spec(color_mode)
if (all(is.na(big_df[[ms$var]]))) { message("no ", color_mode, " values; skipping ", out_pdf); q(save = "no") }

lab_tbl <- big_df %>% distinct(ColKey, ColOrder, Cnt_Full, Chip_Full, Region_Label) %>% arrange(ColOrder)
col_lab <- setNames(mapply(label_check, lab_tbl$Cnt_Full, lab_tbl$Chip_Full, lab_tbl$Region_Label), lab_tbl$ColKey)

big_df$Status    <- factor(big_df$Status, levels = STATUS_LEVELS)
big_df$ColKey    <- factor(big_df$ColKey, levels = lab_tbl$ColKey)
big_df$Epitope   <- factor(big_df$Epitope, levels = unique(big_df$Epitope[order(big_df$ColOrder)]))
big_df$Cell_Line <- factor(big_df$Cell_Line, levels = unique(big_df$Cell_Line))

viol <- build_violin(big_df, facet_wrap(~ ColKey, nrow = 1, labeller = labeller(ColKey = col_lab)),
                     element_markdown(size = fs(8), lineheight = 1.15, halign = 0.5, margin = margin(b = 6)),
                     ms$var, ms$y_label)
volc <- build_volcano(big_df, facet_wrap(~ ColKey, nrow = 1, scales = "free"), ms$var, ms$scale)
body <- viol / volc + plot_layout(heights = c(0.8, 1))

final <- wrap_plots(body, side_legend_plot(), ncol = 2, widths = c(1, 0.20), guides = "collect") +
  plot_annotation(
    subtitle = sprintf("K562 CUT&Tag vs ChIP-seq over ChromHMM ground-truth bins  |  FDR < %s, |log2FC| > %s",
                       FDR_THRESH, FC_THRESH),
    caption  = caption_for(color_mode),
    theme = theme(plot.subtitle = element_text(size = fs(14), face = "italic"),
                  plot.caption  = element_text(size = fs(8), face = "italic", hjust = 0,
                                               colour = "grey30", lineheight = 1.25))) &
  theme(legend.position = "right",
        legend.text  = element_text(size = fs(15), face = "bold"),
        legend.title = element_text(size = fs(16), face = "bold"),
        legend.key.height = unit(1.8, "cm"))

write_reference_table(big_df, sub("\\.pdf$", "_zscore_reference.tsv", out_pdf))
ggsave(out_pdf, final, width = 9 + nrow(lab_tbl) * 4.6, height = 16,
       device = "pdf", bg = "white", limitsize = FALSE)
message("wrote ", out_pdf)
RMAIN

# --- SUPP figure: violin + volcano + MA, one block per feature ---------------
cat > "$RS_DIR/fig_supp_body.R" <<'RSUPP'
args <- commandArgs(trailingOnly = TRUE)
manifest <- args[1]; out_pdf <- args[2]
# args[3] is a comma-separated feature list; one three-row block is drawn per
# feature, stacked in the order given.
color_modes <- strsplit(args[3], ",", fixed = TRUE)[[1]]
color_modes <- trimws(color_modes[nzchar(trimws(color_modes))])
library(ggh4x)

# The manifest is loaded ONCE and reused by every feature block. Rendering the
# features as separate PDFs re-read and re-joined every DESeq2 and feature
# table per feature, which was the dominant cost of this script.
big_df <- load_manifest(manifest)
if (nrow(big_df) == 0) { message("empty manifest; skipping ", out_pdf); q(save = "no") }

# Epitopes in first-appearance order, so a split file names exactly what it
# holds ("H3K36me3 + H3K27me3: all ...") instead of still claiming "All".
epi_present <- unique(as.character(big_df$Epitope))
scope_txt <- if (length(epi_present) <= 3) paste0(paste(epi_present, collapse = " + "), ": all") else "All"

lab_tbl <- big_df %>% distinct(ColKey, ColOrder, Epitope, Cell_Line, Cnt_Full, Chip_Full, Region_Label) %>% arrange(ColOrder)
# label_check_supp, not label_check: the supplementary strips carry the large
# coloured "Author Year" header above the grey detail lines.
col_lab <- setNames(mapply(label_check_supp, lab_tbl$Cnt_Full, lab_tbl$Chip_Full, lab_tbl$Region_Label), lab_tbl$ColKey)

epi_levels <- unique(lab_tbl$Epitope)
ep_lab <- setNames(paste0("<span style='font-size:", fs(19), "pt'>**", epi_levels, "**</span>"), epi_levels)
cl_levels <- unique(lab_tbl$Cell_Line)
cl_lab <- setNames(paste0("<span style='font-size:", fs(14), "pt'>**", cl_levels, "**</span>"), cl_levels)

big_df$Status    <- factor(big_df$Status, levels = STATUS_LEVELS)
big_df$Epitope   <- factor(big_df$Epitope, levels = epi_levels)
big_df$Cell_Line <- factor(big_df$Cell_Line, levels = cl_levels)
big_df$ColKey    <- factor(big_df$ColKey, levels = lab_tbl$ColKey)

viol_facet <- facet_nested_wrap(vars(Epitope, Cell_Line, ColKey), nrow = 1,
  nest_line = element_line(colour = "grey15", linewidth = 0.7), solo_line = TRUE,
  labeller = labeller(Epitope = ep_lab, Cell_Line = cl_lab, ColKey = col_lab),
  strip = strip_nested(background_x = elem_list_rect(fill = NA, colour = NA)))
sub_facet <- facet_nested_wrap(vars(Epitope, Cell_Line, ColKey), nrow = 1,
                               scales = "free", nest_line = element_blank())

# The header is four lines rather than three, so the strip needs a little more
# line spacing and bottom margin than the main figure's.
strip_el <- element_markdown(size = fs(8), lineheight = 1.25, halign = 0.5, margin = margin(b = 6))

# The side legend holds a fixed block of text, so a fixed 12% share becomes
# unreadably narrow once a figure is split down to two or three panels.
legend_w <- min(0.35, max(0.12, 2.2 / nrow(lab_tbl)))

# One self-contained block per feature: violin, volcano, MA, plus its own
# colour bar and the shared enrichment key. Guides are collected INSIDE the
# block rather than at the top of the page, so each block travels with the
# colour bar that explains it and a block can be cut out on its own.
build_feature_block <- function(mode) {
  ms <- metric_spec(mode)
  if (all(is.na(big_df[[ms$var]]))) {
    message("  no ", mode, " values; block skipped")
    return(NULL)
  }
  viol <- build_violin(big_df, viol_facet, strip_el, ms$var, ms$y_label) +
    labs(title = ms$y_label) +
    theme(plot.title = element_text(size = fs(24), face = "bold", hjust = 0,
                                    margin = margin(b = 10)))
  volc <- build_volcano(big_df, sub_facet, ms$var, ms$scale)
  ma   <- build_ma(big_df, sub_facet, ms$var, ms$scale)
  body <- viol / volc / ma + plot_layout(heights = c(0.9, 1, 1))
  wrap_plots(body, side_legend_plot(), ncol = 2, widths = c(1, legend_w), guides = "collect")
}

blocks <- Filter(Negate(is.null), lapply(color_modes, build_feature_block))
if (length(blocks) == 0) { message("no plottable features; skipping ", out_pdf); q(save = "no") }

feat_txt <- paste(vapply(color_modes, function(m) metric_spec(m)$y_label, character(1)), collapse = ", ")

final <- wrap_plots(blocks, ncol = 1) +
  plot_annotation(
    subtitle = sprintf("%s CUT&Tag vs ChIP-seq comparisons over ChromHMM ground-truth bins  |  colour features, top to bottom: %s  |  FDR < %s, |log2FC| > %s",
                       scope_txt, feat_txt, FDR_THRESH, FC_THRESH),
    # The caption text is the same for every feature; the first one only sets
    # which of the GC and accessibility paragraphs is printed first.
    caption  = caption_for(color_modes[1]),
    theme = theme(plot.subtitle = element_text(size = fs(16), face = "italic"),
                  plot.caption  = element_text(size = fs(8), face = "italic", hjust = 0,
                                               colour = "grey30", lineheight = 1.25))) &
  theme(legend.position = "right",
        legend.text  = element_text(size = fs(15), face = "bold"),
        legend.title = element_text(size = fs(16), face = "bold"),
        legend.key.height = unit(1.8, "cm"))

write_reference_table(big_df, sub("\\.pdf$", "_zscore_reference.tsv", out_pdf))
# Height scales with the number of feature blocks; the extra 5 in is the
# subtitle and the multi-line caption, which are drawn once per page.
ggsave(out_pdf, final, width = 9 + nrow(lab_tbl) * 4.4,
       height = 5 + length(blocks) * 23,
       device = "pdf", bg = "white", limitsize = FALSE)
message("wrote ", out_pdf, "  (", length(blocks), " feature blocks x 3 rows, ",
        nrow(lab_tbl), " comparisons)")
RSUPP

# =============================================================================
# RENDER
# =============================================================================
cat "$RS_DIR/fig_common.R" "$RS_DIR/fig_main_body.R" > "$RS_DIR/run_main_figure.R"
cat "$RS_DIR/fig_common.R" "$RS_DIR/fig_supp_body.R" > "$RS_DIR/run_supp_figure.R"

export PARAMS_FILE="$PARAMS" FONT_SCALE="$FONT_SCALE"

# One MAIN PDF per feature.
for MODE in "${COLOR_MODES[@]}"; do
    echo "=== MAIN figure coloured by ${MODE} ==="
    Rscript "$RS_DIR/run_main_figure.R" "$MAIN_MANIFEST" \
        "$FIG_DIR/MAIN_K562_violin_volcano_${MODE}.pdf" "$MODE"
done

# One SUPP PDF per epitope group, each holding every feature.
MODE_LIST=$(IFS=','; echo "${COLOR_MODES[*]}")
echo "=== SUPP figures (features per file: ${MODE_LIST}) ==="
while IFS=$'\t' read -r KEY MANIFEST_PART; do
    [[ -n "$KEY" ]] || continue
    echo "  -> SUPP ${KEY}"
    Rscript "$RS_DIR/run_supp_figure.R" "$MANIFEST_PART" \
        "$FIG_DIR/SUPP_${KEY}_violin_volcano_MA.pdf" "$MODE_LIST"
done < <(split_manifest "$SUPP_MANIFEST" "$FIG_DIR/manifests_split")

echo "========================================================="
echo " Figures: $FIG_DIR"
echo "   MAIN : MAIN_K562_violin_volcano_<feature>.pdf"
echo "   SUPP : SUPP_{H3K4me3,H3K27ac,broadMarks}_violin_volcano_MA.pdf"
echo "========================================================="
