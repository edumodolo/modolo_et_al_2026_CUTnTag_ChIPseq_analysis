#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# FEATURE IMPORTANCE: does the CUT&Tag vs ChIP-seq difference track GC content
# or chromatin accessibility, and in which direction?
#
# For every comparison the per-region log2FC is correlated against each feature
# and the SIGNED Pearson coefficient is reported.
#
# Positive r means the feature increases with log2FC, and log2FC is
# log2(CUT&Tag / ChIP-seq), so positive r = the feature is higher in
# CUT&Tag-enriched regions and negative r = higher in ChIP-seq-enriched regions.
#
# Reads manifest_supp.tsv and manifest_main.tsv from 01_differential_analysis.sh
# and produces:
#   feature_importance_direction_main.pdf   the main-figure comparisons
#   feature_importance_direction_supp.pdf   every comparison
#   feature_importance_correlations.tsv     n, r, R^2, p, p_adj for every cell
#
# ONE MULTIPLE-TESTING CORRECTION FOR BOTH FIGURES. Every correlation in the
# whole analysis is computed in a single R session and Benjamini-Hochberg is
# applied ONCE across all of them. The main figure is then drawn as a subset of
# those already-adjusted values.
#
# This matters because the main comparisons also appear in the supplementary
# figure. Correcting each figure separately adjusted those shared correlations
# against a different number of tests in each, so one correlation could carry
# two different adjusted p-values depending on which panel it was read from --
# a discrepancy a reader comparing the two figures would notice and could not
# explain. Correcting once removes it: a given correlation now has exactly one
# adjusted p-value, and the captions can say so.
#
# The y-axis limit is likewise derived once from the full set, so bar heights
# are directly comparable between the main and supplementary figures.
#
# REGION GEOMETRY: the binned sets only -- promoter marks in 500 bp promoter
# bins, broad marks in 3 kb bins. This is the geometry used for the
# violin/volcano/MA figures, so the correlation figure is measured over the same
# regions as the differential results it is describing.
#
# Usage:
#   conda activate chrom_diff_figures
#   bash 04_feature_importance.sh
# =============================================================================

# --------------------------- USER CONFIG -------------------------------------
WORK_DIR="${WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output}"
DIFF_DIR="$WORK_DIR/differential"
PARAMS="$DIFF_DIR/analysis_params.tsv"
OUT_DIR="${OUT_DIR:-$WORK_DIR/feature_importance_v2}"
RS_DIR="$WORK_DIR/Rscripts"
mkdir -p "$OUT_DIR" "$RS_DIR"

MAIN_MANIFEST="$DIFF_DIR/manifest_main.tsv"
SUPP_MANIFEST="$DIFF_DIR/manifest_supp.tsv"

# Print the value and its BH-adjusted p-value above every bar. With tens of
# thousands of regions per comparison almost every correlation is significant,
# so these labels are there to show the EFFECT SIZE and to make the handful of
# genuinely non-significant cases visible -- not to establish significance.
SHOW_BAR_STATS="${SHOW_BAR_STATS:-1}"

# 1.5x the scale used by the other figure scripts.
FONT_SCALE="${FONT_SCALE:-2.25}"

# --- Axis limits -------------------------------------------------------------
# Derived from the data: each end of the axis is the observed extreme in that
# direction plus 10% headroom, computed once over ALL comparisons so that the
# main and supplementary figures share one scale and a bar is the same height
# in both. There is no fixed-limit option: a hand-set limit either wastes panel
# space or clips a bar, and both figures already share a scale without one.
# --------------------------- PREFLIGHT ---------------------------------------
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found" >&2; exit 1; }
[[ -s "$PARAMS" ]] || { echo "ERROR: missing $PARAMS -- run 01_differential_analysis.sh first" >&2; exit 1; }
[[ -s "$SUPP_MANIFEST" ]] || { echo "ERROR: missing $SUPP_MANIFEST -- run 01_differential_analysis.sh first" >&2; exit 1; }
[[ -s "$MAIN_MANIFEST" ]] || echo "WARNING: $MAIN_MANIFEST not found -- only the supplementary figure will be drawn" >&2

echo "========================================================="
echo " Feature importance (direction only, binned regions)"
echo " Differential input : $DIFF_DIR"
echo " Output             : $OUT_DIR"
echo "========================================================="

# =============================================================================
# R
# =============================================================================
cat > "$RS_DIR/feature_importance.R" <<'RFI'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(ggh4x); library(ggtext); library(scales)
})
options(warn = -1)

a <- commandArgs(trailingOnly = TRUE)
supp_manifest <- a[1]; main_manifest <- a[2]
out_pdf_supp  <- a[3]; out_pdf_main  <- a[4]
out_tsv       <- a[5]
region_note   <- a[6]; params_file <- a[7]
show_stats    <- as.integer(a[8])

FONT_SCALE <- as.numeric(Sys.getenv("FONT_SCALE", "2.25"))
fs <- function(x) x * FONT_SCALE

prm <- suppressMessages(read_tsv(params_file, col_types = cols(.default = col_character())))
P <- setNames(as.list(prm$value), prm$key)

# Colours match the tracks these features are drawn as in the control figures,
# so a feature keeps one identity across the whole manuscript.
FEATURE_LEVELS <- c("GC content", "ATAC-seq", "DNase-seq")
FEATURE_COLORS <- c("GC content" = "#5D5D5D",   # rgb(93,93,93)
                    "ATAC-seq"   = "#ED00E6",   # rgb(237,0,230)
                    "DNase-seq"  = "#02C49B")   # rgb(2,196,155)

# --- minimal loader ----------------------------------------------------------
# No Z-scores are needed anywhere in this script: Pearson r is invariant to
# linear rescaling, so standardising the features would give identical
# coefficients. The raw quantifications are used directly.
read_feature <- function(path, value_name) {
  if (!file.exists(path)) return(data.frame())
  df <- suppressMessages(read_tsv(path, comment = "#",
                                  col_names = c("Chr", "Start", "End", value_name),
                                  show_col_types = FALSE))
  df$Start <- as.integer(df$Start); df$End <- as.integer(df$End)
  df[[value_name]] <- suppressWarnings(as.numeric(df[[value_name]]))
  df
}

load_one <- function(row) {
  if (!file.exists(row$diff_path)) { message("  missing: ", row$diff_path); return(data.frame()) }
  d <- suppressMessages(read_tsv(row$diff_path, show_col_types = FALSE))
  d$Start <- as.integer(d$Start); d$End <- as.integer(d$End)

  gc <- read_feature(row$gc_path, "GC_Raw")
  if (nrow(gc) == 0) return(data.frame())
  # Unit detection once per file, from the file maximum, never per row.
  mx <- suppressWarnings(max(gc$GC_Raw, na.rm = TRUE))
  gc$GC_Percent <- gc$GC_Raw * (if (is.finite(mx) && mx <= 1) 100 else 1)
  d <- inner_join(d, gc[, c("Chr", "Start", "End", "GC_Percent")], by = c("Chr", "Start", "End"))

  for (sp in list(c("atac_path", "ATAC_log2"), c("dnase_path", "DNase_log2"))) {
    s <- read_feature(row[[sp[1]]], "raw")
    if (nrow(s) > 0) {
      s[[sp[2]]] <- log2(pmax(s$raw, 0) + 1)
      d <- left_join(d, s[, c("Chr", "Start", "End", sp[2])], by = c("Chr", "Start", "End"))
    } else d[[sp[2]]] <- NA_real_
  }

  d %>%
    filter(is.finite(log2FC)) %>%
    mutate(epitope = row$epitope, cell_line = row$cell_line, region_set = row$region_set,
           comparison = row$comparison, cnt_label = row$cnt_label, chip_label = row$chip_label,
           col_order = row$.order)
}

# EVERY comparison is loaded here, from the supplementary manifest. The main
# figure is a subset drawn later; it is never loaded or tested separately.
m <- suppressMessages(read_tsv(supp_manifest, show_col_types = FALSE))
if (nrow(m) == 0) { message("empty manifest"); q(save = "no") }
m$.order <- seq_len(nrow(m))
big <- bind_rows(lapply(seq_len(nrow(m)), function(i) load_one(m[i, ])))
if (nrow(big) == 0) { message("no data loaded"); q(save = "no") }

# --- statistics --------------------------------------------------------------
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10 || sd(x[ok]) == 0 || sd(y[ok]) == 0)
    return(list(r = NA_real_, r2 = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok]))
  list(r = unname(ct$estimate), r2 = unname(ct$estimate)^2, p = ct$p.value, n = sum(ok))
}

keys <- big %>% distinct(epitope, cell_line, region_set, comparison, cnt_label, chip_label, col_order)
rows <- list()
for (i in seq_len(nrow(keys))) {
  d <- big %>% filter(comparison == keys$comparison[i], region_set == keys$region_set[i])
  one <- function(feature, col) {
    s <- safe_cor(d$log2FC, d[[col]])
    data.frame(keys[i, ], Feature = feature, n_bins = s$n,
               pearson_r = s$r, r_squared = s$r2, p_value = s$p, stringsAsFactors = FALSE)
  }
  rows[[length(rows) + 1]] <- one("GC content", "GC_Percent")
  rows[[length(rows) + 1]] <- one("ATAC-seq",  "ATAC_log2")
  rows[[length(rows) + 1]] <- one("DNase-seq", "DNase_log2")
}
stats_df <- bind_rows(rows) %>%
  mutate(Feature = factor(Feature, levels = FEATURE_LEVELS))

# --- ONE Benjamini-Hochberg correction, across every correlation -------------
# p-values are the standard test of H0: rho = 0 from cor.test(). The adjustment
# spans the COMPLETE set of comparisons, not the contents of one figure, so the
# adjusted value for a correlation is the same number wherever it is displayed.
# Note the sample size: each test uses tens of thousands of regions, so the test
# has enormous power and detects correlations far too small to matter.
stats_df$p_adj <- p.adjust(stats_df$p_value, method = "BH")
N_TESTS <- sum(!is.na(stats_df$p_value))
message("  ", nrow(stats_df), " correlations across ", nrow(keys), " comparisons; ",
        "BH applied once over ", N_TESTS, " tests")

# --- which correlations appear in the main figure ----------------------------
# Membership is by comparison name, so the main figure is literally a row subset
# of the table written below and cannot diverge from it.
main_cmp <- character(0)
if (nzchar(main_manifest) && file.exists(main_manifest)) {
  mm <- suppressMessages(read_tsv(main_manifest, show_col_types = FALSE))
  if (nrow(mm) > 0) main_cmp <- unique(mm$comparison)
}
stats_df$in_main_figure <- stats_df$comparison %in% main_cmp
message("  main figure: ", length(unique(stats_df$comparison[stats_df$in_main_figure])),
        " of ", nrow(keys), " comparisons")

fmt_p <- function(p) {
  ifelse(is.na(p), "p n/a",
  ifelse(p <= 0,   paste0("p < ", formatC(.Machine$double.xmin, format = "e", digits = 0)),
  ifelse(p < 0.001, paste0("p = ", formatC(p, format = "e", digits = 1)),
                    paste0("p = ", formatC(p, format = "f", digits = 3)))))
}
stats_df <- stats_df %>%
  mutate(lab_r = ifelse(is.na(pearson_r), "r n/a", sprintf("r = %.2f", pearson_r)),
         lab_p = fmt_p(p_adj),
         # Labels alternate between the far end of the bar and the zero
         # baseline, so neighbouring labels in a cluster never sit at the same
         # height. Parity of the feature order gives a stable pattern:
         # GC outer, ATAC inner, DNase outer.
         lab_outer = as.integer(Feature) %% 2 == 1)

# The coefficient is the number being read, so it is larger and bold and sits
# nearer the bar; the p-value is secondary and sits beyond it.
SIZE_R <- fs(2.25); SIZE_P <- fs(1.8)
vj_near <- function(up) ifelse(up, -0.25, 1.15)    # first line, closest to the anchor
vj_far  <- function(up) ifelse(up, -1.55, 2.45)    # second line, beyond it

# r_squared is kept in the table even though no figure plots it: it is simply
# r^2 and reviewers ask for it, but reporting it alongside the sign avoids the
# ambiguity that a variance-explained figure has on its own.
write_tsv(stats_df %>% arrange(col_order, Feature) %>%
            select(epitope, cell_line, region_set, comparison, cnt_label, chip_label,
                   Feature, n_bins, pearson_r, r_squared, p_value, p_adj, in_main_figure),
          out_tsv)
message("  wrote ", out_tsv)

# --- axis limits, computed once and shared by both figures -------------------
# Each end is derived from the data it has to hold, rather than the negative end
# mirroring the positive one. A symmetric axis wastes half the panel whenever
# the correlations lean one way -- with an all-positive set the bars occupied
# about 30% of the panel height, the rest being empty space below zero that no
# bar could ever reach.
#
# Zero always stays inside the range: the bars rise from the zero baseline, and
# the "inner" labels are deliberately drawn just across it, so an all-positive
# set still needs room below zero and an all-negative set above it.
obs_hi <- suppressWarnings(max(stats_df$pearson_r, na.rm = TRUE))
obs_lo <- suppressWarnings(min(stats_df$pearson_r, na.rm = TRUE))
if (!is.finite(obs_hi)) obs_hi <- 0
if (!is.finite(obs_lo)) obs_lo <- 0
r_hi <- if (obs_hi > 0) obs_hi * 1.10 else 0      # 10% headroom above the data
r_lo <- if (obs_lo < 0) obs_lo * 1.10 else 0      # 10% headroom below the data

# Labels sit beyond the bar ends and across the baseline, so a further fixed
# allowance is added to BOTH ends. It is a fraction of the plotted span rather
# than a multiplier on each limit, because a multiplier adds nothing to an end
# that sits at zero -- which is exactly the end the inner labels overflow.
span <- r_hi - r_lo
if (!is.finite(span) || span <= 0) span <- 1
lab_pad <- if (show_stats == 1) span * 0.18 else 0
Y_LIM <- c(r_lo - lab_pad, r_hi + lab_pad)
message(sprintf("  axis: [%.3f, %.3f] from observed r [%.3f, %.3f] + 10%% headroom%s; shared by both figures",
                Y_LIM[1], Y_LIM[2], obs_lo, obs_hi,
                if (show_stats == 1) " and label allowance" else ""))

# --- axis labels -------------------------------------------------------------
# Inside an epitope x cell-line facet the "CnT_H3K4me3_K562_" prefix is the same
# for every bar, so it is dropped; what remains is the study identifier, which
# is what distinguishes the comparisons. Both sides are named because the same
# CUT&Tag sample appears against more than one ChIP-seq source.
short <- function(x) sub("^[^_]+_[^_]+_[^_]+_", "", x)
stats_df <- stats_df %>%
  mutate(x_lab = paste0(short(cnt_label), "\nvs ", short(chip_label)))
lev <- stats_df %>% distinct(col_order, x_lab) %>% arrange(col_order) %>% pull(x_lab)
stats_df$x_lab <- factor(stats_df$x_lab, levels = unique(lev))
stats_df$epitope <- factor(stats_df$epitope,
                           levels = c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3"))

base_theme <- theme_minimal(base_size = fs(9)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = fs(6)),
        axis.text.y = element_text(colour = "black", size = fs(9)),
        axis.title = element_text(face = "bold", size = fs(8), lineheight = 1.15),
        plot.title = element_text(face = "bold", size = fs(12)),
        plot.subtitle = element_text(face = "italic", size = fs(9), colour = "grey30"),
        plot.caption = element_text(face = "italic", size = fs(6), hjust = 0,
                                    colour = "grey30", lineheight = 1.3),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_markdown(colour = "black"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = fs(9), face = "bold"))

# --- one renderer, called twice ----------------------------------------------
# The main and supplementary figures differ ONLY in which rows are passed in.
# Every number on them -- coefficient, adjusted p-value, axis limit -- was
# computed before this point, over the full set.
render_direction <- function(d, out_pdf, scope_txt) {
  d <- d %>% filter(!is.na(pearson_r))
  if (nrow(d) == 0) { message("  no plottable correlations for ", out_pdf); return(invisible(NULL)) }
  # Unused levels would otherwise leave empty x slots and empty facet strips.
  d <- d %>% mutate(x_lab = droplevels(x_lab), epitope = droplevels(epitope))

  epi_levels <- levels(d$epitope)
  ep_lab <- setNames(paste0("<span style='font-size:", fs(11), "pt'>**", epi_levels, "**</span>"),
                     epi_levels)
  facet_layer <- facet_nested(. ~ epitope + cell_line, scales = "free_x", space = "free_x",
                              nest_line = element_line(colour = "grey15", linewidth = 0.7),
                              labeller = labeller(epitope = ep_lab))

  n_cmp  <- length(unique(d$comparison))
  pdf_w  <- max(16, n_cmp * 1.7 + 6)

  # Bars run in both directions from zero. An "outer" label goes past the bar
  # end; an "inner" label goes just across the baseline on the opposite side,
  # which is always free because that dodge slot's own bar is on the other side.
  d <- d %>% mutate(lab_anchor = ifelse(lab_outer, pearson_r, 0),
                    lab_up     = ifelse(lab_outer, pearson_r >= 0, pearson_r < 0))

  p <- ggplot(d, aes(x = x_lab, y = pearson_r, fill = Feature)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.8), colour = "black", linewidth = 0.3, width = 0.7) +
    { if (show_stats == 1)
        geom_text(aes(y = lab_anchor, label = lab_r, vjust = vj_near(lab_up)),
                  position = position_dodge(width = 0.8), size = SIZE_R,
                  fontface = "bold", colour = "grey10") } +
    { if (show_stats == 1)
        geom_text(aes(y = lab_anchor, label = lab_p, vjust = vj_far(lab_up)),
                  position = position_dodge(width = 0.8), size = SIZE_P, colour = "grey30") } +
    scale_fill_manual(values = FEATURE_COLORS, drop = FALSE) +
    coord_cartesian(ylim = Y_LIM) +
    facet_layer +
    # The quantity lives in the title, so the axis can stay short. The sign
    # convention moves to the subtitle rather than being dropped: without it a
    # reader cannot tell which method a positive bar favours.
    labs(title = "Correlation between genomic feature and log2(CUT&Tag / ChIP-seq)",
         subtitle = sprintf(paste0("Positive r: feature higher in CUT&Tag-enriched regions   |   ",
                                   "negative r: feature higher in ChIP-seq-enriched regions   |   %s   |   %s"),
                            scope_txt, region_note),
         y = "Pearson r", x = NULL,
         caption = paste0(
           "log2FC is log2(CUT&Tag / ChIP-seq), so POSITIVE r means the feature is higher in regions enriched for CUT&Tag and NEGATIVE r\n",
           "means it is higher in regions enriched for ChIP-seq. Correlations use every tested region, not only those called differential.\n",
           "Bar labels give r (bold) and its p-value from cor.test (H0: rho = 0), Benjamini-Hochberg adjusted ONCE across all ", N_TESTS,
           " correlations\n",
           "in the complete comparison set, so a given correlation carries the same adjusted p-value in the main and supplementary figures.\n",
           "With tens of thousands of regions per test almost every correlation is significant, so the p-value is a floor on chance rather\n",
           "than evidence of a large effect; the size and sign of r are what carry the result. The three features are correlated with one\n",
           "another across the genome, so their coefficients are marginal associations rather than independent contributions, and Pearson r\n",
           "captures only the linear component of each relationship. The y-axis limit is shared by both figures, so a bar is the same height\n",
           "in each.")) +
    base_theme

  ggsave(out_pdf, p, width = pdf_w, height = 11, device = "pdf", bg = "white", limitsize = FALSE)
  message("  wrote ", out_pdf, "  (", n_cmp, " comparisons, ", nrow(d), " bars)")
}

render_direction(stats_df %>% filter(in_main_figure), out_pdf_main,
                 "one representative comparison per mark")
render_direction(stats_df, out_pdf_supp, "all comparisons")
RFI

# =============================================================================
# RUN
# =============================================================================
export FONT_SCALE

REGION_NOTE="promoter marks in ${PROM_BIN_SIZE:-500} bp bins, broad marks in 3 kb bins"

n_supp=$(( $(wc -l < "$SUPP_MANIFEST") - 1 ))
echo "--- ${n_supp} comparisons ---"

Rscript "$RS_DIR/feature_importance.R" \
    "$SUPP_MANIFEST" "$MAIN_MANIFEST" \
    "$OUT_DIR/feature_importance_direction_supp.pdf" \
    "$OUT_DIR/feature_importance_direction_main.pdf" \
    "$OUT_DIR/feature_importance_correlations.tsv" \
    "$REGION_NOTE" "$PARAMS" "$SHOW_BAR_STATS"

echo "========================================================="
echo " Done. Output: $OUT_DIR"
echo "   main figure : feature_importance_direction_main.pdf"
echo "   supp figure : feature_importance_direction_supp.pdf"
echo "   table       : feature_importance_correlations.tsv"
echo "========================================================="