#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CUT&Tag vs ENCODE ChIP-seq: GC content of the peaks that each method finds
# on its own.
#
# THE QUESTION. Abbasova et al. profiled H3K27ac and H3K27me3 in K562 by
# CUT&Tag. ENCODE profiled the same marks in the same cell line by ChIP-seq.
# Splitting the two peak sets into three groups --
#
#     CUT&Tag unique   peaks called by CUT&Tag with no ChIP-seq peak there
#     Shared           peaks called by both
#     ChIP-seq unique  peaks called by ChIP-seq with no CUT&Tag peak there
#
# -- lets us ask whether sequence composition distinguishes the regions each
# method finds alone.
#
# WHAT THIS SCRIPT DOES, END TO END. Nothing needs to be downloaded by hand.
#
#     PART 1  Reference signal track : hg19 GC content bigWig
#     PART 2  Peak sets              : CUT&Tag from CyVerse, ChIP-seq from ENCODE
#     PART 3  Overlap                : the three subsets above, per sample
#     PART 4  Signal extraction      : mean GC per peak
#     PART 5  Statistics and figure
#
# EVERY STEP IS IDEMPOTENT. Files that already exist are not re-downloaded or
# recomputed, so the script can be re-run safely after a failure or after you
# change only the plotting section.
#
# GENOME BUILD: hg19 throughout. The CUT&Tag peaks, the ENCODE peaks and the GC
# track are all hg19. Do not mix in an hg38 file.
#
# Usage:
#   conda activate <env with bedtools + deeptools + R>
#   bash cnt_vs_encode_gc.sh
# =============================================================================

# =============================================================================
# USER CONFIG
# =============================================================================
BASE="${BASE:-/home/emodolo/gpfs/2026_modolo_et_al/abbasova_reanalysis/output}"

PEAKS_DIR="$BASE/abbasova_peaks"                 # CUT&Tag peaks (downloaded here)
ENCODE_DIR="$BASE/encode_peaks"                  # ENCODE ChIP-seq peaks
REF_DIR="$BASE/reference_signal"                 # GC bigWig
OUT="$BASE/GC_content_analysis/cnt_vs_chip_gc"

# `export` is required: the R block at the bottom is a separate process and
# reads OUT through Sys.getenv(). A plain shell variable would not reach it.
export OUT

THREADS="${THREADS:-16}"

# --- Stage toggles -----------------------------------------------------------
# Set a stage to 0 to skip it. Useful when only the figure needs redrawing:
# DO_DOWNLOAD=0 DO_OVERLAP=0 bash cnt_vs_encode_gc.sh
DO_DOWNLOAD="${DO_DOWNLOAD:-1}"    # PART 1 + PART 2
DO_OVERLAP="${DO_OVERLAP:-1}"      # PART 3 + PART 4
DO_FIGURE="${DO_FIGURE:-1}"        # PART 5

# --- Reference signal track --------------------------------------------------
# UCSC ships GC content as a variable-step wiggle, not a bigWig, so it has to
# be converted once with the UCSC tools. The result is 5 bp resolution.
GC_WIG_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wigVarStep.gz"
GC_WIG="$REF_DIR/hg19.gc5Base.wigVarStep.gz"
GC_CHROM_SIZES="$REF_DIR/hg19.chrom.sizes"
GC_BW="$REF_DIR/hg19.gc5Base.bw"

# --- ENCODE ChIP-seq peak sets -----------------------------------------------
# Downloaded gzipped, decompressed, and given descriptive names.
ENCODE_H3K27AC_URL="https://www.encodeproject.org/files/ENCFF044JNJ/@@download/ENCFF044JNJ.bed.gz"
ENCODE_H3K27AC="$ENCODE_DIR/chip_H3K27ac_ENCODE_2017_hg19_replicated_narrowPeaks_ENCFF044JNJ.bed"
ENCODE_H3K27ME3_URL="https://www.encodeproject.org/files/ENCFF001SZF/@@download/ENCFF001SZF.bed.gz"
ENCODE_H3K27ME3="$ENCODE_DIR/chip_H3K27me3_ENCODE_2010_hg19_replicated_broadPeaks_ENCFF001SZF.bed"

# --- CUT&Tag peak sets (Abbasova et al., hosted on CyVerse) ------------------
# The server returns a browsable HTML index. The script reads that index and
# pulls out the filenames itself, so you do not have to know them in advance.
# %26 is the URL encoding of "&" in "H3K27_CUT&Tag_Benchmark".
CYVERSE_BASE="https://data.cyverse.org/dav-anon/iplant/home/paulinaurbana/H3K27_CUT%26Tag_Benchmark/"

# Only files whose name ends in one of these are wanted: the SEACR and MACS2
# peak calls. Everything else in that directory is ignored.
#
# broadPeak is in this list because MACS2 was run in broad mode for H3K27me3,
# writing "..._q1e-5_peaks.broadPeak". Without it those four files were never
# downloaded and the H3K27me3 MACS2 rows were silently missing from the
# analysis, while the SEACR H3K27me3 rows were present -- which is what makes
# that kind of omission easy to miss.
PEAK_SUFFIX_RE='_peaks\.(stringent\.bed|narrowPeak|broadPeak)'

# Sample-ID to readable-prefix map. The downloaded filenames start with an
# opaque sequencing ID (IGF128036); this turns them into names that state the
# antibody and replicate, which the analysis then parses.
#
# IMPORTANT: the analysis extracts sample_id as fields 3 and 4 of the final
# filename split on "_" (e.g. CnT_H3K27ac_ab177_R1_IGF... -> "ab177_R1"), so
# every prefix must keep the shape  CnT_<epitope>_<antibody>_<replicate>.
declare -A ID_MAP=(
    ["IGF128036"]="CnT_H3K27ac_ab177_R1"
    ["IGF128037"]="CnT_H3K27ac_ab177_R2"
    ["IGF128038"]="CnT_H3K27me3_ab9733_R1"
    ["IGF128039"]="CnT_H3K27me3_ab9733_R2"
    ["IGF128040"]="CnT_H3K27ac_ab4729_R1"
    ["IGF128041"]="CnT_H3K27ac_ab4729_R2"
    ["IGF128042"]="CnT_H3K27ac_diag_R1"
    ["IGF128043"]="CnT_H3K27ac_diag_R2"
)

# --- Output files ------------------------------------------------------------
COUNTS="$OUT/peak_counts.tsv"
TMP="$OUT/_tmp"
mkdir -p "$PEAKS_DIR" "$ENCODE_DIR" "$REF_DIR" \
         "$OUT/intersections" "$OUT/gc_matrices" "$OUT/plots" "$TMP"

# =============================================================================
# PART 0: PREFLIGHT
# =============================================================================
echo "========================================================="
echo " CUT&Tag vs ENCODE: GC content of method-unique peaks"
echo " Peaks     : $PEAKS_DIR"
echo " Reference : $REF_DIR"
echo " Output    : $OUT"
echo "========================================================="

miss=0
need_tool() {   # tool name, human explanation
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: $1 not in PATH ($2)" >&2; miss=1; }
}
need_tool Rscript  "R, for the statistics and figure"
if [[ "$DO_DOWNLOAD" == "1" ]]; then
    need_tool wget    "downloads"
    need_tool gunzip  "decompressing ENCODE bed.gz"
    # The UCSC converters are only needed the first time, to build the GC bigWig.
    if [[ ! -s "$GC_BW" ]]; then
        need_tool fetchChromSizes "UCSC tool, only needed to build the GC bigWig once"
        need_tool wigToBigWig     "UCSC tool, only needed to build the GC bigWig once"
    fi
fi
if [[ "$DO_OVERLAP" == "1" ]]; then
    need_tool bedtools            "peak intersections"
    need_tool multiBigwigSummary  "deepTools, for the per-peak GC mean"
fi
[[ $miss -eq 0 ]] || { echo "Preflight failed." >&2; exit 1; }

# Downloads to a temporary name first and only moves it into place on success,
# so an interrupted download can never leave a truncated file that later steps
# would treat as valid because it merely exists.
fetch() {   # url, destination
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then
        echo "    cached: $(basename "$dest")"
        return 0
    fi
    echo "    downloading: $(basename "$dest")"
    wget -q --show-progress -c -O "${dest}.part" "$url" || {
        echo "    ERROR: download failed: $url" >&2; rm -f "${dest}.part"; return 1; }
    mv "${dest}.part" "$dest"
}

# =============================================================================
# PART 1: REFERENCE SIGNAL TRACK
# =============================================================================
if [[ "$DO_DOWNLOAD" == "1" ]]; then
echo "--- PART 1: GC content track ---"

if [[ -s "$GC_BW" ]]; then
    echo "    cached: $(basename "$GC_BW")"
else
    fetch "$GC_WIG_URL" "$GC_WIG"
    if [[ ! -s "$GC_CHROM_SIZES" ]]; then
        echo "    fetching hg19 chromosome sizes"
        fetchChromSizes hg19 > "${GC_CHROM_SIZES}.part" 2>/dev/null
        mv "${GC_CHROM_SIZES}.part" "$GC_CHROM_SIZES"
    fi
    echo "    converting wiggle to bigWig (several minutes, memory-hungry)"
    wigToBigWig "$GC_WIG" "$GC_CHROM_SIZES" "${GC_BW}.part"
    mv "${GC_BW}.part" "$GC_BW"
    echo "    built: $(basename "$GC_BW")"
fi

# =============================================================================
# PART 2: PEAK SETS
# =============================================================================
echo "--- PART 2: peak sets ---"

# --- ENCODE ChIP-seq ---------------------------------------------------------
get_encode_bed() {   # url, final bed path
    local url="$1" bed="$2" gz
    if [[ -s "$bed" ]]; then echo "    cached: $(basename "$bed")"; return 0; fi
    gz="${bed}.gz"
    fetch "$url" "$gz"
    # -c writes to stdout and leaves the .gz in place, so a re-run still has
    # the original archive to fall back on.
    gunzip -c "$gz" > "${bed}.part"
    mv "${bed}.part" "$bed"
    echo "    decompressed: $(basename "$bed")"
}
get_encode_bed "$ENCODE_H3K27AC_URL"  "$ENCODE_H3K27AC"
get_encode_bed "$ENCODE_H3K27ME3_URL" "$ENCODE_H3K27ME3"

# --- CUT&Tag: rename anything already downloaded by hand ---------------------
# Handles the case where peaks were fetched previously and still carry their
# raw IGF names. Files already renamed start with "CnT_" and are left alone.
shopt -s nullglob
for f in "$PEAKS_DIR"/IGF*; do
    [[ -f "$f" ]] || continue
    bn=$(basename "$f")
    sid="${bn%%_*}"                       # strip from the first "_" onward
    if [[ -n "${ID_MAP[$sid]:-}" ]]; then
        mv "$f" "$PEAKS_DIR/${ID_MAP[$sid]}_${bn}"
        echo "    renamed existing: $bn -> ${ID_MAP[$sid]}_${bn}"
    else
        echo "    WARNING: no prefix mapping for '$sid' (file: $bn)" >&2
    fi
done
shopt -u nullglob

# --- CUT&Tag: read the CyVerse index and download what is missing ------------
INDEX="$TMP/cyverse_index.html"
if [[ ! -s "$INDEX" ]]; then
    echo "    reading CyVerse directory listing"
    wget -q -O "${INDEX}.part" "$CYVERSE_BASE" || {
        echo "ERROR: could not read the CyVerse listing at:" >&2
        echo "       $CYVERSE_BASE" >&2
        echo "       Open it in a browser to confirm it is still up." >&2
        rm -f "${INDEX}.part"; exit 1; }
    mv "${INDEX}.part" "$INDEX"
fi

n_dl=0
for sid in "${!ID_MAP[@]}"; do
    prefix="${ID_MAP[$sid]}"
    # Pull every filename in the listing that starts with this sample ID and
    # ends in a wanted peak suffix. The character class includes "-" so that
    # MACS2 names carrying a threshold ("..._macs2_q1e-5_peaks.broadPeak") are
    # matched in full. `|| true` keeps a no-match grep (exit 1) from killing
    # the script under `set -e`.
    mapfile -t hits < <(grep -oE "${sid}[A-Za-z0-9._-]*${PEAK_SUFFIX_RE}" "$INDEX" | sort -u || true)
    if [[ ${#hits[@]} -eq 0 ]]; then
        echo "    WARNING: no peak files found in the listing for $sid" >&2
        continue
    fi
    for fname in "${hits[@]}"; do
        dest="$PEAKS_DIR/${prefix}_${fname}"
        [[ -s "$dest" ]] && continue
        # An `if` rather than `fetch ... && n_dl=...`: under `set -e` a failing
        # `&&` chain would abort the whole run, so one dead link would throw
        # away the other downloads. This warns and carries on instead.
        if fetch "${CYVERSE_BASE}${fname}" "$dest"; then
            n_dl=$((n_dl + 1))
        else
            echo "    WARNING: skipping $fname" >&2
        fi
    done
done
echo "    CUT&Tag peak files newly downloaded: $n_dl"
fi   # end DO_DOWNLOAD

# The analysis cannot start without these, whether they came from this run or
# an earlier one.
for f in "$GC_BW" "$ENCODE_H3K27AC" "$ENCODE_H3K27ME3"; do
    [[ -s "$f" ]] || { echo "ERROR: required file missing: $f" >&2; exit 1; }
done

# =============================================================================
# PART 3 + 4: OVERLAP SUBSETS, AND PER-PEAK GC
# =============================================================================
if [[ "$DO_OVERLAP" == "1" ]]; then
echo "--- PART 3 + 4: overlaps and GC extraction ---"

printf "caller\tepitope\tsample_id\tn_cnt_total\tn_chip_total\tn_cnt_unique\tn_chip_unique\tn_cnt_shared\tn_chip_shared\n" > "$COUNTS"

shopt -s nullglob
for cnt_file in "$PEAKS_DIR"/*.narrowPeak "$PEAKS_DIR"/*.broadPeak "$PEAKS_DIR"/*.stringent.bed; do
    bname=$(basename "$cnt_file")

    # Route each CUT&Tag file to the matching ENCODE ChIP-seq peak set.
    if   [[ "$bname" == *H3K27ac*  ]]; then epitope="H3K27ac";  encode_file="$ENCODE_H3K27AC"
    elif [[ "$bname" == *H3K27me3* ]]; then epitope="H3K27me3"; encode_file="$ENCODE_H3K27ME3"
    else continue
    fi

    # narrowPeak (H3K27ac) and broadPeak (H3K27me3) are both MACS2 output; the
    # mode differs, the caller does not.
    case "$bname" in
      *.narrowPeak|*.broadPeak) caller="macs2";;
      *.stringent.bed)          caller="seacr";;
      *) continue;;
    esac

    # Fields 3 and 4 of the renamed file: antibody and replicate, e.g. ab177_R1
    sample_id=$(echo "$bname" | awk -F'_' '{print $3"_"$4}')

    tag="${caller}__${epitope}__${sample_id}"
    prefix="$OUT/intersections/$tag"
    echo "  >>> $tag"

    # Reduce both files to plain, sorted, de-duplicated 3-column BED. Dropping
    # the score and summit columns means the two sides are treated identically
    # regardless of which caller produced them.
    awk -v OFS='\t' '$1 ~ /^chr/ && $3>$2 {print $1,$2,$3}' "$cnt_file"    | sort -k1,1 -k2,2n -u > "${prefix}__cnt_all.bed"
    awk -v OFS='\t' '$1 ~ /^chr/ && $3>$2 {print $1,$2,$3}' "$encode_file" | sort -k1,1 -k2,2n -u > "${prefix}__chip_all.bed"

    # -v keeps intervals with NO overlap; -u keeps each interval once if it has
    # at least one overlap. "shared" is defined on ENCODE boundaries, so the
    # shared regions are a subset of the ChIP-seq peaks.
    bedtools intersect -a "${prefix}__cnt_all.bed"  -b "${prefix}__chip_all.bed" -v > "${prefix}__cnt_unique.bed"
    bedtools intersect -a "${prefix}__chip_all.bed" -b "${prefix}__cnt_all.bed"  -v > "${prefix}__chip_unique.bed"
    bedtools intersect -a "${prefix}__chip_all.bed" -b "${prefix}__cnt_all.bed"  -u > "${prefix}__shared.bed"

    n_cnt_total=$(   wc -l <"${prefix}__cnt_all.bed")
    n_chip_total=$(  wc -l <"${prefix}__chip_all.bed")
    n_cnt_unique=$(  wc -l <"${prefix}__cnt_unique.bed")
    n_chip_unique=$( wc -l <"${prefix}__chip_unique.bed")
    n_chip_shared=$( wc -l <"${prefix}__shared.bed")
    n_cnt_shared=$(  bedtools intersect -a "${prefix}__cnt_all.bed" -b "${prefix}__chip_all.bed" -u | wc -l)

    printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" \
      "$caller" "$epitope" "$sample_id" \
      "$n_cnt_total" "$n_chip_total" "$n_cnt_unique" "$n_chip_unique" "$n_cnt_shared" "$n_chip_shared" \
      >> "$COUNTS"

    # An existing .tab is reused, because this is the slow step. That cache is
    # keyed only on the filename, so if the underlying peaks change (a new
    # download, a different ENCODE file) the old values would be reused:
    #   rm -rf "$OUT/gc_matrices"   to force a full recompute.
    for subset in cnt_unique chip_unique shared; do
      sub_bed="${prefix}__${subset}.bed"
      out_tab="$OUT/gc_matrices/${tag}__${subset}.tab"
      if [[ -s "$out_tab" ]]; then
        continue
      elif [[ -s "$sub_bed" ]]; then
        multiBigwigSummary BED-file \
          --bwfiles "$GC_BW" \
          --labels GC \
          --BED "$sub_bed" \
          -o "$OUT/gc_matrices/${tag}__${subset}.npz" \
          --outRawCounts "$out_tab" \
          -p "$THREADS" > /dev/null
      else
        # An empty subset still needs a placeholder file so the R side does not
        # trip over a missing path.
        printf "#chr\tstart\tend\tGC\n" > "$out_tab"
      fi
    done
done
shopt -u nullglob
echo "  overlaps and GC extraction done."
fi   # end DO_OVERLAP

[[ -s "$COUNTS" ]] || { echo "ERROR: $COUNTS is missing or empty; run with DO_OVERLAP=1" >&2; exit 1; }

# =============================================================================
# PART 5: STATISTICS AND FIGURE
# =============================================================================
if [[ "$DO_FIGURE" != "1" ]]; then
    echo "DO_FIGURE=0, stopping before the figure."
    exit 0
fi

echo "--- PART 5: statistics and figure ---"

Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

OUT    <- Sys.getenv("OUT")
gc_dir <- file.path(OUT, "gc_matrices")

counts <- read_tsv(file.path(OUT, "peak_counts.tsv"), show_col_types = FALSE)

# --- Read the per-peak GC tables --------------------------------------------
# Filenames encode the metadata, separated by "__":
#   caller __ epitope __ sample_id __ subset .tab
gc_files <- list.files(gc_dir, pattern = "\\.tab$", full.names = TRUE)

read_gc <- function(f) {
  bn    <- str_remove(basename(f), "\\.tab$")
  parts <- str_split(bn, "__", simplify = TRUE)
  tab <- tryCatch(
    read.table(f, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  colnames(tab) <- c("chr", "start", "end", "GC")
  tibble(
    caller    = parts[1, 1],
    epitope   = parts[1, 2],
    sample_id = parts[1, 3],
    subset    = parts[1, 4],
    # deepTools writes "nan" for intervals with no data, which makes the whole
    # column character. Blanket coercion back to numeric handles that.
    GC = suppressWarnings(as.numeric(tab$GC))
  )
}
gc_df <- map_dfr(gc_files, read_gc) %>% filter(!is.na(GC))
stopifnot(nrow(gc_df) > 0)

# UCSC gc5Base is a percentage, but guard against a fractional 0-1 track being
# swapped in later: rescale only if the observed maximum says it is a fraction.
if (max(gc_df$GC, na.rm = TRUE) <= 1) gc_df$GC <- gc_df$GC * 100

# --- Sample ordering: H3K27ac first, then H3K27me3 --------------------------
label_lookup <- counts %>%
  distinct(epitope, sample_id) %>%
  arrange(factor(epitope, levels = c("H3K27ac", "H3K27me3")), sample_id) %>%
  mutate(label = paste0(sample_id, " (", epitope, ")"))

sample_order <- label_lookup$label
gc_df <- gc_df %>% left_join(label_lookup, by = c("epitope", "sample_id"))
# rev() because ggplot draws the first factor level at the BOTTOM of a
# horizontal plot; reversing puts the first sample at the top.
gc_df$label <- factor(gc_df$label, levels = rev(sample_order))

subset_names <- c(cnt_unique = "CUT&Tag unique", shared = "Shared", chip_unique = "ChIP-seq unique")
subset_cols  <- c("CUT&Tag unique"  = "#ff3b3b",
                  "Shared"          = "#9e9e9e",
                  "ChIP-seq unique" = "#3936ff")
gc_df$subset_label <- factor(subset_names[gc_df$subset],
                             levels = c("CUT&Tag unique", "Shared", "ChIP-seq unique"))

# --- Per-subset descriptive statistics --------------------------------------
subset_summary <- gc_df %>%
  group_by(caller, epitope, sample_id, label, subset, subset_label) %>%
  summarise(n         = n(),
            mean_gc   = mean(GC),
            median_gc = median(GC),
            sd_gc     = sd(GC),
            .groups   = "drop")
write_tsv(subset_summary, file.path(OUT, "gc_subset_summary.tsv"))

# --- Pairwise Wilcoxon rank-sum tests ---------------------------------------
# Rank-based, so it does not assume normality.
stat_one <- function(d) {
  pairs <- list(c("cnt_unique", "chip_unique"),
                c("cnt_unique", "shared"),
                c("shared",     "chip_unique"))
  map_dfr(pairs, function(p) {
    x <- d$GC[d$subset == p[1]]
    y <- d$GC[d$subset == p[2]]
    out <- tibble(group1 = p[1], group2 = p[2],
                  n1 = length(x), n2 = length(y),
                  median1 = if (length(x)) median(x) else NA_real_,
                  median2 = if (length(y)) median(y) else NA_real_,
                  p_value = NA_real_)
    if (length(x) >= 3 && length(y) >= 3) {
      out$p_value <- suppressWarnings(wilcox.test(x, y))$p.value
    }
    out
  })
}

# Benjamini-Hochberg correction is applied WITHIN each caller, i.e. per panel
# of the figure, not across the whole figure.
stats <- gc_df %>%
  group_by(caller, epitope, sample_id, label) %>%
  group_modify(~ stat_one(.x)) %>%
  ungroup() %>%
  group_by(caller) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Numeric p-values rather than asterisks: readers can see how far below the
# threshold a comparison sits instead of everything collapsing to "***".
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "n.d.",
    p == 0     ~ "p < 1e-300",             # underflow at very large n
    p >= 0.001 ~ sprintf("p = %.3f", p),
    TRUE       ~ sprintf("p = %.1e", p)
  )
}
stats$p_label <- fmt_p(stats$p_adj)
write_tsv(stats, file.path(OUT, "gc_pairwise_stats.tsv"))

# Only the CUT&Tag-unique vs ChIP-unique contrast is drawn on the figure; the
# other two pairs stay in the stats table.
anno <- stats %>%
  filter(group1 == "cnt_unique", group2 == "chip_unique") %>%
  mutate(label = factor(as.character(label), levels = rev(sample_order)),
         y_num = as.numeric(label))

# --- Plot constants ----------------------------------------------------------
BASE_SIZE  <- 20
DODGE_W    <- 0.85
N_SUBSETS  <- 3
# Half the distance between the outermost two dodged violins: the bracket has
# to span from the CUT&Tag-unique violin to the ChIP-unique violin.
DODGE_HALF <- DODGE_W/2 - DODGE_W/(2 * N_SUBSETS)
MEDIAN_COL <- "#1b7837"   # green, for the median annotations

common_theme <- theme_minimal(base_size = BASE_SIZE) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title         = element_text(face = "bold", size = 20, hjust = 0))

build_violin <- function(data_gc, data_anno, caller_use, title, show_y_axis = FALSE) {
  d  <- data_gc   %>% filter(caller == caller_use)
  an <- data_anno %>% filter(caller == caller_use)
  if (nrow(d) == 0) return(patchwork::plot_spacer())

  # position_dodge places subset level i at:
  #     y_num + (i - (N+1)/2) * DODGE_W/N
  # Reproducing that formula by hand is what lets the median labels and the
  # significance brackets line up with the violins they belong to.
  medians <- d %>%
    group_by(label, subset_label) %>%
    summarise(med = median(GC, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_num = as.numeric(label) +
             (as.integer(subset_label) - (N_SUBSETS + 1)/2) * (DODGE_W / N_SUBSETS))

  # All x positions are expressed as fractions of the data range, so the layout
  # does not need retuning if the GC range shifts.
  rng  <- range(d$GC, na.rm = TRUE)
  span <- diff(rng); if (!is.finite(span) || span == 0) span <- 1

  bracket_x <- rng[2] + span * 0.05
  tick_len  <- span * 0.02
  p_x       <- bracket_x + span * 0.03   # numeric p-value, right of the bracket
  med_x     <- rng[2] + span * 0.42      # green median column, right of the p-values
  xlim_hi   <- rng[2] + span * 0.72      # room for both annotation columns
  xlim_lo   <- rng[1] - span * 0.03

  p <- ggplot(d, aes(x = GC, y = label, fill = subset_label)) +
    geom_violin(orientation = "y", scale = "width",
                position = position_dodge(width = DODGE_W),
                width = 0.8, alpha = 0.85, colour = "grey30", linewidth = 0.3) +
    # The explicit group is required: without it the boxplots collapse onto a
    # single position instead of following the violins' dodge.
    geom_boxplot(aes(group = interaction(label, subset_label)),
                 orientation = "y", width = 0.12, fill = "white",
                 colour = "grey25", outlier.shape = NA,
                 position = position_dodge(width = DODGE_W)) +
    # Significance bracket: vertical spine plus two inward ticks.
    geom_segment(data = an, aes(y = y_num - DODGE_HALF, yend = y_num + DODGE_HALF),
                 x = bracket_x, xend = bracket_x,
                 inherit.aes = FALSE, linewidth = 0.5, colour = "grey20") +
    geom_segment(data = an, aes(y = y_num - DODGE_HALF, yend = y_num - DODGE_HALF),
                 x = bracket_x, xend = bracket_x - tick_len,
                 inherit.aes = FALSE, linewidth = 0.5, colour = "grey20") +
    geom_segment(data = an, aes(y = y_num + DODGE_HALF, yend = y_num + DODGE_HALF),
                 x = bracket_x, xend = bracket_x - tick_len,
                 inherit.aes = FALSE, linewidth = 0.5, colour = "grey20") +
    geom_text(data = an, aes(y = y_num, label = p_label),
              x = p_x, hjust = 0, vjust = 0.5, size = 4.4,
              inherit.aes = FALSE) +
    # One green median per violin, on the violin's own dodge offset.
    geom_text(data = medians, aes(y = y_num, label = sprintf("%.1f", med)),
              x = med_x, hjust = 0, vjust = 0.5, size = 4.8,
              colour = MEDIAN_COL, fontface = "bold",
              inherit.aes = FALSE) +
    annotate("text", x = med_x, y = max(as.numeric(d$label)) + 0.72,
             label = "median GC (%)", hjust = 0, vjust = 0.5,
             size = 4.2, colour = MEDIAN_COL, fontface = "italic") +
    scale_fill_manual(values = subset_cols, name = "Peak subset") +
    coord_cartesian(xlim = c(xlim_lo, xlim_hi), clip = "off") +
    labs(title = title, x = "Per-peak mean GC content (%)", y = NULL) +
    common_theme

  if (!show_y_axis) {
    p <- p + theme(axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank(),
                   plot.margin  = margin(8, 25, 8, 8))
  } else {
    p <- p + theme(plot.margin = margin(8, 25, 8, 8))
  }
  p
}

# =========================================================================
# FIGURE: all replicates, one panel per caller
# =========================================================================
macs_gc <- build_violin(gc_df, anno, "macs2", "MACS2: GC content",  show_y_axis = TRUE)
seac_gc <- build_violin(gc_df, anno, "seacr", "SEACR: GC content",  show_y_axis = FALSE)

fig <- macs_gc + seac_gc +
  plot_layout(ncol = 2, widths = c(1.25, 1), guides = "collect") +
  plot_annotation(
    title    = "CUT&Tag vs ENCODE ChIP-seq: GC content of method-unique peaks",
    subtitle = paste0(
      "MACS2 (left)  |  SEACR (right).  Shared subset uses ENCODE peak boundaries.  ",
      "Bracketed contrast: CUT&Tag-unique vs ChIP-seq-unique, Wilcoxon rank-sum, ",
      "Benjamini-Hochberg adjusted within each panel.  Green values: median per subset."),
    theme = theme(plot.title    = element_text(face = "bold", size = 25),
                  plot.subtitle = element_text(size = 14, colour = "grey30"))
  ) & theme(legend.position = "top",
            legend.text  = element_text(size = 20),
            legend.title = element_text(size = 20))

ggsave(file.path(OUT, "plots", "cnt_vs_chip_gc_figure.pdf"),
       fig, width = 20, height = 10, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(OUT, "plots", "cnt_vs_chip_gc_figure.png"),
       fig, width = 20, height = 10, dpi = 200, limitsize = FALSE)

cat("\n---- DONE ----\n",
    "counts table : ", file.path(OUT, "peak_counts.tsv"), "\n",
    "stats table  : ", file.path(OUT, "gc_pairwise_stats.tsv"), "\n",
    "summary table: ", file.path(OUT, "gc_subset_summary.tsv"), "\n",
    "figure (pdf) : ", file.path(OUT, "plots/cnt_vs_chip_gc_figure.pdf"), "\n",
    "figure (png) : ", file.path(OUT, "plots/cnt_vs_chip_gc_figure.png"), "\n",
    sep = "")
RSCRIPT

echo "========================================================="
echo " Done."
echo "   Peaks     : $PEAKS_DIR"
echo "   Reference : $REF_DIR"
echo "   Tables    : $OUT"
echo "   Figures   : $OUT/plots"
echo "========================================================="