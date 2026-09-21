#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# DEPTH-MATCHED CUT&Tag COMPARISON -- K562 broad marks
#
# Compares many CUT&Tag protocol variants against one another over the same
# differential region groups used elsewhere in the manuscript. Because the
# comparison is CUT&Tag vs CUT&Tag, every CUT&Tag library is first downsampled
# to a common fragment count; without that, apparent differences between
# protocols would largely be differences in sequencing depth.
#
#   1. downsample every CUT&Tag BAM to TARGET_FRAGMENTS fragments
#   2. build a raw bigWig from each downsampled BAM, exactly as the processing
#      pipeline builds its per-library tracks
#   3. one combined heatmap: H3K36me3 gene-body bins above H3K27me3 Polycomb
#      bins, six region groups
#   4. one combined violin figure, same six groups, all fifteen pairwise tests
#
# ChIP-seq is NOT downsampled. It is shown for reference at its native depth,
# and is not part of the depth-matched comparison.
#
# PAIRED-END DOWNSAMPLING. samtools --subsample decides per READ NAME, so both
# mates of a pair are kept or dropped together. Subsampling by the fraction
# 1e6/fragments therefore yields ~1e6 fragments (~2e6 reads) with no orphaned
# mates. The fragment count is the read-1 count (-f 64), not half the read
# count, so it stays correct even if the BAM is not perfectly paired.
#
# Usage:
#   conda activate chrom_diff_figures     # needs samtools + deepTools + R
#   bash 06_downsampled_broadmark_figures.sh
# =============================================================================

# --- STAGE TOGGLES -----------------------------------------------------------
DO_DOWNSAMPLE=1
DO_BIGWIGS=1
COMPUTE_MATRICES=1     # this analysis has its own matrix; nothing to reuse
PLOT_HEATMAP=1
PLOT_VIOLINS=1

CELL="K562"

# --------------------------- USER CONFIG -------------------------------------
DIFF_WORK_DIR="${DIFF_WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output}"
DIFF_DIR="$DIFF_WORK_DIR/differential"
REGIONS_DIR="$DIFF_WORK_DIR/regions"
PARAMS="$DIFF_DIR/analysis_params.tsv"

WORK_DIR="${WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output/downsampled_broadmarks}"
DS_BAM_DIR="$WORK_DIR/downsampled_bam"
DS_BW_DIR="$WORK_DIR/downsampled_bigwig"
REG_DIR="$WORK_DIR/region_groups"
MAT_DIR="$WORK_DIR/matrices"
HM_DIR="$WORK_DIR/heatmap"
VI_DIR="$WORK_DIR/violin"
RS_DIR="$WORK_DIR/Rscripts"
TMP="$WORK_DIR/_tmp"
mkdir -p "$DS_BAM_DIR" "$DS_BW_DIR" "$REG_DIR" "$MAT_DIR" "$HM_DIR" "$VI_DIR" "$RS_DIR" "$TMP"

MASTER_TSV="${MASTER_TSV:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/summary/Master_Sample_Metrics.tsv}"
GC_BW="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/GC_content/output/bigwigs/hg38_gc5Base.bw"
REF_ROOT="/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/reference_datasets"
DNASE_BW="${REF_ROOT}/DNase-seq/output/bigwigs/DNase_K562_ENCFF972GVB.bigWig"
ATAC_BW="/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/alignment/bigwig/ATAC_NA_K562_ID173_merged_CPM.bw"

# --- Downsampling -------------------------------------------------------------
TARGET_FRAGMENTS="${TARGET_FRAGMENTS:-1000000}"
DOWNSAMPLE_SEED="${DOWNSAMPLE_SEED:-42}"     # fixed, so a re-run selects the same fragments
DS_TAG="dwnsmpl1M"

# --- bigWig settings, matched to the processing pipeline ----------------------
# Paired-end libraries are extended to their true fragment span and the track is
# left RAW: after depth matching, raw coverage is what makes the samples
# comparable, and any normalisation here would undo the downsampling.
BIGWIG_BIN_SIZE="${BIGWIG_BIN_SIZE:-1}"

# --- Figure geometry ----------------------------------------------------------
BINSIZE=25
BROAD_WINDOW=1500
THREADS="${THREADS:-24}"
FONT_SCALE="${FONT_SCALE:-1.5}"
VIOLIN_FONT_SCALE="${VIOLIN_FONT_SCALE:-2.25}"
HEATMAP_COLOR_CAP_FRACTION=0.03
HEATMAP_WIDTH_MODE="fixed"
HEATMAP_WIDTH_CM="${HEATMAP_WIDTH_CM:-4}"
HEATMAP_HEIGHT_CM="${HEATMAP_HEIGHT_CM:-16}"

# --- Which differential comparison defines each mark's region groups ---------
declare -A CNT_GROUP=(
  [H3K36me3]="CnT_H3K36me3_K562_Wu_2025_rep1_2025_rep2_2025_rep3_2025_rep4_2025_rep5"
  [H3K27me3]="CnT_H3K27me3_K562_Abbasova_2025_rep1_2025_rep2"
)
declare -A CHIP_GROUP=(
  [H3K36me3]="chip_H3K36me3_K562_Bernstein_2011_rep2_2022_rep3"
  [H3K27me3]="chip_H3K27me3_K562_Bernstein_2011_rep1_2022_rep3"
)
declare -A EPI_REGION_SET=( [H3K36me3]=genebody_bins [H3K27me3]=polycomb_bins )
declare -A SET_DESC=(
  [genebody_bins]="active gene-body bins"
  [polycomb_bins]="Polycomb-repressed bins"
)

# --- Track list, in figure order ---------------------------------------------
# display label | plot_label | mode
#   ref        a bigWig path given directly, shown as-is
#   full       a sample from the master metrics table at its native depth
#   downsample a sample downsampled to TARGET_FRAGMENTS, then re-tracked
# The order here is the column order of the heatmap and the panel order of the
# violin figure; nothing downstream re-sorts it.
TRACKS=(
  "GC content|${GC_BW}|ref"
  "DNase-seq|${DNASE_BW}|ref"
  "ATAC-seq|${ATAC_BW}|ref"

  "H3K36me3 ChIP Bernstein 2022 rep3|chip_H3K36me3_K562_Bernstein_2022_rep3|full"
  "H3K36me3 CnT Wu 2025 rep4|CnT_H3K36me3_K562_Wu_2025_rep4|downsample"
  "H3K36me3 CnT KayaOkur 2020 NNPC|CnT_H3K36me3_K562_KayaOkur_2020_NNPC|downsample"
  "H3K36me3 CnT KayaOkur 2020 NCSD|CnT_H3K36me3_K562_KayaOkur_2020_NCSD|downsample"
  "H3K36me3 CnT KayaOkur 2020 FNPC|CnT_H3K36me3_K562_KayaOkur_2020_FNPC|downsample"
  "H3K36me3 CnT KayaOkur 2020 FCSD|CnT_H3K36me3_K562_KayaOkur_2020_FCSD|downsample"

  "H3K27me3 ChIP Bernstein 2022 rep3|chip_H3K27me3_K562_Bernstein_2022_rep3|full"
  "H3K27me3 CnT Abbasova 2025 rep2|CnT_H3K27me3_K562_Abbasova_2025_rep2|downsample"
  "H3K27me3 CnT KayaOkur 2019 rep1|CnT_H3K27me3_K562_KayaOkur_2019_rep1|downsample"
  "H3K27me3 CnT KayaOkur 2020 NNPC|CnT_H3K27me3_K562_KayaOkur_2020_NNPC|downsample"
  "H3K27me3 CnT KayaOkur 2020 NCSD|CnT_H3K27me3_K562_KayaOkur_2020_NCSD|downsample"
  "H3K27me3 CnT KayaOkur 2020 FNPC|CnT_H3K27me3_K562_KayaOkur_2020_FNPC|downsample"
  "H3K27me3 CnT KayaOkur 2020 FCSD|CnT_H3K27me3_K562_KayaOkur_2020_FCSD|downsample"
  "H3K27me3 CnT KayaOkur 2020 HSLT|CnT_H3K27me3_K562_KayaOkur_2020_HSLT|downsample"
  "H3K27me3 CnT KayaOkur 2020 LSLT|CnT_H3K27me3_K562_KayaOkur_2020_LSLT|downsample"
)

# --- Heatmap colour scales ----------------------------------------------------
declare -A COLOR_OF=(
  ["GC content"]="#000004,#51127c,#b73779,#fc8961,#fcfdbf"
  ["DNase-seq"]="white,#02C49B"
  ["ATAC-seq"]="white,#C2185B"
)
# Histone-mark columns are green; anything else without a colour is reported
# rather than silently defaulting.
color_for() {
    local l="$1"
    if [[ -n "${COLOR_OF[$l]:-}" ]]; then echo "${COLOR_OF[$l]}"; return 0; fi
    if [[ "$l" == H3K* ]]; then echo "white,green"; return 0; fi
    echo "     WARNING: no colour configured for column '${l}' -- drawn in grey" >&2
    echo "white,grey"
}

# =============================================================================
# PREFLIGHT AND LOOKUPS
# =============================================================================
echo "========================================================="
echo " Depth-matched CUT&Tag comparison -- ${CELL} broad marks"
echo " Target depth : ${TARGET_FRAGMENTS} fragments per CUT&Tag library"
echo " Output       : $WORK_DIR"
echo "========================================================="

miss=0
for t in samtools bamCoverage computeMatrix plotHeatmap multiBigwigSummary Rscript awk sort; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: missing tool: $t" >&2; miss=1; }
done
for f in "$MASTER_TSV" "$PARAMS" "$GC_BW" "$DNASE_BW" "$ATAC_BW"; do
    [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; miss=1; }
done
[[ $miss -eq 0 ]] || { echo "Preflight failed." >&2; exit 1; }

FC_THRESH=$(awk -F'\t' '$1=="fc_thresh"{print $2}'  "$PARAMS")
FDR_THRESH=$(awk -F'\t' '$1=="fdr_thresh"{print $2}' "$PARAMS")
echo "  thresholds from analysis_params.tsv: FDR < ${FDR_THRESH}, |log2FC| > ${FC_THRESH}"

col_index() {
    awk -F'\t' -v want="$1" 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$i); if($i==want){print i; exit}}}' "$MASTER_TSV"
}
COL_LABEL=$(col_index "plot_label")
COL_BAM=$(col_index "location_final_bam")
COL_BW=$(col_index "location_final_bigwig")
for v in COL_LABEL COL_BAM COL_BW; do
    [[ -n "${!v}" ]] || { echo "ERROR: ${v} not found in $MASTER_TSV" >&2; exit 1; }
done
lookup_by_label() {
    awk -F'\t' -v s="$1" -v n="$COL_LABEL" -v c="$2" \
        'NR>1{gsub(/\r/,"",$n); if($n==s){gsub(/\r/,"",$c); print $c; exit}}' "$MASTER_TSV"
}
bam_of()    { lookup_by_label "$1" "$COL_BAM"; }
bigwig_of() { lookup_by_label "$1" "$COL_BW"; }

# Each awk reads a FILE and the result is captured per iteration; piping the
# loop into `head -1` would SIGPIPE the later awks under pipefail.
diff_rows_for() {   # epitope, region_set, cnt_group, chip_group
    local epi="$1" rs="$2" cg="$3" chg="$4" m
    for m in "$DIFF_DIR/manifest_supp.tsv" "$DIFF_DIR/manifest_promoter_full.tsv"; do
        [[ -s "$m" ]] || continue
        awk -F'\t' -v c="$CELL" -v e="$epi" -v r="$rs" -v g="$cg" -v h="$chg" \
            'NR>1 && $2==c && $1==e && $3==r && $5==g && $6==h {print $11}' "$m"
    done | sort -u
}

# =============================================================================
# PART 1: DOWNSAMPLE THE CUT&Tag LIBRARIES
# =============================================================================
DS_METRICS="$WORK_DIR/downsampling_metrics.tsv"
printf 'plot_label\tsource_bam\tfragments_before\ttarget_fragments\tfraction_used\tfragments_after\tdownsampled\tdownsampled_bam\tdownsampled_bigwig\n' > "$DS_METRICS"

# Fragment count = reads carrying the read-1 flag, with secondary and
# supplementary alignments excluded. Counting reads and halving would be wrong
# for any BAM that is not perfectly paired.
count_fragments() { samtools view -c -f 64 -F 0x900 -@ "$THREADS" "$1"; }

downsample_one() {   # plot_label -> echoes the downsampled BAM path
    local lab="$1" src out frag frac after
    out="$DS_BAM_DIR/${lab}_${DS_TAG}.bam"
    src=$(bam_of "$lab")
    if [[ -z "$src" ]]; then
        echo "     ERROR: plot_label '${lab}' not found in $(basename "$MASTER_TSV")" >&2; echo ""; return 0
    fi
    if [[ ! -s "$src" ]]; then
        echo "     ERROR: BAM missing for '${lab}': ${src}" >&2; echo ""; return 0
    fi
    if [[ -s "$out" && -s "${out}.bai" ]]; then
        echo "     cached: $(basename "$out")" >&2; echo "$out"; return 0
    fi

    frag=$(count_fragments "$src")
    if [[ "$frag" -le "$TARGET_FRAGMENTS" ]]; then
        # Not enough depth to match the others. Linking rather than subsampling
        # keeps the library in the figure, but it is NOT depth-matched and the
        # metrics table records that.
        echo "     WARNING: ${lab} has ${frag} fragments, below the ${TARGET_FRAGMENTS} target." >&2
        echo "              Kept at native depth -- NOT comparable to the downsampled libraries." >&2
        ln -sf "$src" "$out"
        [[ -s "${src}.bai" ]] && ln -sf "${src}.bai" "${out}.bai" || samtools index -@ "$THREADS" "$out"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lab" "$src" "$frag" "$TARGET_FRAGMENTS" \
            "1.000000" "$frag" "no" "$out" "$DS_BW_DIR/${lab}_${DS_TAG}.bw" >> "$DS_METRICS"
        echo "$out"; return 0
    fi

    frac=$(awk -v t="$TARGET_FRAGMENTS" -v f="$frag" 'BEGIN{ printf "%.6f", t / f }')
    echo "     ${lab}: ${frag} fragments -> fraction ${frac}" >&2
    # --subsample selects by read name, so mates are kept or dropped together
    # and no fragment is ever split. A fixed seed makes the run reproducible.
    samtools view -b -@ "$THREADS" \
        --subsample "$frac" --subsample-seed "$DOWNSAMPLE_SEED" \
        -o "$out" "$src"
    samtools index -@ "$THREADS" "$out"
    after=$(count_fragments "$out")
    echo "     ${lab}: ${after} fragments after downsampling" >&2
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$lab" "$src" "$frag" "$TARGET_FRAGMENTS" \
        "$frac" "$after" "yes" "$out" "$DS_BW_DIR/${lab}_${DS_TAG}.bw" >> "$DS_METRICS"
    echo "$out"
}

# =============================================================================
# PART 2: BIGWIGS FROM THE DOWNSAMPLED BAMs
# =============================================================================
# Same call the processing pipeline uses for a per-library paired-end track:
# --extendReads with no argument, so each fragment contributes its true
# sequenced span, at binSize 1 and with no normalisation.
bigwig_from_bam() {   # bam, out_bw
    local bam="$1" out="$2"
    [[ -s "$out" ]] && { echo "     cached: $(basename "$out")" >&2; return 0; }
    echo "     bamCoverage -> $(basename "$out")" >&2
    bamCoverage -b "$bam" --extendReads --binSize "$BIGWIG_BIN_SIZE" \
        -p "$THREADS" -o "$out" >/dev/null 2>&1 || {
        echo "     ERROR: bamCoverage failed for $(basename "$bam")" >&2; return 1; }
}

# =============================================================================
# PART 3: ASSEMBLE THE TRACK LIST
# =============================================================================
declare -A TRACK_BW=()
BW=(); LBL=()
echo "--- Preparing tracks ---"
for rec in "${TRACKS[@]}"; do
    IFS='|' read -r DISP SRC MODE <<< "$rec"
    path=""
    case "$MODE" in
        ref)
            path="$SRC" ;;
        full)
            path=$(bigwig_of "$SRC")
            [[ -n "$path" ]] || echo "     ERROR: plot_label '${SRC}' not found in $(basename "$MASTER_TSV")" >&2 ;;
        downsample)
            if [[ "$DO_DOWNSAMPLE" == "1" ]]; then
                dsbam=$(downsample_one "$SRC")
            else
                dsbam="$DS_BAM_DIR/${SRC}_${DS_TAG}.bam"
            fi
            if [[ -n "$dsbam" && -s "$dsbam" ]]; then
                path="$DS_BW_DIR/${SRC}_${DS_TAG}.bw"
                if [[ "$DO_BIGWIGS" == "1" ]]; then
                    bigwig_from_bam "$dsbam" "$path" || path=""
                fi
            fi ;;
        *)  echo "     ERROR: unknown track mode '${MODE}' for ${DISP}" >&2 ;;
    esac
    if [[ -n "$path" && -s "$path" ]]; then
        BW+=("$path"); LBL+=("$DISP"); TRACK_BW["$DISP"]="$path"
    else
        echo "     skipping ${DISP}: no usable bigWig" >&2
    fi
done
[[ ${#BW[@]} -gt 0 ]] || { echo "ERROR: no tracks available." >&2; exit 1; }
echo "  ${#LBL[@]} tracks: $(printf '%s; ' "${LBL[@]}")"

# =============================================================================
# PART 4: REGION GROUPS
# =============================================================================
declare -A REGION_COMPARISON=()
for EPI in H3K36me3 H3K27me3; do
    RS="${EPI_REGION_SET[$EPI]}"
    mapfile -t rows < <(diff_rows_for "$EPI" "$RS" "${CNT_GROUP[$EPI]}" "${CHIP_GROUP[$EPI]}")
    if [[ ${#rows[@]} -eq 0 ]]; then
        echo "ERROR: no differential result for ${CELL} ${EPI} ${RS}" >&2
        echo "       CUT&Tag group : ${CNT_GROUP[$EPI]}" >&2
        echo "       ChIP-seq group: ${CHIP_GROUP[$EPI]}" >&2
        exit 1
    fi
    if [[ ${#rows[@]} -gt 1 ]]; then
        echo "ERROR: ${EPI} ${RS} matched ${#rows[@]} differential results; expected one." >&2
        printf '         %s\n' "${rows[@]}" >&2; exit 1
    fi
    DIFF="${rows[0]}"
    REGION_COMPARISON[$EPI]="${CNT_GROUP[$EPI]} vs ${CHIP_GROUP[$EPI]}"
    echo "--- ${EPI} regions (${SET_DESC[$RS]}) ---"
    echo "  from: ${REGION_COMPARISON[$EPI]}"

    PRE="$REG_DIR/${CELL}_${EPI}"
    awk -F'\t' -v fc="$FC_THRESH" -v fdr="$FDR_THRESH" -v pre="$PRE" '
      NR==1 { for(i=1;i<=NF;i++){ if($i=="Chr")c=i; if($i=="Start")s=i; if($i=="End")e=i;
                                  if($i=="BinID")b=i; if($i=="log2FC")l=i; if($i=="padj")p=i }
              next }
      { st="nondiff"
        if ($p+0 < fdr && $l+0 >  fc) st="cnt"
        else if ($p+0 < fdr && $l+0 < -fc) st="chip"
        print $c "\t" $s "\t" $e "\t" $b "\t0\t." > (pre "_" st ".bed") }' "$DIFF"
    for g in chip nondiff cnt; do
        [[ -s "${PRE}_${g}.bed" ]] || : > "${PRE}_${g}.bed"
        sort -k1,1 -k2,2n -o "${PRE}_${g}.bed" "${PRE}_${g}.bed"
        printf "    %-9s %8d bins\n" "$g" "$(wc -l < "${PRE}_${g}.bed")"
    done
done

K36="$REG_DIR/${CELL}_H3K36me3"; K27="$REG_DIR/${CELL}_H3K27me3"

# =============================================================================
# PART 5: HEATMAP
# =============================================================================
# Per-column limits, as in the other figure scripts: each column gets its own
# colour ceiling and its own profile axis, so a track with a wide dynamic range
# cannot flatten the others.
compute_matrix_limits() {   # matrix.gz
    Rscript - "$1" "$HEATMAP_COLOR_CAP_FRACTION" <<'EOF'
suppressPackageStartupMessages({ library(jsonlite); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); mgz <- a[1]; capf <- as.numeric(a[2])
if (!is.finite(capf) || capf < 0) capf <- 0
upper_q <- 1 - capf; lower_q <- capf
con <- gzfile(mgz, "r"); h1 <- readLines(con, n = 1); close(con)
hdr <- fromJSON(sub("^@", "", h1)); sb <- hdr$sample_boundaries; gb <- hdr$group_boundaries
dt <- fread(cmd = paste("zcat", shQuote(mgz)), skip = 1, header = FALSE, sep = "\t",
            na.strings = c("nan", "-nan", "NA", ""))
valm <- as.matrix(dt[, 7:ncol(dt)]); valm[is.na(valm)] <- 0
zmin <- c(); zmax <- c(); ymin <- c(); ymax <- c()
for (i in 1:(length(sb) - 1)) {
  m <- valm[, (sb[i] + 1):sb[i + 1], drop = FALSE]
  zx <- as.numeric(quantile(m, upper_q, na.rm = TRUE)); if (!is.finite(zx) || zx <= 0) zx <- 1
  zmax <- c(zmax, round(zx, 3))
  if (i == 1) {
    zn <- as.numeric(quantile(m, lower_q, na.rm = TRUE))
    if (!is.finite(zn) || zn < 0) zn <- 0
    if (zn >= zx) zn <- 0
  } else zn <- 0
  zmin <- c(zmin, round(zn, 3))
  gmat <- c()
  for (g in 1:(length(gb) - 1)) {
    gm <- m[(gb[g] + 1):gb[g + 1], , drop = FALSE]
    gmat <- c(gmat, colMeans(gm, na.rm = TRUE))
  }
  lo <- min(gmat, na.rm = TRUE); hi <- max(gmat, na.rm = TRUE)
  rng <- hi - lo; if (!is.finite(rng) || rng == 0) rng <- 1
  ymin <- c(ymin, round(lo - 0.1 * rng, 3)); ymax <- c(ymax, round(hi + 0.1 * rng, 3))
}
cat("ZMIN:", paste(zmin, collapse = " "), "\n")
cat("ZMAX:", paste(zmax, collapse = " "), "\n")
cat("YMIN:", paste(ymin, collapse = " "), "\n")
cat("YMAX:", paste(ymax, collapse = " "), "\n")
EOF
}

# The title names the source bigWig of every column, so the figure records
# which file produced each track and which of them were downsampled.
build_title() {   # header line
    local hdr="$1" i t
    [[ ${#HM_LABELS[@]} -gt 0 ]] || { printf '%s' "$hdr"; return 0; }
    t="${hdr}"$'\n'"columns, in order (dataset = bigWig):"
    for i in "${!HM_LABELS[@]}"; do
        t+=$'\n'"$((i + 1)). ${HM_LABELS[$i]}  =  $(basename "${TRACK_BW[${HM_LABELS[$i]}]:-n/a}")"
    done
    printf '%s' "$t"
}

HM_BED=("${K36}_nondiff.bed" "${K36}_cnt.bed" "${K36}_chip.bed"
        "${K27}_nondiff.bed" "${K27}_cnt.bed" "${K27}_chip.bed")
HM_GROUP_DISP=(
    "H3K36me3 non-differential (gene body)" "H3K36me3 CUT&Tag-enriched (gene body)"
    "H3K36me3 ChIP-seq-enriched (gene body)" "H3K27me3 non-differential (Polycomb)"
    "H3K27me3 CUT&Tag-enriched (Polycomb)"  "H3K27me3 ChIP-seq-enriched (Polycomb)")

MAT="$MAT_DIR/${CELL}_broadMarks_${DS_TAG}_matrix.gz"
SIG="${MAT}.signature"
SIG_NOW=$(printf '%s\n' "${LBL[@]}" "--" "${HM_BED[@]}")
if [[ -s "$MAT" && -s "$SIG" && "$(cat "$SIG")" == "$SIG_NOW" ]]; then
    echo "--- Matrix cached ---"
elif [[ "$COMPUTE_MATRICES" == "1" ]]; then
    echo "--- computeMatrix (${#BW[@]} tracks, ${#HM_BED[@]} region groups) ---"
    computeMatrix reference-point --referencePoint center \
        -S "${BW[@]}" --samplesLabel "${LBL[@]}" -R "${HM_BED[@]}" \
        -b "$BROAD_WINDOW" -a "$BROAD_WINDOW" --binSize "$BINSIZE" --missingDataAsZero \
        -p "$THREADS" -o "$MAT" >/dev/null 2>&1 && printf '%s\n' "$SIG_NOW" > "$SIG"
fi

if [[ "$PLOT_HEATMAP" == "1" && -s "$MAT" ]]; then
    HM_LABELS=("${LBL[@]}")
    CLIST=(); for l in "${LBL[@]}"; do CLIST+=("$(color_for "$l")"); done
    LIMITS=$(compute_matrix_limits "$MAT")
    ZMIN=$(echo "$LIMITS" | awk '/^ZMIN:/ { $1=""; print $0 }')
    ZMAX=$(echo "$LIMITS" | awk '/^ZMAX:/ { $1=""; print $0 }')
    YMIN=$(echo "$LIMITS" | awk '/^YMIN:/ { $1=""; print $0 }')
    YMAX=$(echo "$LIMITS" | awk '/^YMAX:/ { $1=""; print $0 }')
    OUT_HM="$HM_DIR/${CELL}_broadMarks_${DS_TAG}_heatmap.pdf"
    echo "--- plotHeatmap (${#LBL[@]} columns @ ${HEATMAP_WIDTH_CM} cm) ---"
    # shellcheck disable=SC2086
    plotHeatmap -m "$MAT" -o "$OUT_HM" --plotType lines \
        --plotTitle "$(build_title "${CELL} broad marks, CUT&Tag depth-matched to ${TARGET_FRAGMENTS} fragments
H3K36me3 regions: ${REGION_COMPARISON[H3K36me3]}
H3K27me3 regions: ${REGION_COMPARISON[H3K27me3]}")" \
        --sortRegions descend --sortUsing mean --sortUsingSamples 1 \
        --colorList "${CLIST[@]}" --zMin $ZMIN --zMax $ZMAX --yMin $YMIN --yMax $YMAX \
        --regionsLabel "${HM_GROUP_DISP[@]}" --refPointLabel "bin center" \
        --legendLocation center-right \
        --heatmapHeight "$HEATMAP_HEIGHT_CM" --heatmapWidth "$HEATMAP_WIDTH_CM" \
        >/dev/null 2>&1 || echo "  plotHeatmap failed" >&2
    echo "  wrote $OUT_HM"
fi

# =============================================================================
# PART 6: VIOLINS
# =============================================================================
cat > "$RS_DIR/violin.R" <<'RVIOLIN'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(scales)
})
options(warn = -1)

a <- commandArgs(trailingOnly = TRUE)
comb <- a[1]; out_pdf <- a[2]; stats_tsv <- a[3]
title_txt <- a[4]; subtitle_txt <- a[5]; caption_txt <- a[6]
gc_col <- a[7]
GROUP_LEVELS <- strsplit(a[8], "|", fixed = TRUE)[[1]]
GROUP_COLORS <- setNames(strsplit(a[9], "|", fixed = TRUE)[[1]], GROUP_LEVELS)

FONT_SCALE <- as.numeric(Sys.getenv("VIOLIN_FONT_SCALE", "2.25"))
fs <- function(x) x * FONT_SCALE
MEDIAN_COL <- "#1a9d3b"
LEGEND_LEVELS <- GROUP_LEVELS[nzchar(trimws(GROUP_LEVELS))]

df <- suppressMessages(read_tsv(comb, show_col_types = FALSE, na = c("nan", "-nan", "NA", "")))
if (nrow(df) == 0) { message("empty summary table"); q(save = "no") }

long <- df %>%
  pivot_longer(cols = -c(chr, start, end, Group), names_to = "panel", values_to = "Signal") %>%
  mutate(Signal = suppressWarnings(as.numeric(Signal)),
         Signal = ifelse(is.na(Signal), 0, Signal),
         Group  = factor(Group, levels = GROUP_LEVELS))

# GC is a percentage, not an assay readout, so it is plotted as mean % GC with
# no log transform. Unit detection happens once, from the column maximum.
gc_max <- suppressWarnings(max(long$Signal[long$panel == gc_col], na.rm = TRUE))
gc_mult <- if (is.finite(gc_max) && gc_max <= 1) 100 else 1
long <- long %>%
  mutate(is_pct = panel == gc_col,
         Value  = ifelse(is_pct, Signal * gc_mult, log2(Signal + 1)))

panel_levels <- names(df)[!(names(df) %in% c("chr", "start", "end", "Group"))]
long$panel <- factor(long$panel, levels = panel_levels)

# Every pair of region groups within a panel, two-sided Wilcoxon rank-sum, with
# Benjamini-Hochberg applied WITHIN each panel so a panel's labels never depend
# on which other panels share the figure.
pairs <- list()
for (p in levels(long$panel)) {
  d <- long %>% filter(panel == p)
  grps <- levels(droplevels(d$Group))
  if (length(grps) < 2) next
  for (cb in combn(grps, 2, simplify = FALSE)) {
    x <- d$Value[d$Group == cb[1]]; y <- d$Value[d$Group == cb[2]]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) next
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    pairs[[length(pairs) + 1]] <- data.frame(
      panel = p, group1 = cb[1], group2 = cb[2], n1 = length(x), n2 = length(y),
      median1 = median(x), median2 = median(y), p = wt$p.value, stringsAsFactors = FALSE)
  }
}
pair_tests <- if (length(pairs) > 0) bind_rows(pairs) else
  data.frame(panel = character(), group1 = character(), group2 = character(),
             n1 = integer(), n2 = integer(), median1 = double(), median2 = double(),
             p = double(), stringsAsFactors = FALSE)
if (nrow(pair_tests) > 0) {
  pair_tests <- pair_tests %>% group_by(panel) %>%
    mutate(p.adj = p.adjust(p, method = "BH")) %>% ungroup()
  write_tsv(pair_tests, stats_tsv)
  message("  wrote ", nrow(pair_tests), " pairwise tests -> ", stats_tsv)
} else { pair_tests$p.adj <- double() }

# Numeric adjusted p-values, never stars. Wilcoxon uses the normal
# approximation at these sample sizes and underflows to exactly 0, so a floored
# value is reported as a bound.
P_FLOOR <- .Machine$double.xmin
fmt_p <- function(p) ifelse(is.na(p), "n/a",
                    ifelse(p <= 0, paste0("p < ", formatC(P_FLOOR, format = "e", digits = 0)),
                    ifelse(p >= 0.001, paste0("p = ", formatC(p, format = "f", digits = 3)),
                                       paste0("p = ", formatC(p, format = "e", digits = 1)))))

meds <- long %>% group_by(panel, Group) %>%
  summarise(m = median(Value, na.rm = TRUE), is_pct = first(is_pct), .groups = "drop") %>%
  mutate(lbl = ifelse(is_pct, sprintf("%.1f%%", m), sprintf("%.2f", m)))
pan <- long %>% group_by(panel) %>%
  summarise(top = max(Value, na.rm = TRUE), bot = min(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(rng = ifelse(top - bot <= 0, 1, top - bot))
meds <- left_join(meds, pan, by = "panel") %>% mutate(y_med = top + rng * 0.06)

xpos <- setNames(seq_along(GROUP_LEVELS), GROUP_LEVELS)
br0 <- pair_tests %>% filter(!is.na(p.adj))
# Fifteen brackets per panel at the default spacing would push the violins into
# the bottom third, so the gap tightens once there are more than eight tiers.
n_tiers <- if (nrow(br0) > 0) max(table(br0$panel)) else 0
tier_gap <- if (n_tiers > 8) 0.055 else 0.085
br <- br0 %>%
  left_join(pan, by = "panel") %>%
  mutate(x1 = unname(xpos[group1]), x2 = unname(xpos[group2]), span = abs(x2 - x1)) %>%
  group_by(panel) %>% arrange(span, x1, .by_group = TRUE) %>%
  mutate(tier = row_number(),
         y = top + rng * (0.16 + tier_gap * (tier - 1)),
         label = fmt_p(p.adj)) %>% ungroup()

p <- ggplot(long, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(alpha = 0.7, colour = "black", width = 0.9, scale = "width") +
  geom_boxplot(width = 0.15, colour = "black", outlier.shape = NA) +
  geom_text(data = meds, aes(x = Group, y = y_med, label = lbl), inherit.aes = FALSE,
            colour = MEDIAN_COL, fontface = "bold.italic", size = fs(3.2)) +
  { if (nrow(br) > 0) geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE, linewidth = 0.6, colour = "grey15") } +
  { if (nrow(br) > 0) geom_text(data = br, aes(x = (x1 + x2) / 2, y = y, label = label),
      inherit.aes = FALSE, size = fs(2.3), vjust = -0.35, fontface = "bold", colour = "grey15") } +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE, name = "Region group",
                    breaks = LEGEND_LEVELS, na.value = "white") +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "log2(mean signal + 1)   |   GC panel: mean GC content (%)",
       title = title_txt, subtitle = subtitle_txt, caption = caption_txt) +
  theme_minimal(base_size = fs(10)) +
  theme(strip.text = element_text(size = fs(9), face = "bold", colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        panel.spacing = unit(0.7, "cm"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = fs(7), colour = "black"),
        axis.text.y = element_text(size = fs(9), colour = "black"),
        axis.title.y = element_text(size = fs(10), face = "bold"),
        plot.title = element_text(size = fs(13), face = "bold"),
        plot.subtitle = element_text(size = fs(9), face = "italic"),
        plot.caption = element_text(size = fs(6.5), hjust = 0, colour = "grey30", lineheight = 1.3),
        legend.position = "bottom",
        legend.title = element_text(size = fs(10), face = "bold"),
        legend.text = element_text(size = fs(9)))

n_pan <- length(panel_levels); ncols <- 3; nrows <- ceiling(n_pan / ncols)
ggsave(out_pdf, p, width = ncols * (4.0 + 0.62 * length(GROUP_LEVELS)),
       height = nrows * (5.6 + 0.42 * n_tiers) + 4.5,
       device = "pdf", bg = "white", limitsize = FALSE)
message("wrote ", out_pdf, " (", n_tiers, " bracket tiers per panel)")
RVIOLIN

export VIOLIN_FONT_SCALE

if [[ "$PLOT_VIOLINS" == "1" ]]; then
    echo "--- Violins ---"
    # Same six groups as the heatmap, in the ChIP -> CUT&Tag order used by every
    # other violin figure, with a blank level separating the two marks.
    VI_BED=("${K36}_chip.bed" "${K36}_nondiff.bed" "${K36}_cnt.bed" ""
            "${K27}_chip.bed" "${K27}_nondiff.bed" "${K27}_cnt.bed")
    VI_DISP=("H3K36me3 ChIP-seq Enriched" "H3K36me3 Non-differential" "H3K36me3 CUT&Tag Enriched" " "
             "H3K27me3 ChIP-seq Enriched" "H3K27me3 Non-differential" "H3K27me3 CUT&Tag Enriched")

    COMBINED="$VI_DIR/${CELL}_broadMarks_${DS_TAG}_signal.tsv"
    COLSIG="$VI_DIR/${CELL}_broadMarks_${DS_TAG}_columns.txt"
    { printf 'chr\tstart\tend'; for l in "${LBL[@]}"; do printf '\t%s' "$l"; done
      printf '\tGroup\n'; } > "$COMBINED"

    # A cached .tab holds one column per bigWig it was made from, so the column
    # set is fingerprinted; reusing a stale one would pair old numbers with new
    # column names.
    COLNOW=$(printf '%s\n' "${BW[@]}")
    if [[ ! -f "$COLSIG" || "$(cat "$COLSIG")" != "$COLNOW" ]]; then
        echo "  track list changed -- clearing cached violin summaries" >&2
        rm -f "$VI_DIR/${CELL}_broadMarks_${DS_TAG}_sig_"*.tab
        printf '%s\n' "$COLNOW" > "$COLSIG"
    fi

    for i in "${!VI_BED[@]}"; do
        bed="${VI_BED[$i]}"
        [[ -n "$bed" && -s "$bed" ]] || continue          # spacer level
        tab="$VI_DIR/${CELL}_broadMarks_${DS_TAG}_sig_${i}.tab"
        if [[ ! -s "$tab" || "$bed" -nt "$tab" ]]; then
            multiBigwigSummary BED-file -b "${BW[@]}" --BED "$bed" \
                -o "${tab%.tab}.npz" --outRawCounts "$tab" -p "$THREADS" >/dev/null 2>&1
            rm -f "${tab%.tab}.npz"
        fi
        tail -n +2 "$tab" | awk -v g="${VI_DISP[$i]}" 'BEGIN{OFS="\t"}{print $0, g}' >> "$COMBINED"
    done

    Rscript "$RS_DIR/violin.R" "$COMBINED" \
        "$VI_DIR/${CELL}_broadMarks_${DS_TAG}_violins.pdf" \
        "$VI_DIR/${CELL}_broadMarks_${DS_TAG}_pairwise_stats.tsv" \
        "${CELL} broad marks: CUT&Tag protocols at matched depth over differential region groups" \
        "CUT&Tag downsampled to ${TARGET_FRAGMENTS} fragments; ChIP-seq shown at native depth  |  raw log2(signal + 1)  |  status: FDR < ${FDR_THRESH}, |log2FC| > ${FC_THRESH}" \
        "Left block: active gene-body bins split by H3K36me3 differential status. Right block: Polycomb-repressed bins split by H3K27me3
differential status. Because the two marks are mutually exclusive, each block also serves as the other's negative control.
Every CUT&Tag library was downsampled to the same fragment count before its track was built, so panel heights ARE comparable between CUT&Tag
samples. The ChIP-seq panels are at native depth and are shown for reference only; their heights are not comparable to the CUT&Tag panels.
Green numbers are group MEDIANS. Brackets carry Benjamini-Hochberg adjusted p-values from two-sided Wilcoxon rank-sum tests over all fifteen
pairs, corrected WITHIN each panel. With thousands of regions per group nearly every comparison is significant, so read the effect size
rather than the p-value. See downsampling_metrics.tsv for the depth of every library before and after subsampling." \
        "GC content" \
        "H3K36me3 ChIP-seq Enriched|H3K36me3 Non-differential|H3K36me3 CUT&Tag Enriched| |H3K27me3 ChIP-seq Enriched|H3K27me3 Non-differential|H3K27me3 CUT&Tag Enriched" \
        "#3936ff|grey70|#ff3b3b|white|#3936ff|grey70|#ff3b3b"
fi

echo "========================================================="
echo " Done."
echo "   Downsampled BAMs : $DS_BAM_DIR"
echo "   Downsampled bigWigs : $DS_BW_DIR"
echo "   Depth metrics    : $DS_METRICS"
echo "   Region groups    : $REG_DIR"
echo "   Heatmap          : $HM_DIR"
echo "   Violins          : $VI_DIR"
echo "========================================================="