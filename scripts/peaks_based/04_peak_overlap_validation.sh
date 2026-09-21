#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# VALIDATION: do method-enriched bins carry real signal in their own method?
#
# The differential analysis says a bin is CUT&Tag-enriched or ChIP-seq-enriched
# relative to the other method. That alone cannot distinguish two possibilities:
#
#   (a) the bin carries genuine, method-specific enrichment, or
#   (b) the bin is background noise in one method that DESeq2 read as a
#       difference because the other method was near zero there.
#
# Peaks called independently within each method separate the two. If ChIP-seq-
# enriched bins are real ChIP-seq signal, they should overlap ChIP-seq peaks at
# a high rate -- peaks called from the ChIP-seq data alone, with no reference to
# CUT&Tag. Noise would not.
#
#   PART 1  Peak calling with MACS2
#   PART 2  Merged peak sets: peaks from both replicates merged together
#   PART 3  Region sets, split by differential status
#   PART 4  Overlap counting
#   PART 5  Bar figure
#
# SIX FIGURE PANELS PER CELL LINE. The promoter marks are evaluated twice, on
# both region geometries used elsewhere in the manuscript:
#
#   1. H3K4me3   active promoter, 500 bp bins      <- the geometry used for the
#   2. H3K27ac   active promoter, 500 bp bins         violin/volcano/MA figures
#   3. H3K4me3   whole active promoter, 5 kb
#   4. H3K27ac   whole active promoter, 5 kb
#   5. H3K27me3  Polycomb repressed, 3 kb bins
#   6. H3K36me3  active gene body, 3 kb bins
#
# Peaks are called ONCE per cell line x epitope and reused across both promoter
# geometries: the peak set does not depend on which regions it is counted over.
#
# ONE CALLER FOR BOTH METHODS. MACS2 is run on ChIP-seq and on CUT&Tag with the
# same statistical settings, so the peak-calling algorithm cannot itself explain
# any difference between the two peak sets. The only deliberate difference is
# the input format: ChIP-seq libraries may be single-end and are read as -f BAM
# with a fixed --extsize 200, while CUT&Tag is always paired-end and is read as
# -f BAMPE so MACS2 uses the true fragment ends.
#
# HOW TO READ THE RESULT. The informative axis is WITHIN one peak set, across
# the three bin classes: does the ChIP-seq peak set cover ChIP-seq-enriched bins
# far more often than CUT&Tag-enriched bins, and vice versa? Absolute
# percentages are NOT comparable between panels. A 500 bp bin has less chance of
# touching a peak than a 5 kb window purely because it is shorter, so the
# 500 bp panels sit lower than the whole-promoter panels throughout; that offset
# is geometry, not biology. Each panel carries its own size-matched negative
# control for exactly this reason.
#
# Usage:
#   conda activate peak_calling_analysis
#   bash 03_peak_overlap_validation.sh
# =============================================================================

# --------------------------- USER CONFIG -------------------------------------
# Output of 01_differential_analysis.sh
DIFF_WORK_DIR="${DIFF_WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output}"
DIFF_DIR="$DIFF_WORK_DIR/differential"
REGIONS_DIR="$DIFF_WORK_DIR/regions"     # negative-control BEDs come from here
PARAMS="$DIFF_DIR/analysis_params.tsv"

WORK_DIR="${WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output/peak_overlap_validation_MACS2_wPromBins}"
PEAKS_DIR="$WORK_DIR/peaks"
MERGED_DIR="$WORK_DIR/peaks_merged_replicates"
BINS_DIR="$WORK_DIR/bins_by_status"
DERIVED_DIR="$WORK_DIR/regions_derived"   # 500 bp repressed-promoter bins
TABLE_DIR="$WORK_DIR/tables"
FIG_DIR="$WORK_DIR/figures"
RS_DIR="$WORK_DIR/Rscripts"
TMP="$WORK_DIR/_tmp"
mkdir -p "$PEAKS_DIR" "$MERGED_DIR" "$BINS_DIR" "$DERIVED_DIR" "$TABLE_DIR" "$FIG_DIR" "$RS_DIR" "$TMP"

MASTER_TSV="${MASTER_TSV:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/summary/Master_Sample_Metrics.tsv}"
BLACKLIST="$HOME/gpfs/genomes/bedfiles/hg38-blacklist.v2.bed"
GENOME="hg38"

# Must match PROM_BIN_SIZE in 01_differential_analysis.sh: the 500 bp repressed
# promoter bins built below are the size-matched control for the 500 bp active
# promoter bins that 01 tested, and a mismatch would make the control bar
# incomparable to the bars beside it.
PROM_BIN_SIZE=500

# --- Stage toggles -----------------------------------------------------------
DO_MACS2=1
DO_OVERLAP=1
DO_FIGURE=1

# --- Overlap definition ------------------------------------------------------
# Minimum fraction of a BIN that must be covered by a peak for it to count as
# overlapping. 0 = any overlap of at least 1 bp (bedtools default). Raise it
# (e.g. 0.2) if you want to require substantive coverage rather than a clipped
# edge.
OVERLAP_MIN_FRAC="${OVERLAP_MIN_FRAC:-0}"

# --- MACS2 parameters --------------------------------------------------------
# No input control is available for these libraries, so --nolambda disables the
# local lambda model. Duplicates were already removed upstream, so --keep-dup
# all keeps every remaining fragment. Everything except the input format and
# read extension is identical between the two methods.
CHIP_NARROW_ARGS="-f BAM -g hs -q 1e-5 --nolambda --keep-dup all --nomodel --extsize 200"
CHIP_BROAD_ARGS="$CHIP_NARROW_ARGS --broad --broad-cutoff 0.1"
CNT_NARROW_ARGS="-f BAMPE -g hs -q 1e-5 --nolambda --keep-dup all --nomodel"
CNT_BROAD_ARGS="$CNT_NARROW_ARGS --broad --broad-cutoff 0.1"

# --- Per-epitope settings ----------------------------------------------------
# Narrow mode for the punctate promoter marks, broad for the domain marks.
declare -A EPI_PEAK_MODE=(
  [H3K4me3]=narrow [H3K27ac]=narrow [H3K27me3]=broad [H3K36me3]=broad
)

# Which region sets each epitope is evaluated on, space separated. The promoter
# marks get BOTH geometries, so each produces two figure panels; the broad
# marks have one geometry and produce one panel each.
declare -A EPI_REGION_SETS=(
  [H3K4me3]="promoter_bins promoter_full"
  [H3K27ac]="promoter_bins promoter_full"
  [H3K36me3]="genebody_bins"
  [H3K27me3]="polycomb_bins"
)

# NEGATIVE CONTROL, keyed on the REGION SET rather than the epitope, so each
# panel's control is the same size as the bins beside it. Without a control, a
# peak set that covered most of the genome would produce high overlap
# everywhere and still look specific, because every bar would be drawn from
# regions the mark genuinely occupies. The control answers the separate
# question of whether the peaks are confined to the chromatin class they should
# be.
#   promoter bins  -> 500 bp bins of Polycomb-repressed, untranscribed promoters
#   whole promoter -> whole 5 kb Polycomb-repressed promoters
#   gene body      -> Polycomb domains (devoid of a transcription mark)
#   Polycomb       -> active gene bodies (devoid of a repressive mark)
declare -A CONTROL_FOR_SET=(
  [promoter_bins]=promoter_repressed_bins
  [promoter_full]=promoter_repressed
  [genebody_bins]=polycomb_bins
  [polycomb_bins]=genebody_bins
)

# Human-readable names used in the tables and in the figure's region line.
declare -A SET_DESC=(
  [promoter_bins]="active promoters, 500 bp bins"
  [promoter_full]="whole active promoters, 5 kb"
  [promoter_repressed_bins]="Polycomb-repressed promoters, 500 bp bins"
  [promoter_repressed]="whole Polycomb-repressed promoters, 5 kb"
  [genebody_bins]="active gene-body bins, 3 kb"
  [polycomb_bins]="Polycomb-repressed bins, 3 kb"
)

# --- One representative comparison per cell line x epitope -------------------
# cell | epitope | CUT&Tag group | ChIP-seq group | CnT peak reps | ChIP peak reps
# The groups locate the differential result; the replicate lists are the two
# libraries peaks are called on. They need not be the whole group: K562
# H3K36me3 was tested with five CUT&Tag replicates but peaks are called on two.
PEAK_JOBS=(
  "K562|H3K4me3|CnT_H3K4me3_K562_KayaOkur_2020_rep1_2020_rep2|chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2|CnT_H3K4me3_K562_KayaOkur_2020_rep1,CnT_H3K4me3_K562_KayaOkur_2020_rep2|chip_H3K4me3_K562_Bernstein_2017_rep1,chip_H3K4me3_K562_Bernstein_2017_rep2"
  "K562|H3K27ac|CnT_H3K27ac_K562_Abbasova_2025_ab177r1_2025_ab177r2|chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3|CnT_H3K27ac_K562_Abbasova_2025_ab177r1,CnT_H3K27ac_K562_Abbasova_2025_ab177r2|chip_H3K27ac_K562_Bernstein_2011_rep2,chip_H3K27ac_K562_Bernstein_2022_rep3"
  "K562|H3K27me3|CnT_H3K27me3_K562_Abbasova_2025_rep1_2025_rep2|chip_H3K27me3_K562_Bernstein_2011_rep1_2022_rep3|CnT_H3K27me3_K562_Abbasova_2025_rep1,CnT_H3K27me3_K562_Abbasova_2025_rep2|chip_H3K27me3_K562_Bernstein_2011_rep1,chip_H3K27me3_K562_Bernstein_2022_rep3"
  "K562|H3K36me3|CnT_H3K36me3_K562_Wu_2025_rep1_2025_rep2_2025_rep3_2025_rep4_2025_rep5|chip_H3K36me3_K562_Bernstein_2011_rep2_2022_rep3|CnT_H3K36me3_K562_Wu_2025_rep4,CnT_H3K36me3_K562_Wu_2025_rep5|chip_H3K36me3_K562_Bernstein_2011_rep2,chip_H3K36me3_K562_Bernstein_2022_rep3"
  "MCF7|H3K4me3|CnT_H3K4me3_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K4me3_MCF7_Bernstein_2017_rep1_2017_rep2|CnT_H3K4me3_MCF7_Tian_2023_rep1,CnT_H3K4me3_MCF7_Tian_2023_rep2|chip_H3K4me3_MCF7_Bernstein_2017_rep1,chip_H3K4me3_MCF7_Bernstein_2017_rep2"
  "MCF7|H3K27ac|CnT_H3K27ac_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K27ac_MCF7_Bernstein_2017_rep1_2017_rep2|CnT_H3K27ac_MCF7_Tian_2023_rep1,CnT_H3K27ac_MCF7_Tian_2023_rep2|chip_H3K27ac_MCF7_Bernstein_2017_rep1,chip_H3K27ac_MCF7_Bernstein_2017_rep2"
  "MCF7|H3K27me3|CnT_H3K27me3_MCF7_Tian_2023_rep1_2023_rep2_2023_rep3_2023_rep4|chip_H3K27me3_MCF7_Bernstein_2017_rep1_2017_rep2|CnT_H3K27me3_MCF7_Tian_2023_rep3,CnT_H3K27me3_MCF7_Tian_2023_rep4|chip_H3K27me3_MCF7_Bernstein_2017_rep1,chip_H3K27me3_MCF7_Bernstein_2017_rep2"
  "MCF7|H3K36me3|CnT_H3K36me3_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K36me3_MCF7_Bernstein_2021_rep1_2021_rep2|CnT_H3K36me3_MCF7_Tian_2023_rep1,CnT_H3K36me3_MCF7_Tian_2023_rep2|chip_H3K36me3_MCF7_Bernstein_2021_rep1,chip_H3K36me3_MCF7_Bernstein_2021_rep2"
)

# =============================================================================
# PART 0: PREFLIGHT AND LOOKUPS
# =============================================================================
echo "========================================================="
echo " Peak-overlap validation (MACS2)"
echo " Differential input : $DIFF_DIR"
echo " Output             : $WORK_DIR"
echo "========================================================="

miss=0
for t in bedtools awk sort Rscript; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: missing tool: $t" >&2; miss=1; }
done
[[ "$DO_MACS2" == "1" ]] && { command -v macs2 >/dev/null 2>&1 || { echo "ERROR: macs2 not found" >&2; miss=1; }; }
for f in "$MASTER_TSV" "$BLACKLIST" "$PARAMS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; miss=1; }
done
[[ $miss -eq 0 ]] || { echo "Preflight failed." >&2; exit 1; }

# Thresholds come from the differential analysis, never redeclared here.
FC_THRESH=$(awk -F'\t' '$1=="fc_thresh"{print $2}'  "$PARAMS")
FDR_THRESH=$(awk -F'\t' '$1=="fdr_thresh"{print $2}' "$PARAMS")
[[ -n "$FC_THRESH" && -n "$FDR_THRESH" ]] || { echo "ERROR: could not read thresholds from $PARAMS" >&2; exit 1; }
echo "  thresholds from analysis_params.tsv: FDR < ${FDR_THRESH}, |log2FC| > ${FC_THRESH}"

# The promoter bin size is recorded by 01; warn loudly if it disagrees with the
# size used to build the control bins here.
PARAM_BIN=$(awk -F'\t' '$1=="prom_bin_size"{print $2}' "$PARAMS")
if [[ -n "$PARAM_BIN" && "$PARAM_BIN" != "$PROM_BIN_SIZE" ]]; then
    echo "  WARNING: PROM_BIN_SIZE here is ${PROM_BIN_SIZE} but 01 used ${PARAM_BIN}." >&2
    echo "           The 500 bp control bins would not be size-matched. Fix before trusting the figure." >&2
fi

col_index() {
    awk -F'\t' -v want="$1" 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$i); if($i==want){print i; exit}}}' "$MASTER_TSV"
}
COL_LABEL=$(col_index "plot_label")
COL_BAM=$(col_index "location_final_bam")
for v in COL_LABEL COL_BAM; do
    [[ -n "${!v}" ]] || { echo "ERROR: ${v} not found in $MASTER_TSV" >&2; exit 1; }
done
lookup_by_label() {
    awk -F'\t' -v s="$1" -v n="$COL_LABEL" -v c="$2" \
        'NR>1{gsub(/\r/,"",$n); if($n==s){gsub(/\r/,"",$c); print $c; exit}}' "$MASTER_TSV"
}
bam_of() { lookup_by_label "$1" "$COL_BAM"; }

# Finds the DESeq2 result for a comparison by searching the manifests, so a
# path change in 01_differential_analysis.sh cannot silently break this script.
#
# BOTH sides of the comparison are matched. Five of the eight jobs here use a
# CUT&Tag group that 01 tested against two independent ChIP-seq sources, so
# matching on the CUT&Tag group alone returns whichever row happens to come
# first in COMPARISONS. That currently gives the right answer for all five, but
# only by luck: reordering COMPARISONS in 01 would silently pair these bars
# with a different ChIP-seq dataset than the one the peaks were called on.
diff_path_for() {   # cell, epitope, region_set, cnt_group, chip_group
    local cell="$1" epi="$2" rs="$3" cg="$4" chg="$5" m hit
    # Each awk reads a FILE and the result is captured per iteration. Piping the
    # whole loop into `head -1` would work only by luck: head exits after the
    # first line and the remaining awks then die of SIGPIPE, which pipefail
    # turns into a fatal error partway through the run.
    for m in "$DIFF_DIR/manifest_promoter_full.tsv" "$DIFF_DIR/manifest_supp.tsv" "$DIFF_DIR/manifest_main.tsv"; do
        [[ -s "$m" ]] || continue
        hit=$(awk -F'\t' -v c="$cell" -v e="$epi" -v r="$rs" -v g="$cg" -v h="$chg" \
              'NR>1 && $2==c && $1==e && $3==r && $5==g && $6==h {print $11; exit}' "$m")
        [[ -n "$hit" ]] && { echo "$hit"; return 0; }
    done
    echo ""
}

# --- Peak cleanup: main chromosomes, blacklist removed, merged ---------------
# Applied identically to both methods' output so the two peak sets differ only
# by the data that produced them.
BL_CLEAN="$TMP/blacklist_clean.bed"
[[ -s "$BL_CLEAN" ]] || cut -f1-3 "$BLACKLIST" | tr -d '\r' \
    | awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' | sort -k1,1 -k2,2n | bedtools merge -i - > "$BL_CLEAN"

clean_peaks() {   # raw peak file, out bed
    local IN="$1" OUT="$2"
    [[ -s "$IN" ]] || { : > "$OUT"; return 0; }
    cut -f1-3 "$IN" | tr -d '\r' \
      | awk -F'\t' 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|[XY])$/ && $2 ~ /^[0-9]+$/ && $3 > $2' \
      | sort -k1,1 -k2,2n | bedtools merge -i - \
      | bedtools intersect -a - -b "$BL_CLEAN" -v > "$OUT"
}

# --- 500 bp repressed-promoter bins, the control for the promoter_bins panel -
# 01 writes the repressed promoters as whole 5 kb windows only. Tiling them
# here reproduces exactly what 01's tile_regions does to the ACTIVE promoters:
# makewindows, keep only tiles of the full bin size, drop blacklist, dedupe.
# Built once and cached, since it depends on nothing that changes per run.
ensure_repressed_bins() {   # cell -> echoes path, or empty
    local CL="$1"
    local SRC="$REGIONS_DIR/${CL}_promoter_repressed.bed"
    local OUT="$DERIVED_DIR/${CL}_promoter_repressed_bins.bed"
    if [[ -s "$OUT" ]]; then echo "$OUT"; return 0; fi
    if [[ ! -s "$SRC" ]]; then
        echo "      WARNING: repressed promoters not found: $SRC" >&2; echo ""; return 0
    fi
    bedtools makewindows -b "$SRC" -w "$PROM_BIN_SIZE" -i srcwinnum \
      | awk -F'\t' -v w="$PROM_BIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' \
      | sort -k1,1 -k2,2n \
      | bedtools intersect -a - -b "$BL_CLEAN" -v \
      | awk -F'\t' '!seen[$1"\t"$2"\t"$3]++' > "$OUT"
    echo "      built ${PROM_BIN_SIZE} bp repressed-promoter bins: $(wc -l < "$OUT") from $(wc -l < "$SRC") promoters" >&2
    echo "$OUT"
}

# Resolves a control region set to a BED path, building the derived one on
# demand and taking every other one straight from 01's output.
control_bed_for() {   # cell, control set name -> echoes path, or empty
    local CL="$1" CS="$2" p
    if [[ "$CS" == "promoter_repressed_bins" ]]; then
        ensure_repressed_bins "$CL"
    else
        p="$REGIONS_DIR/${CL}_${CS}.bed"
        [[ -s "$p" ]] && echo "$p" || echo ""
    fi
}

# =============================================================================
# PART 1: PEAK CALLING
# =============================================================================
call_macs2() {   # plot_label, epitope, method(chip|cnt) -> writes cleaned bed
    local lab="$1" epi="$2" meth="$3"
    local out="$PEAKS_DIR/${lab}_macs2.bed"
    [[ -s "$out" ]] && { echo "      cached: $(basename "$out")"; return 0; }
    local bam args mode ext raw
    bam=$(bam_of "$lab"); mode="${EPI_PEAK_MODE[$epi]}"
    [[ -s "$bam" ]] || { echo "      ERROR: BAM missing for ${lab}" >&2; return 1; }

    if [[ "$meth" == "cnt" ]]; then
        [[ "$mode" == "broad" ]] && args="$CNT_BROAD_ARGS" || args="$CNT_NARROW_ARGS"
    else
        [[ "$mode" == "broad" ]] && args="$CHIP_BROAD_ARGS" || args="$CHIP_NARROW_ARGS"
        # -f BAMPE is not used for ChIP because these libraries may be single-end;
        # the fixed --extsize 200 applies to both cases.
    fi
    [[ "$mode" == "broad" ]] && ext="broadPeak" || ext="narrowPeak"

    echo "      MACS2 (${mode}): ${lab}"
    # shellcheck disable=SC2086
    macs2 callpeak -t "$bam" -n "${lab}_macs2" $args --outdir "$PEAKS_DIR" \
        > "$PEAKS_DIR/${lab}_macs2.log" 2>&1 || {
        echo "      ERROR: macs2 failed for ${lab}; see ${PEAKS_DIR}/${lab}_macs2.log" >&2; return 1; }
    raw="$PEAKS_DIR/${lab}_macs2_peaks.${ext}"
    [[ -s "$raw" ]] || { echo "      ERROR: no ${ext} produced for ${lab}" >&2; return 1; }
    clean_peaks "$raw" "$out"
}

# =============================================================================
# PART 2: MERGED PEAK SETS
# =============================================================================
# A peak enters the master set if it is present in EITHER replicate. Taking the
# union of the two replicates keeps all potential signal regions.
merged_peaks() {   # out bed, rep1 bed, rep2 bed
    local OUT="$1" A="$2" B="$3"
    [[ -s "$OUT" ]] && { echo "      cached: $(basename "$OUT")"; return 0; }
    if [[ ! -s "$A" || ! -s "$B" ]]; then
        echo "      WARNING: missing replicate peak file; wrote empty $(basename "$OUT")" >&2
        : > "$OUT"; return 0
    fi
    cat "$A" "$B" \
      | sort -k1,1 -k2,2n | bedtools merge -i - > "$OUT"
    printf "      merged peaks: %-46s %8d (from %d and %d)\n" "$(basename "$OUT")" \
        "$(wc -l < "$OUT")" "$(wc -l < "$A")" "$(wc -l < "$B")"
}

# =============================================================================
# MAIN LOOP: CALL PEAKS
# =============================================================================
# Once per cell line x epitope. The peak set is a property of the libraries, not
# of the regions it is later counted over, so both promoter geometries share it.
for job in "${PEAK_JOBS[@]}"; do
    IFS='|' read -r CELL EPI CNT_G CHIP_G CNT_REPS CHIP_REPS <<< "$job"
    echo "--- ${CELL} ${EPI} (${EPI_PEAK_MODE[$EPI]}) ---"
    IFS=',' read -ra cnt_reps  <<< "$CNT_REPS"
    IFS=',' read -ra chip_reps <<< "$CHIP_REPS"

    for lab in "${chip_reps[@]}"; do
        [[ "$DO_MACS2" == "1" ]] && call_macs2 "$lab" "$EPI" chip || true
    done
    for lab in "${cnt_reps[@]}"; do
        [[ "$DO_MACS2" == "1" ]] && call_macs2 "$lab" "$EPI" cnt || true
    done

    merged_peaks "$MERGED_DIR/${CELL}_${EPI}_ChIPseq_MACS2.bed" \
        "$PEAKS_DIR/${chip_reps[0]}_macs2.bed" "$PEAKS_DIR/${chip_reps[1]}_macs2.bed"
    merged_peaks "$MERGED_DIR/${CELL}_${EPI}_CUTandTag_MACS2.bed" \
        "$PEAKS_DIR/${cnt_reps[0]}_macs2.bed" "$PEAKS_DIR/${cnt_reps[1]}_macs2.bed"
done

# =============================================================================
# PART 3 + 4: REGIONS BY STATUS, AND OVERLAP COUNTING
# =============================================================================
OVERLAP_TSV="$TABLE_DIR/peak_overlap_by_status.tsv"

if [[ "$DO_OVERLAP" == "1" ]]; then
printf 'cell_line\tepitope\tregion_role\tregion_set\tregion_desc\tcomparison\tbin_class\tn_bins\tmethod\tpeak_caller\tpeak_set\tn_peaks\tn_bins_overlapping\tpct_bins_overlapping\n' > "$OVERLAP_TSV"

# bedtools -f requires a positive fraction, so it is only added when asked for.
FRAC_ARG=()
awk -v f="$OVERLAP_MIN_FRAC" 'BEGIN{exit !(f+0 > 0)}' && FRAC_ARG=( -f "$OVERLAP_MIN_FRAC" )

# Counts how many of a BED's intervals overlap a peak set, and appends a row.
# PANEL_SET is the region set that names the figure panel, so the control row
# is filed under the panel it belongs to rather than under its own region set.
count_and_record() {   # cell, epitope, role, panel_set, region_set, comparison, bin_class, bins_bed
    local CELL="$1" EPI="$2" ROLE="$3" PANEL_SET="$4" RS="$5" CMP="$6" CLASS="$7" BINS="$8"
    local N_BINS N_PEAKS N_OV PCT ps PS_LABEL PS_METHOD PS_CALLER PS_FILE
    N_BINS=$(wc -l < "$BINS")
    for ps in "${PEAK_SETS[@]}"; do
        IFS='|' read -r PS_LABEL PS_METHOD PS_CALLER PS_FILE <<< "$ps"
        N_PEAKS=0; [[ -s "$PS_FILE" ]] && N_PEAKS=$(wc -l < "$PS_FILE")
        N_OV=0
        if [[ "$N_BINS" -gt 0 && "$N_PEAKS" -gt 0 ]]; then
            N_OV=$(bedtools intersect -a "$BINS" -b "$PS_FILE" -u "${FRAC_ARG[@]+"${FRAC_ARG[@]}"}" | wc -l)
        fi
        PCT=$(awk -v a="$N_OV" -v b="$N_BINS" 'BEGIN{ if(b>0) printf "%.2f", 100*a/b; else print "NA" }')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$CELL" "$EPI" "$ROLE" "$PANEL_SET" "${SET_DESC[$RS]:-$RS}" "$CMP" "$CLASS" "$N_BINS" \
            "$PS_METHOD" "$PS_CALLER" "$PS_LABEL" "$N_PEAKS" "$N_OV" "$PCT" >> "$OVERLAP_TSV"
    done
    printf "        %-26s %-42s %8d regions\n" "$CLASS" "${SET_DESC[$RS]:-$RS}" "$N_BINS"
}

declare -A CLASS_DISP=(
  [chip_enriched]="ChIP-seq enriched" [cnt_enriched]="CUT&Tag enriched"
  [nondifferential]="Non-differential"
)

echo "--- Overlap counting ---"
for job in "${PEAK_JOBS[@]}"; do
    IFS='|' read -r CELL EPI CNT_G CHIP_G CNT_REPS CHIP_REPS <<< "$job"
    COMPARISON="${CNT_G}_VS_${CHIP_G}"
    echo "  ${CELL} ${EPI}"

    # peak set label | method | caller | file
    # The ChIP-seq entry comes first so the figure can shade it as the left block.
    # Called once per cell x epitope and reused for every region set below.
    PEAK_SETS=(
      "ChIP-seq MACS2|ChIP-seq|MACS2|$MERGED_DIR/${CELL}_${EPI}_ChIPseq_MACS2.bed"
      "CUT&Tag MACS2|CUT&Tag|MACS2|$MERGED_DIR/${CELL}_${EPI}_CUTandTag_MACS2.bed"
    )

    for RS in ${EPI_REGION_SETS[$EPI]}; do
        CONTROL="${CONTROL_FOR_SET[$RS]}"
        DIFF=$(diff_path_for "$CELL" "$EPI" "$RS" "$CNT_G" "$CHIP_G")
        if [[ -z "$DIFF" || ! -s "$DIFF" ]]; then
            echo "    WARNING: no differential result for ${CELL} ${EPI} ${RS}" >&2
            echo "             CUT&Tag group : ${CNT_G}" >&2
            echo "             ChIP-seq group: ${CHIP_G}" >&2
            continue
        fi
        echo "    region set: ${RS}  (${SET_DESC[$RS]:-$RS})"

        # --- primary set, split by differential status ----------------------
        PRE="$BINS_DIR/${CELL}_${EPI}_${RS}"
        awk -F'\t' -v fc="$FC_THRESH" -v fdr="$FDR_THRESH" -v pre="$PRE" '
          NR==1 { for(i=1;i<=NF;i++){ if($i=="Chr")c=i; if($i=="Start")s=i; if($i=="End")e=i;
                                      if($i=="BinID")b=i; if($i=="log2FC")l=i; if($i=="padj")p=i }
                  if(!c||!s||!e||!l||!p){ print "ERROR: expected columns not found in " FILENAME > "/dev/stderr"; exit 2 }
                  next }
          {
            st = "nondifferential"
            if ($p+0 < fdr && $l+0 >  fc) st = "cnt_enriched"
            else if ($p+0 < fdr && $l+0 < -fc) st = "chip_enriched"
            print $c "\t" $s "\t" $e "\t" $b > (pre "_" st ".bed")
          }' "$DIFF"

        for st in chip_enriched cnt_enriched nondifferential; do
            f="${PRE}_${st}.bed"; [[ -s "$f" ]] || : > "$f"
            sort -k1,1 -k2,2n -o "$f" "$f"
            count_and_record "$CELL" "$EPI" "primary" "$RS" "$RS" "$COMPARISON" "${CLASS_DISP[$st]}" "$f"
        done

        # --- negative control set, taken whole ------------------------------
        # Not split by status: these regions were never part of this comparison,
        # so there is no differential result to split them by. Every interval in
        # the set contributes to one bar. It is filed under the PANEL's region
        # set so the control bar lands in the panel it controls for.
        CTRL_BED=$(control_bed_for "$CELL" "$CONTROL")
        if [[ -n "$CTRL_BED" && -s "$CTRL_BED" ]]; then
            CTRL_SORTED="$BINS_DIR/${CELL}_${EPI}_${RS}_control_${CONTROL}.bed"
            cut -f1-4 "$CTRL_BED" | sort -k1,1 -k2,2n > "$CTRL_SORTED"
            count_and_record "$CELL" "$EPI" "control" "$RS" "$CONTROL" "$COMPARISON" "Negative control" "$CTRL_SORTED"
        else
            echo "        WARNING: control regions unavailable for ${CONTROL}" >&2
        fi
    done
done
echo "  wrote $OVERLAP_TSV"
fi

# =============================================================================
# PART 5: FIGURE
# =============================================================================
cat > "$RS_DIR/plot_peak_overlap.R" <<'RPLOT'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales)
})
options(warn = -1)

a <- commandArgs(trailingOnly = TRUE)
in_tsv <- a[1]; out_pdf <- a[2]; fc <- a[3]; fdr <- a[4]

FONT_SCALE <- 1.5
fs <- function(x) x * FONT_SCALE

df <- suppressMessages(read_tsv(in_tsv, show_col_types = FALSE))
if (nrow(df) == 0) { message("empty overlap table"); q(save = "no") }

# Two peak sets, both MACS2. The x-axis labels name the method rather than the
# caller, since the caller is the same on both sides and is stated in the caption.
PEAK_LEVELS <- c("ChIP-seq MACS2", "CUT&Tag MACS2")
PEAK_SHORT  <- c("ChIP-seq MACS2" = "ChIP-seq\npeaks", "CUT&Tag MACS2" = "CUT&Tag\npeaks")
CLASS_LEVELS <- c("ChIP-seq enriched", "Non-differential", "CUT&Tag enriched", "Negative control")
CLASS_COLORS <- c("ChIP-seq enriched" = "#3936ff", "Non-differential" = "grey70",
                  "CUT&Tag enriched"  = "#ff3b3b", "Negative control" = "grey25")

# --- PANEL ORDER ------------------------------------------------------------
# Stated explicitly as an ordered epitope x region-set table rather than
# recovered from the data: the two promoter marks each appear twice, so
# faceting on epitope alone would collapse the two geometries on top of each
# other, and anything that infers the order from the table risks falling back
# to alphabetical.
PANEL_ORDER <- data.frame(
  epitope    = c("H3K4me3", "H3K27ac", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3"),
  region_set = c("promoter_bins", "promoter_bins", "promoter_full", "promoter_full",
                 "polycomb_bins", "genebody_bins"),
  stringsAsFactors = FALSE)

RS_SHORT <- c(promoter_bins = "active promoter, 500 bp bins",
              promoter_full = "whole active promoter, 5 kb",
              genebody_bins = "active gene body, 3 kb bins",
              polycomb_bins = "Polycomb repressed, 3 kb bins")
panel_name <- function(epi, rs) {
  short <- unname(RS_SHORT[rs]); short[is.na(short)] <- rs[is.na(short)]
  paste0(epi, "\n", short)
}
PANEL_ORDER$panel <- panel_name(PANEL_ORDER$epitope, PANEL_ORDER$region_set)

df <- df %>%
  filter(!is.na(pct_bins_overlapping), peak_set %in% PEAK_LEVELS) %>%
  mutate(peak_set  = factor(peak_set,  levels = PEAK_LEVELS),
         bin_class = factor(bin_class, levels = CLASS_LEVELS),
         panel     = panel_name(epitope, region_set))

lev   <- PANEL_ORDER$panel[PANEL_ORDER$panel %in% unique(df$panel)]
extra <- setdiff(unique(df$panel), lev)
if (length(extra) > 0) {
  message("  WARNING: panels not named in PANEL_ORDER, appended last: ",
          paste(gsub("\n", " / ", extra), collapse = "; "))
  lev <- c(lev, extra)
}
df$panel <- factor(df$panel, levels = lev)
if (any(is.na(df$panel))) stop("panel assignment produced NA")
message("  panel order: ",
        paste(seq_along(lev), gsub("\n", " / ", lev), sep = ". ", collapse = "  |  "))

# --- layout ------------------------------------------------------------------
# x = peak set, fill = bin class: the four bars that must be compared sit
# adjacent on a shared baseline inside each cluster.
n_chip <- sum(grepl("^ChIP", PEAK_LEVELS))
shade  <- data.frame(xmin = 0.5, xmax = n_chip + 0.5, ymin = -Inf, ymax = Inf)

# Dashed reference at each peak set's negative-control rate: the level a bar
# would sit at if the peaks were indifferent to chromatin class. Per panel,
# since each panel has its own size-matched control.
base_df <- df %>% filter(bin_class == "Negative control") %>%
  mutate(xnum = as.numeric(peak_set)) %>%
  select(cell_line, panel, xnum, y = pct_bins_overlapping)

# Names the region set behind each bar, since the primary and control bars in a
# panel are drawn from different regions.
n_txt <- df %>% distinct(cell_line, panel, bin_class, region_desc, n_bins) %>%
  arrange(cell_line, panel, bin_class) %>%
  group_by(cell_line, panel) %>%
  summarise(lbl = paste0(
      "regions: ", paste(sprintf("%s %s", bin_class, format(n_bins, big.mark = ",", trim = TRUE)),
                         collapse = " | "),
      "\nsets: ", paste(unique(region_desc), collapse = "  vs  ")),
    .groups = "drop")

dodge <- position_dodge(width = 0.85)

p <- ggplot(df, aes(x = peak_set, y = pct_bins_overlapping, fill = bin_class)) +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey80", alpha = 0.30) +
  geom_col(position = dodge, width = 0.78, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct_bins_overlapping)),
            position = dodge, vjust = -0.35, size = fs(2.4), fontface = "bold") +
  geom_segment(data = base_df, aes(x = xnum - 0.44, xend = xnum + 0.44, y = y, yend = y),
               inherit.aes = FALSE, linetype = "22", linewidth = 0.5, colour = "grey20") +
  geom_text(data = n_txt, aes(x = 0.5, y = 112, label = lbl), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = fs(1.9), colour = "grey35", lineheight = 1.1) +
  scale_fill_manual(values = CLASS_COLORS, name = "Region class", drop = FALSE) +
  scale_x_discrete(labels = PEAK_SHORT) +
  scale_y_continuous(limits = c(-4, 125), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  facet_grid(cell_line ~ panel) +
  labs(
    x = NULL, y = "Regions overlapping a peak (%)",
    title = "Are method-enriched regions real signal within their own method, and are the peaks chromatin-class specific?",
    subtitle = sprintf(paste0("Within each peak set, compare the four adjacent bars  |  MACS2 peaks called per method and merged across both replicates ",
                              "  |  status: FDR < %s, |log2FC| > %s"), fdr, fc),
    caption = paste0(
      "Each cluster is ONE peak set. The first three bars are the percentage of ChIP-seq-enriched, non-differential and CUT&Tag-enriched\n",
      "regions overlapping it; the fourth is a NEGATIVE CONTROL region set where the mark should be largely absent, SIZE-MATCHED to the\n",
      "panel it sits in (500 bp bins of Polycomb-repressed promoters for the 500 bp promoter panels, whole 5 kb repressed promoters for\n",
      "the whole-promoter panels, Polycomb bins for H3K36me3, active gene bodies for H3K27me3). The dashed line marks that control rate.\n",
      "Peaks were called on the same libraries used in the differential test but WITHOUT reference to the other method, and ONE peak set\n",
      "per cell line x mark serves both promoter geometries.\n",
      "Two things are being asked at once. Real signal: a method's peaks should cover its own enriched regions far more than the other\n",
      "method's. Specificity: all three primary bars should sit well above the control, or the peak set is calling too much of the genome\n",
      "for any of the differences between the first three bars to mean much.\n",
      "READ WITHIN A PANEL, AND WITHIN A CLUSTER. The two leftmost panel pairs are the SAME marks over different region geometries: a\n",
      "500 bp bin has less chance of touching a peak than a 5 kb window purely because it is shorter, so the 500 bp panels sit lower\n",
      "throughout and that offset is geometry, not biology. MACS2 was used for both methods with matched settings (no input control,\n",
      "--nolambda, q < 1e-5), so the caller is not a confounder, but library depth differs and ChIP-seq was called as -f BAM with\n",
      "--extsize 200 against -f BAMPE for CUT&Tag. Shaded column is the ChIP-seq-derived peak set.")) +
  theme_bw(base_size = fs(11)) +
  theme(
    strip.text = element_text(size = fs(9), face = "bold", colour = "black", lineheight = 1.1),
    strip.background = element_rect(fill = "grey92", colour = "black"),
    axis.text.x = element_text(size = fs(8.5), colour = "black", lineheight = 0.95),
    axis.text.y = element_text(size = fs(10), colour = "black"),
    axis.title.y = element_text(size = fs(12), face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    plot.title = element_text(size = fs(13), face = "bold"),
    plot.subtitle = element_text(size = fs(9.5), face = "italic"),
    plot.caption = element_text(size = fs(6.5), hjust = 0, colour = "grey30", lineheight = 1.3),
    legend.position = "bottom",
    legend.text = element_text(size = fs(10)),
    legend.title = element_text(size = fs(11), face = "bold"))

# Panel width scales with the number of peak sets plotted, so two clusters give
# a narrower figure rather than one padded with whitespace.
n_panel <- length(levels(df$panel)); n_cell <- length(unique(df$cell_line))
n_ps    <- length(PEAK_LEVELS)
panel_w <- 1.6 * n_ps + 0.8
ggsave(out_pdf, p, width = panel_w * n_panel + 2, height = 5.0 * n_cell + 5.0,
       device = "pdf", bg = "white", limitsize = FALSE)
message("wrote ", out_pdf, "  (", n_panel, " panels x ", n_cell, " cell lines)")
RPLOT

if [[ "$DO_FIGURE" == "1" ]]; then
    echo "--- Figure ---"
    Rscript "$RS_DIR/plot_peak_overlap.R" "$OVERLAP_TSV" \
        "$FIG_DIR/peak_overlap_groupedByPeakSet.pdf" "$FC_THRESH" "$FDR_THRESH"
fi

echo "========================================================="
echo " Done."
echo "   Per-replicate peaks : $PEAKS_DIR"
echo "   Merged peaks        : $MERGED_DIR"
echo "   Derived regions     : $DERIVED_DIR  (500 bp repressed-promoter bins)"
echo "   Regions by status   : $BINS_DIR"
echo "   Overlap table       : $OVERLAP_TSV"
echo "   Figure              : $FIG_DIR"
echo "========================================================="