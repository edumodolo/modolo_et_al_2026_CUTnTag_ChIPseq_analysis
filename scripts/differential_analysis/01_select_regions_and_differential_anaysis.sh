#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# GROUND-TRUTH REGIONS, QUANTIFICATION AND DIFFERENTIAL ANALYSIS
# CUT&Tag vs ChIP-seq over ChromHMM-defined region sets
#
#   PART 1  Region sets from ChromHMM 18-state annotations
#   PART 2  Per-bin features: GC content, ATAC-seq, DNase-seq
#   PART 3  Sample groups and comparisons
#   PART 4  Fragment-level counting with deepTools multiBamSummary
#   PART 5  Differential analysis with DESeq2
#   PART 6  Pi-score ranked summary, one per comparison
#   PART 7  UCSC custom tracks coloured by differential status
#   PART 8  Manifests + parameter file consumed by 02_make_figures.sh
#
# METHODS SUMMARY
#   Regions. Four region sets are built per cell line from the ENCODE ChromHMM
#   18-state model. Active promoter = +/-2.5 kb around the TSS of active genes,
#   used both as whole 5 kb windows and tiled into 500 bp bins. Active gene body
#   = ChromHMM Tx/TxWk gene bodies in 3 kb bins with TSS +/-3 kb excluded.
#   Polycomb repressed = ChromHMM ReprPC/ReprPCWk domains in 3 kb bins. A gene
#   counts as active when >=20% of its body is Tx/TxWk and its TSS +/-1 kb
#   overlaps a ChromHMM TSS state. ENCODE blacklist regions are subtracted and
#   analysis is restricted to chr1-22, X, Y.
#
#   Quantification. Fragments are counted per bin directly from the filtered,
#   deduplicated BAMs with deepTools multiBamSummary
#
#   Differential analysis. DESeq2 with a ~condition design, ChIP-seq as the
#   reference level, tested as contrast (CnT vs ChIP). Positive log2FC is
#   therefore CUT&Tag enrichment and negative log2FC is ChIP-seq enrichment,
#   with no downstream sign flipping anywhere in the pipeline.
#
# NAMING
#   Samples are referred to throughout by their plot_label from the processing
#   pipeline's metrics table (e.g. chip_H3K4me3_K562_Bernstein_2017_rep1), not
#   by accession-based sample names, so the comparisons read as biology.
#   A sample GROUP name carries the shared prefix plus the year and replicate
#   tag of every member, e.g.
#       chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2
#   and every output file is named <CnT group>_VS_<ChIP group>_<region set>.
#
# Usage:
#   conda activate <env with bedtools deeptools R/DESeq2 [bedToBigBed]>
#   bash 01_differential_analysis.sh
#   bash 02_make_figures.sh          # reads this script's output
# =============================================================================

# --------------------------- USER CONFIG -------------------------------------
WORK_DIR="${WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output}"

REGIONS_DIR="$WORK_DIR/regions"
TRACKS_DIR="$WORK_DIR/tracks"
FEATURES_DIR="$WORK_DIR/features"
COUNTS_DIR="$WORK_DIR/counts"
DIFF_DIR="$WORK_DIR/differential"
SUMMARY_DIR="$WORK_DIR/differential_summaries"
RS_DIR="$WORK_DIR/Rscripts"
TMP="$WORK_DIR/_tmp"
mkdir -p "$REGIONS_DIR" "$TRACKS_DIR" "$FEATURES_DIR" "$COUNTS_DIR" "$DIFF_DIR" \
         "$SUMMARY_DIR" "$RS_DIR" "$TMP"

# Master metrics table written by the unified processing pipeline. Supplies the
# final BAM path and library type for every sample, looked up by plot_label.
MASTER_TSV="${MASTER_TSV:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/summary/Master_Sample_Metrics.tsv}"

GENE_BODIES="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/hg38_genomic_annotations/gene_body/canonical_geneBodies.hg38.bed"
BLACKLIST="$HOME/gpfs/genomes/bedfiles/hg38-blacklist.v2.bed"
CHROM_SIZES="/home/emodolo/gpfs/Homer/data/genomes/hg38/chrom.sizes"
GC_BW="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/GC_content/output/bigwigs/hg38_gc5Base.bw"

CHROMHMM_DIR="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/chromHMM_annotations/output"
declare -A CHROMHMM_BED=(
  [K562]="${CHROMHMM_DIR}/K562_hg38_ChromHMM_18state_ENCFF963KIA.bed"
  [MCF7]="${CHROMHMM_DIR}/MCF7_hg38_ChromHMM_18state_ENCFF985EWD.bed"
)

# Merged biological-replicate ATAC-seq tracks from the processing pipeline.
ATAC_BW_DIR="${ATAC_BW_DIR:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/alignment/bigwig}"
declare -A ATAC_BW=(
  [K562]="${ATAC_BW_DIR}/ATAC_NA_K562_ID173_merged_CPM.bw"
  # <<< FILL IN >>> MCF7 has two merged ATAC groups; keep one, delete the other.
  [MCF7]="${ATAC_BW_DIR}/ATAC_NA_MCF7_ID262_merged_CPM.bw"     # Tian 2023
  # [MCF7]="${ATAC_BW_DIR}/ATAC_NA_MCF7_ID940_merged_CPM.bw"   # Guan 2019
)

DNASE_BW_DIR="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/DNase-seq/output/bigwigs"
declare -A DNASE_BW=(
  [K562]="${DNASE_BW_DIR}/DNase_K562_ENCFF972GVB.bigWig"
  [MCF7]="${DNASE_BW_DIR}/DNase_MCF7_ENCFF038NFC.bigWig"
)

WEB_DIR="/home/emodolo/web/2026_modolo_et_al/ground_truth_regions_final_diff"
WEB_URL="http://homer.ucsd.edu/emodolo/2026_modolo_et_al/ground_truth_regions_final_diff"

THREADS="${THREADS:-24}"

# --- Differential thresholds: ONE pair of values for the whole manuscript -----
FC_THRESH=1.0
FDR_THRESH=0.05

# --- Region set geometry -----------------------------------------------------
PROM_UP=2500              # active promoter window, upstream of TSS
PROM_DN=2500              # active promoter window, downstream of TSS
PROM_FULL_SIZE=$((PROM_UP + PROM_DN))
PROM_BIN_SIZE=500         # promoter tiling for the binned promoter set
BROAD_BIN_SIZE=3000       # gene body and Polycomb tiling

# --- Active-gene criteria ----------------------------------------------------
TSS_PROMOTER_WINDOW=1000  # TSS +/- this must overlap a ChromHMM TSS state
MIN_BODY_TX_FRAC=0.20     # >= this fraction of the gene body must be Tx/TxWk
TSS_EXCLUDE_WINDOW=3000   # TSS +/- this is removed from the gene body set

# --- Pi-score summary --------------------------------------------------------
# Browser context padding added around each bin in the UCSC_Location column,
# and how many representative regions to tag per status group.
PAD_PROMOTER=1000
PAD_BROAD=10000
TARGET_REPS=2

# Counting mode. 0 = count a fragment in every bin it overlaps (deepTools
# default). 1 = --centerReads, which collapses each fragment towards its centre
# before counting and reduces the edge effect that longer ChIP-seq fragments
# have on 500 bp bins. Changing this changes every count, so re-run from clean.
CENTER_READS="${CENTER_READS:-0}"

CELL_LINES=(K562 MCF7)

# --- Region set definitions --------------------------------------------------
# Each epitope is tested on its primary region set. Promoter marks are ALSO
# tested on whole promoters; those results feed the heatmaps and browser-region
# selection and are deliberately kept out of the violin/volcano/MA figures.
declare -A EPI_REGION_SET=(
  [H3K4me3]=promoter_bins  [H3K27ac]=promoter_bins
  [H3K36me3]=genebody_bins [H3K27me3]=polycomb_bins
)
declare -A EPI_EXTRA_SET=( [H3K4me3]=promoter_full [H3K27ac]=promoter_full )

declare -A SET_WIDTH=(
  [promoter_bins]=$PROM_BIN_SIZE  [promoter_full]=$PROM_FULL_SIZE
  [genebody_bins]=$BROAD_BIN_SIZE [polycomb_bins]=$BROAD_BIN_SIZE
  [promoter_repressed]=$PROM_FULL_SIZE
)
declare -A SET_PAD=(
  [promoter_bins]=$PAD_PROMOTER   [promoter_full]=$PAD_PROMOTER
  [genebody_bins]=$PAD_BROAD      [polycomb_bins]=$PAD_BROAD
)
# promoter_repressed is a NEGATIVE CONTROL set: promoters of genes that are
# Polycomb-repressed and not transcribed. No epitope is tested on it, so it
# never enters EPI_REGION_SET and produces no differential result. It exists so
# that the peak-overlap validation and the control figures can ask
# what a method's peaks and signal look like where the mark should be absent.
ALL_REGION_SETS=(promoter_bins promoter_full genebody_bins polycomb_bins promoter_repressed)

# TSS +/- this must overlap ReprPC for a gene to count as repressed.
REPRESSED_TSS_WINDOW=1000

# =============================================================================
# PART 0: PREFLIGHT
# =============================================================================
echo "========================================================="
echo " Differential analysis: CUT&Tag vs ChIP-seq"
echo " Output: $WORK_DIR"
echo "========================================================="

HAVE_BB=1
command -v bedToBigBed >/dev/null 2>&1 || { echo "WARNING: bedToBigBed not in PATH -- skipping .bb tracks." >&2; HAVE_BB=0; }

miss=0
for t in bedtools multiBamSummary multiBigwigSummary Rscript awk sort; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: missing tool: $t" >&2; miss=1; }
done
for f in "$MASTER_TSV" "$GENE_BODIES" "$BLACKLIST" "$CHROM_SIZES" "$GC_BW"; do
    [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; miss=1; }
done
for CL in "${CELL_LINES[@]}"; do
    [[ -s "${CHROMHMM_BED[$CL]}" ]] || { echo "ERROR: missing ChromHMM bed for ${CL}: ${CHROMHMM_BED[$CL]}" >&2; miss=1; }
    [[ -s "${ATAC_BW[$CL]:-}"  ]] || echo "WARNING: ATAC bigWig missing for ${CL} -- ATAC figures will be skipped." >&2
    [[ -s "${DNASE_BW[$CL]:-}" ]] || echo "WARNING: DNase bigWig missing for ${CL} -- DNase figures will be skipped." >&2
done
[[ $miss -eq 0 ]] || { echo "Preflight failed." >&2; exit 1; }

# --- Master table lookups, by plot_label, resolved by header name ------------
col_index() {
    awk -F'\t' -v want="$1" 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$i); if($i==want){print i; exit}}}' "$MASTER_TSV"
}
COL_LABEL=$(col_index "plot_label")
COL_BAM=$(col_index "location_final_bam")
COL_ENDS=$(col_index "read_ends")
for v in COL_LABEL COL_BAM COL_ENDS; do
    [[ -n "${!v}" ]] || { echo "ERROR: ${v} not found in header of $MASTER_TSV" >&2; exit 1; }
done

lookup_by_label() {   # plot_label, column index -> value
    awk -F'\t' -v s="$1" -v n="$COL_LABEL" -v c="$2" \
        'NR>1{gsub(/\r/,"",$n); if($n==s){gsub(/\r/,"",$c); print $c; exit}}' "$MASTER_TSV"
}
bam_of()  { lookup_by_label "$1" "$COL_BAM"; }
ends_of() { lookup_by_label "$1" "$COL_ENDS"; }

# Display label: the group's shared identity without the replicate enumeration,
# used for figure strip labels. Keeps antibody suffixes (ab4729r1 -> ab4729) so
# that different antibodies from one study stay distinguishable.
short_label() { echo "$1" | sed -E 's/(_rep[0-9]+|r[0-9]+)$//'; }

# =============================================================================
# PART 1: REGION SETS
# =============================================================================
mainchr() { awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' "$@"; }
dedupe()  { awk -F'\t' '!seen[$1"\t"$2"\t"$3]++'; }

get_states() {   # ChromHMM bed, space-separated state names -> merged intervals
    local CHMM="$1" STATES="$2"
    awk -F'\t' -v st="$STATES" 'BEGIN{split(st,a," "); for(i in a) s[a[i]]=""} $4 in s {print $1"\t"$2"\t"$3}' "$CHMM" \
        | mainchr | sort -k1,1 -k2,2n | bedtools merge -i -
}

tile_regions() {   # named bed, bin size, blacklist, out
    bedtools makewindows -b "$1" -w "$2" -i srcwinnum \
      | awk -F'\t' -v w="$2" 'BEGIN{OFS="\t"} ($3-$2)==w' \
      | sort -k1,1 -k2,2n | bedtools intersect -a - -b "$3" -v | dedupe > "$4"
}

for CL in "${CELL_LINES[@]}"; do
    echo "--- Region sets: $CL ---"
    CHMM="${CHROMHMM_BED[$CL]}"

    BL="$TMP/${CL}_blacklist.bed"
    cut -f1-3 "$BLACKLIST" | tr -d '\r' | mainchr | sort -k1,1 -k2,2n | bedtools merge -i - > "$BL"

    PROM_ANCHOR="$TMP/${CL}_promoter_anchor.bed"
    get_states "$CHMM" "TssA TssFlnk TssFlnkU TssFlnkD" \
        | bedtools intersect -a - -b "$BL" -v | sort -k1,1 -k2,2n > "$PROM_ANCHOR"

    GENES="$TMP/${CL}_genes.bed"
    cut -f1-6 "$GENE_BODIES" | tr -d '\r' | mainchr | sort -k1,1 -k2,2n \
        | bedtools intersect -a - -b "$BL" -v | sort -k1,1 -k2,2n > "$GENES"

    ACT_BLK="$TMP/${CL}_txBlocks.bed"
    get_states "$CHMM" "Tx TxWk" > "$ACT_BLK"
    ACT_BODY_NAMES="$TMP/${CL}_actBody.names"
    bedtools coverage -a "$GENES" -b "$ACT_BLK" \
      | awk -F'\t' -v f="$MIN_BODY_TX_FRAC" '($NF+0) >= f {print $4}' | sort -u > "$ACT_BODY_NAMES"

    TSSWIN="$TMP/${CL}_tssWindows.bed"
    awk -F'\t' -v w="$TSS_PROMOTER_WINDOW" 'BEGIN{OFS="\t"}
        { chr=$1; s=$2; e=$3; name=$4; str=$6;
          if(str=="-") tss=e; else tss=s;
          ws=tss-w; if(ws<0)ws=0; we=tss+w; if(we<=ws) we=ws+1;
          print chr, ws, we, name, ".", str }' "$GENES" | sort -k1,1 -k2,2n > "$TSSWIN"
    PROMACT_NAMES="$TMP/${CL}_promActive.names"
    bedtools intersect -a "$TSSWIN" -b "$PROM_ANCHOR" -u | cut -f4 | sort -u > "$PROMACT_NAMES"

    ACTIVE_NAMES="$TMP/${CL}_active.names"
    comm -12 "$ACT_BODY_NAMES" "$PROMACT_NAMES" > "$ACTIVE_NAMES"
    ACT_GENES="$TMP/${CL}_activeGenes.bed"
    awk -F'\t' 'NR==FNR{k[$1]=1; next} ($4 in k)' "$ACTIVE_NAMES" "$GENES" | sort -k1,1 -k2,2n > "$ACT_GENES"

    # ---- Active promoter: whole windows, then the same windows tiled --------
    PROM_FULL_RAW="$TMP/${CL}_promoter_full_raw.bed"
    awk -F'\t' -v up="$PROM_UP" -v dn="$PROM_DN" '
        NR==FNR { size[$1]=$2; next }
        { chr=$1; s=$2; e=$3; str=$6;
          if (str=="-") { tss=e; ws=tss-dn; we=tss+up } else { tss=s; ws=tss-up; we=tss+dn }
          if (ws<0) ws=0;
          if ((chr in size) && we>size[chr]) we=size[chr];
          if (we>ws) print chr"\t"ws"\t"we }' "$CHROM_SIZES" "$ACT_GENES" \
      | awk -F'\t' -v w="$PROM_FULL_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' \
      | sort -k1,1 -k2,2n -u | bedtools intersect -a - -b "$BL" -v | dedupe > "$PROM_FULL_RAW"

    awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1,$2,$3, cl"_promoter_"NR}' "$PROM_FULL_RAW" \
        > "$REGIONS_DIR/${CL}_promoter_full.bed"
    tile_regions "$REGIONS_DIR/${CL}_promoter_full.bed" "$PROM_BIN_SIZE" "$BL" \
        "$REGIONS_DIR/${CL}_promoter_bins.bed"

    # ---- Active gene body, TSS-proximal signal excluded ---------------------
    ACT_TSS_EXCL="$TMP/${CL}_activeTSS_excl.bed"
    awk -F'\t' 'NR==FNR{k[$1]=1; next} ($4 in k)' "$PROMACT_NAMES" "$GENES" \
      | awk -F'\t' -v w="$TSS_EXCLUDE_WINDOW" 'BEGIN{OFS="\t"}
          { chr=$1; s=$2; e=$3; str=$6;
            if(str=="-") tss=e; else tss=s;
            ws=tss-w; if(ws<0)ws=0; we=tss+w; if(we<=ws) we=ws+1;
            print chr, ws, we }' \
      | sort -k1,1 -k2,2n | bedtools merge -i - > "$ACT_TSS_EXCL"

    ACT_NAMED="$TMP/${CL}_genebody_named.bed"
    sort -k1,1 -k2,2n "$ACT_GENES" | bedtools merge -i - \
      | awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1,$2,$3, cl"_geneBody_"NR}' > "$ACT_NAMED"
    tile_regions "$ACT_NAMED" "$BROAD_BIN_SIZE" "$BL" "$TMP/${CL}_genebody_tiled.bed"
    bedtools intersect -a "$TMP/${CL}_genebody_tiled.bed" -b "$ACT_TSS_EXCL" -v \
        > "$REGIONS_DIR/${CL}_genebody_bins.bed"

    # ---- Polycomb repressed -------------------------------------------------
    # Bivalent states are kept only where they sit inside a core ReprPC domain,
    # so isolated bivalent promoters do not enter the repressed set.
    REPR_ALL="$TMP/${CL}_repr_all.bed";  get_states "$CHMM" "ReprPC ReprPCWk TssBiv EnhBiv" > "$REPR_ALL"
    REPR_CORE="$TMP/${CL}_repr_core.bed"; get_states "$CHMM" "ReprPC ReprPCWk" > "$REPR_CORE"
    REPR_NAMED="$TMP/${CL}_repr_named.bed"
    bedtools intersect -a "$REPR_ALL" -b "$REPR_CORE" -u | sort -k1,1 -k2,2n \
      | awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1,$2,$3, cl"_polycomb_"NR}' > "$REPR_NAMED"
    tile_regions "$REPR_NAMED" "$BROAD_BIN_SIZE" "$BL" "$REGIONS_DIR/${CL}_polycomb_bins.bed"

    # ---- Repressed promoters: the negative control --------------------------
    # A gene qualifies when its TSS +/- REPRESSED_TSS_WINDOW overlaps a ReprPC
    # domain, does NOT overlap any active-promoter state, and its body carries
    # no Tx/TxWk. All three conditions are needed: ReprPC proximity alone would
    # admit bivalent promoters that are also active, and dropping the Tx filter
    # would admit transcribed genes whose promoter merely sits near a domain.
    REPR_ONLY="$TMP/${CL}_state_ReprPConly.bed"; get_states "$CHMM" "ReprPC" > "$REPR_ONLY"
    R_NAMES="$TMP/${CL}_repressed.names"
    bedtools intersect -a "$TSSWIN" -b "$REPR_ONLY" -u | cut -f4 | sort -u > "$R_NAMES"
    TX_NAMES="$TMP/${CL}_txBody.names"
    bedtools intersect -a "$GENES" -b "$ACT_BLK" -u | cut -f4 | sort -u > "$TX_NAMES"
    NOT_ACT="$TMP/${CL}_repressed_notActive.names"
    comm -23 "$R_NAMES" "$PROMACT_NAMES" > "$NOT_ACT"
    REPRESSED_NAMES="$TMP/${CL}_repressedFinal.names"
    comm -23 "$NOT_ACT" "$TX_NAMES" > "$REPRESSED_NAMES"

    REPRESSED_GENES="$TMP/${CL}_repressedGenes.bed"
    awk -F'\t' 'NR==FNR{k[$1]=1; next} ($4 in k)' "$REPRESSED_NAMES" "$GENES" \
        | sort -k1,1 -k2,2n > "$REPRESSED_GENES"

    awk -F'\t' -v up="$PROM_UP" -v dn="$PROM_DN" '
        NR==FNR { size[$1]=$2; next }
        { chr=$1; s=$2; e=$3; str=$6;
          if (str=="-") { tss=e; ws=tss-dn; we=tss+up } else { tss=s; ws=tss-up; we=tss+dn }
          if (ws<0) ws=0;
          if ((chr in size) && we>size[chr]) we=size[chr];
          if (we>ws) print chr"\t"ws"\t"we }' "$CHROM_SIZES" "$REPRESSED_GENES" \
      | awk -F'\t' -v w="$PROM_FULL_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' \
      | sort -k1,1 -k2,2n -u | bedtools intersect -a - -b "$BL" -v | dedupe \
      | awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1,$2,$3, cl"_repressedPromoter_"NR}' \
      > "$REGIONS_DIR/${CL}_promoter_repressed.bed"

    for RS in "${ALL_REGION_SETS[@]}"; do
        printf "    %-16s %8d intervals of %s bp\n" "$RS" \
            "$(wc -l < "$REGIONS_DIR/${CL}_${RS}.bed")" "${SET_WIDTH[$RS]}"
    done
done

# =============================================================================
# PART 2: PER-BIN FEATURES (GC, ATAC-seq, DNase-seq)
# =============================================================================
# One table per cell line x region set x feature. Because the table is defined
# by the region set and not by the comparison, every comparison over the same
# bins reads identical feature values, and the Z-score references built in the
# figure script are automatically restricted to one bin size.
summarise_bigwig() {   # bigwig, region bed, out tab
    local BW="$1" BED="$2" OUT="$3"
    [[ -s "$BW" ]] || return 0
    [[ -s "$OUT" ]] && return 0
    local NPZ="${OUT%.tab}.npz"
    multiBigwigSummary BED-file -b "$BW" --BED "$BED" -o "$NPZ" \
        --outRawCounts "$OUT" -p "$THREADS" >/dev/null 2>&1
    rm -f "$NPZ"
}

echo "--- Per-bin features ---"
for CL in "${CELL_LINES[@]}"; do
    for RS in "${ALL_REGION_SETS[@]}"; do
        BED="$REGIONS_DIR/${CL}_${RS}.bed"
        echo "  ${CL} ${RS}"
        summarise_bigwig "$GC_BW"              "$BED" "$FEATURES_DIR/gc_${CL}_${RS}.tab"
        summarise_bigwig "${ATAC_BW[$CL]:-}"   "$BED" "$FEATURES_DIR/atac_${CL}_${RS}.tab"
        summarise_bigwig "${DNASE_BW[$CL]:-}"  "$BED" "$FEATURES_DIR/dnase_${CL}_${RS}.tab"
    done
done

# =============================================================================
# PART 3: SAMPLE SAMPLE_GROUPS AND COMPARISONS
# =============================================================================
# Samples are listed by plot_label. A group's key is its full descriptive name:
# the shared prefix plus the year and replicate tag of every member. That key
# is what appears in output filenames, so a result file states exactly which
# replicates produced it.
# NOTE: this array must NOT be called GROUPS. That name is a bash built-in
# special variable holding the user's group IDs, and `declare -A GROUPS` fails
# with "cannot convert indexed to associative array", leaving it empty.
declare -A SAMPLE_GROUPS=(
  # --- ChIP-seq ---
  [chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3]="chip_H3K27ac_K562_Bernstein_2011_rep2 chip_H3K27ac_K562_Bernstein_2022_rep3"
  [chip_H3K27ac_K562_Tak_2020_rep1_2020_rep2]="chip_H3K27ac_K562_Tak_2020_rep1 chip_H3K27ac_K562_Tak_2020_rep2"
  [chip_H3K27ac_MCF7_Bernstein_2017_rep1_2017_rep2]="chip_H3K27ac_MCF7_Bernstein_2017_rep1 chip_H3K27ac_MCF7_Bernstein_2017_rep2"
  [chip_H3K27me3_K562_Bernstein_2011_rep1_2022_rep3]="chip_H3K27me3_K562_Bernstein_2011_rep1 chip_H3K27me3_K562_Bernstein_2022_rep3"
  [chip_H3K27me3_K562_Dou_2023_rep1_2023_rep2]="chip_H3K27me3_K562_Dou_2023_rep1 chip_H3K27me3_K562_Dou_2023_rep2"
  [chip_H3K27me3_MCF7_Bernstein_2017_rep1_2017_rep2]="chip_H3K27me3_MCF7_Bernstein_2017_rep1 chip_H3K27me3_MCF7_Bernstein_2017_rep2"
  [chip_H3K27me3_MCF7_Farnham_2012_rep1_2012_rep2]="chip_H3K27me3_MCF7_Farnham_2012_rep1 chip_H3K27me3_MCF7_Farnham_2012_rep2"
  [chip_H3K36me3_K562_Bernstein_2011_rep2_2022_rep3]="chip_H3K36me3_K562_Bernstein_2011_rep2 chip_H3K36me3_K562_Bernstein_2022_rep3"
  [chip_H3K36me3_MCF7_Bernstein_2021_rep1_2021_rep2]="chip_H3K36me3_MCF7_Bernstein_2021_rep1 chip_H3K36me3_MCF7_Bernstein_2021_rep2"
  [chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2]="chip_H3K4me3_K562_Bernstein_2017_rep1 chip_H3K4me3_K562_Bernstein_2017_rep2"
  [chip_H3K4me3_K562_Stamatoyannopoulous_2012_rep1_2012_rep2]="chip_H3K4me3_K562_Stamatoyannopoulous_2012_rep1 chip_H3K4me3_K562_Stamatoyannopoulous_2012_rep2"
  [chip_H3K4me3_MCF7_Bernstein_2017_rep1_2017_rep2]="chip_H3K4me3_MCF7_Bernstein_2017_rep1 chip_H3K4me3_MCF7_Bernstein_2017_rep2"
  [chip_H3K4me3_MCF7_Stamatoyannopoulous_2012_rep1_2012_rep2]="chip_H3K4me3_MCF7_Stamatoyannopoulous_2012_rep1 chip_H3K4me3_MCF7_Stamatoyannopoulous_2012_rep2"

  # --- CUT&Tag ---
  [CnT_H3K27ac_K562_Abbasova_2025_ab177r1_2025_ab177r2]="CnT_H3K27ac_K562_Abbasova_2025_ab177r1 CnT_H3K27ac_K562_Abbasova_2025_ab177r2"
  [CnT_H3K27ac_K562_Abbasova_2025_ab4729r1_2025_ab4729r2]="CnT_H3K27ac_K562_Abbasova_2025_ab4729r1 CnT_H3K27ac_K562_Abbasova_2025_ab4729r2"
  [CnT_H3K27ac_K562_Abbasova_2025_abdiagr1_2025_abdiagr2]="CnT_H3K27ac_K562_Abbasova_2025_abdiagr1 CnT_H3K27ac_K562_Abbasova_2025_abdiagr2"
  [CnT_H3K27ac_K562_KayaOkur_2019_rep1_2019_rep2]="CnT_H3K27ac_K562_KayaOkur_2019_rep1 CnT_H3K27ac_K562_KayaOkur_2019_rep2"
  [CnT_H3K27ac_MCF7_Fischer_2025_rep1_2025_rep2]="CnT_H3K27ac_MCF7_Fischer_2025_rep1 CnT_H3K27ac_MCF7_Fischer_2025_rep2"
  [CnT_H3K27ac_MCF7_Tian_2023_rep1_2023_rep2]="CnT_H3K27ac_MCF7_Tian_2023_rep1 CnT_H3K27ac_MCF7_Tian_2023_rep2"
  [CnT_H3K27me3_K562_Abbasova_2025_rep1_2025_rep2]="CnT_H3K27me3_K562_Abbasova_2025_rep1 CnT_H3K27me3_K562_Abbasova_2025_rep2"
  [CnT_H3K27me3_K562_KayaOkur_2019_rep1_2019_rep2]="CnT_H3K27me3_K562_KayaOkur_2019_rep1 CnT_H3K27me3_K562_KayaOkur_2019_rep2"
  [CnT_H3K27me3_MCF7_Tian_2023_rep1_2023_rep2_2023_rep3_2023_rep4]="CnT_H3K27me3_MCF7_Tian_2023_rep1 CnT_H3K27me3_MCF7_Tian_2023_rep2 CnT_H3K27me3_MCF7_Tian_2023_rep3 CnT_H3K27me3_MCF7_Tian_2023_rep4"
  [CnT_H3K36me3_K562_Wu_2025_rep1_2025_rep2_2025_rep3_2025_rep4_2025_rep5]="CnT_H3K36me3_K562_Wu_2025_rep1 CnT_H3K36me3_K562_Wu_2025_rep2 CnT_H3K36me3_K562_Wu_2025_rep3 CnT_H3K36me3_K562_Wu_2025_rep4 CnT_H3K36me3_K562_Wu_2025_rep5"
  [CnT_H3K36me3_MCF7_Tian_2023_rep1_2023_rep2]="CnT_H3K36me3_MCF7_Tian_2023_rep1 CnT_H3K36me3_MCF7_Tian_2023_rep2"
  [CnT_H3K4me3_K562_KayaOkur_2019_rep1_2019_rep2]="CnT_H3K4me3_K562_KayaOkur_2019_rep1 CnT_H3K4me3_K562_KayaOkur_2019_rep2"
  [CnT_H3K4me3_K562_KayaOkur_2020_rep1_2020_rep2]="CnT_H3K4me3_K562_KayaOkur_2020_rep1 CnT_H3K4me3_K562_KayaOkur_2020_rep2"
  [CnT_H3K4me3_MCF7_Tian_2023_rep1_2023_rep2]="CnT_H3K4me3_MCF7_Tian_2023_rep1 CnT_H3K4me3_MCF7_Tian_2023_rep2"
)

# epitope | cell | CUT&Tag group | ChIP-seq group | figure set
#   main = appears in both the main and supplementary figures
#   supp = supplementary figure only
#
# Each CUT&Tag sample is compared against every independent ChIP-seq source
# available for that epitope and cell line, so a difference that depends on one
# ChIP-seq dataset is visible as such. Where both sources resolve to the SAME
# ChIP-seq group the comparison appears once, not twice: K562 and MCF7 H3K36me3
# and MCF7 H3K27ac have only one ChIP-seq source between them.
COMPARISONS=(
  "H3K4me3|K562|CnT_H3K4me3_K562_KayaOkur_2019_rep1_2019_rep2|chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2|supp"
  "H3K4me3|K562|CnT_H3K4me3_K562_KayaOkur_2020_rep1_2020_rep2|chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2|main"
  "H3K4me3|K562|CnT_H3K4me3_K562_KayaOkur_2019_rep1_2019_rep2|chip_H3K4me3_K562_Stamatoyannopoulous_2012_rep1_2012_rep2|supp"
  "H3K4me3|K562|CnT_H3K4me3_K562_KayaOkur_2020_rep1_2020_rep2|chip_H3K4me3_K562_Stamatoyannopoulous_2012_rep1_2012_rep2|supp"
  "H3K4me3|MCF7|CnT_H3K4me3_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K4me3_MCF7_Bernstein_2017_rep1_2017_rep2|supp"
  "H3K4me3|MCF7|CnT_H3K4me3_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K4me3_MCF7_Stamatoyannopoulous_2012_rep1_2012_rep2|supp"

  "H3K27ac|K562|CnT_H3K27ac_K562_KayaOkur_2019_rep1_2019_rep2|chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3|supp"
  "H3K27ac|K562|CnT_H3K27ac_K562_Abbasova_2025_ab4729r1_2025_ab4729r2|chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3|supp"
  "H3K27ac|K562|CnT_H3K27ac_K562_Abbasova_2025_ab177r1_2025_ab177r2|chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3|main"
  "H3K27ac|K562|CnT_H3K27ac_K562_Abbasova_2025_abdiagr1_2025_abdiagr2|chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3|supp"
  "H3K27ac|K562|CnT_H3K27ac_K562_KayaOkur_2019_rep1_2019_rep2|chip_H3K27ac_K562_Tak_2020_rep1_2020_rep2|supp"
  "H3K27ac|K562|CnT_H3K27ac_K562_Abbasova_2025_ab4729r1_2025_ab4729r2|chip_H3K27ac_K562_Tak_2020_rep1_2020_rep2|supp"
  "H3K27ac|K562|CnT_H3K27ac_K562_Abbasova_2025_ab177r1_2025_ab177r2|chip_H3K27ac_K562_Tak_2020_rep1_2020_rep2|supp"
  "H3K27ac|K562|CnT_H3K27ac_K562_Abbasova_2025_abdiagr1_2025_abdiagr2|chip_H3K27ac_K562_Tak_2020_rep1_2020_rep2|supp"
  "H3K27ac|MCF7|CnT_H3K27ac_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K27ac_MCF7_Bernstein_2017_rep1_2017_rep2|supp"
  "H3K27ac|MCF7|CnT_H3K27ac_MCF7_Fischer_2025_rep1_2025_rep2|chip_H3K27ac_MCF7_Bernstein_2017_rep1_2017_rep2|supp"

  "H3K27me3|K562|CnT_H3K27me3_K562_KayaOkur_2019_rep1_2019_rep2|chip_H3K27me3_K562_Bernstein_2011_rep1_2022_rep3|supp"
  "H3K27me3|K562|CnT_H3K27me3_K562_Abbasova_2025_rep1_2025_rep2|chip_H3K27me3_K562_Bernstein_2011_rep1_2022_rep3|main"
  "H3K27me3|K562|CnT_H3K27me3_K562_KayaOkur_2019_rep1_2019_rep2|chip_H3K27me3_K562_Dou_2023_rep1_2023_rep2|supp"
  "H3K27me3|K562|CnT_H3K27me3_K562_Abbasova_2025_rep1_2025_rep2|chip_H3K27me3_K562_Dou_2023_rep1_2023_rep2|supp"
  "H3K27me3|MCF7|CnT_H3K27me3_MCF7_Tian_2023_rep1_2023_rep2_2023_rep3_2023_rep4|chip_H3K27me3_MCF7_Bernstein_2017_rep1_2017_rep2|supp"
  "H3K27me3|MCF7|CnT_H3K27me3_MCF7_Tian_2023_rep1_2023_rep2_2023_rep3_2023_rep4|chip_H3K27me3_MCF7_Farnham_2012_rep1_2012_rep2|supp"

  "H3K36me3|K562|CnT_H3K36me3_K562_Wu_2025_rep1_2025_rep2_2025_rep3_2025_rep4_2025_rep5|chip_H3K36me3_K562_Bernstein_2011_rep2_2022_rep3|main"
  "H3K36me3|MCF7|CnT_H3K36me3_MCF7_Tian_2023_rep1_2023_rep2|chip_H3K36me3_MCF7_Bernstein_2021_rep1_2021_rep2|supp"
)

# --- Validate every label before any compute happens -------------------------
# A plot_label that is absent from the master table would otherwise surface much
# later as an empty BAM path or a missing count column.
echo "--- Validating sample groups ---"
bad=0
for g in "${!SAMPLE_GROUPS[@]}"; do
    for lab in ${SAMPLE_GROUPS[$g]}; do
        e=$(ends_of "$lab")
        b=$(bam_of "$lab")
        if [[ -z "$e" || -z "$b" ]]; then
            echo "  ERROR: plot_label not found in ${MASTER_TSV}: '${lab}' (group ${g})" >&2; bad=1
        elif [[ ! -s "$b" ]]; then
            echo "  ERROR: BAM missing for '${lab}': ${b}" >&2; bad=1
        fi
    done
done
for rec in "${COMPARISONS[@]}"; do
    IFS='|' read -r EPI CELL CNT_G CHIP_G FIGSET <<< "$rec"
    [[ -n "${SAMPLE_GROUPS[$CNT_G]:-}"  ]] || { echo "  ERROR: undefined CUT&Tag group: ${CNT_G}" >&2; bad=1; }
    [[ -n "${SAMPLE_GROUPS[$CHIP_G]:-}" ]] || { echo "  ERROR: undefined ChIP-seq group: ${CHIP_G}" >&2; bad=1; }
done
[[ $bad -eq 0 ]] || { echo "Fix the sample groups before continuing." >&2; exit 1; }
echo "  ${#SAMPLE_GROUPS[@]} groups, ${#COMPARISONS[@]} comparisons validated"

# --- Expand every comparison into the region sets it must be run on ----------
declare -a JOBS=()
declare -A NEEDED=()    # "cell|region_set|read_ends" -> space-separated labels
for rec in "${COMPARISONS[@]}"; do
    IFS='|' read -r EPI CELL CNT_G CHIP_G FIGSET <<< "$rec"
    sets=("${EPI_REGION_SET[$EPI]}")
    [[ -n "${EPI_EXTRA_SET[$EPI]:-}" ]] && sets+=("${EPI_EXTRA_SET[$EPI]}")
    for RS in "${sets[@]}"; do
        JOBS+=("${EPI}|${CELL}|${CNT_G}|${CHIP_G}|${FIGSET}|${RS}")
        for lab in ${SAMPLE_GROUPS[$CHIP_G]} ${SAMPLE_GROUPS[$CNT_G]}; do
            ends=$(ends_of "$lab")
            key="${CELL}|${RS}|${ends}"
            [[ " ${NEEDED[$key]:-} " == *" $lab "* ]] || NEEDED[$key]+="$lab "
        done
    done
done
echo "  ${#COMPARISONS[@]} comparisons -> ${#JOBS[@]} region-set jobs"

# =============================================================================
# PART 4: FRAGMENT COUNTING (deepTools multiBamSummary)
# =============================================================================
# Counting is batched by cell line x region set x library type, so each sample
# is counted once over a given region set no matter how many comparisons reuse
# it. Single-end batches pass an explicit 200 bp fragment estimate, since a
# single-end read never measured one. The two library types must be separate
# calls because --extendReads is one global setting per invocation.
#
# deepTools counts BOTH mates of a pair even with --extendReads, so
# --samFlagInclude 64 (first mate only) is required to count each fragment
# once. 

echo "--- Counting fragments per bin ---"
for key in "${!NEEDED[@]}"; do
    IFS='|' read -r CELL RS ENDS <<< "$key"
    OUT="$COUNTS_DIR/counts_${CELL}_${RS}_${ENDS}.tab"
    [[ -s "$OUT" ]] && { echo "  cached: $(basename "$OUT")"; continue; }

    read -ra labels <<< "${NEEDED[$key]}"
    bams=()
    for lab in "${labels[@]}"; do
        b=$(bam_of "$lab")
        [[ -s "${b}.bai" || -s "${b%.bam}.bai" ]] || echo "  WARNING: no index for ${b}" >&2
        bams+=("$b")
    done

    if [[ "$ENDS" == "Single" ]]; then
        opts=( --extendReads 200 )
    else
        opts=( --extendReads --samFlagInclude 64 )
    fi
    [[ "$CENTER_READS" == "1" ]] && opts+=( --centerReads )

    echo "  ${CELL} ${RS} ${ENDS}: ${#bams[@]} samples"
    NPZ="${OUT%.tab}.npz"
    multiBamSummary BED-file --BED "$REGIONS_DIR/${CELL}_${RS}.bed" \
        -b "${bams[@]}" --labels "${labels[@]}" \
        "${opts[@]}" -p "$THREADS" -o "$NPZ" --outRawCounts "$OUT" >/dev/null 2>&1
    rm -f "$NPZ"
    [[ -s "$OUT" ]] || { echo "ERROR: multiBamSummary produced no counts: $OUT" >&2; exit 1; }
done

# =============================================================================
# PART 5: DESeq2
# =============================================================================
cat > "$RS_DIR/run_deseq2.R" <<'RDESEQ'
#!/usr/bin/env Rscript
# Usage: run_deseq2.R <region_bed> <samples_tsv> <out_tsv> <counts_tab> [...]
#
# samples_tsv: two columns, plot_label and condition (ChIP or CnT), no header.
# The count matrix is assembled by name, never by column position, so the order
# of samples in the file and the order of columns in the count tables cannot
# disagree. ChIP is the reference level and the tested contrast is CnT vs ChIP,
# so positive log2FoldChange means CUT&Tag enrichment.
suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr) })

args <- commandArgs(trailingOnly = TRUE)
region_bed <- args[1]; samples_tsv <- args[2]; out_tsv <- args[3]
count_files <- args[-(1:3)]

# deepTools writes its header as a comment line of quoted names.
read_counts_wide <- function(path) {
  hdr <- readLines(path, n = 1)
  nm  <- gsub("^'|'$", "", trimws(strsplit(sub("^#", "", hdr), "\t")[[1]]))
  df  <- suppressMessages(read_tsv(path, comment = "#", col_names = nm, show_col_types = FALSE))
  names(df)[1:3] <- c("Chr", "Start", "End")
  df$Start <- as.integer(df$Start); df$End <- as.integer(df$End)
  df
}

regions <- suppressMessages(read_tsv(region_bed, col_names = FALSE, show_col_types = FALSE))
regions <- regions[, 1:4]; colnames(regions) <- c("Chr", "Start", "End", "BinID")
regions$Start <- as.integer(regions$Start); regions$End <- as.integer(regions$End)

samples <- suppressMessages(read_tsv(samples_tsv, col_names = c("sample", "condition"), show_col_types = FALSE))

mat <- regions
for (f in count_files) {
  cw <- read_counts_wide(f)
  keep <- intersect(samples$sample, colnames(cw))
  if (length(keep) == 0) next
  mat <- dplyr::left_join(mat, cw[, c("Chr", "Start", "End", keep)], by = c("Chr", "Start", "End"))
}

missing <- setdiff(samples$sample, colnames(mat))
if (length(missing) > 0) stop("no counts found for: ", paste(missing, collapse = ", "))

counts <- as.data.frame(mat[, samples$sample, drop = FALSE])
counts[] <- lapply(counts, function(x) { x <- as.numeric(x); x[is.na(x)] <- 0; as.integer(round(x)) })
rownames(counts) <- mat$BinID

# Bins with no coverage anywhere carry no information and only cost power in
# the multiple-testing correction.
keep_rows <- rowSums(counts) > 0
message("  bins: ", nrow(counts), " total, ", sum(keep_rows), " with non-zero coverage")
counts <- counts[keep_rows, , drop = FALSE]
meta   <- mat[keep_rows, c("BinID", "Chr", "Start", "End")]

colData <- data.frame(condition = factor(samples$condition, levels = c("ChIP", "CnT")))
rownames(colData) <- samples$sample
stopifnot(identical(rownames(colData), colnames(counts)))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds, fitType = "local", quiet = TRUE)
res <- results(dds, contrast = c("condition", "CnT", "ChIP"))

res$padj[is.na(res$padj)] <- 1
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0

out <- data.frame(meta,
                  baseMean = res$baseMean,
                  log2FC   = res$log2FoldChange,
                  lfcSE    = res$lfcSE,
                  pvalue   = res$pvalue,
                  padj     = res$padj,
                  check.names = FALSE)

norm <- as.data.frame(counts(dds, normalized = TRUE))
colnames(norm) <- paste0("norm.", colnames(norm))
write_tsv(cbind(out, norm), out_tsv)
RDESEQ

# =============================================================================
# PART 6: PI-SCORE RANKED SUMMARY
# =============================================================================
cat > "$RS_DIR/rank_by_pi_score.R" <<'RPI'
#!/usr/bin/env Rscript
# Usage: rank_by_pi_score.R <diff_tsv> <out_csv> <fc> <fdr> <pad_bp> <n_reps>
#
# Ranks every bin by a signed pi-score,
#     pi = log2FC * -log10(padj)
# which rewards bins that are both strongly and confidently shifted -- the
# corners of the volcano plot -- rather than bins that are merely significant
# or merely large in magnitude.
#
# Rows are sorted by pi-score ASCENDING, so the strongest ChIP-seq-enriched
# bins are at the top of the file and the strongest CUT&Tag-enriched bins at
# the bottom, with non-differential bins in the middle.
suppressPackageStartupMessages({ library(readr); library(dplyr) })
options(scipen = 999)

a <- commandArgs(trailingOnly = TRUE)
diff_tsv <- a[1]; out_csv <- a[2]
fc <- as.numeric(a[3]); fdr <- as.numeric(a[4])
pad <- as.numeric(a[5]); n_reps <- as.numeric(a[6])

df <- suppressMessages(read_tsv(diff_tsv, show_col_types = FALSE))
if (nrow(df) == 0) { message("  empty differential table, skipping"); quit(save = "no") }

# padj can underflow to exactly 0; flooring keeps -log10 finite so those bins
# stay rankable instead of becoming Inf.
PADJ_FLOOR <- 1e-300

out <- df %>%
  filter(!is.na(log2FC), !is.na(padj)) %>%
  mutate(
    pi_score = log2FC * -log10(pmax(padj, PADJ_FLOOR)),
    Status = case_when(
      log2FC >  fc & padj < fdr ~ "CUT&Tag Enriched",
      log2FC < -fc & padj < fdr ~ "ChIP-seq Enriched",
      TRUE ~ "Non-differential"),
    # Start is 0-based half-open; UCSC positions are 1-based inclusive.
    UCSC_Location = paste0(Chr, ":", pmax(1, Start + 1 - pad), "-", End + pad)) %>%
  group_by(Status) %>%
  arrange(pi_score, .by_group = TRUE) %>%
  mutate(
    rn = row_number(),
    dist_to_med = abs(rn - median(seq_len(n()))),
    Representative_Region = case_when(
      Status == "ChIP-seq Enriched" & rn <= n_reps ~ "Top ChIP-seq enriched",
      Status == "CUT&Tag Enriched"  &
        rank(desc(pi_score), ties.method = "first") <= n_reps ~ "Top CUT&Tag enriched",
      Status == "Non-differential"  &
        rank(dist_to_med, ties.method = "first") <= n_reps ~ "Middle non-differential",
      TRUE ~ "")) %>%
  ungroup() %>%
  arrange(pi_score) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, BinID, UCSC_Location, Chr, Start, End,
         baseMean, log2FC, padj, pi_score, Status, Representative_Region)

write_csv(out, out_csv)
message(sprintf("  ranked %d bins (ChIP-enriched %d, CUT&Tag-enriched %d, non-differential %d)",
                nrow(out),
                sum(out$Status == "ChIP-seq Enriched"),
                sum(out$Status == "CUT&Tag Enriched"),
                sum(out$Status == "Non-differential")))
RPI

# =============================================================================
# PART 7: UCSC TRACKS
# =============================================================================
# Positive log2FC = CUT&Tag enriched = red; negative = ChIP-seq enriched = blue.
# Region coordinates are already 0-based half-open, straight from the BED.
generate_colored_track() {
    local DIFF="$1" PREFIX="$2" DESC="$3"
    local BB="$TRACKS_DIR/${PREFIX}.bb" B9="$TMP/${PREFIX}_bed9.bed"
    awk -F'\t' -v fc="$FC_THRESH" -v fdr="$FDR_THRESH" '
    NR==1 { for(i=1;i<=NF;i++){ if($i=="Chr")c_chr=i; if($i=="Start")c_s=i; if($i=="End")c_e=i;
                                if($i=="BinID")c_id=i; if($i=="log2FC")c_fc=i; if($i=="padj")c_p=i } next }
    {
      p=$c_p+0; l=$c_fc+0
      if (p < fdr && l >  fc) col="255,59,59";       # CUT&Tag enriched
      else if (p < fdr && l < -fc) col="57,54,255";  # ChIP-seq enriched
      else col="178,178,178"
      print $c_chr"\t"$c_s"\t"$c_e"\t"$c_id"\t0\t.\t"$c_s"\t"$c_e"\t"col
    }' "$DIFF" | sort -k1,1 -k2,2n > "$B9"

    [[ "$HAVE_BB" == "1" ]] || return 0
    bedToBigBed -type=bed9 "$B9" "$CHROM_SIZES" "$BB" 2>/dev/null || true
    [[ -s "$BB" ]] || { echo "    (bedToBigBed produced no track for ${PREFIX})"; return 0; }
    mkdir -p "$WEB_DIR"; cp -f "$BB" "$WEB_DIR/"; chmod 644 "$WEB_DIR/$(basename "$BB")" 2>/dev/null || true
    grep -q "name=\"${PREFIX}\"" "$TRACK_TXT" || \
      echo "track type=bigBed name=\"${PREFIX}\" description=\"${DESC}\" bigDataUrl=\"${WEB_URL}/$(basename "$BB")\" visibility=pack itemRgb=\"On\"" >> "$TRACK_TXT"
}

TRACK_TXT="$TRACKS_DIR/UCSC_Custom_Tracks.txt"
: > "$TRACK_TXT"

# =============================================================================
# PART 8: RUN EVERY COMPARISON
# =============================================================================
MAIN_MANIFEST="$DIFF_DIR/manifest_main.tsv"
SUPP_MANIFEST="$DIFF_DIR/manifest_supp.tsv"
FULL_MANIFEST="$DIFF_DIR/manifest_promoter_full.tsv"
MHDR='epitope\tcell_line\tregion_set\tcomparison\tcnt_group\tchip_group\tcnt_label\tchip_label\tn_chip\tn_cnt\tdiff_path\tsummary_path\tgc_path\tatac_path\tdnase_path\n'
printf "$MHDR" > "$MAIN_MANIFEST"; printf "$MHDR" > "$SUPP_MANIFEST"; printf "$MHDR" > "$FULL_MANIFEST"

echo "--- Differential analysis ---"
for job in "${JOBS[@]}"; do
    IFS='|' read -r EPI CELL CNT_G CHIP_G FIGSET RS <<< "$job"
    COMPARISON="${CNT_G}_VS_${CHIP_G}"
    PREFIX="${COMPARISON}_${RS}"
    DIFF="$DIFF_DIR/${PREFIX}_DiffBins.tsv"
    SUMMARY="$SUMMARY_DIR/${PREFIX}_DiffBins_PiRanked.csv"
    REGION="$REGIONS_DIR/${CELL}_${RS}.bed"

    if [[ ! -s "$DIFF" ]]; then
        echo "  -> ${PREFIX}"
        SAMP="$TMP/${PREFIX}_samples.tsv"; : > "$SAMP"
        for lab in ${SAMPLE_GROUPS[$CHIP_G]}; do printf '%s\tChIP\n' "$lab" >> "$SAMP"; done
        for lab in ${SAMPLE_GROUPS[$CNT_G]};  do printf '%s\tCnT\n'  "$lab" >> "$SAMP"; done
        cfiles=("$COUNTS_DIR/counts_${CELL}_${RS}_Paired.tab" "$COUNTS_DIR/counts_${CELL}_${RS}_Single.tab")
        present=(); for f in "${cfiles[@]}"; do [[ -s "$f" ]] && present+=("$f"); done
        Rscript "$RS_DIR/run_deseq2.R" "$REGION" "$SAMP" "$DIFF" "${present[@]}"
    else
        echo "  cached: ${PREFIX}"
    fi

    # Ranked summary, regenerated whenever it is missing or older than the diff.
    if [[ -s "$DIFF" && ( ! -s "$SUMMARY" || "$DIFF" -nt "$SUMMARY" ) ]]; then
        Rscript "$RS_DIR/rank_by_pi_score.R" "$DIFF" "$SUMMARY" \
            "$FC_THRESH" "$FDR_THRESH" "${SET_PAD[$RS]}" "$TARGET_REPS"
    fi

    generate_colored_track "$DIFF" "$PREFIX" \
        "${CELL} ${EPI} ${RS}: ${CNT_G} vs ${CHIP_G} (red=CUT&Tag, blue=ChIP-seq)"

    CNT_LBL=$(short_label "$(echo "${SAMPLE_GROUPS[$CNT_G]}"  | awk '{print $1}')")
    CHIP_LBL=$(short_label "$(echo "${SAMPLE_GROUPS[$CHIP_G]}" | awk '{print $1}')")
    N_CHIP=$(echo "${SAMPLE_GROUPS[$CHIP_G]}" | wc -w); N_CNT=$(echo "${SAMPLE_GROUPS[$CNT_G]}" | wc -w)
    ROW=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
        "$EPI" "$CELL" "$RS" "$COMPARISON" "$CNT_G" "$CHIP_G" "$CNT_LBL" "$CHIP_LBL" \
        "$N_CHIP" "$N_CNT" "$DIFF" "$SUMMARY" \
        "$FEATURES_DIR/gc_${CELL}_${RS}.tab" \
        "$FEATURES_DIR/atac_${CELL}_${RS}.tab" \
        "$FEATURES_DIR/dnase_${CELL}_${RS}.tab")

    # Whole-promoter results are kept for heatmaps and browser-region selection
    # and are deliberately excluded from the figure manifests.
    if [[ "$RS" == "promoter_full" ]]; then
        echo "$ROW" >> "$FULL_MANIFEST"
    else
        echo "$ROW" >> "$SUPP_MANIFEST"
        [[ "$FIGSET" == "main" ]] && echo "$ROW" >> "$MAIN_MANIFEST"
    fi
done

# --- Parameter file: single source of truth for the figure script ------------
PARAMS="$DIFF_DIR/analysis_params.tsv"
{
  printf 'key\tvalue\n'
  printf 'fc_thresh\t%s\n'          "$FC_THRESH"
  printf 'fdr_thresh\t%s\n'         "$FDR_THRESH"
  printf 'prom_up\t%s\n'            "$PROM_UP"
  printf 'prom_dn\t%s\n'            "$PROM_DN"
  printf 'prom_bin_size\t%s\n'      "$PROM_BIN_SIZE"
  printf 'broad_bin_size\t%s\n'     "$BROAD_BIN_SIZE"
  printf 'tss_exclude_window\t%s\n' "$TSS_EXCLUDE_WINDOW"
  printf 'min_body_tx_frac\t%s\n'   "$MIN_BODY_TX_FRAC"
  printf 'center_reads\t%s\n'       "$CENTER_READS"
  printf 'pe_counting\textendReads+samFlagInclude64\n'
  printf 'se_counting\textendReads200\n'
  printf 'se_fragment_length\t200\n'
  printf 'pi_target_reps\t%s\n'     "$TARGET_REPS"
} > "$PARAMS"

echo "========================================================="
echo " Done."
echo "   Regions     : $REGIONS_DIR"
echo "   Features    : $FEATURES_DIR"
echo "   Counts      : $COUNTS_DIR"
echo "   DESeq2      : $DIFF_DIR"
echo "   Pi-ranked   : $SUMMARY_DIR"
echo "   Tracks      : $TRACK_TXT"
echo "   Manifests   : $(basename "$MAIN_MANIFEST"), $(basename "$SUPP_MANIFEST"), $(basename "$FULL_MANIFEST")"
echo "========================================================="