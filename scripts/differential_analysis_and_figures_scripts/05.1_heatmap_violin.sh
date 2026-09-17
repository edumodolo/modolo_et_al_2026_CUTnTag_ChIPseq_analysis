#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONTROL / VALIDATION FIGURES
#
# For each cell line x histone mark, splits the differential regions into
# ChIP-seq-enriched, non-differential and CUT&Tag-enriched groups, adds a
# negative-control region set where the mark should be absent, and asks what a
# panel of independent reference datasets looks like over each group.
#
#   (A) HEATMAPS   per cell line x mark. K562 gets a MAIN figure (a short list
#                  of reference tracks) and a SUPPLEMENTARY figure (everything);
#                  MCF7 gets the supplementary only.
#   (B) VIOLINS    one per cell line x mark; every region group side by side,
#                  RAW log2(signal + 1), numeric BH-adjusted p-values on every
#                  pairwise bracket, green numbers are MEDIANS, and the legend
#                  carries the bin count of each region group.
#
# PER-COLUMN HEATMAP SCALING. Every heatmap column gets its OWN colour scale and
# its OWN profile y-axis, both derived from that column's data. Without this a
# single shared scale is set by whichever dataset has the widest dynamic range,
# and every other column is flattened into it -- the figure then says more about
# the loudest track than about the regions.
#
# No background normalisation anywhere: every violin is raw log2(signal + 1)
# and every region group appears on the axis, so a panel is read by comparing
# groups within it. Absolute heights are NOT comparable between panels.
#
# VIOLIN PANEL ORDER is the order tracks are added to the track list, which is
# what gets written as the column order of the summary table and read back as
# the panel factor levels. For the broad-mark figures that gives, at four
# panels per row:
#   K562 (16 panels, 4 x 4)  GC, DNase, ATAC, POLR2A, POLR2A pS5, POLR2A pS2,
#                            Pol II pS5 CnT, Pol II pS2 CnT, Pol II pS2S5 CnT,
#                            PRO-seq Core, PRO-seq Dastidar, EZH2,
#                            H3K36me3 ChIP, H3K36me3 CnT,
#                            H3K27me3 ChIP, H3K27me3 CnT
#   MCF7 (9 panels, 4/4/1)   the same list with every track MCF7 lacks dropped
# To change the order, reorder the add_track calls in build_tracks_broad or
# build_tracks_promoter; nothing downstream hard-codes a panel list.
#
# Region sets come from 01_differential_analysis.sh; differential results are
# located through its manifests, so the sample groups here cannot drift from
# the ones actually tested.
#
# Usage:
#   conda activate chrom_diff_figures
#   bash 05_control_validation_figures.sh
# =============================================================================

# --- STAGE TOGGLES -----------------------------------------------------------
# 0 = reuse the matrices already on disk and go straight to plotting. The
# heatmap column labels are then read FROM the matrix rather than from the
# config, so the colours and the title always describe what is actually in the
# file even if the config has since changed. Set to 1 to (re)build.
COMPUTE_MATRICES=0
# Violins only: this run remakes the violin figures and nothing else.
PLOT_HEATMAPS=0
PLOT_VIOLINS=1

CELL_LINES=(K562 MCF7)
EPITOPES=(H3K4me3 H3K27ac H3K27me3 H3K36me3)

# Which cell lines get the short MAIN heatmap in addition to the supplementary.
MAIN_FIGURE_CELLS="K562"

# --------------------------- USER CONFIG -------------------------------------
DIFF_WORK_DIR="${DIFF_WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output}"
DIFF_DIR="$DIFF_WORK_DIR/differential"
REGIONS_DIR="$DIFF_WORK_DIR/regions"
PARAMS="$DIFF_DIR/analysis_params.tsv"

WORK_DIR="${WORK_DIR:-/home/emodolo/gpfs/2026_modolo_et_al/differential_analysis/final_diff_output/heatmaps_violin_meta}"
REG_DIR="$WORK_DIR/region_groups"
MAT_DIR="$WORK_DIR/matrices"
HM_DIR="$WORK_DIR/heatmap"
VI_DIR="$WORK_DIR/violin"
RS_DIR="$WORK_DIR/Rscripts"
TMP="$WORK_DIR/_tmp"
mkdir -p "$REG_DIR" "$MAT_DIR" "$HM_DIR" "$VI_DIR" "$RS_DIR" "$TMP"

MASTER_TSV="${MASTER_TSV:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/summary/Master_Sample_Metrics.tsv}"
GENE_BODIES="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/hg38_genomic_annotations/gene_body/canonical_geneBodies.hg38.bed"
CHROM_SIZES="/home/emodolo/gpfs/Homer/data/genomes/hg38/chrom.sizes"
GC_BW="/home/emodolo/gpfs/2026_modolo_et_al/reference_datasets/GC_content/output/bigwigs/hg38_gc5Base.bw"

REF_ROOT="/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/reference_datasets"

BINSIZE=25
THREADS="${THREADS:-24}"
BROAD_WINDOW=1500       # +/- around the 3 kb bin centre for the broad marks
PROM_WINDOW=2500        # +/- around the TSS for the promoter marks
VIOLIN_FONT_SCALE="${VIOLIN_FONT_SCALE:-2.25}"
HEATMAP_COLOR_CAP_FRACTION=0.03   # clip the top (and GC's bottom) 3% for colour scaling

# Panels per row in the violin figures. 16 K562 broad panels give 4 x 4; the 9
# MCF7 panels give 4 / 4 / 1.
VIOLIN_NCOL="${VIOLIN_NCOL:-4}"

# --- Heatmap geometry ---------------------------------------------------------
# deepTools sizes a heatmap by the width of ONE column, so a figure with fewer
# columns is narrower overall even though each column is the same thickness.
#   fixed  every column is HEATMAP_WIDTH_CM wide. Main and supplementary columns
#          match, and the main figure is simply a narrower page.
#   fill   the column width is chosen so the whole figure lands near
#          HEATMAP_TARGET_TOTAL_CM, clamped to the min/max below. Main and
#          supplementary then come out a similar overall size, at the cost of
#          columns no longer being the same width between the two.
HEATMAP_WIDTH_MODE="${HEATMAP_WIDTH_MODE:-fixed}"
HEATMAP_WIDTH_CM="${HEATMAP_WIDTH_CM:-4}"          # used when mode = fixed
HEATMAP_TARGET_TOTAL_CM="${HEATMAP_TARGET_TOTAL_CM:-64}"
HEATMAP_WIDTH_MIN_CM="${HEATMAP_WIDTH_MIN_CM:-4}"
HEATMAP_WIDTH_MAX_CM="${HEATMAP_WIDTH_MAX_CM:-9}"
HEATMAP_HEIGHT_CM="${HEATMAP_HEIGHT_CM:-16}"

# --- Reference bigWigs, per cell line ----------------------------------------
# Leave a value empty when a dataset does not exist for that cell line; the
# track is then dropped from that cell line's figures rather than plotted blank.
# This is what turns the 16-panel K562 broad figure into the 9-panel MCF7 one
# without a second track list.
declare -A ATAC_BW=(
  [K562]="/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/alignment/bigwig/ATAC_NA_K562_ID173_merged_CPM.bw"
  [MCF7]="/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19/alignment/bigwig/ATAC_NA_MCF7_ID262_merged_CPM.bw"
)
declare -A DNASE_BW=(
  [K562]="${REF_ROOT}/DNase-seq/output/bigwigs/DNase_K562_ENCFF972GVB.bigWig"
  [MCF7]="${REF_ROOT}/DNase-seq/output/bigwigs/DNase_MCF7_ENCFF038NFC.bigWig"
)
# EZH2 is K562 only: the MCF7 dataset is too noisy to interpret and was dropped.
declare -A EZH2_BW=(
  [K562]="${REF_ROOT}/EZH2_ChIP-seq/output/bigwigs/chip_EZH2_K562_ENCFF163LOW.bigWig"
  [MCF7]=""
)
declare -A EP300_BW=(
  [K562]="${REF_ROOT}/EP300_ChIP-seq/output/bigwigs/chip_EP300_K562_ENCFF325DSL.bigWig"
  [MCF7]="${REF_ROOT}/EP300_ChIP-seq/output/bigwigs/chip_EP300_MCF7_ENCFF708NMR.bigWig"
)
declare -A MNASE_BW=(
  [K562]="${REF_ROOT}/MNase-seq/output/MNase-seq_K562_ENCFF000VNN_hg19.hg38lift.hg38.bw"
  [MCF7]=""
)
# RNA Pol II ChIP-seq: K562 has total plus the two phospho-forms; MCF7 total only.
declare -A POLR2A_BW=(
  [K562]="${REF_ROOT}/RNA_Poll_ChIP-seq/output/bigwigs/chip_POLR2A_K562_ENCFF914WIS.bigWig"
  [MCF7]="${REF_ROOT}/RNA_Poll_ChIP-seq/output/bigwigs/chip_POLR2A_MCF7_ENCFF827YIP.bigWig"
)
declare -A POLR2A_S5_BW=(
  [K562]="${REF_ROOT}/RNA_Poll_ChIP-seq/output/bigwigs/chip_POLR2AphosphoS5_K562_ENCFF677XKP.bigWig"
  [MCF7]=""
)
declare -A POLR2A_S2_BW=(
  [K562]="${REF_ROOT}/RNA_Poll_ChIP-seq/output/bigwigs/chip_POLR2AphosphoS2_K562_ENCFF957YRD.bigWig"
  [MCF7]=""
)

# RNA Pol II CUT&Tag, given as plot_labels and resolved through the master
# metrics table like the epitope tracks. K562 only.
declare -A POLII_CNT_S5=( [K562]="CnT_PolSer5P_K562_KayaOkur_2019_rep1"  [MCF7]="" )
declare -A POLII_CNT_S2=( [K562]="CnT_PolSer2P_K562_KayaOkur_2019_rep1"  [MCF7]="" )
declare -A POLII_CNT_S25=( [K562]="CnT_PolSer25P_K562_KayaOkur_2019_rep1" [MCF7]="" )

# --- PRO-seq: stranded, merged on the fly ------------------------------------
PRO_DIR="${REF_ROOT}/PRO-seq/output/hg38_bigwigs"
declare -A PRO_DAST_PLUS=(
  [K562]="${PRO_DIR}/PRO-seq_Dastidar_2023_GSM6383641_K562_NHS_for_hg38.bw"
  [MCF7]="${PRO_DIR}/PRO-seq_Dastidar_2023_GSM6383644_MCF7_NHS_rep1_for_hg38.bw"
)
declare -A PRO_DAST_MINUS=(
  [K562]="${PRO_DIR}/PRO-seq_Dastidar_2023_GSM6383641_K562_NHS_rev_hg38.bw"
  [MCF7]="${PRO_DIR}/PRO-seq_Dastidar_2023_GSM6383644_MCF7_NHS_rep1_rev_hg38.bw"
)
declare -A PRO_CORE_PLUS=(
  [K562]="${PRO_DIR}/PRO-seq_Core_2014_GSM1480327_K562_plus_hg38.bw"
  [MCF7]=""
)
declare -A PRO_CORE_MINUS=(
  [K562]="${PRO_DIR}/PRO-seq_Core_2014_GSM1480327_K562_minus_hg38.bw"
  [MCF7]=""
)

# --- Which sample supplies each epitope's ChIP-seq and CUT&Tag track ----------
declare -A CHIP_TRACK=(
  [K562,H3K4me3]="chip_H3K4me3_K562_Bernstein_2017_rep1"
  [K562,H3K27ac]="chip_H3K27ac_K562_Bernstein_2022_rep3"
  [K562,H3K27me3]="chip_H3K27me3_K562_Bernstein_2022_rep3"
  [K562,H3K36me3]="chip_H3K36me3_K562_Bernstein_2022_rep3"
  [MCF7,H3K4me3]="chip_H3K4me3_MCF7_Bernstein_2017_rep1"
  [MCF7,H3K27ac]="chip_H3K27ac_MCF7_Bernstein_2017_rep1"
  [MCF7,H3K27me3]="chip_H3K27me3_MCF7_Bernstein_2017_rep2"
  [MCF7,H3K36me3]="chip_H3K36me3_MCF7_Bernstein_2021_rep1"
)
declare -A CNT_TRACK=(
  [K562,H3K4me3]="CnT_H3K4me3_K562_KayaOkur_2020_rep2"
  [K562,H3K27ac]="CnT_H3K27ac_K562_Abbasova_2025_ab177r1"
  [K562,H3K27me3]="CnT_H3K27me3_K562_Abbasova_2025_rep2"
  [K562,H3K36me3]="CnT_H3K36me3_K562_Wu_2025_rep4"
  [MCF7,H3K4me3]="CnT_H3K4me3_MCF7_Tian_2023_rep2"
  [MCF7,H3K27ac]="CnT_H3K27ac_MCF7_Tian_2023_rep2"
  [MCF7,H3K27me3]="CnT_H3K27me3_MCF7_Tian_2023_rep3"
  [MCF7,H3K36me3]="CnT_H3K36me3_MCF7_Tian_2023_rep1"
)

# --- Which differential comparison defines the region groups ------------------
declare -A CNT_GROUP=(
  [K562,H3K4me3]="CnT_H3K4me3_K562_KayaOkur_2020_rep1_2020_rep2"
  [K562,H3K27ac]="CnT_H3K27ac_K562_Abbasova_2025_ab177r1_2025_ab177r2"
  [K562,H3K27me3]="CnT_H3K27me3_K562_Abbasova_2025_rep1_2025_rep2"
  [K562,H3K36me3]="CnT_H3K36me3_K562_Wu_2025_rep1_2025_rep2_2025_rep3_2025_rep4_2025_rep5"
  [MCF7,H3K4me3]="CnT_H3K4me3_MCF7_Tian_2023_rep1_2023_rep2"
  [MCF7,H3K27ac]="CnT_H3K27ac_MCF7_Tian_2023_rep1_2023_rep2"
  [MCF7,H3K27me3]="CnT_H3K27me3_MCF7_Tian_2023_rep1_2023_rep2_2023_rep3_2023_rep4"
  [MCF7,H3K36me3]="CnT_H3K36me3_MCF7_Tian_2023_rep1_2023_rep2"
)

# The ChIP-seq side of that comparison. This is NOT optional: every CUT&Tag
# group is tested against more than one independent ChIP-seq source, so the
# CUT&Tag group alone does not identify a single differential result. Naming
# both sides is what makes the region groups reproducible.
declare -A CHIP_GROUP=(
  [K562,H3K4me3]="chip_H3K4me3_K562_Bernstein_2017_rep1_2017_rep2"
  [K562,H3K27ac]="chip_H3K27ac_K562_Bernstein_2011_rep2_2022_rep3"
  [K562,H3K27me3]="chip_H3K27me3_K562_Bernstein_2011_rep1_2022_rep3"
  [K562,H3K36me3]="chip_H3K36me3_K562_Bernstein_2011_rep2_2022_rep3"
  [MCF7,H3K4me3]="chip_H3K4me3_MCF7_Bernstein_2017_rep1_2017_rep2"
  [MCF7,H3K27ac]="chip_H3K27ac_MCF7_Bernstein_2017_rep1_2017_rep2"
  [MCF7,H3K27me3]="chip_H3K27me3_MCF7_Bernstein_2017_rep1_2017_rep2"
  [MCF7,H3K36me3]="chip_H3K36me3_MCF7_Bernstein_2021_rep1_2021_rep2"
)

# --- Region sets, matching 03_peak_overlap_validation.sh ----------------------
declare -A EPI_PRIMARY_SET=(
  [H3K4me3]=promoter_full  [H3K27ac]=promoter_full
  [H3K36me3]=genebody_bins [H3K27me3]=polycomb_bins
)
declare -A EPI_CONTROL_SET=(
  [H3K4me3]=promoter_repressed [H3K27ac]=promoter_repressed
  [H3K36me3]=polycomb_bins     [H3K27me3]=genebody_bins
)
declare -A SET_DESC=(
  [promoter_full]="active promoters"
  [promoter_repressed]="Polycomb-repressed promoters"
  [genebody_bins]="active gene-body bins"
  [polycomb_bins]="Polycomb-repressed bins"
)
PROMOTER_MARKS="H3K4me3 H3K27ac"
# The broad marks are shown alongside each other, as in the original figures.
declare -A OTHER_BROAD=( [H3K36me3]=H3K27me3 [H3K27me3]=H3K36me3 )

# --- Heatmap colour scales, one per data type --------------------------------
# Carried over from the original script so a data type keeps one identity
# across every figure. Every Pol II form, ChIP-seq or CUT&Tag, is the same
# purple; DNase-seq is the new teal.
declare -A COLOR_OF=(
  ["GC content"]="#000004,#51127c,#b73779,#fc8961,#fcfdbf"
  ["DNase-seq"]="white,#02C49B"
  ["ATAC-seq"]="white,#C2185B"
  ["MNase-seq"]="white,#8c510a"
  ["POLR2A ChIP-seq"]="white,purple"
  ["POLR2A pS5 ChIP-seq"]="white,purple"
  ["POLR2A pS2 ChIP-seq"]="white,purple"
  ["Pol II pS5 CUT&Tag"]="white,purple"
  ["Pol II pS2 CUT&Tag"]="white,purple"
  ["Pol II pS2S5 CUT&Tag"]="white,purple"
  ["PRO-seq (Core 2014)"]="white,orange"
  ["PRO-seq (Dastidar 2023)"]="white,orange"
  ["EP300 ChIP-seq"]="white,blue"
  ["EZH2 ChIP-seq"]="white,red"
)

# Histone-mark columns are named after their epitope and are always green, so
# they are matched by pattern rather than listed. Anything else that falls
# through is a LABEL THAT DOES NOT MATCH ITS COLOUR KEY -- the failure that
# silently painted EP300 and EZH2 green -- so it is reported instead of
# quietly defaulting.
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
echo " Control / validation figures"
echo " Differential input : $DIFF_DIR"
echo " Output             : $WORK_DIR"
echo "========================================================="

miss=0
for t in bedtools computeMatrix plotHeatmap multiBigwigSummary bigwigCompare Rscript awk sort; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: missing tool: $t" >&2; miss=1; }
done
command -v computeMatrixOperations >/dev/null 2>&1 || \
    echo "  NOTE: computeMatrixOperations not found; MAIN figures will be skipped" >&2
for f in "$MASTER_TSV" "$GENE_BODIES" "$CHROM_SIZES" "$GC_BW" "$PARAMS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing file: $f" >&2; miss=1; }
done
[[ $miss -eq 0 ]] || { echo "Preflight failed." >&2; exit 1; }

FC_THRESH=$(awk -F'\t' '$1=="fc_thresh"{print $2}'  "$PARAMS")
FDR_THRESH=$(awk -F'\t' '$1=="fdr_thresh"{print $2}' "$PARAMS")
echo "  thresholds from analysis_params.tsv: FDR < ${FDR_THRESH}, |log2FC| > ${FC_THRESH}"

col_index() {
    awk -F'\t' -v want="$1" 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$i); if($i==want){print i; exit}}}' "$MASTER_TSV"
}
COL_LABEL=$(col_index "plot_label"); COL_BW=$(col_index "location_final_bigwig")
[[ -n "$COL_LABEL" && -n "$COL_BW" ]] || { echo "ERROR: plot_label/location_final_bigwig not in $MASTER_TSV" >&2; exit 1; }
bigwig_of() {
    [[ -n "${1:-}" ]] || { echo ""; return 0; }
    awk -F'\t' -v s="$1" -v n="$COL_LABEL" -v c="$COL_BW" \
        'NR>1{gsub(/\r/,"",$n); if($n==s){gsub(/\r/,"",$c); print $c; exit}}' "$MASTER_TSV"
}

# Returns every manifest row matching the request, as "chip_group<TAB>diff_path".
# Each awk reads a FILE and the result is captured per iteration; piping the
# loop into `head -1` would SIGPIPE the later awks under pipefail.
diff_rows_for() {   # cell, epitope, region_set, cnt_group, chip_group ("" = any)
    local cell="$1" epi="$2" rs="$3" cg="$4" chg="${5:-}" m
    for m in "$DIFF_DIR/manifest_supp.tsv" "$DIFF_DIR/manifest_promoter_full.tsv"; do
        [[ -s "$m" ]] || continue
        awk -F'\t' -v c="$cell" -v e="$epi" -v r="$rs" -v g="$cg" -v h="$chg" \
            'NR>1 && $2==c && $1==e && $3==r && $5==g && (h == "" || $6==h) {print $6 "\t" $11}' "$m"
    done | sort -u
}

# Resolves one cell x epitope to exactly one differential result, setting
# DIFF_PATH and DIFF_CHIP. Ambiguity is a hard error rather than a silent pick:
# with several ChIP-seq sources per CUT&Tag group, taking "the first match"
# would make the region groups depend on the order comparisons happen to be
# listed in 01_differential_analysis.sh.
DIFF_PATH=""; DIFF_CHIP=""
resolve_diff() {   # cell, epitope, region_set
    local cell="$1" epi="$2" rs="$3"
    local cg="${CNT_GROUP[$cell,$epi]:-}" chg="${CHIP_GROUP[$cell,$epi]:-}"
    DIFF_PATH=""; DIFF_CHIP=""
    [[ -n "$cg" ]] || { echo "  ERROR: no CNT_GROUP entry for ${cell},${epi}" >&2; return 1; }
    local rows n
    mapfile -t rows < <(diff_rows_for "$cell" "$epi" "$rs" "$cg" "$chg")
    n=${#rows[@]}
    if [[ $n -eq 0 ]]; then
        echo "  ERROR: no differential result for ${cell} ${epi} ${rs}" >&2
        echo "         CUT&Tag group : ${cg}" >&2
        echo "         ChIP-seq group: ${chg:-<unset>}" >&2
        mapfile -t rows < <(diff_rows_for "$cell" "$epi" "$rs" "$cg" "")
        if [[ ${#rows[@]} -gt 0 ]]; then
            echo "         that CUT&Tag group IS present, but against:" >&2
            printf '           %s\n' "${rows[@]%%$'\t'*}" >&2
        fi
        return 1
    fi
    if [[ $n -gt 1 ]]; then
        echo "  ERROR: ${cell} ${epi} ${rs} matches ${n} differential results:" >&2
        printf '           %s\n' "${rows[@]%%$'\t'*}" >&2
        echo "         set CHIP_GROUP[${cell},${epi}] to the one this figure should use." >&2
        return 1
    fi
    DIFF_CHIP="${rows[0]%%$'\t'*}"; DIFF_PATH="${rows[0]#*$'\t'}"
    [[ -s "$DIFF_PATH" ]] || { echo "  ERROR: differential file missing: $DIFF_PATH" >&2; return 1; }
    return 0
}

is_promoter_mark() { [[ " $PROMOTER_MARKS " == *" $1 "* ]]; }
is_main_cell()     { [[ " $MAIN_FIGURE_CELLS " == *" $1 "* ]]; }

# --- PRO-seq strand merge -----------------------------------------------------
# The minus strand is stored as negative values in some datasets and positive in
# others, so the operation is chosen from the file's own minimum rather than
# assumed; adding a negative track would cancel real signal. Cached on disk.
merge_proseq() {   # plus, minus, out -> echoes the merged path, or empty
    local P="$1" M="$2" OUT="$3"
    [[ -n "$P" && -n "$M" ]] || { echo ""; return 0; }
    if [[ ! -s "$P" || ! -s "$M" ]]; then
        echo "  WARNING: PRO-seq strand missing, dataset skipped: $P / $M" >&2; echo ""; return 0
    fi
    if [[ ! -s "$OUT" ]]; then
        local IS_NEG OP
        IS_NEG=$(python -c "import pyBigWig; bw=pyBigWig.open('$M'); print(1 if bw.header()['minVal'] < 0 else 0)")
        OP="add"; [[ "$IS_NEG" == "1" ]] && OP="subtract"
        echo "  -> merging PRO-seq strands (${OP}) -> $(basename "$OUT")" >&2
        bigwigCompare -b1 "$P" -b2 "$M" --operation "$OP" --binSize "$BINSIZE" \
            -p "$THREADS" -o "$OUT" >/dev/null 2>&1
    else
        echo "  PRO-seq merged bigWig cached: $(basename "$OUT")" >&2
    fi
    [[ -s "$OUT" ]] && echo "$OUT" || echo ""
}

# =============================================================================
# REGION GROUPS
# =============================================================================
mainchr() { awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' "$@"; }

# Recovers strand and the TSS point for promoter windows. The differential
# regions carry no strand, so they are matched back to the stranded all-gene
# promoter windows of the same geometry by exact reciprocal overlap.
ALL_PROM_STRANDED="$TMP/all_genes_promoter_${PROM_WINDOW}.bed"
build_all_prom_stranded() {
    [[ -s "$ALL_PROM_STRANDED" ]] && return 0
    echo "  -> building stranded all-gene promoter windows" >&2
    local w=$((PROM_WINDOW * 2))
    cut -f1-6 "$GENE_BODIES" | tr -d '\r' | mainchr | sort -k1,1 -k2,2n \
      | awk -F'\t' -v up="$PROM_WINDOW" -v dn="$PROM_WINDOW" '
          NR==FNR{size[$1]=$2; next}
          { chr=$1;s=$2;e=$3;name=$4;str=$6;
            if(str=="-"){tss=e;ws=tss-dn;we=tss+up}else{tss=s;ws=tss-up;we=tss+dn}
            if(ws<0)ws=0; if((chr in size)&&we>size[chr])we=size[chr];
            if(we>ws)print chr"\t"ws"\t"we"\t"name"\t.\t"str }' "$CHROM_SIZES" - \
      | awk -F'\t' -v w="$w" 'BEGIN{OFS="\t"}($3-$2)==w' \
      | sort -k1,1 -k2,2n | awk -F'\t' '!seen[$1"\t"$2"\t"$3]++' > "$ALL_PROM_STRANDED"
}

window_to_tss() {   # window bed, out
    bedtools intersect -a "$1" -b "$ALL_PROM_STRANDED" -wa -wb -f 1.0 -r \
      | awk -F'\t' -v up="$PROM_WINDOW" 'BEGIN{OFS="\t"}
          { chr=$1;s=$2;e=$3;id=$4;str=$(NF);
            if(str=="-") tss=e-up; else tss=s+up;
            if(tss<0)tss=0; k=chr"\t"tss; if(!seen[k]++) print chr,tss,tss+1,id,0,str }' \
      | sort -k1,1 -k2,2n > "$2"
}

# Splits one epitope's primary region set into the three differential classes,
# and writes its negative-control set alongside. Returns non-zero when the
# differential result is missing so the caller can skip that epitope.
declare -A REGION_COMPARISON=()
build_region_groups() {   # cell, epitope
    local CL="$1" EPI="$2"
    local PRIMARY="${EPI_PRIMARY_SET[$EPI]}" CONTROL="${EPI_CONTROL_SET[$EPI]}"
    local DIFF PRE g
    resolve_diff "$CL" "$EPI" "$PRIMARY" || return 1
    DIFF="$DIFF_PATH"
    REGION_COMPARISON["$CL,$EPI"]="${CNT_GROUP[$CL,$EPI]} vs ${DIFF_CHIP}"
    echo "  regions from: ${CNT_GROUP[$CL,$EPI]}  vs  ${DIFF_CHIP}" >&2
    PRE="$REG_DIR/${CL}_${EPI}"
    awk -F'\t' -v fc="$FC_THRESH" -v fdr="$FDR_THRESH" -v pre="$PRE" '
      NR==1 { for(i=1;i<=NF;i++){ if($i=="Chr")c=i; if($i=="Start")s=i; if($i=="End")e=i;
                                  if($i=="BinID")b=i; if($i=="log2FC")l=i; if($i=="padj")p=i }
              next }
      { st="nondiff"
        if ($p+0 < fdr && $l+0 >  fc) st="cnt"
        else if ($p+0 < fdr && $l+0 < -fc) st="chip"
        print $c "\t" $s "\t" $e "\t" $b "\t0\t." > (pre "_" st ".bed") }' "$DIFF"
    cut -f1-4 "$REGIONS_DIR/${CL}_${CONTROL}.bed" \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,0,"."}' > "${PRE}_control.bed"
    for g in chip nondiff cnt control; do
        [[ -s "${PRE}_${g}.bed" ]] || : > "${PRE}_${g}.bed"
        sort -k1,1 -k2,2n -o "${PRE}_${g}.bed" "${PRE}_${g}.bed"
    done
    # Promoter marks are anchored on the TSS, so the windows are reduced to
    # strand-aware single-base points for computeMatrix.
    if is_promoter_mark "$EPI"; then
        for g in chip nondiff cnt control; do
            window_to_tss "${PRE}_${g}.bed" "${PRE}_${g}_TSS.bed"
        done
    fi
    return 0
}

# =============================================================================
# TRACK LISTS
# =============================================================================
# BW/LBL are built in figure order; TRACK_BW keeps label -> bigWig so the figure
# title can name the exact file behind every column. LBL order is also the
# violin panel order, since it becomes the column order of the summary table.
# NOTE ON ARRAY EXPANSIONS. Every "${arr[@]}" below is written as
# ${arr[@]+"${arr[@]}"}. On bash older than 4.4 -- which is what most HPC
# images ship -- expanding an EMPTY array under `set -u` raises "unbound
# variable" and kills the run. That is what happened when MCF7 reached the
# heatmap step with an empty MAIN-figure column list, since MCF7 has no main
# figure. The guarded form expands to nothing when the array is empty and is
# identical otherwise.
declare -A TRACK_BW=()
BW=(); LBL=()
reset_tracks() { BW=(); LBL=(); }
# Adds a track given a bigWig PATH. The three ways a track can fail are
# reported separately: a track that quietly vanishes from a figure because a
# label was mistyped looks identical to one that was never configured.
add_track() {   # bigwig path (may be empty), label
    if [[ -z "${1:-}" ]]; then
        echo "     skipping ${2}: not configured for this cell line" >&2
    elif [[ ! -s "${1}" ]]; then
        echo "     skipping ${2}: bigWig missing on disk -- ${1}" >&2
    else
        BW+=("$1"); LBL+=("$2"); TRACK_BW["$2"]="$1"
    fi
}

# Adds a track given a plot_label to look up in the master metrics table. A
# label that is absent from the table is named explicitly, because that is the
# failure that silently drops a dataset from every figure.
add_track_by_label() {   # plot_label (may be empty), label
    local lab="${1:-}" disp="$2" bw
    if [[ -z "$lab" ]]; then
        echo "     skipping ${disp}: not configured for this cell line" >&2; return 0
    fi
    bw=$(bigwig_of "$lab")
    if [[ -z "$bw" ]]; then
        echo "     skipping ${disp}: plot_label '${lab}' not found in $(basename "$MASTER_TSV")" >&2
        return 0
    fi
    add_track "$bw" "$disp"
}

# Reference blocks shared by both figure types, in the order they appear.
add_accessibility() {   # cell [, "with_mnase"]
    add_track "$GC_BW" "GC content"
    # MNase-seq sits directly after GC content, and only in the promoter
    # figures: nucleosome positioning is what the TSS-anchored view is about.
    [[ "${2:-}" == "with_mnase" ]] && add_track "${MNASE_BW[$1]:-}" "MNase-seq"
    add_track "${DNASE_BW[$1]:-}" "DNase-seq"
    add_track "${ATAC_BW[$1]:-}"  "ATAC-seq"
}
add_polii() {   # cell -- ChIP-seq forms first, then CUT&Tag forms
    add_track "${POLR2A_BW[$1]:-}"                "POLR2A ChIP-seq"
    add_track "${POLR2A_S5_BW[$1]:-}"             "POLR2A pS5 ChIP-seq"
    add_track "${POLR2A_S2_BW[$1]:-}"             "POLR2A pS2 ChIP-seq"
    add_track_by_label "${POLII_CNT_S5[$1]:-}"  "Pol II pS5 CUT&Tag"
    add_track_by_label "${POLII_CNT_S2[$1]:-}"  "Pol II pS2 CUT&Tag"
    add_track_by_label "${POLII_CNT_S25[$1]:-}" "Pol II pS2S5 CUT&Tag"
}
add_proseq() {   # cell
    add_track "$PRO_CORE" "PRO-seq (Core 2014)"
    add_track "$PRO_DAST" "PRO-seq (Dastidar 2023)"
}
add_mark_pair() {   # cell, epitope -- the two assays under comparison
    add_track_by_label "${CHIP_TRACK[$1,$2]:-}" "${2} ChIP-seq"
    add_track_by_label "${CNT_TRACK[$1,$2]:-}"  "${2} CUT&Tag"
}

# Promoter figure: accessibility, Pol II, PRO-seq, EP300, EZH2, the repressive
# counter-example, and the mark under test last.
build_tracks_promoter() {   # cell, epitope
    reset_tracks
    add_accessibility "$1" with_mnase
    add_polii "$1"
    add_proseq "$1"
    add_track "${EP300_BW[$1]:-}" "EP300 ChIP-seq"
    add_track "${EZH2_BW[$1]:-}"  "EZH2 ChIP-seq"
    add_mark_pair "$1" H3K27me3
    add_mark_pair "$1" "$2"
}

# Broad figure: no EP300, no MNase, and both broad marks at the end. This call
# order IS the violin panel order: GC, DNase, ATAC, the six Pol II forms, the
# two PRO-seq datasets, EZH2, then H3K36me3 and H3K27me3 as ChIP/CUT&Tag pairs.
build_tracks_broad() {   # cell
    reset_tracks
    add_accessibility "$1"
    add_polii "$1"
    add_proseq "$1"
    add_track "${EZH2_BW[$1]:-}" "EZH2 ChIP-seq"
    add_mark_pair "$1" H3K36me3
    add_mark_pair "$1" H3K27me3
}

# MAIN figure column lists. Any label a cell line does not have is dropped.
main_labels_promoter() {   # epitope
    printf '%s\n' "GC content" "MNase-seq" "DNase-seq" "ATAC-seq" "PRO-seq (Core 2014)" \
                  "${1} ChIP-seq" "${1} CUT&Tag"
}
# The two Pol II pS2 columns sit between ATAC-seq and PRO-seq deliberately.
# That is where they already sit in the SUPPLEMENTARY matrix (add_polii runs
# before add_proseq), so the main figure keeps the same left-to-right order as
# the supplementary regardless of whether computeMatrixOperations honours the
# order of --samples or the order stored in the matrix. Serine-2 phosphorylated
# Pol II marks elongating polymerase, which is what H3K36me3 tracks, so it also
# belongs next to the accessibility and nascent-transcription block rather than
# beside the histone marks at the end.
main_labels_broad() {
    printf '%s\n' "GC content" "DNase-seq" "ATAC-seq" \
                  "POLR2A pS2 ChIP-seq" "Pol II pS2 CUT&Tag" \
                  "PRO-seq (Core 2014)" "EZH2 ChIP-seq" \
                  "H3K36me3 ChIP-seq" "H3K36me3 CUT&Tag" "H3K27me3 ChIP-seq" "H3K27me3 CUT&Tag"
}

# =============================================================================
# HEATMAP RENDERING
# =============================================================================
# Per-column limits, exactly as the original script derived them:
#   zMax  the (1 - cap) quantile of THAT column's values
#   zMin  0, except the first column (GC content), which never approaches zero
#         and so takes its own lower quantile instead
#   yMin/yMax  that column's group-mean profile, padded by 10% of its range
compute_matrix_limits() {   # matrix.gz -> ZMIN/ZMAX/YMIN/YMAX lines on stdout
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

# The title carries the source bigWig of every column, in column order, so a
# figure states exactly which file produced each heatmap rather than leaving it
# to be reconstructed from a config file later.
build_title() {   # header line; reads HM_LABELS and TRACK_BW
    local hdr="$1" i t
    [[ ${#HM_LABELS[@]} -gt 0 ]] || { printf '%s' "$hdr"; return 0; }
    t="${hdr}"$'\n'"columns, in order (dataset = bigWig):"
    for i in "${!HM_LABELS[@]}"; do
        t+=$'\n'"$((i + 1)). ${HM_LABELS[$i]}  =  $(basename "${TRACK_BW[${HM_LABELS[$i]}]:-n/a}")"
    done
    printf '%s' "$t"
}

plot_heatmap() {   # matrix.gz, out_pdf, header text, reflabel
    local MAT="$1" OUT="$2" HDR="$3" REFLAB="$4"
    [[ -s "$MAT" ]] || { echo "     no matrix for $(basename "$OUT")" >&2; return 0; }
    [[ ${#HM_LABELS[@]} -gt 0 ]] || { echo "     no columns for $(basename "$OUT")" >&2; return 0; }
    local CLIST=() l
    for l in ${HM_LABELS[@]+"${HM_LABELS[@]}"}; do CLIST+=("$(color_for "$l")"); done

    local LIMITS ZMIN ZMAX YMIN YMAX
    LIMITS=$(compute_matrix_limits "$MAT")
    ZMIN=$(echo "$LIMITS" | awk '/^ZMIN:/ { $1=""; print $0 }')
    ZMAX=$(echo "$LIMITS" | awk '/^ZMAX:/ { $1=""; print $0 }')
    YMIN=$(echo "$LIMITS" | awk '/^YMIN:/ { $1=""; print $0 }')
    YMAX=$(echo "$LIMITS" | awk '/^YMAX:/ { $1=""; print $0 }')

    local WCM="$HEATMAP_WIDTH_CM"
    if [[ "$HEATMAP_WIDTH_MODE" == "fill" ]]; then
        WCM=$(awk -v t="$HEATMAP_TARGET_TOTAL_CM" -v n="${#HM_LABELS[@]}" \
                  -v lo="$HEATMAP_WIDTH_MIN_CM" -v hi="$HEATMAP_WIDTH_MAX_CM" \
                  'BEGIN{ w = t / n; if (w < lo) w = lo; if (w > hi) w = hi; printf "%.2f", w }')
    fi
    echo "  -> plotHeatmap $(basename "$OUT")  (${#HM_LABELS[@]} columns @ ${WCM} cm, per-column scales)" >&2
    # shellcheck disable=SC2086
    plotHeatmap -m "$MAT" -o "$OUT" --plotType lines \
        --plotTitle "$(build_title "$HDR")" \
        --sortRegions descend --sortUsing mean --sortUsingSamples 1 \
        --colorList ${CLIST[@]+"${CLIST[@]}"} --zMin $ZMIN --zMax $ZMAX --yMin $YMIN --yMax $YMAX \
        --regionsLabel ${HM_GROUP_DISP[@]+"${HM_GROUP_DISP[@]}"} --refPointLabel "$REFLAB" \
        --legendLocation center-right --heatmapHeight "$HEATMAP_HEIGHT_CM" --heatmapWidth "$WCM" \
        >/dev/null 2>&1 || echo "     plotHeatmap failed for $(basename "$OUT")" >&2
}

# Subsets an existing matrix down to a list of columns. computeMatrixOperations
# does this without touching the bigWigs, so the MAIN figure costs nothing
# beyond the full matrix that was already built.
#
# The cache is keyed on the COLUMN LIST, not just on file timestamps. Checking
# only that the subset is newer than the full matrix was wrong: editing
# main_labels_broad() changes neither file, so the stale subset would be reused
# and the requested columns would never appear.
subset_matrix() {   # in_matrix, out_matrix, labels...
    local IN="$1" OUT="$2"; shift 2
    command -v computeMatrixOperations >/dev/null 2>&1 || return 1
    local SIG="${OUT}.signature" SIG_NOW
    SIG_NOW=$(printf '%s\n' "$@")
    if [[ -s "$OUT" && "$OUT" -nt "$IN" && -s "$SIG" && "$(cat "$SIG")" == "$SIG_NOW" ]]; then
        return 0
    fi
    [[ -s "$OUT" ]] && echo "     main-figure column list changed -- rebuilding subset matrix" >&2
    computeMatrixOperations subset -m "$IN" --samples "$@" -o "$OUT" >/dev/null 2>&1 || return 1
    [[ -s "$OUT" ]] || return 1
    printf '%s\n' "$SIG_NOW" > "$SIG"
}

# Builds (or reuses) a matrix. The signature records the exact columns and
# region files it was built from: renaming a column, adding a track or
# reordering the region groups all change what is stored inside the matrix, so
# the signature is what decides whether the cached file can be reused.
ensure_matrix() {   # matrix path, ref_point, before/after, beds...
    local MAT="$1" REFP="$2" FLANK="$3"; shift 3
    [[ ${#BW[@]} -gt 0 ]] || { echo "  ERROR: no tracks to build $(basename "$MAT")" >&2; return 1; }
    local BEDS=("$@") SIG="${MAT}.signature" SIG_NOW
    SIG_NOW=$(printf '%s\n' ${LBL[@]+"${LBL[@]}"} "--" ${BEDS[@]+"${BEDS[@]}"})
    if [[ -s "$MAT" && -s "$SIG" && "$(cat "$SIG")" == "$SIG_NOW" ]]; then
        echo "  matrix cached: $(basename "$MAT")" >&2; return 0
    fi
    if [[ "$COMPUTE_MATRICES" != "1" ]]; then
        if [[ -s "$MAT" ]]; then
            echo "  reusing existing matrix (COMPUTE_MATRICES=0): $(basename "$MAT")" >&2
            return 0
        fi
        echo "  WARNING: $(basename "$MAT") does not exist and COMPUTE_MATRICES=0 -- skipped" >&2
        return 1
    fi
    [[ -s "$MAT" ]] && echo "  matrix contents changed (columns or region order) -- recomputing" >&2
    echo "  -> computeMatrix (${REFP}, +/-${FLANK}, ${#BW[@]} tracks, ${#BEDS[@]} region groups)" >&2
    computeMatrix reference-point --referencePoint "$REFP" \
        -S ${BW[@]+"${BW[@]}"} --samplesLabel ${LBL[@]+"${LBL[@]}"} -R ${BEDS[@]+"${BEDS[@]}"} \
        -b "$FLANK" -a "$FLANK" --binSize "$BINSIZE" --missingDataAsZero \
        -p "$THREADS" -o "$MAT" >/dev/null 2>&1 && printf '%s\n' "$SIG_NOW" > "$SIG"
}

# Reads the column labels out of a matrix header. This is the source of truth
# for anything drawn from that matrix: with COMPUTE_MATRICES=0 the file on disk
# may predate a config change, and taking labels from the config instead would
# mis-assign every colour and every name in the title.
matrix_labels() {   # matrix.gz
    python - "$1" <<'PYEOF'
import gzip, json, sys
with gzip.open(sys.argv[1], "rt") as f:
    print("\n".join(json.loads(f.readline()[1:])["sample_labels"]))
PYEOF
}

# Draws the supplementary figure and, for the main cell lines, a subset figure.
render_heatmap_pair() {   # matrix, out_prefix, header, reflabel, main_labels...
    local MAT="$1" PREFIX="$2" HDR="$3" REFLAB="$4"; shift 4
    local WANT=("$@") KEEP=() MISSING=() l have found
    mapfile -t HM_LABELS < <(matrix_labels "$MAT")
    if [[ "$(printf '%s\n' ${HM_LABELS[@]+"${HM_LABELS[@]}"})" != "$(printf '%s\n' ${LBL[@]+"${LBL[@]}"})" ]]; then
        echo "     NOTE: matrix columns differ from the configured track list;" >&2
        echo "           plotting what is in the matrix. Set COMPUTE_MATRICES=1 to rebuild." >&2
    fi
    plot_heatmap "$MAT" "${PREFIX}_supplementary.pdf" "$HDR" "$REFLAB"

    [[ ${#WANT[@]} -gt 0 ]] || return 0
    for l in ${WANT[@]+"${WANT[@]}"}; do
        found=0
        for have in ${HM_LABELS[@]+"${HM_LABELS[@]}"}; do [[ "$l" == "$have" ]] && { found=1; break; }; done
        if [[ $found -eq 1 ]]; then KEEP+=("$l"); else MISSING+=("$l"); fi
    done
    # A requested column that is not in the matrix is reported by name. Silently
    # dropping it is how a main figure ends up missing a track that was added to
    # the label list but never made it into the matrix.
    if [[ ${#MISSING[@]} -gt 0 ]]; then
        echo "     NOTE: requested MAIN columns absent from $(basename "$MAT"):" >&2
        printf '           %s\n' ${MISSING[@]+"${MISSING[@]}"} >&2
        echo "           set COMPUTE_MATRICES=1 to rebuild the matrix with them." >&2
    fi
    [[ ${#KEEP[@]} -gt 0 ]] || return 0
    local MAT_MAIN="${MAT%.gz}_main.gz"
    if subset_matrix "$MAT" "$MAT_MAIN" ${KEEP[@]+"${KEEP[@]}"}; then
        mapfile -t HM_LABELS < <(matrix_labels "$MAT_MAIN")
        plot_heatmap "$MAT_MAIN" "${PREFIX}_main.pdf" "$HDR" "$REFLAB"
    else
        echo "     could not build the MAIN subset matrix -- supplementary only" >&2
    fi
}

# =============================================================================
# R: VIOLINS
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
# Levels, colours and block membership are supplied by the caller, pipe
# separated and aligned, so one script serves both the four-group promoter
# figure and the six-group combined broad figure.
GROUP_LEVELS <- strsplit(a[8], "|", fixed = TRUE)[[1]]
GROUP_COLORS <- setNames(strsplit(a[9],  "|", fixed = TRUE)[[1]], GROUP_LEVELS)
GROUP_BLOCK  <- setNames(strsplit(a[10], "|", fixed = TRUE)[[1]], GROUP_LEVELS)
# "all" brackets every pair; "within_block" restricts brackets to pairs inside
# one block. The block assignment is recorded in the stats table either way.
BRACKET_SCOPE <- if (length(a) >= 11 && nzchar(a[11])) a[11] else "all"
# Panel order, pipe separated, in the order the tracks were added to the track
# list -- the same order that becomes the heatmap's column order. Passed
# explicitly rather than inferred from the table's column names: anything that
# recovers the order by reading it back out of the data is one silent coercion
# away from ggplot falling back to alphabetical.
PANEL_LEVELS <- if (length(a) >= 12 && nzchar(a[12]))
  strsplit(a[12], "|", fixed = TRUE)[[1]] else character(0)

FONT_SCALE <- as.numeric(Sys.getenv("VIOLIN_FONT_SCALE", "2.25"))
fs <- function(x) x * FONT_SCALE
NCOL_PANELS <- as.integer(Sys.getenv("VIOLIN_NCOL", "4"))
if (!is.finite(NCOL_PANELS) || NCOL_PANELS < 1) NCOL_PANELS <- 4

MEDIAN_COL <- "#1a9d3b"
# A level with no rows renders as an empty slot, which is how the two marks are
# separated on the axis of the combined broad figure.
SPACER <- GROUP_LEVELS[!nzchar(trimws(GROUP_LEVELS))]
LEGEND_LEVELS <- GROUP_LEVELS[nzchar(trimws(GROUP_LEVELS))]

df <- suppressMessages(read_tsv(comb, show_col_types = FALSE, na = c("nan", "-nan", "NA", "")))
if (nrow(df) == 0) { message("empty summary table"); q(save = "no") }

# --- region group sizes -----------------------------------------------------
# One row of the summary table is one BIN, and every panel is drawn over the
# same bins, so n is a property of the region group and not of the panel. It is
# therefore reported once in the legend rather than repeated on all sixteen
# panels. Counted here, before the pivot, where one row is still one bin.
grp_n  <- df %>%
  mutate(Group = factor(Group, levels = GROUP_LEVELS)) %>%
  count(Group, .drop = FALSE)
n_of <- setNames(grp_n$n, as.character(grp_n$Group))

legend_labels <- vapply(LEGEND_LEVELS, function(g) {
  k <- unname(n_of[g])
  if (is.na(k)) g else sprintf("%s  (n = %s bins)", g, format(k, big.mark = ",", trim = TRUE))
}, character(1))
names(legend_labels) <- LEGEND_LEVELS
for (g in LEGEND_LEVELS) message("  ", legend_labels[[g]])

# The same counts go to disk beside the stats table, so the numbers quoted in
# the manuscript text can be taken from a file rather than read off a figure.
write_tsv(grp_n %>% filter(nzchar(trimws(as.character(Group)))) %>% rename(n_bins = n),
          sub("_pairwise_stats\\.tsv$", "_group_sizes.tsv", stats_tsv))

long <- df %>%
  pivot_longer(cols = -c(chr, start, end, Group), names_to = "panel", values_to = "Signal") %>%
  mutate(Signal = suppressWarnings(as.numeric(Signal)),
         Signal = ifelse(is.na(Signal), 0, Signal),
         Group  = factor(Group, levels = GROUP_LEVELS))

# GC is a percentage, not an assay readout: it is plotted as mean % GC with no
# log transform. Unit detection happens once, from the column maximum.
gc_max <- suppressWarnings(max(long$Signal[long$panel == gc_col], na.rm = TRUE))
gc_mult <- if (is.finite(gc_max) && gc_max <= 1) 100 else 1
long <- long %>%
  mutate(is_pct = panel == gc_col,
         Value  = ifelse(is_pct, Signal * gc_mult, log2(Signal + 1)))

# --- PANEL ORDER ------------------------------------------------------------
# The order is taken from PANEL_LEVELS, supplied by the caller, and only
# filtered down to the columns actually present. Any column the caller did not
# name is appended and reported rather than dropped, so a track that appears in
# the data but not in the order is visible instead of vanishing.
data_panels <- names(df)[!(names(df) %in% c("chr", "start", "end", "Group"))]
if (length(PANEL_LEVELS) == 0) PANEL_LEVELS <- data_panels
panel_levels <- PANEL_LEVELS[PANEL_LEVELS %in% data_panels]
extra <- setdiff(data_panels, panel_levels)
if (length(extra) > 0) {
  message("  WARNING: columns absent from the supplied panel order, appended last: ",
          paste(extra, collapse = ", "))
  panel_levels <- c(panel_levels, extra)
}
long$panel <- factor(long$panel, levels = panel_levels)
# A single unmatched name would turn its rows into NA and silently move that
# dataset into a stray "NA" facet, so the mismatch is a hard error.
if (any(is.na(long$panel))) {
  stop("panel names in the table do not match the supplied panel order; got: ",
       paste(setdiff(unique(as.character(df$Group)), panel_levels), collapse = ", "))
}
message("  panel order (", length(panel_levels), " panels): ",
        paste(seq_along(panel_levels), panel_levels, sep = ". ", collapse = "  |  "))

# --- pairwise statistics -----------------------------------------------------
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
      panel = p, group1 = cb[1], group2 = cb[2],
      n1 = length(x), n2 = length(y),
      median1 = median(x), median2 = median(y),
      unit = if (d$is_pct[1]) "percent" else "log2",
      p = wt$p.value, stringsAsFactors = FALSE)
  }
}
# bind_rows(list()) returns a 0x0 tibble with no columns, which would break the
# joins below, so the empty case is given the right shape explicitly.
pair_tests <- if (length(pairs) > 0) bind_rows(pairs) else
  data.frame(panel = character(), group1 = character(), group2 = character(),
             n1 = integer(), n2 = integer(), median1 = double(), median2 = double(),
             unit = character(), p = double(), stringsAsFactors = FALSE)
if (nrow(pair_tests) > 0) {
  pair_tests <- pair_tests %>% group_by(panel) %>%
    mutate(p.adj = p.adjust(p, method = "BH")) %>% ungroup() %>%
    # Every pair of groups present in a panel is tested and BH-corrected within
    # that panel. Whether a pair also gets a bracket drawn is BRACKET_SCOPE's
    # decision; same_block is recorded either way so cross-block pairs stay
    # identifiable in the output table.
    mutate(same_block = unname(GROUP_BLOCK[group1] == GROUP_BLOCK[group2]),
           bracket_shown = if (BRACKET_SCOPE == "within_block") same_block else TRUE)
  write_tsv(pair_tests, stats_tsv)
  message("  wrote ", nrow(pair_tests), " pairwise tests (",
          sum(pair_tests$bracket_shown), " bracketed) -> ", stats_tsv)
  # pair_tests$panel was built from a character loop variable. Putting it back
  # on the same factor levels as the plot data keeps the brackets attached to
  # the facets they belong to, and keeps the join with `pan` factor-to-factor
  # instead of silently coercing both sides to character.
  pair_tests$panel <- factor(pair_tests$panel, levels = panel_levels)
} else {
  pair_tests$p.adj <- double(); pair_tests$bracket_shown <- logical()
  message("  no testable pairs")
}

# Numeric adjusted p-values, never stars. Wilcoxon uses the normal
# approximation at these sample sizes and underflows to exactly 0 for very
# large samples, so a floored value is reported as a bound.
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
br0 <- pair_tests %>% filter(!is.na(p.adj), bracket_shown)
# Six groups give fifteen brackets per panel. At the default spacing that stack
# would push the violins into the bottom third of the panel, so the gap tightens
# once there are more than eight tiers. Brackets are ordered shortest span
# first, so the stack rises from within-block pairs to the widest cross-block
# ones and the two kinds stay visually separable.
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
  # The legend doubles as the bin-count key: every entry carries the number of
  # regions behind that violin, so the figure states its own group sizes.
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE, name = "Region group (bins per group)",
                    breaks = LEGEND_LEVELS, labels = legend_labels, na.value = "white") +
  guides(fill = guide_legend(nrow = if (length(LEGEND_LEVELS) > 4) 2 else 1,
                             byrow = TRUE, override.aes = list(alpha = 1))) +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ panel, scales = "free_y", ncol = NCOL_PANELS) +
  labs(x = NULL, y = "log2(mean signal + 1)   |   GC panel: mean GC content (%)",
       title = title_txt, subtitle = subtitle_txt, caption = caption_txt) +
  theme_minimal(base_size = fs(10)) +
  theme(strip.text = element_text(size = fs(10), face = "bold", colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        panel.spacing = unit(0.7, "cm"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = fs(7.5), colour = "black"),
        axis.text.y = element_text(size = fs(9), colour = "black"),
        axis.title.y = element_text(size = fs(10), face = "bold"),
        plot.title = element_text(size = fs(13), face = "bold"),
        plot.subtitle = element_text(size = fs(9), face = "italic"),
        plot.caption = element_text(size = fs(6.5), hjust = 0, colour = "grey30", lineheight = 1.3),
        legend.position = "bottom",
        legend.title = element_text(size = fs(10), face = "bold"),
        legend.text = element_text(size = fs(9)))

n_pan <- length(panel_levels); nrows <- ceiling(n_pan / NCOL_PANELS)
ggsave(out_pdf, p, width = NCOL_PANELS * (4.0 + 0.62 * length(GROUP_LEVELS)),
       height = nrows * (5.6 + 0.42 * n_tiers) + 4.5,
       device = "pdf", bg = "white", limitsize = FALSE)
message("  ", n_pan, " panels in ", nrows, " rows of ", NCOL_PANELS,
        "  |  ", n_tiers, " bracket tiers per panel")
message("wrote ", out_pdf)
RVIOLIN

export VIOLIN_FONT_SCALE VIOLIN_NCOL

# =============================================================================
# VIOLIN HELPER
# =============================================================================
# Summarises each region group over the current BW/LBL track list and renders
# one figure. VI_BED / VI_DISP hold the groups in axis order; a display name of
# a single space marks a spacer column and has no BED behind it.
run_violin() {   # tag, out_pdf, stats_tsv, title, subtitle, caption, levels, colours, blocks
    local TAG="$1" OUT="$2" STATS="$3" TITLE="$4" SUB="$5" CAP="$6" LEV="$7" COLS="$8" BLKS="$9"
    local BRACKET_SCOPE="${BRACKET_SCOPE:-all}"
    [[ ${#BW[@]} -gt 0 ]] || { echo "     no tracks -- violin ${TAG} skipped" >&2; return 0; }
    local COMBINED="$VI_DIR/${TAG}_signal.tsv" COLSIG="$VI_DIR/${TAG}_columns.txt" COLNOW i bed tab l

    # The header is written in LBL order, and that column order becomes the
    # panel order in the figure.
    { printf 'chr\tstart\tend'; for l in ${LBL[@]+"${LBL[@]}"}; do printf '\t%s' "$l"; done
      printf '\tGroup\n'; } > "$COMBINED"

    # A cached .tab holds one column per bigWig it was made from. If the track
    # list changes, reusing it would silently pair old numbers with new column
    # names, so the column set is fingerprinted.
    COLNOW=$(printf '%s\n' ${BW[@]+"${BW[@]}"})
    if [[ ! -f "$COLSIG" || "$(cat "$COLSIG")" != "$COLNOW" ]]; then
        echo "     track list changed -- clearing cached violin summaries" >&2
        rm -f "$VI_DIR/${TAG}_sig_"*.tab
        printf '%s\n' "$COLNOW" > "$COLSIG"
    fi

    [[ ${#VI_BED[@]} -gt 0 ]] || { echo "     no region groups -- ${TAG} skipped" >&2; return 0; }
    for i in "${!VI_BED[@]}"; do
        bed="${VI_BED[$i]}"
        [[ -n "$bed" && -s "$bed" ]] || continue      # spacer, or an empty group
        tab="$VI_DIR/${TAG}_sig_${i}.tab"
        if [[ ! -s "$tab" || "$bed" -nt "$tab" ]]; then
            multiBigwigSummary BED-file -b ${BW[@]+"${BW[@]}"} --BED "$bed" \
                -o "${tab%.tab}.npz" --outRawCounts "$tab" -p "$THREADS" >/dev/null 2>&1
            rm -f "${tab%.tab}.npz"
        fi
        tail -n +2 "$tab" | awk -v g="${VI_DISP[$i]}" 'BEGIN{OFS="\t"}{print $0, g}' >> "$COMBINED"
    done

    # The panel order handed to R is the track list itself, pipe separated, so
    # the violin panels and the heatmap columns are driven by the same array
    # and cannot drift apart. Track labels never contain "|".
    local PANELS
    PANELS=$(IFS='|'; printf '%s' "${LBL[*]}")

    echo "  -> violins: ${TAG}" >&2
    echo "     panel order: $(printf '%s; ' ${LBL[@]+"${LBL[@]}"})" >&2
    Rscript "$RS_DIR/violin.R" "$COMBINED" "$OUT" "$STATS" \
        "$TITLE" "$SUB" "$CAP" "GC content" "$LEV" "$COLS" "$BLKS" "$BRACKET_SCOPE" "$PANELS"
}

# =============================================================================
# MAIN LOOP
# =============================================================================
build_all_prom_stranded

for CL in "${CELL_LINES[@]}"; do
    echo "#########  ${CL}  #########"
    PRO_DAST=$(merge_proseq "${PRO_DAST_PLUS[$CL]:-}" "${PRO_DAST_MINUS[$CL]:-}" \
                 "$MAT_DIR/PROseq_Dastidar_${CL}_bothStrands.bw")
    PRO_CORE=$(merge_proseq "${PRO_CORE_PLUS[$CL]:-}" "${PRO_CORE_MINUS[$CL]:-}" \
                 "$MAT_DIR/PROseq_Core_${CL}_bothStrands.bw")

    # --- region groups for every epitope, once --------------------------------
    HAVE=()
    for EPI in "${EPITOPES[@]}"; do
        if build_region_groups "$CL" "$EPI"; then
            HAVE+=("$EPI")
        else
            echo "  WARNING: no differential result for ${CL} ${EPI} -- skipped" >&2
        fi
    done

    for EPI in ${HAVE[@]+"${HAVE[@]}"}; do
        PRIMARY="${EPI_PRIMARY_SET[$EPI]}"; CONTROL="${EPI_CONTROL_SET[$EPI]}"
        PRE="$REG_DIR/${CL}_${EPI}"
        echo "--- ${CL} ${EPI} (${SET_DESC[$PRIMARY]} vs ${SET_DESC[$CONTROL]}) ---"

        # The track list is built ONCE here and used by both the violins and the
        # supplementary heatmap, so the two figures cannot show different
        # datasets. Building it separately per figure made that agreement a
        # coincidence rather than a guarantee.
        if is_promoter_mark "$EPI"; then build_tracks_promoter "$CL" "$EPI"
        else                             build_tracks_broad "$CL"; fi
        if [[ ${#BW[@]} -eq 0 ]]; then
            echo "     no tracks available -- ${CL} ${EPI} skipped" >&2; continue
        fi
        echo "     ${#LBL[@]} datasets: $(printf '%s; ' ${LBL[@]+"${LBL[@]}"})" >&2

        # --- violins: promoter marks only ------------------------------------
        # The broad marks share a single combined figure, built after this loop.
        # Violins keep the ChIP -> CUT&Tag left-to-right order used by every
        # other figure in the manuscript, and keep the per-epitope negative
        # control.
        if [[ "$PLOT_VIOLINS" == "1" ]] && is_promoter_mark "$EPI"; then
            VI_BED=("${PRE}_chip.bed" "${PRE}_nondiff.bed" "${PRE}_cnt.bed" "${PRE}_control.bed")
            VI_DISP=("ChIP-seq Enriched" "Non-differential" "CUT&Tag Enriched" "Negative Control")
            run_violin "${CL}_${EPI}" \
                "$VI_DIR/${CL}_${EPI}_violins.pdf" "$VI_DIR/${CL}_${EPI}_pairwise_stats.tsv" \
                "${CL} ${EPI}: reference signal over differential region groups" \
                "raw log2(signal + 1), no background normalisation  |  regions from ${REGION_COMPARISON[$CL,$EPI]}  |  status: FDR < ${FDR_THRESH}, |log2FC| > ${FC_THRESH}" \
                "Groups: ${SET_DESC[$PRIMARY]} split by differential status, plus a negative control of ${SET_DESC[$CONTROL]}.
The legend gives the number of bins behind each violin; the same counts are written to *_group_sizes.tsv.
Green numbers are group MEDIANS. Brackets carry Benjamini-Hochberg adjusted p-values from two-sided Wilcoxon rank-sum tests, corrected
WITHIN each panel so a panel's labels do not depend on which other panels share the figure. With thousands of regions per group nearly
every comparison is significant, so read the effect size rather than the p-value; the CUT&Tag-enriched vs ChIP-seq-enriched pair in
particular reflects the DESeq2 contrast that DEFINED those groups, whereas comparisons against the non-differential and negative-control
groups are independent. Panels are raw signal and are NOT comparable to one another: height reflects each track's depth and normalisation." \
                "ChIP-seq Enriched|Non-differential|CUT&Tag Enriched|Negative Control" \
                "#3936ff|grey70|#ff3b3b|grey25" \
                "all|all|all|all"
        fi

        # --- promoter matrix and heatmap -------------------------------------
        if [[ "$PLOT_HEATMAPS" == "1" ]] && is_promoter_mark "$EPI"; then
            HM_BED=("${PRE}_nondiff_TSS.bed" "${PRE}_cnt_TSS.bed" "${PRE}_chip_TSS.bed" "${PRE}_control_TSS.bed")
            HM_GROUP_DISP=("${EPI} non-differential" "${EPI} CUT&Tag-enriched" \
                           "${EPI} ChIP-seq-enriched" "Polycomb-repressed promoters (control)")
            MAT="$MAT_DIR/${CL}_${EPI}_matrix.gz"
            if ensure_matrix "$MAT" TSS "$PROM_WINDOW" ${HM_BED[@]+"${HM_BED[@]}"}; then
                mapfile -t MAIN_WANT < <(main_labels_promoter "$EPI")
                WANT=(); is_main_cell "$CL" && WANT=(${MAIN_WANT[@]+"${MAIN_WANT[@]}"})
                render_heatmap_pair "$MAT" "$HM_DIR/${CL}_${EPI}_heatmap" \
                    "${CL} ${EPI}: active promoters by differential status, vs Polycomb-repressed promoters
regions from: ${REGION_COMPARISON[$CL,$EPI]}" \
                    "TSS" ${WANT[@]+"${WANT[@]}"}
            fi
        fi
    done

    # --- broad marks: ONE combined heatmap --------------------------------
    # H3K36me3 and H3K27me3 are mutually exclusive, so each mark's bins act as
    # the other's negative control and both belong in a single figure. The six
    # region groups run H3K36me3 non-diff / CUT&Tag / ChIP-seq over gene-body
    # bins, then the same three for H3K27me3 over Polycomb bins.
    if [[ "$PLOT_HEATMAPS" == "1" ]]; then
        K36="$REG_DIR/${CL}_H3K36me3"; K27="$REG_DIR/${CL}_H3K27me3"
        if [[ -s "${K36}_nondiff.bed" && -s "${K27}_nondiff.bed" ]]; then
            echo "--- ${CL} broad marks (combined H3K36me3 + H3K27me3) ---"
            build_tracks_broad "$CL"
            HM_BED=("${K36}_nondiff.bed" "${K36}_cnt.bed" "${K36}_chip.bed"
                    "${K27}_nondiff.bed" "${K27}_cnt.bed" "${K27}_chip.bed")
            HM_GROUP_DISP=(
                "H3K36me3 non-differential (gene body)" "H3K36me3 CUT&Tag-enriched (gene body)"
                "H3K36me3 ChIP-seq-enriched (gene body)" "H3K27me3 non-differential (Polycomb)"
                "H3K27me3 CUT&Tag-enriched (Polycomb)"  "H3K27me3 ChIP-seq-enriched (Polycomb)")
            MAT_B="$MAT_DIR/${CL}_broadMarks_matrix.gz"
            if ensure_matrix "$MAT_B" center "$BROAD_WINDOW" ${HM_BED[@]+"${HM_BED[@]}"}; then
                mapfile -t MAIN_WANT < <(main_labels_broad)
                WANT=(); is_main_cell "$CL" && WANT=(${MAIN_WANT[@]+"${MAIN_WANT[@]}"})
                render_heatmap_pair "$MAT_B" "$HM_DIR/${CL}_broadMarks_heatmap" \
                    "${CL} broad marks: H3K36me3 gene-body bins and H3K27me3 Polycomb bins by differential status
regions from: ${REGION_COMPARISON[$CL,H3K36me3]}  |  ${REGION_COMPARISON[$CL,H3K27me3]}" \
                    "bin center" ${WANT[@]+"${WANT[@]}"}
            fi
        else
            echo "  WARNING: broad-mark region groups missing for ${CL} -- combined heatmap skipped" >&2
        fi
    fi

    # --- broad marks: ONE combined violin figure ----------------------------
    # Same six region groups and the same track list as the combined heatmap,
    # so the two figures show the same data. A blank level between the two
    # marks separates the gene-body block from the Polycomb block on the axis.
    if [[ "$PLOT_VIOLINS" == "1" ]]; then
        K36="$REG_DIR/${CL}_H3K36me3"; K27="$REG_DIR/${CL}_H3K27me3"
        if [[ -s "${K36}_nondiff.bed" && -s "${K27}_nondiff.bed" ]]; then
            echo "--- ${CL} broad marks (combined H3K36me3 + H3K27me3) violins ---"
            build_tracks_broad "$CL"
            VI_BED=("${K36}_chip.bed" "${K36}_nondiff.bed" "${K36}_cnt.bed" ""
                    "${K27}_chip.bed" "${K27}_nondiff.bed" "${K27}_cnt.bed")
            VI_DISP=("H3K36me3 ChIP-seq Enriched" "H3K36me3 Non-differential" "H3K36me3 CUT&Tag Enriched" " "
                     "H3K27me3 ChIP-seq Enriched" "H3K27me3 Non-differential" "H3K27me3 CUT&Tag Enriched")
            run_violin "${CL}_broadMarks" \
                "$VI_DIR/${CL}_broadMarks_violins.pdf" "$VI_DIR/${CL}_broadMarks_pairwise_stats.tsv" \
                "${CL} broad marks: reference signal over differential region groups" \
                "raw log2(signal + 1), no background normalisation  |  H3K36me3 over active gene-body bins, H3K27me3 over Polycomb-repressed bins  |  status: FDR < ${FDR_THRESH}, |log2FC| > ${FC_THRESH}" \
                "Left block: active gene-body bins split by H3K36me3 differential status. Right block: Polycomb-repressed bins split by H3K27me3
differential status. Because the two marks are mutually exclusive, each block also serves as the other's negative control.
The legend gives the number of bins behind each violin; the same counts are written to *_group_sizes.tsv.
Green numbers are group MEDIANS. Brackets carry Benjamini-Hochberg adjusted p-values from two-sided Wilcoxon rank-sum tests over all fifteen
pairs of region groups, corrected WITHIN each panel so a panel's labels do not depend on which other panels share the figure. Cross-block
pairs compare a different mark over a different region type, so they measure how far apart the two chromatin classes sit rather than a
difference between methods; the within-block pairs are the method comparison. With thousands of regions per group nearly every comparison is
significant, so read the effect size rather than the p-value. Panels are raw signal and are NOT comparable to one another: height reflects
each track's depth and normalisation. The same_block column of *_pairwise_stats.tsv marks which pairs are which." \
                "H3K36me3 ChIP-seq Enriched|H3K36me3 Non-differential|H3K36me3 CUT&Tag Enriched| |H3K27me3 ChIP-seq Enriched|H3K27me3 Non-differential|H3K27me3 CUT&Tag Enriched" \
                "#3936ff|grey70|#ff3b3b|white|#3936ff|grey70|#ff3b3b" \
                "H3K36me3|H3K36me3|H3K36me3|spacer|H3K27me3|H3K27me3|H3K27me3"
        fi
    fi
done

echo "========================================================="
echo " Done."
echo "   Region groups : $REG_DIR"
echo "   Matrices      : $MAT_DIR"
echo "   Heatmaps      : $HM_DIR"
echo "                   <cell>_<mark>_heatmap_{main,supplementary}.pdf   (promoter marks)"
echo "                   <cell>_broadMarks_heatmap_{main,supplementary}.pdf"
echo "   Violins       : $VI_DIR"
echo "                   <cell>_<mark>_violins.pdf, <cell>_broadMarks_violins.pdf"
echo "                   + *_pairwise_stats.tsv, *_group_sizes.tsv"
echo "========================================================="