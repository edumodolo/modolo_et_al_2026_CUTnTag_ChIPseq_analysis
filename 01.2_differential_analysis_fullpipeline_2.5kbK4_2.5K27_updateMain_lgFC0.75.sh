#!/usr/bin/env bash
set -euo pipefail

# =========================================================================
# UNIFIED GROUND TRUTH REGIONS & DIFFERENTIAL ANALYSIS PIPELINE
# -------------------------------------------------------------------------
# 1. Selects regions using ChromHMM (+/- 1kb TSS for H3K4me3, +/- 2.5kb for H3K27ac).
# 2. Runs differential analysis using NATIVE DESeq2. 
#    * Skips calculation if outputs already exist.
# 3. Generates UCSC Custom Tracks colored by differential status.
# 4. Generates all visualization plots (harmonized Left=ChIP, Right=CnT).
# 5. Calculates Variance Explained (R-squared) feature importance.
# =========================================================================

# --------------------------- USER CONFIG ---------------------------------
WORK_DIR="/home/emodolo/gpfs/modolo_et_al/try26_ground_truth_analysis/output/differential_analysis_nativeDESeq"
REGIONS_DIR="$WORK_DIR/regions"
TRACKS_DIR="$REGIONS_DIR/tracks"
DIFF_DIR="$WORK_DIR/differential"
RS_DIR="$WORK_DIR/Rscripts"
FILT_DIR="$DIFF_DIR/regions_filtered"
PLOTS_DIR="$WORK_DIR/feature_importance"
TMP="$WORK_DIR/_tmp"

mkdir -p "$REGIONS_DIR" "$TRACKS_DIR" "$DIFF_DIR" "$RS_DIR" "$FILT_DIR" "$PLOTS_DIR" "$TMP"

CSV_FILE="/gpfs/data01/gorenlab/emodolo/modolo_et_al/master_sample_file_hub/Master_sample_metrics_correctedCnTBWs.tsv"
GENE_BODIES="/gpfs/data01/gorenlab/emodolo/modolo_et_al/genome_data/ucsc_hg38/canonical_geneBodies.hg38.bed"
BLACKLIST="$HOME/gpfs/genomes/bedfiles/hg38-blacklist.v2.bed"
CHROM_SIZES="/home/emodolo/gpfs/Homer/data/genomes/hg38/chrom.sizes"

WEB_DIR="/home/emodolo/web/modolo_et_al/try26_ground_truth_regions"
WEB_URL="http://homer.ucsd.edu/emodolo/modolo_et_al/try26_ground_truth_regions"

GENOME="hg38"
GC_5bp_BW="/home/emodolo/gpfs/modolo_et_al/GC_content_bias_analysis/genome_wide_gc_bigwigs/gc5Base.bw"
THREADS=24

# --- Analysis Thresholds ---
FC_THRESH=0.75
FDR_THRESH=0.05

# --- Region Size Cutoffs ---
# H3K4me3 (Narrower, +/- 1kb)
H3K4ME3_UP=2500
H3K4ME3_DN=2500
H3K4ME3_WIN_SIZE=$((H3K4ME3_UP + H3K4ME3_DN))

# H3K27ac (Broader, +/- 2.5kb)
H3K27AC_UP=2500
H3K27AC_DN=2500
H3K27AC_WIN_SIZE=$((H3K27AC_UP + H3K27AC_DN))

# H3K27me3/H3K36me3 Bin Size
BIN_SIZE=3000

# --- Active-gene criteria ---
TSS_PROMOTER_WINDOW=1000
MIN_BODY_TX_FRAC=0.20
TSS_EXCLUDE_WINDOW=3000

CELL_LINES=(K562 MCF7)

declare -A CHROMHMM_BED
CHROMHMM_BED[K562]="/gpfs/data01/gorenlab/emodolo/modolo_et_al/try19_comprehensive_pipeline/output_files/repressed_regions_check/K562_hg38_ChromHMM_18state_ENCFF963KIA.bed"
CHROMHMM_BED[MCF7]="/gpfs/data01/gorenlab/emodolo/modolo_et_al/try19_comprehensive_pipeline/output_files/repressed_regions_check/MCF7_hg38_ChromHMM_18state_ENCFF985EWD.bed"

declare -A ATAC_BW
ATAC_BW[K562]="/home/emodolo/gpfs/modolo_et_al/merged_atac-seq/output/ATAC_SRR_9234_9235_cpm.bw"
ATAC_BW[MCF7]="/gpfs/data01/gorenlab/emodolo/modolo_et_al/try13_unified_cnt_chip_processing/output_files/alignment/bigwig/CnT_ATAC_merg_SRR_2779_2778_MCF7_ID262_ATAC_cpm.bw"

# =========================================================================
# PART 1: REGION GENERATION 
# =========================================================================
HAVE_BB=1
command -v bedToBigBed >/dev/null 2>&1 || { echo "⚠️  bedToBigBed not found in PATH -- skipping .bb/track creation." >&2; HAVE_BB=0; }

[[ -s "$BLACKLIST" ]]    || { echo "❌ ERROR: blacklist not found: $BLACKLIST" >&2; exit 1; }
[[ -s "$GENE_BODIES" ]]  || { echo "❌ ERROR: gene bodies not found: $GENE_BODIES" >&2; exit 1; }
[[ -s "$CHROM_SIZES" ]]  || { echo "❌ ERROR: chrom.sizes not found: $CHROM_SIZES" >&2; exit 1; }

mainchr() { awk -F'\t' '$1 ~ /^chr([0-9]+|[XY])$/' "$@"; }

get_states() {
    local CHMM="$1" STATES="$2"
    awk -F'\t' -v st="$STATES" 'BEGIN{split(st, a, " "); for(i in a) s[a[i]]=""} $4 in s {print $1"\t"$2"\t"$3}' "$CHMM" \
        | mainchr | sort -k1,1 -k2,2n | bedtools merge -i -
}

OUT_TRACK_TXT="$TRACKS_DIR/UCSC_Custom_Tracks.txt"
> "$OUT_TRACK_TXT"

publish_track() {
    [[ "$HAVE_BB" == "1" ]] || return 0
    local BB="$1" NAME="$2" DESC="$3"
    mkdir -p "$WEB_DIR"
    cp -f "$BB" "$WEB_DIR/"
    chmod 644 "$WEB_DIR/$(basename "$BB")" 2>/dev/null || true
    
    # Check if this track name is already in the file before appending
    if ! grep -q "name=\"${NAME}\"" "$OUT_TRACK_TXT"; then
        echo "track type=bigBed name=\"${NAME}\" description=\"${DESC}\" bigDataUrl=\"${WEB_URL}/$(basename "$BB")\" visibility=pack itemRgb=\"On\"" >> "$OUT_TRACK_TXT"
        echo "" >> "$OUT_TRACK_TXT"
    fi
}

for CL in "${CELL_LINES[@]}"; do
    echo "========================================================="
    echo " Generating Regions for $CL"
    echo "========================================================="
    CHMM="${CHROMHMM_BED[$CL]}"
    
    BL="$TMP/${CL}_blacklist.bed"
    cut -f1-3 "$BLACKLIST" | tr -d '\r' | mainchr | sort -k1,1 -k2,2n | bedtools merge -i - > "$BL"

    PROM_ANCHOR="$TMP/${CL}_activePromoter_anchor.bed"
    get_states "$CHMM" "TssA TssFlnk TssFlnkU TssFlnkD" | bedtools intersect -a - -b "$BL" -v | sort -k1,1 -k2,2n > "$PROM_ANCHOR"

    GENES_ALL="$TMP/${CL}_genes6.bed"
    cut -f1-6 "$GENE_BODIES" | tr -d '\r' | mainchr | sort -k1,1 -k2,2n > "$GENES_ALL"
    GENES="$TMP/${CL}_genes6_noBL.bed"
    bedtools intersect -a "$GENES_ALL" -b "$BL" -v | sort -k1,1 -k2,2n > "$GENES"

    ACT_BLK="$TMP/${CL}_actBody.bed"
    get_states "$CHMM" "Tx TxWk" > "$ACT_BLK"

    ACT_BODY_NAMES="$TMP/${CL}_actBody.names"
    bedtools coverage -a "$GENES" -b "$ACT_BLK" | awk -F'\t' -v f="$MIN_BODY_TX_FRAC" '($NF+0) >= f {print $4}' | sort -u > "$ACT_BODY_NAMES"

    TSSWIN="$TMP/${CL}_tssWindows.bed"
    awk -F'\t' -v w="$TSS_PROMOTER_WINDOW" 'BEGIN{OFS="\t"}
        { chr=$1; s=$2; e=$3; name=$4; str=$6;
          if(str=="-") tss=e; else tss=s;
          ws=tss-w; if(ws<0)ws=0; we=tss+w; if(we<=ws) we=ws+1;
          print chr, ws, we, name, ".", str }' "$GENES" | sort -k1,1 -k2,2n > "$TSSWIN"

    PROMACT_NAMES="$TMP/${CL}_promActive.names"
    bedtools intersect -a "$TSSWIN" -b "$PROM_ANCHOR" -u | cut -f4 | sort -u > "$PROMACT_NAMES"

    ACTIVE_NAMES="$TMP/${CL}_ACTIVE.names"
    comm -12 "$ACT_BODY_NAMES" "$PROMACT_NAMES" > "$ACTIVE_NAMES"

    ACT_GENES="$TMP/${CL}_active_selectedGenes.bed"
    awk -F'\t' 'NR==FNR{k[$1]=1; next} ($4 in k)' "$ACTIVE_NAMES" "$GENES" | sort -k1,1 -k2,2n > "$ACT_GENES"

    # --- H3K4me3 ---
    H3K4ME3_CLEAN="$TMP/${CL}_H3K4me3_promoter_clean.bed"
    awk -F'\t' -v up="$H3K4ME3_UP" -v dn="$H3K4ME3_DN" '
        NR==FNR { size[$1]=$2; next }
        {
          chr=$1; s=$2; e=$3; str=$6;
          if (str=="-") { tss=e; ws=tss-dn; we=tss+up }
          else          { tss=s; ws=tss-up; we=tss+dn }
          if (ws<0) ws=0;
          if ((chr in size) && we>size[chr]) we=size[chr];
          if (we>ws) print chr"\t"ws"\t"we
        }' "$CHROM_SIZES" "$ACT_GENES" \
      | awk -F'\t' -v w="$H3K4ME3_WIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' \
      | sort -k1,1 -k2,2n -u | bedtools intersect -a - -b "$BL" -v > "$H3K4ME3_CLEAN"

    awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1, $2, $3, cl"_H3K4me3_"(++n)}' "$H3K4ME3_CLEAN" > "$REGIONS_DIR/${CL}_H3K4me3_regions.bed"

    # --- H3K27ac ---
    H3K27AC_CLEAN="$TMP/${CL}_H3K27ac_promoter_clean.bed"
    awk -F'\t' -v up="$H3K27AC_UP" -v dn="$H3K27AC_DN" '
        NR==FNR { size[$1]=$2; next }
        {
          chr=$1; s=$2; e=$3; str=$6;
          if (str=="-") { tss=e; ws=tss-dn; we=tss+up }
          else          { tss=s; ws=tss-up; we=tss+dn }
          if (ws<0) ws=0;
          if ((chr in size) && we>size[chr]) we=size[chr];
          if (we>ws) print chr"\t"ws"\t"we
        }' "$CHROM_SIZES" "$ACT_GENES" \
      | awk -F'\t' -v w="$H3K27AC_WIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' \
      | sort -k1,1 -k2,2n -u | bedtools intersect -a - -b "$BL" -v > "$H3K27AC_CLEAN"

    awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1, $2, $3, cl"_H3K27ac_"(++n)}' "$H3K27AC_CLEAN" > "$REGIONS_DIR/${CL}_H3K27ac_regions.bed"


    # --- H3K36me3 ---
    ACT_TSS_EXCL="$TMP/${CL}_activeTSS_excl.bed"
    awk -F'\t' 'NR==FNR{k[$1]=1; next} ($4 in k)' "$PROMACT_NAMES" "$GENES" \
        | awk -F'\t' -v w="$TSS_EXCLUDE_WINDOW" 'BEGIN{OFS="\t"}
            { chr=$1; s=$2; e=$3; str=$6;
              if(str=="-") tss=e; else tss=s;
              ws=tss-w; if(ws<0)ws=0; we=tss+w; if(we<=ws) we=ws+1;
              print chr, ws, we }' \
        | sort -k1,1 -k2,2n | bedtools merge -i - > "$ACT_TSS_EXCL"

    ACT_REGIONS="$TMP/${CL}_active_geneRegions.bed"
    sort -k1,1 -k2,2n "$ACT_GENES" | bedtools merge -i - -c 4 -o distinct > "$ACT_REGIONS"
    ACT_NAMED="$TMP/${CL}_act_named.bed"
    awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1, $2, $3, cl"_H3K36me3_act_"NR}' "$ACT_REGIONS" > "$ACT_NAMED"

    bedtools makewindows -b "$ACT_NAMED" -w "$BIN_SIZE" -i srcwinnum | \
        awk -F'\t' -v w="$BIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' | \
        sort -k1,1 -k2,2n | bedtools intersect -a - -b "$BL" -v | bedtools intersect -a - -b "$ACT_TSS_EXCL" -v > "$REGIONS_DIR/${CL}_H3K36me3_regions.bed"
        
    # --- H3K27me3 ---
    REPR_ALL="$TMP/${CL}_repr_all.bed"
    get_states "$CHMM" "ReprPC ReprPCWk TssBiv EnhBiv" > "$REPR_ALL"
    REPR_CORE="$TMP/${CL}_repr_core.bed"
    get_states "$CHMM" "ReprPC ReprPCWk" > "$REPR_CORE"
    REPR_FILTERED="$TMP/${CL}_repr_filtered.bed"
    bedtools intersect -a "$REPR_ALL" -b "$REPR_CORE" -u | sort -k1,1 -k2,2n > "$REPR_FILTERED"
    
    REPR_NAMED="$TMP/${CL}_repr_named.bed"
    awk -F'\t' -v cl="$CL" 'BEGIN{OFS="\t"} {print $1, $2, $3, cl"_H3K27me3_rep_"NR}' "$REPR_FILTERED" > "$REPR_NAMED"
    
    bedtools makewindows -b "$REPR_NAMED" -w "$BIN_SIZE" -i srcwinnum | \
        awk -F'\t' -v w="$BIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2)==w' | \
        sort -k1,1 -k2,2n | bedtools intersect -a - -b "$BL" -v > "$REGIONS_DIR/${CL}_H3K27me3_regions.bed"
done


# =========================================================================
# PART 2: DIFFERENTIAL ANALYSIS (NATIVE DESeq2) HELPERS
# =========================================================================

cat > "$RS_DIR/run_deseq2_matrix.R" <<'RDESEQ'
#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr) })

args <- commandArgs(trailingOnly=TRUE)
raw_homer <- args[1]; out_diff  <- args[2]
n_chip <- as.numeric(args[3]); n_cnt  <- as.numeric(args[4])

# Suppress readr parsing warnings
df <- suppressWarnings(suppressMessages(read_tsv(raw_homer, comment="#", show_col_types=FALSE, name_repair="minimal")))
if(nrow(df) == 0) { message("Empty regions"); quit(save="no") }

meta <- df[, 1:19]
counts <- df[, 20:ncol(df)]

# Clean horrible Homer headers
clean_names <- sapply(colnames(counts), function(x) {
  prefix <- strsplit(x, " Tag Count")[[1]][1]
  basename(prefix)
})
colnames(counts) <- clean_names

# Force strict integers, clean NAs, suppress round() type warnings
counts <- as.data.frame(lapply(counts, function(x) as.integer(round(as.numeric(x)))))
counts[is.na(counts)] <- 0

condition <- factor(c(rep("ChIP", n_chip), rep("CnT", n_cnt)))
condition <- relevel(condition, ref="ChIP")

colData <- data.frame(condition = condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# Suppress dispersion trend warning by explicitly calling fitType='local'
dds <- DESeq(dds, fitType="local", quiet=TRUE)
res <- results(dds)

res$padj[is.na(res$padj)] <- 1
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0

out <- meta
out$"Log2 Fold Change" <- res$log2FoldChange
out$"adj. p-value" <- res$padj

norm_counts <- counts(dds, normalized=TRUE)
colnames(norm_counts) <- paste0(colnames(norm_counts), " Tag Count")
out <- cbind(out, norm_counts)

write_tsv(out, out_diff)
RDESEQ
chmod +x "$RS_DIR/run_deseq2_matrix.R"

get_tagdirs() {
    local tdirs=""
    for s in "$@"; do
        local tagdir
        tagdir=$(awk -F'\t' -v s="$s" 'NR>1 { gsub(/\r/,"",$2); if ($2==s) print $15 }' "$CSV_FILE" | head -n1)
        tdirs="$tdirs $tagdir"
    done
    echo "$tdirs"
}

get_group_plot_label() {
    local s="$1" raw
    raw=$(awk -F'\t' -v s="$s" 'NR>1 { gsub(/\r/,"",$2); if ($2==s) print $3 }' "$CSV_FILE" | head -n1)
    if [[ -z "$raw" ]]; then echo "$s"; else echo "${raw%_rep*}"; fi
}

filter_regions() {
    local EPI="$1" IN="$2" OUT="$3"
    case "$EPI" in
        H3K4me3)           awk -F'\t' -v w="$H3K4ME3_WIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2) == w' "$IN" > "$OUT" ;;
        H3K27ac)           awk -F'\t' -v w="$H3K27AC_WIN_SIZE" 'BEGIN{OFS="\t"} ($3-$2) == w' "$IN" > "$OUT" ;;
        H3K27me3|H3K36me3) awk -F'\t' -v w="$BIN_SIZE"         'BEGIN{OFS="\t"} ($3-$2) == w' "$IN" > "$OUT" ;;
        *)                 cp -f "$IN" "$OUT" ;;
    esac
}

compute_gc() {
    local REGION="$1" GC_OUT="$2"
    if [[ ! -s "$GC_OUT" ]]; then
        local NPZ="${GC_OUT%.tab}.npz"
        multiBigwigSummary BED-file -b "$GC_5bp_BW" --BED "$REGION" -o "$NPZ" --outRawCounts "$GC_OUT" -p "$THREADS" 2>/dev/null
        rm -f "$NPZ"
    fi
}

compute_atac() {
    local REGION="$1" ATAC_OUT="$2" CELL="$3" BW="${ATAC_BW[$CELL]:-}"
    if [[ -n "$BW" && -s "$BW" && ! -s "$ATAC_OUT" ]]; then
        local NPZ="${ATAC_OUT%.tab}.npz"
        multiBigwigSummary BED-file -b "$BW" --BED "$REGION" -o "$NPZ" --outRawCounts "$ATAC_OUT" -p "$THREADS" 2>/dev/null
        rm -f "$NPZ"
    fi
}

run_diff() {
    local DIFF_OUT="$1" REGION="$2" CHIP_DIRS="$3" CNT_DIRS="$4"
    if [[ ! -s "$DIFF_OUT" ]]; then
        local n_chip=$(echo $CHIP_DIRS | wc -w)
        local n_cnt=$(echo $CNT_DIRS | wc -w)
        local RAW_OUT="${DIFF_OUT%.txt}.raw.txt"
        
        echo "  -> Running DESeq2 -> $(basename "$DIFF_OUT")" >&2
        annotatePeaks.pl "$REGION" "$GENOME" -d $CHIP_DIRS $CNT_DIRS -noadj > "$RAW_OUT" 2>/dev/null
        Rscript "$RS_DIR/run_deseq2_matrix.R" "$RAW_OUT" "$DIFF_OUT" "$n_chip" "$n_cnt"
    else
        echo "  ⏩ Using cached diff: $(basename "$DIFF_OUT")" >&2
    fi
}

generate_colored_track() {
    local DIFF="$1" PREFIX="$2" DESC="$3"
    local BB_OUT="$TRACKS_DIR/${PREFIX}.bb"
    local B9="$TMP/${PREFIX}_bed9.bed"
    
    awk -F'\t' -v fc="$FC_THRESH" -v fdr="$FDR_THRESH" '
    NR==1 {
        for(i=1;i<=NF;i++) {
            if($i=="Log2 Fold Change") c_fc=i;
            if($i=="adj. p-value") c_fdr=i;
        }
        next
    }
    {
        pval=$c_fdr
        logfc=$c_fc
        
        # HARMONIZED ALIGNMENT: log2(CnT/ChIP) 
        # Positive = CnT Enriched (Red) | Negative = ChIP Enriched (Blue)
        if (pval != "NA" && pval < fdr && logfc > fc) color="255,59,59"      # Red
        else if (pval != "NA" && pval < fdr && logfc < -fc) color="57,54,255" # Blue
        else color="178,178,178" # Grey
        
        name=$1; if(name=="") name="region_"NR
        print $2 "\t" $3-1 "\t" $4 "\t" name "\t0\t.\t" $3-1 "\t" $4 "\t" color
    }' "$DIFF" | sort -k1,1 -k2,2n > "$B9"

    if [[ "$HAVE_BB" == "1" ]]; then
        bedToBigBed -type=bed9 "$B9" "$CHROM_SIZES" "$BB_OUT" 2>/dev/null
        publish_track "$BB_OUT" "$PREFIX" "$DESC"
    fi
}

add_comparison() {
    local MANIFEST="$1" EPI="$2" CELL="$3" COMPTAG="$4" CHIP_SAMPLES="$5" CNT_SAMPLES="$6"
    local RAW_REGION="$REGIONS_DIR/${CELL}_${EPI}_regions.bed"
    local REGION="$FILT_DIR/${CELL}_${EPI}_regions.binfilt.bed"
    
    filter_regions "$EPI" "$RAW_REGION" "$REGION"
    
    local GC="$DIFF_DIR/gc_${CELL}_${EPI}.tab"
    compute_gc "$REGION" "$GC"
    
    local ATAC="$DIFF_DIR/atac_${CELL}_${EPI}.tab"
    compute_atac "$REGION" "$ATAC" "$CELL"
    
    local DIFF="$DIFF_DIR/${CELL}_${EPI}_${COMPTAG}_DiffPeaks.txt"
    run_diff "$DIFF" "$REGION" "$(get_tagdirs $CHIP_SAMPLES)" "$(get_tagdirs $CNT_SAMPLES)"
    
    generate_colored_track "$DIFF" "${CELL}_${EPI}_${COMPTAG}" "$CELL $EPI Diff vs $COMPTAG (Red=CnT, Blue=ChIP)"
    
    local CNT_LBL CHIP_LBL first_cnt first_chip
    first_cnt=$(echo $CNT_SAMPLES | awk '{print $1}')
    first_chip=$(echo $CHIP_SAMPLES | awk '{print $1}')
    CNT_LBL=$(get_group_plot_label "$first_cnt")
    CHIP_LBL=$(get_group_plot_label "$first_chip")
    
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$EPI" "$CELL" "$CNT_LBL" "$CHIP_LBL" "$DIFF" "$GC" "$ATAC" >> "$MANIFEST"
}

# =========================================================================
# PART 3: R PLOTTING GENERATORS
# =========================================================================
cat > "$RS_DIR/common.R" <<'RCOMMON'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(patchwork); library(ggrastr); library(ggtext); library(scales)
})
options(warn = -1)

FC_THRESH   <- 0.75
FDR_THRESH  <- 0.05
CHIP_COL    <- "#3936ff"
CNT_COL     <- "#ff3b3b"
NONDIFF_COL <- "grey70"
MEAN_COL    <- "#1a9d3b"

# HARMONIZED ALIGNMENT: ChIP on the Left (Level 1), CnT on the Right (Level 3)
STATUS_LEVELS <- c("ChIP-seq Enriched", "Non-differential", "CUT&Tag Enriched")
STATUS_COLORS <- c("ChIP-seq Enriched" = CHIP_COL, "Non-differential" = NONDIFF_COL, "CUT&Tag Enriched" = CNT_COL)

prep_gc <- function(file) {
  if (!file.exists(file)) return(data.frame())
  df <- suppressMessages(read_tsv(file, show_col_types = FALSE))
  colnames(df) <- c("Chr", "Start", "End", "GC_Raw")
  df %>% mutate(GC_Percent = ifelse(GC_Raw <= 1.0, GC_Raw * 100, GC_Raw))
}

prep_atac <- function(file) {
  if (!file.exists(file)) return(data.frame())
  df <- suppressMessages(read_tsv(file, show_col_types = FALSE))
  colnames(df) <- c("Chr", "Start", "End", "ATAC_Raw")
  df %>% mutate(ATAC_log2 = log2(pmax(ATAC_Raw, 0) + 1))
}

load_comparison <- function(diff_path, gc_path, atac_path, epitope, cell, cnt_label, chip_label, col_order) {
  if (!file.exists(diff_path)) return(data.frame())
  diff <- suppressMessages(read_tsv(diff_path, show_col_types = FALSE, name_repair = "unique", comment = ""))
  colnames(diff)[1:4] <- c("PeakID", "Chr", "Start", "End")
  diff$Start <- diff$Start - 1
  fc_i <- which(grepl("Log2 Fold Change", colnames(diff), ignore.case = TRUE))[1]
  fd_i <- which(grepl("adj. p-value",     colnames(diff), ignore.case = TRUE))[1]
  colnames(diff)[fc_i] <- "log2FC"
  colnames(diff)[fd_i] <- "FDR"
  cnt_cols <- grep("Tag Count", colnames(diff), ignore.case = TRUE)
  diff$BaseMean <- rowMeans(diff[, cnt_cols, drop = FALSE], na.rm = TRUE)

  gc <- prep_gc(gc_path)
  df <- inner_join(diff, gc, by = c("Chr", "Start", "End")) %>%
    filter(!is.na(GC_Percent), !is.na(log2FC), !is.na(FDR)) %>%
    mutate(
      Status = case_when(
        log2FC >  FC_THRESH & FDR < FDR_THRESH ~ "CUT&Tag Enriched",
        log2FC < -FC_THRESH & FDR < FDR_THRESH ~ "ChIP-seq Enriched",
        TRUE ~ "Non-differential"),
      nLog10_FDR    = pmin(ifelse(FDR == 0, 300, -log10(FDR)), 300),
      BaseMean_Log2 = log2(BaseMean + 1),
      Epitope   = epitope, Cell_Line = cell, ColKey = paste0(epitope, " | ", cell, " | ", cnt_label),
      ColOrder  = col_order, Cnt_Full = cnt_label, Chip_Full = chip_label)

  atac <- prep_atac(atac_path)
  if (nrow(atac) > 0) df <- left_join(df, atac[, c("Chr", "Start", "End", "ATAC_log2")], by = c("Chr", "Start", "End")) else df$ATAC_log2 <- NA_real_

  df %>% select(PeakID, Chr, Start, End, GC_Percent, ATAC_log2, log2FC, FDR, nLog10_FDR, BaseMean_Log2, Status, Epitope, Cell_Line, ColKey, ColOrder, Cnt_Full, Chip_Full)
}

load_manifest <- function(path) {
  m <- suppressMessages(read_tsv(path, show_col_types = FALSE))
  out <- vector("list", nrow(m))
  for (i in seq_len(nrow(m))) out[[i]] <- load_comparison(m$diff_path[i], m$gc_path[i], m$atac_path[i], m$epitope[i], m$cell_line[i], m$cnt_label[i], m$chip_label[i], i)
  bind_rows(out)
}

format_p <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.0001) return("****"); if (p < 0.001) return("***")
  if (p < 0.01)  return("**");    if (p < 0.05)  return("*")
  "ns"
}

make_annotations <- function(big_df, mean_y, b_y1, b_y2, metric_col) {
  stat_df <- big_df %>% group_by(ColKey, Status, .drop = FALSE) %>% summarise(n = n(), mean_val = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>% mutate(y_pos_mean = mean_y, mean_label = ifelse(n > 0, sprintf("%.2f", mean_val), ""))
  b_df <- data.frame()
  for (k in unique(as.character(big_df$ColKey))) {
    sub <- big_df %>% filter(as.character(ColKey) == k)
    # HARMONIZED ALIGNMENT: Map factors strictly left-to-right to build the significance bars
    g1 <- sub[[metric_col]][sub$Status == "ChIP-seq Enriched"]
    g2 <- sub[[metric_col]][sub$Status == "Non-differential"]
    g3 <- sub[[metric_col]][sub$Status == "CUT&Tag Enriched"]
    padj <- p.adjust(c(
      tryCatch(wilcox.test(g1, g2)$p.value, error = function(e) NA),
      tryCatch(wilcox.test(g2, g3)$p.value, error = function(e) NA),
      tryCatch(wilcox.test(g1, g3)$p.value, error = function(e) NA)), method = "BH")
    b_df <- rbind(b_df,
      data.frame(x = 1, xend = 2, y = b_y1, label = format_p(padj[1]), ColKey = k),
      data.frame(x = 2, xend = 3, y = b_y1, label = format_p(padj[2]), ColKey = k),
      data.frame(x = 1, xend = 3, y = b_y2, label = format_p(padj[3]), ColKey = k))
  }
  km <- big_df %>% distinct(ColKey, Epitope, Cell_Line)
  b_df$ColKey <- factor(b_df$ColKey, levels = levels(big_df$ColKey))
  b_df <- dplyr::left_join(b_df, km, by = "ColKey")
  stat_df <- dplyr::left_join(stat_df, km, by = "ColKey")
  list(stat_df = stat_df, b_df = b_df)
}

sub_theme <- function() theme_minimal(base_size = 18) + theme(strip.text = element_blank(), strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), axis.title = element_text(size = 22, face = "bold"), axis.text  = element_text(size = 18, colour = "black"), panel.spacing = unit(0.9, "cm"))

# HARMONIZED ALIGNMENT: Visual layout matches x-axis mapping (ChIP top/left, CnT bottom/right)
label_check <- function(cf, pf) paste0("<span style='font-size:8pt;color:grey35'>ChIP: ", pf, "</span><br><span style='font-size:8pt;color:grey35'>CnT: ", cf, "</span>")

side_legend_plot <- function() {
  d1 <- data.frame(x = 0, y = 0.80, label = paste0("<span style='font-size:16pt'>**Enriched in**</span><br><br><span style='font-size:18pt;color:", CHIP_COL, "'>**ChIP-seq**</span><br><span style='font-size:14pt'>vs</span><br><span style='font-size:18pt;color:", CNT_COL, "'>**CUT&Tag**</span>"))
  d2 <- data.frame(x = 0.16, y = 0.34, label = paste0("<span style='font-size:16pt;color:", MEAN_COL, "'>**mean values**</span>"))
  ggplot() + theme_void() + xlim(0, 1) + ylim(0, 1) +
    geom_richtext(data = d1, aes(x, y, label = label), hjust = 0, vjust = 1, fill = NA, label.color = NA) +
    annotate("point", x = 0.06, y = 0.34, colour = MEAN_COL, size = 6) +
    geom_richtext(data = d2, aes(x, y, label = label), hjust = 0, vjust = 0.5, fill = NA, label.color = NA)
}

build_violin <- function(big_df, facet_layer, strip_element, metric_col, y_label) {
  max_val <- max(big_df[[metric_col]], na.rm = TRUE); min_val <- min(big_df[[metric_col]], na.rm = TRUE)
  range_val <- ifelse(max_val - min_val == 0, 1, max_val - min_val)
  mean_y <- max_val + (range_val * 0.08); b_y1 <- max_val + (range_val * 0.20); b_y2 <- max_val + (range_val * 0.32); y_top <- max_val + (range_val * 0.42)
  ann <- make_annotations(big_df, mean_y, b_y1, b_y2, metric_col)
  ggplot(big_df, aes(x = Status, y = .data[[metric_col]], fill = Status)) +
    geom_violin(alpha = 0.7, colour = "black", width = 0.85) + geom_boxplot(width = 0.18, colour = "black", outlier.shape = NA) +
    geom_segment(data = ann$b_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
    geom_segment(data = ann$b_df, aes(x = x, xend = x, y = y - (range_val * 0.03), yend = y), inherit.aes = FALSE) +
    geom_segment(data = ann$b_df, aes(x = xend, xend = xend, y = y - (range_val * 0.03), yend = y), inherit.aes = FALSE) +
    geom_text(data = ann$b_df, aes(x = (x + xend) / 2, y = y + (range_val * 0.02), label = label), inherit.aes = FALSE, size = 6.6, fontface = "bold", vjust = 0) +
    geom_text(data = dplyr::filter(ann$stat_df, n > 0), aes(x = Status, y = y_pos_mean, label = mean_label), inherit.aes = FALSE, size = 5.5, fontface = "bold.italic", colour = MEAN_COL) +
    scale_fill_manual(values = STATUS_COLORS, drop = FALSE, guide = "none") + scale_x_discrete(drop = FALSE) + facet_layer + coord_cartesian(ylim = c(min_val - (range_val * 0.05), y_top)) +
    labs(y = y_label, x = NULL) + theme_minimal(base_size = 18) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 22, colour = "black"), axis.title.y = element_text(size = 24, face = "bold"), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9), panel.spacing = unit(0.9, "cm"), strip.background = element_blank(), strip.text = strip_element)
}

color_spec <- function(mode, big_df) {
  if (toupper(mode) == "ATAC") list(var = "ATAC_Z", scale = scale_colour_viridis_c(name = "ATAC\n(Z-score)", option = "viridis", limits = c(-3, 3), oob = scales::squish))
  else if (toupper(mode) == "COMPOSITE") list(var = "Composite_Z", scale = scale_colour_viridis_c(name = "GC+ATAC\n(Z-score)", option = "plasma", limits = c(-3, 3), oob = scales::squish))
  else list(var = "GC_Percent", scale = scale_colour_viridis_c(name = "GC %", option = "magma", limits = c(30, 75), oob = scales::squish))
}

build_volcano <- function(big_df, facet_layer, color_var, color_scale) {
  km <- big_df %>% distinct(ColKey, Epitope, Cell_Line)
  ann <- big_df %>% group_by(ColKey) %>% summarise(n_chip = sum(Status == "ChIP-seq Enriched", na.rm = TRUE), n_cnt = sum(Status == "CUT&Tag Enriched",  na.rm = TRUE), .groups = "drop") %>% mutate(lbl_chip = paste0("n=", n_chip), lbl_cnt = paste0("n=", n_cnt)) %>% dplyr::left_join(km, by = "ColKey")
  ggplot(big_df, aes(x = log2FC, y = nLog10_FDR, colour = .data[[color_var]])) +
    geom_point_rast(alpha = 0.5, size = 0.5, raster.dpi = 300) + color_scale +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed", colour = "grey50") + geom_hline(yintercept = -log10(FDR_THRESH), linetype = "dashed", colour = "grey50") +
    geom_text(data = ann, aes(x =  Inf, y = Inf, label = lbl_cnt), hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold", colour = CNT_COL, inherit.aes = FALSE) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = lbl_chip),  hjust = -0.1, vjust = 1.5, size = 6, fontface = "bold", colour = CHIP_COL, inherit.aes = FALSE) +
    facet_layer + labs(x = "log2(CUT&Tag / ChIP-seq)", y = "-log10(FDR)") + sub_theme()
}

build_ma <- function(big_df, facet_layer, color_var, color_scale) {
  ggplot(big_df, aes(x = BaseMean_Log2, y = log2FC, colour = .data[[color_var]])) +
    geom_point_rast(alpha = 0.5, size = 0.5, raster.dpi = 300) + color_scale +
    geom_hline(yintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed", colour = "grey50") + geom_hline(yintercept = 0, colour = "black") +
    facet_layer + labs(x = "Avg Signal (log2)", y = "log2(CUT&Tag / ChIP-seq)") + sub_theme()
}
CAPTION <- "Significance:  ns (p \u2265 0.05)   * (p < 0.05)   ** (p < 0.01)   *** (p < 0.001)   **** (p < 0.0001)"
RCOMMON

cat > "$RS_DIR/main_body.R" <<'RMAIN'
args <- commandArgs(trailingOnly = TRUE); manifest <- args[1]; out_pdf <- args[2]; color_mode <- ifelse(length(args) >= 3, args[3], "GC")
big_df <- load_manifest(manifest)
if (nrow(big_df) == 0) q(save = "no")
z_score <- function(x) { if(sum(!is.na(x)) < 2) return(rep(0, length(x))); (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }
big_df <- big_df %>% group_by(Cell_Line) %>% mutate(GC_Z = z_score(GC_Percent), ATAC_Z = z_score(ATAC_log2), Raw_Composite = case_when(is.na(ATAC_Z) ~ GC_Z, TRUE ~ GC_Z + ATAC_Z), Composite_Z = z_score(Raw_Composite)) %>% ungroup()
if (toupper(color_mode) == "ATAC") { target_metric <- "ATAC_Z"; target_y_label <- "ATAC-seq (Z-score)" } else if (toupper(color_mode) == "COMPOSITE") { target_metric <- "Composite_Z"; target_y_label <- "ATAC (Z) + GC (Z)" } else { target_metric <- "GC_Percent"; target_y_label <- "GC Content (%)" }
lab_tbl <- big_df %>% distinct(ColKey, ColOrder, Cnt_Full, Chip_Full) %>% arrange(ColOrder)
col_lab <- setNames(mapply(label_check, lab_tbl$Cnt_Full, lab_tbl$Chip_Full), lab_tbl$ColKey)
big_df$Status <- factor(big_df$Status, levels = STATUS_LEVELS); big_df$ColKey <- factor(big_df$ColKey, levels = lab_tbl$ColKey); big_df$Epitope <- factor(big_df$Epitope, levels = unique(big_df$Epitope[order(big_df$ColOrder)])); big_df$Cell_Line <- factor(big_df$Cell_Line, levels = unique(big_df$Cell_Line))
cs <- color_spec(color_mode, big_df)
viol <- build_violin(big_df, facet_wrap(~ ColKey, nrow = 1, labeller = labeller(ColKey = col_lab)), element_markdown(size = 8, lineheight = 1.15, halign = 0.5, margin = margin(b = 6)), target_metric, target_y_label)
volc <- build_volcano(big_df, facet_wrap(~ ColKey, nrow = 1, scales = "free"), cs$var, cs$scale)
body <- viol / volc + plot_layout(heights = c(0.75, 1))
final <- wrap_plots(body, side_legend_plot(), ncol = 2, widths = c(1, 0.20), guides = "collect") + plot_annotation(subtitle = sprintf("K562 CUT&Tag vs ChIP-seq over ChromHMM ground-truth regions  |  FDR < %s, |log2FC| > %s", FDR_THRESH, FC_THRESH), caption = CAPTION, theme = theme(plot.subtitle = element_text(size = 14, face = "italic"), plot.caption = element_text(size = 16, face = "italic", hjust = 1, colour = "grey30"))) & theme(legend.position = "right", legend.text = element_text(size = 17, face = "bold"), legend.title = element_text(size = 19, face = "bold"), legend.key.height = unit(1.8, "cm"))
ggsave(out_pdf, final, width = 6 + nrow(lab_tbl) * 3.6, height = 11.5, device = "pdf", bg = "white", limitsize = FALSE)
RMAIN

cat > "$RS_DIR/supp_body.R" <<'RSUPP'
args <- commandArgs(trailingOnly = TRUE); manifest <- args[1]; out_pdf <- args[2]; color_mode <- ifelse(length(args) >= 3, args[3], "GC")
library(ggh4x)
big_df <- load_manifest(manifest)
if (nrow(big_df) == 0) q(save = "no")
z_score <- function(x) { if(sum(!is.na(x)) < 2) return(rep(0, length(x))); (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }
big_df <- big_df %>% group_by(Cell_Line) %>% mutate(GC_Z = z_score(GC_Percent), ATAC_Z = z_score(ATAC_log2), Raw_Composite = case_when(is.na(ATAC_Z) ~ GC_Z, TRUE ~ GC_Z + ATAC_Z), Composite_Z = z_score(Raw_Composite)) %>% ungroup()
if (toupper(color_mode) == "ATAC") { target_metric <- "ATAC_Z"; target_y_label <- "ATAC-seq (Z-score)" } else if (toupper(color_mode) == "COMPOSITE") { target_metric <- "Composite_Z"; target_y_label <- "ATAC (Z) + GC (Z)" } else { target_metric <- "GC_Percent"; target_y_label <- "GC Content (%)" }
lab_tbl <- big_df %>% distinct(ColKey, ColOrder, Epitope, Cell_Line, Cnt_Full, Chip_Full) %>% arrange(ColOrder)
col_lab <- setNames(mapply(label_check, lab_tbl$Cnt_Full, lab_tbl$Chip_Full), lab_tbl$ColKey)
epi_levels <- unique(lab_tbl$Epitope); ep_lab <- setNames(paste0("<span style='font-size:19pt'>**", epi_levels, "**</span>"), epi_levels); cl_lab <- c(K562 = "<span style='font-size:14pt'>**K562**</span>", MCF7 = "<span style='font-size:14pt'>**MCF7**</span>")
big_df$Status <- factor(big_df$Status, levels = STATUS_LEVELS); big_df$Epitope <- factor(big_df$Epitope, levels = epi_levels); big_df$Cell_Line <- factor(big_df$Cell_Line, levels = c("K562", "MCF7")); big_df$ColKey <- factor(big_df$ColKey, levels = lab_tbl$ColKey)
viol_facet <- facet_nested_wrap(vars(Epitope, Cell_Line, ColKey), nrow = 1, nest_line = element_line(colour = "grey15", linewidth = 0.7), solo_line = TRUE, labeller = labeller(Epitope = ep_lab, Cell_Line = cl_lab, ColKey = col_lab), strip = strip_nested(background_x = elem_list_rect(fill = NA, colour = NA)))
sub_facet <- facet_nested_wrap(vars(Epitope, Cell_Line, ColKey), nrow = 1, scales = "free", nest_line = element_blank())
viol <- build_violin(big_df, viol_facet, element_markdown(size = 8, lineheight = 1.12, halign = 0.5, margin = margin(b = 4)), target_metric, target_y_label)
cs <- color_spec(color_mode, big_df)
volc <- build_volcano(big_df, sub_facet, cs$var, cs$scale); ma <- build_ma(big_df, sub_facet, cs$var, cs$scale)
body <- viol / volc / ma + plot_layout(heights = c(0.85, 1, 1))
final <- wrap_plots(body, side_legend_plot(), ncol = 2, widths = c(1, 0.12), guides = "collect") + plot_annotation(subtitle = sprintf("All CUT&Tag vs ChIP-seq comparisons over ChromHMM ground-truth regions  |  FDR < %s, |log2FC| > %s", FDR_THRESH, FC_THRESH), caption = CAPTION, theme = theme(plot.subtitle = element_text(size = 16, face = "italic"), plot.caption = element_text(size = 17, face = "italic", hjust = 1, colour = "grey30"))) & theme(legend.position = "right", legend.text = element_text(size = 17, face = "bold"), legend.title = element_text(size = 19, face = "bold"), legend.key.height = unit(1.8, "cm"))
ggsave(out_pdf, final, width = 6 + nrow(lab_tbl) * 3.4, height = 19, device = "pdf", bg = "white", limitsize = FALSE)
RSUPP


# =========================================================================
# PART 4: COMPARISON MANIFESTS & EXECUTION
# =========================================================================
echo "========================================================="
echo " Running Comparisons (NATIVE DESeq2 Strategy)"
echo "========================================================="

MAIN_MANIFEST="$DIFF_DIR/manifest_main.tsv"
SUPP_MANIFEST="$DIFF_DIR/manifest_supp.tsv"
printf 'epitope\tcell_line\tcnt_label\tchip_label\tdiff_path\tgc_path\tatac_path\n' > "$MAIN_MANIFEST"
printf 'epitope\tcell_line\tcnt_label\tchip_label\tdiff_path\tgc_path\tatac_path\n' > "$SUPP_MANIFEST"

CHIP_K562_H3K4="chip_H3K4me3_SRR5339104_K562_ID303_BER1 chip_H3K4me3_SRR5339105_K562_ID303_BER2"
CHIP_MCF7_H3K4="chip_H3K4me3_SRR5339666-SRR5339667-SRR5339668_MCF7_ID506_BER1 chip_H3K4me3_SRR5339663-SRR5339664-SRR5339665_MCF7_ID506_BER2"
CHIP_K562_H3K27ac="chip_H3K27ac_SRR10319903_K562_ID190_TAK1 chip_H3K27ac_SRR10319904_K562_ID190_TAK2"
CHIP_MCF7_H3K27ac="chip_H3K27ac_SRR5339250_MCF7_ID352_BER1 chip_H3K27ac_SRR5339251_MCF7_ID352_BER2"
CHIP_K562_H3K27me3="chip_H3K27me3_SRR227389_K562_ID611_BER1 chip_H3K27me3_ENCFF936TRH-ENCFF699LBD-ENCFF014OTE-ENCFF465KWY_K562_IDAKQ_BER3"
CHIP_MCF7_H3K27me3="chip_H3K27me3_SRR5339277_MCF7_ID363_BER1 chip_H3K27me3_SRR5339274-SRR5339275-SRR5339276_MCF7_ID363_BER2"
CHIP_K562_H3K36me3="chip_H3K36me3_SRR227512_K562_ID611_BER2 chip_H3K36me3_ENCFF195ATD-ENCFF237PVB-ENCFF627AZR-ENCFF544IFT_K562_IDAKR_BER3"
CHIP_MCF7_H3K36me3="chip_H3K36me3_SRR14636642-SRR14636643-SRR14636644_MCF7_ID945_BER1 chip_H3K36me3_SRR14636641-SRR14636640_MCF7_ID945_BER2"

CNT_K562_H3K4_KAY1="CnT_H3K4me3_SRR8383516_K562_ID557_KAY1 CnT_H3K4me3_SRR8383517_K562_ID557_KAY2"
CNT_K562_H3K4_KAY2="CnT_H3K4me3_SRR11074249_K562_ID187_KAY1 CnT_H3K4me3_SRR11074250_K562_ID187_KAY2"
CNT_MCF7_H3K4_TIA1="CnT_H3K4me3_SRR18862763_MCF7_ID262_TIA1 CnT_H3K4me3_SRR18862762_MCF7_ID262_TIA2"
CNT_K562_H3K27ac_KAY1="CnT_H3K27ac_SRR8383507_K562_ID557_KAY1 CnT_H3K27ac_SRR8383508_K562_ID557_KAY2"
CNT_K562_H3K27ac_ABB1="CnT_H3K27ac_SRR31972743_K562_ID492_A471 CnT_H3K27ac_SRR31972742_K562_ID492_A472"
CNT_K562_H3K27ac_ABB2="CnT_H3K27ac_SRR31972753_K562_ID492_A171 CnT_H3K27ac_SRR31972744_K562_ID492_A172"
CNT_K562_H3K27ac_ABB3="CnT_H3K27ac_SRR31972741_K562_ID492_D191 CnT_H3K27ac_SRR31972740_K562_ID492_D192"
CNT_MCF7_H3K27ac_TIA1="CnT_H3K27ac_SRR18862771_MCF7_ID262_TIA1 CnT_H3K27ac_SRR18862770_MCF7_ID262_TIA2"
CNT_MCF7_H3K27ac_FIS1="CnT_H3K27ac_SRR30799704_MCF7_ID033_FIS1 CnT_H3K27ac_SRR30799703_MCF7_ID033_FIS2"
CNT_K562_H3K27me3_KAY1="CnT_H3K27me3_SRR8383510_K562_ID557_KAY1 CnT_H3K27me3_SRR8383511_K562_ID557_KAY2"
CNT_K562_H3K27me3_ABB1="CnT_H3K27me3_SRR31972739_K562_ID492_ABB1 CnT_H3K27me3_SRR31972738_K562_ID492_ABB2"
CNT_MCF7_H3K27me3_TIA1="CnT_H3K27me3_SRR18862774_MCF7_ID262_TIA1 CnT_H3K27me3_SRR18862773_MCF7_ID262_TIA2 CnT_H3K27me3_SRR24561243_MCF7_ID262_TIA3 CnT_H3K27me3_SRR24561242_MCF7_ID262_TIA4"
CNT_K562_H3K36me3_WUW1="CnT_H3K36me3_SRR29493117_K562_ID327_WUW1 CnT_H3K36me3_SRR29493116_K562_ID327_WUW2 CnT_H3K36me3_SRR29493115_K562_ID327_WUW3 CnT_H3K36me3_SRR29493114_K562_ID327_WUW4 CnT_H3K36me3_SRR29493113_K562_ID327_WUW5"
CNT_MCF7_H3K36me3_TIA1="CnT_H3K36me3_SRR18862751_MCF7_ID262_TIA1 CnT_H3K36me3_SRR18862750_MCF7_ID262_TIA2"

# ------------------- SUPPLEMENTARY ---------------------
add_comparison "$SUPP_MANIFEST" H3K4me3  K562 vs_KAY1 "$CHIP_K562_H3K4"      "$CNT_K562_H3K4_KAY1"
add_comparison "$SUPP_MANIFEST" H3K4me3  K562 vs_KAY2 "$CHIP_K562_H3K4"      "$CNT_K562_H3K4_KAY2"
add_comparison "$SUPP_MANIFEST" H3K4me3  MCF7 vs_TIA1 "$CHIP_MCF7_H3K4"      "$CNT_MCF7_H3K4_TIA1"
add_comparison "$SUPP_MANIFEST" H3K27ac  K562 vs_KAY1 "$CHIP_K562_H3K27ac"  "$CNT_K562_H3K27ac_KAY1"
add_comparison "$SUPP_MANIFEST" H3K27ac  K562 vs_ABB1 "$CHIP_K562_H3K27ac"  "$CNT_K562_H3K27ac_ABB1"
add_comparison "$SUPP_MANIFEST" H3K27ac  K562 vs_ABB2 "$CHIP_K562_H3K27ac"  "$CNT_K562_H3K27ac_ABB2"
add_comparison "$SUPP_MANIFEST" H3K27ac  K562 vs_ABB3 "$CHIP_K562_H3K27ac"  "$CNT_K562_H3K27ac_ABB3"
add_comparison "$SUPP_MANIFEST" H3K27ac  MCF7 vs_TIA1 "$CHIP_MCF7_H3K27ac"  "$CNT_MCF7_H3K27ac_TIA1"
add_comparison "$SUPP_MANIFEST" H3K27ac  MCF7 vs_FIS1 "$CHIP_MCF7_H3K27ac"  "$CNT_MCF7_H3K27ac_FIS1"
add_comparison "$SUPP_MANIFEST" H3K27me3 K562 vs_KAY1 "$CHIP_K562_H3K27me3" "$CNT_K562_H3K27me3_KAY1"
add_comparison "$SUPP_MANIFEST" H3K27me3 K562 vs_ABB1 "$CHIP_K562_H3K27me3" "$CNT_K562_H3K27me3_ABB1"
add_comparison "$SUPP_MANIFEST" H3K27me3 MCF7 vs_TIA1 "$CHIP_MCF7_H3K27me3" "$CNT_MCF7_H3K27me3_TIA1"
add_comparison "$SUPP_MANIFEST" H3K36me3 K562 vs_WUW1 "$CHIP_K562_H3K36me3" "$CNT_K562_H3K36me3_WUW1"
add_comparison "$SUPP_MANIFEST" H3K36me3 MCF7 vs_TIA1 "$CHIP_MCF7_H3K36me3" "$CNT_MCF7_H3K36me3_TIA1"

# ------------------- MAIN ---------------------
# Swapped vs_KAY1 out for vs_KAY2, and swapped vs_KAY1 out for vs_ABB1
add_comparison "$MAIN_MANIFEST" H3K4me3  K562 vs_KAY2 "$CHIP_K562_H3K4"      "$CNT_K562_H3K4_KAY2"
add_comparison "$MAIN_MANIFEST" H3K27ac  K562 vs_ABB2 "$CHIP_K562_H3K27ac"  "$CNT_K562_H3K27ac_ABB2"
add_comparison "$MAIN_MANIFEST" H3K27me3 K562 vs_ABB1 "$CHIP_K562_H3K27me3" "$CNT_K562_H3K27me3_ABB1"
add_comparison "$MAIN_MANIFEST" H3K36me3 K562 vs_WUW1 "$CHIP_K562_H3K36me3" "$CNT_K562_H3K36me3_WUW1"

# =========================================================================
# PART 5: RENDER FIGURES
# =========================================================================
MAIN_GC_PDF="$DIFF_DIR/MAIN_K562_perEpitope_violin_volcano_GC_raster.pdf"
MAIN_ATAC_PDF="$DIFF_DIR/MAIN_K562_perEpitope_violin_volcano_ATAC_raster.pdf"
MAIN_COMP_PDF="$DIFF_DIR/MAIN_K562_perEpitope_violin_volcano_COMPOSITE_raster.pdf"

SUPP_GC_PDF="$DIFF_DIR/SUPP_allComparisons_violin_volcano_MA_GC_raster.pdf"
SUPP_ATAC_PDF="$DIFF_DIR/SUPP_allComparisons_violin_volcano_MA_ATAC_raster.pdf"
SUPP_COMP_PDF="$DIFF_DIR/SUPP_allComparisons_violin_volcano_MA_COMPOSITE_raster.pdf"

cat "$RS_DIR/common.R" "$RS_DIR/main_body.R" > "$RS_DIR/run_main.R"
cat "$RS_DIR/common.R" "$RS_DIR/supp_body.R" > "$RS_DIR/run_supp.R"

echo "=== Rendering (GC) ==="
Rscript "$RS_DIR/run_main.R" "$MAIN_MANIFEST" "$MAIN_GC_PDF" "GC"
Rscript "$RS_DIR/run_supp.R" "$SUPP_MANIFEST" "$SUPP_GC_PDF" "GC"

echo "=== Rendering (ATAC) ==="
Rscript "$RS_DIR/run_main.R" "$MAIN_MANIFEST" "$MAIN_ATAC_PDF" "ATAC"
Rscript "$RS_DIR/run_supp.R" "$SUPP_MANIFEST" "$SUPP_ATAC_PDF" "ATAC"

echo "=== Rendering (COMPOSITE) ==="
Rscript "$RS_DIR/run_main.R" "$MAIN_MANIFEST" "$MAIN_COMP_PDF" "COMPOSITE"
Rscript "$RS_DIR/run_supp.R" "$SUPP_MANIFEST" "$SUPP_COMP_PDF" "COMPOSITE"

# =========================================================================
# PART 6: FEATURE IMPORTANCE (VARIANCE EXPLAINED)
# =========================================================================
echo "========================================================="
echo "📊 Running Feature Importance (R-Squared) Analysis..."
echo "========================================================="

MAIN_R2_PDF="$PLOTS_DIR/MAIN_Variance_Explained_R2_BarChart.pdf"
SUPP_R2_PDF="$PLOTS_DIR/SUPP_Variance_Explained_R2_BarChart.pdf"

cat > "$RS_DIR/run_variance_explained.R" <<'EOF'
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(ggh4x); library(ggtext)
})
options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
out_pdf <- args[2]

prep_gc <- function(file) {
  df <- suppressMessages(read_tsv(file, show_col_types = FALSE))
  colnames(df) <- c("Chr", "Start", "End", "GC_Raw")
  df %>% mutate(GC_Percent = ifelse(GC_Raw <= 1.0, GC_Raw * 100, GC_Raw))
}

prep_atac <- function(file) {
  if (!file.exists(file)) return(data.frame())
  df <- suppressMessages(read_tsv(file, show_col_types = FALSE))
  colnames(df) <- c("Chr", "Start", "End", "ATAC_Raw")
  df %>% mutate(ATAC_log2 = log2(pmax(ATAC_Raw, 0) + 1))
}

z_score <- function(x) {
    if(sum(!is.na(x)) < 2) return(rep(NA, length(x)))
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

manifest <- suppressMessages(read_tsv(manifest_path, show_col_types = FALSE))
big_df <- data.frame()

for (i in seq_len(nrow(manifest))) {
  r <- manifest[i, ]
  if (!file.exists(r$diff_path)) next
  
  diff <- suppressMessages(read_tsv(r$diff_path, show_col_types = FALSE, name_repair = "unique", comment = ""))
  colnames(diff)[1:4] <- c("PeakID", "Chr", "Start", "End")
  diff$Start <- diff$Start - 1
  
  fc_col <- grep("Log2 Fold Change", colnames(diff), ignore.case = TRUE)[1]
  colnames(diff)[fc_col] <- "log2FC"
  
  gc <- prep_gc(r$gc_path)
  df <- inner_join(diff[, c("Chr", "Start", "End", "log2FC")], gc, by = c("Chr", "Start", "End"))
  
  atac <- prep_atac(r$atac_path)
  if (nrow(atac) > 0) {
    df <- left_join(df, atac[, c("Chr", "Start", "End", "ATAC_log2")], by = c("Chr", "Start", "End"))
  } else {
    df$ATAC_log2 <- NA_real_
  }
  
  df$Epitope <- r$epitope
  df$Cell_Line <- r$cell_line
  df$Cnt_Label <- r$cnt_label
  
  big_df <- bind_rows(big_df, df)
}

stats_df <- big_df %>%
  filter(!is.na(log2FC) & !is.na(GC_Percent)) %>%
  group_by(Cell_Line) %>%
  mutate(
    GC_Z = z_score(GC_Percent),
    ATAC_Z = z_score(ATAC_log2),
    Raw_Comp = ifelse(is.na(ATAC_Z), GC_Z, GC_Z + ATAC_Z),
    Comp_Z = z_score(Raw_Comp)
  ) %>%
  group_by(Epitope, Cell_Line, Cnt_Label) %>%
  summarize(
    R2_GC = cor(log2FC, GC_Z, use="complete.obs")^2,
    R2_ATAC = ifelse(all(is.na(ATAC_Z)), NA, cor(log2FC, ATAC_Z, use="complete.obs")^2),
    R2_Composite = ifelse(all(is.na(Comp_Z)), NA, cor(log2FC, Comp_Z, use="complete.obs")^2),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("R2_"), names_to = "Feature", values_to = "R_Squared") %>%
  mutate(
    Feature = factor(Feature, 
                     levels = c("R2_GC", "R2_ATAC", "R2_Composite"), 
                     labels = c("GC Content", "ATAC-seq", "Composite (GC+ATAC)")),
    Epitope = factor(Epitope, levels = c("H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3"))
  )

stats_df$Variance_Percent <- stats_df$R_Squared * 100

epi_levels <- unique(stats_df$Epitope)
ep_lab <- setNames(paste0("<span style='font-size:16pt'>**", epi_levels, "**</span>"), epi_levels)

p <- ggplot(stats_df %>% filter(!is.na(Variance_Percent)), 
       aes(x = Cnt_Label, y = Variance_Percent, fill = Feature)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color="black", linewidth=0.3, width = 0.7) +
  scale_fill_manual(values = c("GC Content" = "#51127c", "ATAC-seq" = "#21918c", "Composite (GC+ATAC)" = "#fc8961")) +
  facet_nested(. ~ Epitope + Cell_Line, scales = "free_x", space = "free_x",
               nest_line = element_line(colour = "grey15", linewidth = 0.7),
               labeller = labeller(Epitope = ep_lab)) +
  labs(
    title = "Variance in Log2FC (CUT&Tag vs ChIP-seq) Explained by Sequence/Chromatin Features",
    subtitle = "Higher percentage means the feature is a stronger predictor of the technique bias",
    y = "Variance Explained (%)",
    x = "CUT&Tag Sample Comparison"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black", size=10),
    axis.text.y = element_text(color="black", size=14),
    axis.title = element_text(face = "bold", size=16),
    plot.title = element_text(face = "bold", size=18),
    plot.subtitle = element_text(face = "italic", size=14, color="grey30"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_blank(),
    strip.text.x = element_markdown(color="black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=14, face="bold")
  )

num_samples <- length(unique(paste(stats_df$Epitope, stats_df$Cell_Line, stats_df$Cnt_Label)))
out_width <- max(12, num_samples * 0.8)

ggsave(out_pdf, plot = p, width = out_width, height = 8, device = "pdf", bg = "white")
EOF

echo " -> Rendering Main R2 Plot"
Rscript "$RS_DIR/run_variance_explained.R" "$MAIN_MANIFEST" "$MAIN_R2_PDF"

echo " -> Rendering Supp R2 Plot"
Rscript "$RS_DIR/run_variance_explained.R" "$SUPP_MANIFEST" "$SUPP_R2_PDF"

echo "========================================================="
echo "🎉 Done. Pipeline completed successfully."
echo " Output Directories:"
echo "   - Tracks:      $TRACKS_DIR"
echo "   - Diff:        $DIFF_DIR"
echo "   - R2 Analysis: $PLOTS_DIR"
echo "========================================================="