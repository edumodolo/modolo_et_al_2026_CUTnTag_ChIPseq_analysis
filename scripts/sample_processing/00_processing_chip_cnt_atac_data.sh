#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ==============================================================================
# This is the unified chromatin profiling pipeline -- compatible with ATAC-seq, ChIP-seq, CUT&Tag, CUT&RUN
#
# Usage:
#   conda activate <environment with fastp bowtie2 samtools deeptools bedtools
#                   sra-tools pigz homer R and a couple UCSC tools... > 
#   
# I provide a yml file called "chrom_prof_processing_v2.yml" that can be used to generate the environment necessary
# for running this pipeline, to run the pipeline, call "bash unified_chromatin_pipeline.sh" or "./unified_chromatin_pipeline.sh"
#
# The input for this pipeline is one sample sheet (MAPPING_FILE) labeled "rename_files_plot_labels.tab", specify the path to the directory 
# containing this sample sheet in the INPUT_DIR USER CONFIG section. In Supplementary Table 1, the grey highlight identifies the sample sheet 
# used in this manuscript's analysis. 
#
# The assay/method is read from the first field of the "sample name", and every assay-specific parameter (MAPQ, aligner mode, fragment
# cap, track color, replicate merging) is looked up from that assay tag.
# Steps that do not need to differ between assays are shared.
# 
# ─── SAMPLE_NAME and PLOT_LABEL NAMING CONVENTION ──────────────────────────────────────────────────────────────
#
# Naming convention for the "sample_name", which is used for processing of all biological and technical replicates looks like this: 
# "MethodID_Epitope_RunID_CellType_ExperimentID_BiologicalReplicateID"
#
# "sample_name" examples include:
# ATAC-seq sample:
# ATAC_NA_SRR5809234_K562_ID173_LIU1
# CUT&Tag sample:
# CnT_H3K27ac_SRR31972743_K562_ID492_A471
# ChIP-seq sample:
# chip_H3K4me3_SRR5339104_K562_ID303_BER1
#
# A more aethetic name for plotting/selecting samples in downstream analysis is the "plot_label" 
#
# Naming convention looks like this: "MethodID_Epitope_CellType_AuthorLastName_PublicationYear_ReplicateID"
# Important Note: technical replicates with separate "sample_names" get merged under one "plot_label"
#
# "plot_label" examples include:
# ATAC-seq sample:
# ATAC_NA_K562_Liu_2017_rep1
# CUT&Tag sample:
# CnT_H3K27ac_K562_Abbasova_2025_ab4729r1
# ChIP-seq sample:
# chip_H3K4me3_K562_Bernstein_2017_rep1 
#
#
# NOTE about library sequencing structure: Single-end and paired-end libraries are handled in the same pipeline run using autodetection
# (detected by how the fastq files are named after download (i.e. contains "_2" means its a paired end file) and this information is 
# stored and used throughout the whole script. 
#
# General flowthrough of this CUT&Tag/ChIP-seq/ATAC-seq processing pipeline looks like:
#   download -> trim -> align -> filter -> deduplicate -> chrM/blacklist ->
#   bigWigs -> HOMER tag directories + QC -> processing metrics + UCSC tracks
#
# NOTE about how the script resumes on retry: each step (e.g., fastq download, alignment, filtering), writes a completion marker containing a 
# fingerprint of the library list it covered. Adding samples to the sheet invalidates the markers as the fingerprint will be different, 
# and re-running the script will process only the newly added samples, i.e., work already finished is skipped.
#
# After running this pipeline, you will get a Master_Sample_Metrics.tsv file in the summary directory which provides useful metrics such as total fragments 
# left in a sample after filtering, duplication rate, and locations of important files (bigwig, bam) that will be used in downstream analysis. This pipeline
# also automatically produces bigwigs, and accompanying UCSC track upload text, for easy inspection of processed samples on the UCSC Genome Browser.    
#
# ==============================================================================

# ─── USER CONFIG ──────────────────────────────────────────────────────────────
INPUT_DIR="${INPUT_DIR:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/input}"
BASE_DIR="${BASE_DIR:-/gpfs/data01/gorenlab/emodolo/2026_modolo_et_al/ChIP_CnT_ATAC_datasets/output_full_08-19}"

# Sample sheet: 3 tab-separated columns, no header.
#   run_id <TAB> sample_name <TAB> plot_label
# sample_name is six underscore-separated fields:
#   MethodID_Epitope_RunID_CellType_ExperimentID_BiologicalReplicateID
#   e.g. CnT_H3K4me3_SRR8383516_K562_ID557_KAY1
# Rows that differ ONLY in RunID are considered technical replicates and are concatenated.
# Rows that differ in BiologicalReplicateID are biological replicates.
MAPPING_FILE="${MAPPING_FILE:-${INPUT_DIR}/rename_files_plot_labels.tab}"

GENOME="hg38"
GENOME_INDEX="${GENOME_INDEX:-/gpfs/data01/gorenlab/emodolo/genomes/genome_indexes/hg38_bowtie2_index/hg38_bwt2_index}"
BLACKLIST="${BLACKLIST:-/home/emodolo/gpfs/genomes/bedfiles/hg38-blacklist.v2.bed}"

THREADS="${THREADS:-16}"
FREQ_ZOOM_BP="${FREQ_ZOOM_BP:-50}"
# This section allows you to set the "zoomed-in" nucleotide-frequency QC panel created by plotting results from HOMER's checkGC, 
# set the number of bp to extend from either side of the read 5' end. The zoomed-out panel always shows HOMER's full tagFreq range.

BIGWIG_BIN_SIZE="${BIGWIG_BIN_SIZE:-1}"
MERGED_BIN_SIZE="${MERGED_BIN_SIZE:-1}"

WEB_DIR_MASTER="${WEB_DIR_MASTER:-/homer_data/www/html/emodolo/modolo_et_al_2026/full_run_08-19}"
WEB_URL_MASTER="${WEB_URL_MASTER:-http://homer.ucsd.edu/emodolo/modolo_et_al_2026/full_run_08-19}"

# ==============================================================================
# ASSAY-SPECIFIC PARAMETERS
#
# Everything that must differ between assays is declared here, keyed by assay/method,
# so no processing step contains a hard-coded assay name. 
# ==============================================================================

# --- MAPQ threshold ----------------------------------------------------------
declare -A ASSAY_MAPQ=( [ATAC]=30 [CHIP]=20 [CNT]=20 [CNR]=20 )

# --- Maximum fragment length (bowtie2 -X, used for paired-end libraries only) ----------
declare -A ASSAY_MAXFRAG=( [ATAC]=2000 [CHIP]=700 [CNT]=700 [CNR]=700 )

# --- bowtie2 alignment mode --------------------------------------------------
declare -A ASSAY_BT2_PE=(
    [ATAC]="--mm"
    [CHIP]="--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10"
    [CNT]="--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10"
    [CNR]="--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10"
)
declare -A ASSAY_BT2_SE=(
    [ATAC]="--mm"
    [CHIP]="--local --very-sensitive --phred33"
    [CNT]="--local --very-sensitive --phred33"
    [CNR]="--local --very-sensitive --phred33"
)

# --- fastp paired-end adapter detection --------------------------------------
declare -A ASSAY_DETECT_ADAPTER_PE=( [ATAC]=1 [CHIP]=1 [CNT]=1 [CNR]=1 )

# --- Single-end fragment extension for bigWigs -------------------------------
declare -A ASSAY_SE_EXTEND=( [ATAC]=200 [CHIP]=200 [CNT]=200 [CNR]=200 )

# --- Sequencing format requirements ( 0 = can be either paired or single end, 1 = throw an error if single end) ----------------------
declare -A ASSAY_REQUIRE_PE=( [ATAC]=0 [CHIP]=0 [CNT]=1 [CNR]=1 )

# --- Merging biological-replicate tracks into a single bam file and CPM normalized bigwig -------------------------
# ATAC replicates are merged into one CPM-normalised track, because ATAC is used
# here as a reference dataset of open chromatin rather than as the quantity
# being tested. ChIP/CUT&Tag replicates are kept separate: they are the
# comparison datasets, and merging would hide replicate-level variability.
declare -A ASSAY_MERGE_REPS=( [ATAC]=1 [CHIP]=0 [CNT]=0 [CNR]=0 )

# --- UCSC Genome Browser track colors ------------------------------------------------------
declare -A ASSAY_COLOR=(
    [ATAC]="237,0,230"   # pink
    [CHIP]="57,54,255"   # blue
    [CNT]="255,59,59"    # red
    [CNR]="0,201,64"     # green
)

# ─── MAKING DIRECTORIES & SETTING COMPLETION MARKER FILES ────────────────────────────────────────────────────
FASTQ_DIR="${BASE_DIR}/fastq_files"
CLEAN_FASTQ_DIR="${BASE_DIR}/clean_fastq_files"
ALIGNMENT_DIR="${BASE_DIR}/alignment"
BAM_DIR="${ALIGNMENT_DIR}/bam"
ALIGNED_BAM_DIR="${BAM_DIR}/aligned"
FILTERED_BAM_DIR="${BAM_DIR}/filtered"
DEDUP_BAM_DIR="${BAM_DIR}/dedup"
FINAL_BAM_DIR="${BAM_DIR}/final"
TMP_DIR="${ALIGNMENT_DIR}/tmp"
BIGWIG_DIR="${ALIGNMENT_DIR}/bigwig"
HOMER_QC_DIR="${BASE_DIR}/homer_qc"
SUMMARY_DIR="${BASE_DIR}/summary"
LOG_DIR="${SUMMARY_DIR}/logs"

mkdir -p "$FASTQ_DIR" "$CLEAN_FASTQ_DIR" "$ALIGNED_BAM_DIR" "$FILTERED_BAM_DIR" \
         "$DEDUP_BAM_DIR" "$FINAL_BAM_DIR" "$TMP_DIR" "$BIGWIG_DIR" \
         "$HOMER_QC_DIR" "$SUMMARY_DIR" "$LOG_DIR"
mkdir -p "$WEB_DIR_MASTER" 2>/dev/null || true

DOWNLOAD_COMPLETE="${FASTQ_DIR}/download_complete.txt"
CLEAN_COMPLETE="${CLEAN_FASTQ_DIR}/clean_complete.txt"
ALIGN_COMPLETE="${ALIGNED_BAM_DIR}/align_complete.txt"
FILTER_COMPLETE="${FILTERED_BAM_DIR}/filter_complete.txt"
DUP_COMPLETE="${DEDUP_BAM_DIR}/dup_complete.txt"
FINALBAM_COMPLETE="${FINAL_BAM_DIR}/finalbam_complete.txt"
BIGWIG_COMPLETE="${BIGWIG_DIR}/bigwig_complete.txt"
HOMER_QC_COMPLETE="${HOMER_QC_DIR}/qc_complete.txt"

MASTER_METRICS="${SUMMARY_DIR}/Master_Sample_Metrics.tsv"
MERGED_METRICS="${SUMMARY_DIR}/Merged_Replicate_Bigwigs.tsv"
UCSC_TRACKS="${SUMMARY_DIR}/UCSC_tracks.txt"

echo "════════════════════════════════════════════════════════════"
echo " Unified chromatin pipeline  (ATAC / ChIP / CUT&Tag / CUT&RUN)"
echo " Genome : ${GENOME}"
echo " Sheet  : ${MAPPING_FILE}"
echo " Output : ${BASE_DIR}"
echo "════════════════════════════════════════════════════════════"

# ─── PREFLIGHT ────────────────────────────────────────────────────────────────
echo "🔧 Preflight..."
missing=0
for t in fastp bowtie2 samtools bamCoverage bedtools pigz md5sum awk sort Rscript; do
    command -v "$t" >/dev/null 2>&1 || { echo "  ❌ missing: $t"; missing=1; }
done
for t in prefetch fasterq-dump; do
    command -v "$t" >/dev/null 2>&1 || echo "  ⚠️  missing (needed for SRR downloads): $t"
done
command -v makeTagDirectory >/dev/null 2>&1 || echo "  ⚠️  missing (needed for HOMER QC): makeTagDirectory"
[[ -f "${GENOME_INDEX}.1.bt2" ]] || { echo "  ❌ bowtie2 index not found: ${GENOME_INDEX}.1.bt2"; missing=1; }
[[ -f "$BLACKLIST" ]]           || { echo "  ❌ blacklist not found: $BLACKLIST"; missing=1; }
[[ -f "$MAPPING_FILE" ]]        || { echo "  ❌ sample sheet not found: $MAPPING_FILE"; missing=1; }
[[ $missing -eq 0 ]] || { echo "Preflight failed."; exit 1; }
echo "  ✅ required tools present"

# ─── HELPER FUNCTIONS ──────────────────────────────────────────────────────────────────
bam_mapped_reads() { samtools view -@ "$THREADS" -c -F 0x4 "$1" 2>/dev/null || echo "NA"; }
# counts mapped reads from a bam file ( -F 0x4 filters out reads where the SAM bit flag 0x4 is set, indicating that it is an unmapped read)
safe_div() { awk -v a="$1" -v b="$2" 'BEGIN{ if (b+0>0) printf "%.4f", a/b; else printf "NA" }'; }
# divides a/b safely to four decimal places, protecting against division-by-zero errors

# the command below Normalises the assay/method tag written in the sample sheet (i.e., the first field of the sample_name in the MAPPING_FILE) 
# to the key used by the parameter tables above, so ATAC/atac and CnT/cnt/CUTTAG all resolve correctly.
assay_key() {
    local a="${1^^}"
    case "$a" in
        ATAC|ATACSEQ)      echo "ATAC" ;;
        CHIP|CHIPSEQ)      echo "CHIP" ;;
        CNT|CUTTAG|CUTNTAG) echo "CNT" ;;
        CNR|CUTRUN|CUTNRUN) echo "CNR" ;;
        *)                 echo "UNKNOWN" ;;
    esac
}

# Helper commands to make and check completion markers, which store a fingerprint of the library list they covered. 
grps_fingerprint() { printf '%s\n' "${ALL_GROUPS[@]}" | sort | md5sum | cut -d' ' -f1; }
# ^ this command prints each element from the ALL_GROUPS array, sorts it alphabetically, then calculates a 32 character
# MD5 checksum of that sorted list, cutting off the trailing whitespace and dash, so that we only get the MD5 hash string
marker_valid() { [[ -f "$1" ]] && [[ "$(cat "$1" 2>/dev/null)" == "$(grps_fingerprint)" ]]; }
# ^ this command verifies that the marker file ($1) exists and that the hash written inside the file is the same as the fingerprint 
# computed from ALL_GROUPS using the helper function set above, only returns 0 if everything is true 
mark_done() { grps_fingerprint > "$1"; }
# ^ this command computes the fingerprint and writes it into the marker file passed into the command as $1

# ─── 1. PARSE SAMPLE SHEET ────────────────────────────────────────────────────
# sample sheet Rows are collapsed into library groups (grouping technical replicates together), 
# using all "sample_name" information except the RunID field, which is different between technical replicates.
# this parsing also guards the silent failure modes that all look identical downstream 
# (i.e., fewer libraries than expected, missing trailing newline, spaces instead of tabs, CRLF line endings, and unrecognised assay tags). 

echo "🗂️  Parsing sample sheet..."
declare -A GRP_TO_IDS GRP_TO_PLOT GRP_TO_MERGED GRP_TO_ASSAY GRP_PAIRED
declare -A SAMPLE_TO_GRPS ID_TO_DL
declare -a ALL_GROUPS=() ALL_SAMPLES=()
sheet_rows=0; sheet_bad=0; lineno=0
# ^ this sets dictionaries/arrays and variables (such as line number) that will be important to parse the sample sheet

while IFS=$'\t' read -r id rbase plabel || [[ -n "${id:-}" ]]; do
# ^ reads the sample sheet and splits each line on tabs into three variables (id, rbase, plabel).
# The `|| [[ -n "${id:-}" ]]` part handles a file whose last line has NO trailing newline
# character: `read` still fills the variables but returns non-zero, which would normally end
# the loop and silently drop your last sample. Checking that id is non-empty runs that final
# iteration anyway. This is so you can easily paste your sample information in from Excel or Sheets.
    lineno=$(( lineno + 1 ))
    id="${id//$'\r'/}"; rbase="${rbase//$'\r'/}"; plabel="${plabel//$'\r'/}"
    if [[ -z "${id// /}" || "$id" == \#* ]]; then id=""; continue; fi
    if [[ -z "${rbase:-}" ]]; then
        echo "  ❌ line ${lineno}: no tab separator (spaces are not accepted): ${id}"
        sheet_bad=1; id=""; continue
    fi
# ^ this section cleans up the text and validates the formating of each row as the command reads through the sheet. The variable "lineno"
# keeps track of what line the command is on, so if there is a formating issue you know exactly where to look in the sheet. Examples of 
# cleaning is the id="${id//$'\r'/}" section which saves $id as a safe value without the invisible carriage return characters that can 
# mess things up downstream. Additionally, it makes sure to skip empty lines or lines that start with #, and continue safely to the rest 
# of the sheet. 
# If the sheet is space-separated instead of tab-separated, rbase comes back empty, sheet_bad=1
# is set, and the script reports the offending line number rather than failing later with a
# confusing error.

    id="$(echo -e "${id}" | tr -d '[:space:]')"
    rbase="$(echo -e "${rbase}" | tr -d '[:space:]')"
    plabel="$(echo -e "${plabel:-}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
# ^ this section finalizes the text clean up, deleting whitespace characters from the id and rbase variables, while deleting the leading
# and trailing white space for the plabel variable (in case you want to keep some white spaces for aesthetic plotting in the plot_label.

    IFS='_' read -r m e s cl ex br <<< "$rbase"
    if [[ -z "${br:-}" ]]; then
        echo "  ❌ line ${lineno}: sample_name needs 6 underscore fields, got '${rbase}'"
        sheet_bad=1; id=""; continue
    fi
# ^ this section reads each sample_name (rbase variable) and checks that it has the required 6 fields for sample processing, if it has
# fewer than 6 strings, br variable will be empty, it will set the sheet as bad, but instead of crashing right then and there, it 
# continues through to find any other potential errors and reports all in one run. 

    akey="$(assay_key "$m")"
    if [[ "$akey" == "UNKNOWN" ]]; then
        echo "  ❌ line ${lineno}: unrecognised assay tag '${m}' (expected ATAC, chip, CnT or CnR)"
        sheet_bad=1; id=""; continue
    fi
# ^ this checks that the assay/method assigned to the sample is something actually compatible with the script, 
# using the assay_key() normalization helper function designed above. 

    [[ -z "$plabel" ]] && plabel="${m}_${e}_${cl}_${br}"
    sheet_rows=$(( sheet_rows + 1 ))
# if you leave the plot_label column blank on the sample sheet, this automatically makes one from the sample_name field variables. the 
# sheet_rows variable counts the number of successfully parsed sheet rows 

    grp="${m}_${e}_${cl}_${ex}_${br}"          # library = technical reps merged
    sample="${m}_${e}_${cl}_${ex}"             # sample  = biological reps merged
# this constructs two unique identifier variables, one for the group of technical replicates (grp) 
# and one for the group of biological replicates (sample), useful downstream in the processing

    if [[ -z "${GRP_TO_IDS[$grp]:-}" ]]; then
        ALL_GROUPS+=("$grp")
        GRP_TO_PLOT[$grp]="$plabel"
        GRP_TO_ASSAY[$grp]="$akey"
        if [[ -z "${SAMPLE_TO_GRPS[$sample]:-}" ]]; then ALL_SAMPLES+=("$sample"); fi
        SAMPLE_TO_GRPS[$sample]+="$grp "
    fi
# ^ this section adds sample specific variables to the different arrays we set above, checking first if the grp label has been seen already. It
# then adds grp to ALL_GROUPS array, the plable (plot_label) to the GRP_TO_PLOT array, and the normalized assay key $akey to GRP_TO_ASSAY. 
# The last part of the section checks if the broader biological sample identifier is new, 
# and if so, adds it to the ALL_SAMPLES array, and then finally it links
# the biological replicate library (group of technical replicates) "$grp" to its parent biological sample ARRAY, building a space separated 
# list of all the biological replicate libraries that belong together. 

    GRP_TO_IDS[$grp]+="$id "
    ID_TO_DL[$id]=1
    id=""
done < "$MAPPING_FILE"
# ^ this closes the loop and adds the RunID to that specfic technical replicate's list of total IDs "GRP_TO_IDS", 
# this way when we make the master metrics file we can check which RunIDs (i.e. SRR# or ENC File ID) went into each technical replciate
# it also adds the RunID to the ID_TO_DL group which flags this specific ID for downloading using prefetch and fasterq-dump.
# finally, clearing the id variable makes sure no ID leaks into the next loop. 

expected_rows=$(grep -cvE '^[[:space:]]*(#|$)' "$MAPPING_FILE" || true)
[[ "$sheet_bad" -eq 0 ]] || { echo "  ❌ malformed sheet. Check: cat -A '${MAPPING_FILE}'"; exit 1; }
if [[ "$expected_rows" -ne "$sheet_rows" ]]; then
    echo "  ❌ parsed ${sheet_rows} rows but the file has ${expected_rows}."
    echo "     Usual cause: no trailing newline on the last line."
    echo "     Check: cat -A '${MAPPING_FILE}'   (^I = tab, \$ = line end)"
    exit 1
fi
# ^ this step counts the sheet rows that are neither comments nor blank. The regex matches a line that
# is only whitespace followed by either a "#" (comment) or the end of the line (blank line);
# grep -v inverts that, and -c counts what's left.

# The on-disk sample basename keeps every contributing RunID, so a BAM can always be
# traced back to its raw data without consulting the sheet.
for grp in "${ALL_GROUPS[@]}"; do
    read -ra ids <<< "${GRP_TO_IDS[$grp]}"
    joined_ids=$(IFS=-; echo "${ids[*]}")
    IFS='_' read -r m e cl ex br <<< "$grp"
    GRP_TO_MERGED[$grp]="${m}_${e}_${joined_ids}_${cl}_${ex}_${br}"
done
# ^ this step constructs the sample_name that is output in the Master_Sample_Metrics.tsv, that contains all the RunIDs 
# for a technical replicate in the same name, separated by a "-". The array GRP_TO_MERGED now holds the new, highly
# specific biological replicate names, with the merged technical replicate RunIDs. For example, the biological replicate:
# "chip_H3K27me3_MCF7_Bernstein_2017_rep2" is made up of 3 techical replicates, with RunIDs specified in the sample_name:
# "chip_H3K27me3_SRR5339274-SRR5339275-SRR5339276_MCF7_ID363_BER2". this name is used for all downstream BAM and BigWig file naming

echo "  → ${sheet_rows} rows, ${#ALL_GROUPS[@]} libraries, ${#ALL_SAMPLES[@]} sample groups"
for a in ATAC CHIP CNT CNR; do
    n=0
    for grp in "${ALL_GROUPS[@]}"; do [[ "${GRP_TO_ASSAY[$grp]}" == "$a" ]] && n=$(( n + 1 )); done
    [[ $n -gt 0 ]] && printf "     %-5s %3d libraries\n" "$a" "$n"
done
echo ""
# ^ finally this step prints out a summary of the sample sheet parsing, stating the number of rows found, 
# how many merged libraries there were (ALL_GROUPS), and how many groups of biological replicates there are (ALL_SAMPLES)
# and lastly giving you how many libraries belong to each assay/method type. 


# ─── 2. DOWNLOAD FASTQS ───────────────────────────────────────────────────────
# NOTE, if you add a technical replicate (new run ID) to the samples, the marker file fingerprint will 
# remain the same, because its based on the ${ALL_GROUPS[@]} which does not carry RunID information. If you add a new technical
# replicate, make sure to rm the download_complete.txt marker file
if ! marker_valid "$DOWNLOAD_COMPLETE"; then
    echo "📥 Downloading FASTQs..."
    step_failed=0
    for id in "${!ID_TO_DL[@]}"; do
        [[ -f "${FASTQ_DIR}/${id}_1.fastq.gz" ]] && continue
        echo "  → Fetching ${id}..."
        set +e
        if [[ "$id" == ENCFF* ]]; then
            wget --tries=3 --timeout=60 -O "${FASTQ_DIR}/${id}_1.fastq.gz" \
                "https://www.encodeproject.org/files/${id}/@@download/${id}.fastq.gz"
            [[ -s "${FASTQ_DIR}/${id}_1.fastq.gz" ]] || {
                echo "  ❌ ENCODE download failed for ${id}"
                rm -f "${FASTQ_DIR}/${id}_1.fastq.gz"; step_failed=1; set -e; continue; }
# this step downloads fastq files from ENCODE using the ENCFF ID number, with wget. if its not an ENCFF but rather 
# an SRR ID, the download follows the process below 

        elif [[ "$id" == SRR* ]]; then
            sra_log="${LOG_DIR}/${id}_download.log"
# this next step downloads fastq files from SRA database using prefetch, fasterq-dump, and pigz
# prefetch pulls the compressed archive, fasterq-dump extracts it into FASTQ format 
# (splitting paired ends into separate files), and pigz compresses those extracted files utilizing all 
# available CPU threads ($THREADS) for speed

            echo "     prefetch ${id}..."
            prefetch "$id" --max-size u -O "$FASTQ_DIR" 2>&1 | tee "$sra_log"
            prefetch_rc=${PIPESTATUS[0]}
# --max-size u lifts the 20 GB default cap that deep libraries can commonly hit.
            sra_path="${FASTQ_DIR}/${id}/${id}.sra"
            [[ -s "$sra_path" ]] || sra_path="${FASTQ_DIR}/${id}.sra"
            if [[ $prefetch_rc -ne 0 || ! -s "$sra_path" ]]; then
                echo "  ❌ prefetch failed for ${id} (rc=${prefetch_rc}); see ${sra_log}"
                echo "     If this node has no internet, run downloads on a login node."
                rm -rf "${FASTQ_DIR:?}/${id}"; step_failed=1; set -e; continue
            fi
            
            echo "     fasterq-dump ${id}..."
            fasterq-dump "$sra_path" --threads "$THREADS" \
                --outdir "$FASTQ_DIR" --split-files --force 2>&1 | tee -a "$sra_log"
            fqd_rc=${PIPESTATUS[0]}
# note, PIPESTATUS[0] holds the exit code for the first command in the pipe, in this case fasterq-dump

# --split-files writes "_1/_2" for paired runs and a single file for single-end runs, which is what the PE/SE detection below reads.
# --readids/-I is deliberately NOT used: it appends .1/.2 to read names so bowtie2 no longer recognises mates.

            mapfile -t fq_produced < <(
                compgen -G "${FASTQ_DIR}/${id}.fastq"
                compgen -G "${FASTQ_DIR}/${id}_*.fastq"
            )
            if [[ $fqd_rc -ne 0 || ${#fq_produced[@]} -eq 0 ]]; then
                echo "  ❌ fasterq-dump produced no FASTQ for ${id} (rc=${fqd_rc}); see ${sra_log}"
                rm -rf "${FASTQ_DIR:?}/${id}"; step_failed=1; set -e; continue
            fi

            echo "     compressing $(printf '%s ' "${fq_produced[@]##*/}")"
            pigz -p "$THREADS" "${fq_produced[@]}"
            rm -rf "${FASTQ_DIR:?}/${id}"

            if [[ ! -f "${FASTQ_DIR}/${id}_1.fastq.gz" && -f "${FASTQ_DIR}/${id}.fastq.gz" ]]; then
                mv "${FASTQ_DIR}/${id}.fastq.gz" "${FASTQ_DIR}/${id}_1.fastq.gz"
            fi
        else
            echo "  ❌ unrecognised accession type: ${id}"; step_failed=1
        fi
        [[ -f "${FASTQ_DIR}/${id}_1.fastq.gz" ]] || { echo "  ❌ Failed to download ${id}"; step_failed=1; }
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$DOWNLOAD_COMPLETE"
else echo "⏩ Skipping download."; fi

# ─── 3. FASTQ INVENTORY + SINGLE/PAIRED DETECTION ─────────────────────────────
# Library type (single or paired) is decided here ONCE and for all, from whether an "_2" FASTQ exists, and stored
# in GRP_PAIRED. Every later step reads that value, so no two steps can disagree about whether a library is paired 
# in this pipeline. 

echo "📋 FASTQ inventory and PE/SE detection:"
inv_bad=0
for grp in "${ALL_GROUPS[@]}"; do
    read -ra ids <<< "${GRP_TO_IDS[$grp]}"
    grp_type=""
    for id in "${ids[@]}"; do
        f1="${FASTQ_DIR}/${id}_1.fastq.gz"; f2="${FASTQ_DIR}/${id}_2.fastq.gz"
        if [[ -f "$f1" ]]; then
            if [[ -f "$f2" ]]; then this="PE"; else this="SE"; fi
            if [[ -n "$grp_type" && "$grp_type" != "$this" ]]; then
                echo "  ❌ ${grp}: technical replicates mix SE and PE; cannot concatenate."
                inv_bad=1
            fi
            grp_type="$this"
        else
            printf "  ❌ %-42s %-13s MISSING %s\n" "$grp" "$id" "$f1"
            inv_bad=1
        fi
    done
# this chunk detects whether a technical replicate is paired-end or single-end by checking if a _2.fastq.gz file exists alongside its matching 
# _1 file. if this is true the run is marked PE, if not, its marked SE. It also makes sure that technical replicates are either both PE or
# both SE, which is important to know before merging, right now the pipeline cannot handle mixed end technical replicates, as they are not super
# common. 
    
    [[ "$grp_type" == "PE" ]] && GRP_PAIRED[$grp]=1 || GRP_PAIRED[$grp]=0

    akey="${GRP_TO_ASSAY[$grp]}"
    if [[ "${ASSAY_REQUIRE_PE[$akey]}" -eq 1 && "${GRP_PAIRED[$grp]}" -eq 0 ]]; then
        echo "  ❌ ${grp}: ${akey} must be paired-end, but only _1 FASTQs were found."
        inv_bad=1
    fi
    printf "  ✅ %-42s %-5s %s\n" "$grp" "$akey" \
        "$( [[ ${GRP_PAIRED[$grp]} -eq 1 ]] && echo PAIRED || echo SINGLE )"
done
[[ $inv_bad -eq 0 ]] || { echo "  ❌ Missing or inconsistent FASTQs; fix before continuing."; exit 1; }

# Biological replicates that get merged must share a library type, or the merged
# track would combine fragment-based and extension-based coverage.
for s in "${ALL_SAMPLES[@]}"; do
    read -ra grps <<< "${SAMPLE_TO_GRPS[$s]}"
    akey="${GRP_TO_ASSAY[${grps[0]}]}"
    [[ "${ASSAY_MERGE_REPS[$akey]}" -eq 1 ]] || continue
    [[ ${#grps[@]} -lt 2 ]] && continue
    types=""
    for g in "${grps[@]}"; do [[ "${GRP_PAIRED[$g]}" -eq 1 ]] && types+="PE " || types+="SE "; done
    uniq_types=$(echo "$types" | tr ' ' '\n' | grep -v '^$' | sort -u | tr '\n' ' ')
    if [[ $(echo "$uniq_types" | wc -w) -gt 1 ]]; then
        echo "  ❌ sample '${s}' mixes ${uniq_types}across biological replicates; split it."
        exit 1
    fi
done
echo ""

# ─── 4. MERGE TECHNICAL REPLICATES + ADAPTER TRIM (fastp) ─────────────────────
if ! marker_valid "$CLEAN_COMPLETE"; then
    echo "🔗 Merging technical replicates and trimming (fastp)..."
    step_failed=0
    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        akey="${GRP_TO_ASSAY[$grp]}"
        is_paired="${GRP_PAIRED[$grp]}"
        read -ra ids <<< "${GRP_TO_IDS[$grp]}"

        out1="${CLEAN_FASTQ_DIR}/${merged_base}_1.clean.fastq.gz"
        out2="${CLEAN_FASTQ_DIR}/${merged_base}_2.clean.fastq.gz"
        [[ -f "$out1" ]] && continue

        echo "  → ${merged_base} (${akey}, $( [[ $is_paired -eq 1 ]] && echo PE || echo SE ))"
        raw1="${TMP_DIR}/${merged_base}_RAW_1.fastq.gz"
        raw2="${TMP_DIR}/${merged_base}_RAW_2.fastq.gz"
        files1=(); files2=()
        for id in "${ids[@]}"; do
            files1+=("${FASTQ_DIR}/${id}_1.fastq.gz")
            [[ $is_paired -eq 1 ]] && files2+=("${FASTQ_DIR}/${id}_2.fastq.gz")
        done
        cat "${files1[@]}" > "$raw1"
        [[ $is_paired -eq 1 ]] && cat "${files2[@]}" > "$raw2"

        fastp_opts=()
        [[ $is_paired -eq 1 && "${ASSAY_DETECT_ADAPTER_PE[$akey]}" -eq 1 ]] && \
            fastp_opts+=( --detect_adapter_for_pe )

        set +e
        if [[ $is_paired -eq 1 ]]; then
            fastp -i "$raw1" -I "$raw2" -o "$out1" -O "$out2" \
                  ${fastp_opts[@]+"${fastp_opts[@]}"} --thread "$THREADS" \
                  --json "${LOG_DIR}/${merged_base}_fastp.json" \
                  --html "${LOG_DIR}/${merged_base}_fastp.html" \
                  > "${LOG_DIR}/${merged_base}_fastp.log" 2>&1
        else
            fastp -i "$raw1" -o "$out1" --thread "$THREADS" \
                  --json "${LOG_DIR}/${merged_base}_fastp.json" \
                  --html "${LOG_DIR}/${merged_base}_fastp.html" \
                  > "${LOG_DIR}/${merged_base}_fastp.log" 2>&1
        fi
        rc=$?
        if [[ $rc -eq 0 && -f "$out1" ]]; then
            rm -f "$raw1" "$raw2"
        else
            echo "  ❌ fastp failed for ${merged_base}"; step_failed=1
        fi
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$CLEAN_COMPLETE"
else echo "⏩ Skipping trim."; fi

# ─── 5. ALIGNMENT (bowtie2) ───────────────────────────────────────────────────
# The per-assay option strings and -X fragment cap come from the tables at the
# top; -X is paired-end only because single-end reads have no insert size.
if ! marker_valid "$ALIGN_COMPLETE"; then
    echo "🧬 Aligning (bowtie2)..."
    step_failed=0
    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        akey="${GRP_TO_ASSAY[$grp]}"
        is_paired="${GRP_PAIRED[$grp]}"
        outbam="${ALIGNED_BAM_DIR}/${merged_base}.aligned.bam"
        sumfile="${LOG_DIR}/${merged_base}_bowtie2.txt"
        fq1="${CLEAN_FASTQ_DIR}/${merged_base}_1.clean.fastq.gz"
        fq2="${CLEAN_FASTQ_DIR}/${merged_base}_2.clean.fastq.gz"

        if [[ -f "$outbam" ]] && samtools quickcheck "$outbam" 2>/dev/null; then continue; fi
        [[ -f "$fq1" ]] || { echo "  ❌ missing clean FASTQ for ${merged_base}"; step_failed=1; continue; }

        echo "  → ${merged_base} (${akey}, $( [[ $is_paired -eq 1 ]] && echo PE || echo SE ))"
        set +e
        if [[ $is_paired -eq 1 ]]; then
            read -ra bt2_opts <<< "${ASSAY_BT2_PE[$akey]}"
            bowtie2 "${bt2_opts[@]}" -X "${ASSAY_MAXFRAG[$akey]}" -p "$THREADS" \
                    -x "$GENOME_INDEX" -1 "$fq1" -2 "$fq2" 2> "$sumfile" \
              | samtools view -@ "$THREADS" -b -o "$outbam" -
        else
            read -ra bt2_opts <<< "${ASSAY_BT2_SE[$akey]}"
            bowtie2 "${bt2_opts[@]}" -p "$THREADS" \
                    -x "$GENOME_INDEX" -U "$fq1" 2> "$sumfile" \
              | samtools view -@ "$THREADS" -b -o "$outbam" -
        fi
        rc=$?
        if [[ $rc -ne 0 ]] || ! samtools quickcheck "$outbam" 2>/dev/null; then
            echo "  ❌ alignment failed for ${merged_base}"; rm -f "$outbam"; step_failed=1
        fi
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$ALIGN_COMPLETE"
else echo "⏩ Skipping alignment."; fi

# ─── 6. FILTERING (ENCODE-style, shared by all assays) ────────────────────────
# -F excludes a read if ANY listed bit is set; -f keeps it only if the bit IS set.
#
#   1804 = 4 (unmapped) + 8 (mate unmapped) + 256 (secondary)
#          + 512 (fails vendor QC) + 1024 (duplicate)
#   -f 2 = properly paired                                  [paired-end only]
#
# PAIRED-END runs the filter twice, around fixmate. Pairing information lives
# inside each read's own record, so the first pass can delete one mate and leave
# the other behind still claiming to be properly paired with a read that no
# longer exists. fixmate walks the name-sorted file (mates adjacent), repairs
# the mate flags, coordinates and insert size, and clears the proper-pair bit on
# those orphans, which is what lets the second pass remove them. fixmate -m also
# adds the mate-score tag that `samtools markdup` requires.
#
# Because bowtie2 -X already capped the fragment length and -f 2 keeps only
# proper pairs, this step also performs the fragment-size filtering; the
# separate TLEN-based awk filter used in the old ChIP script is redundant.
if ! marker_valid "$FILTER_COMPLETE"; then
    echo "🔍 Filtering (-F 1804, per-assay MAPQ; PE adds -f 2 twice around fixmate)..."
    step_failed=0
    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        akey="${GRP_TO_ASSAY[$grp]}"
        mapq="${ASSAY_MAPQ[$akey]}"
        is_paired="${GRP_PAIRED[$grp]}"
        inbam="${ALIGNED_BAM_DIR}/${merged_base}.aligned.bam"
        filt_bam="${FILTERED_BAM_DIR}/${merged_base}.filtered.bam"

        if [[ -f "$filt_bam" ]] && samtools quickcheck "$filt_bam" 2>/dev/null; then continue; fi
        samtools quickcheck "$inbam" 2>/dev/null \
            || { echo "  ❌ invalid aligned BAM: ${merged_base}"; step_failed=1; continue; }

        echo "  → ${merged_base} (${akey}, MAPQ ${mapq})"
        set +e
        if [[ $is_paired -eq 1 ]]; then
            samtools view -@ "$THREADS" -F 1804 -f 2 -q "$mapq" -u "$inbam" \
              | samtools sort -@ "$THREADS" -n -T "${TMP_DIR}/${merged_base}_ns" \
                              -o "${TMP_DIR}/${merged_base}_ns.bam" -
            samtools fixmate -@ "$THREADS" -r -m \
                "${TMP_DIR}/${merged_base}_ns.bam" "${TMP_DIR}/${merged_base}_fix.bam"
            samtools view -@ "$THREADS" -F 1804 -f 2 -u "${TMP_DIR}/${merged_base}_fix.bam" \
              | samtools sort -@ "$THREADS" -T "${TMP_DIR}/${merged_base}_cs" -o "$filt_bam" -
            rm -f "${TMP_DIR}/${merged_base}_ns.bam" "${TMP_DIR}/${merged_base}_fix.bam"
        else
            samtools view -@ "$THREADS" -F 1804 -q "$mapq" -u "$inbam" \
              | samtools sort -@ "$THREADS" -T "${TMP_DIR}/${merged_base}_cs" -o "$filt_bam" -
        fi
        rc=$?
        if [[ $rc -ne 0 ]] || ! samtools quickcheck "$filt_bam" 2>/dev/null; then
            echo "  ❌ filtering failed for ${merged_base}"; rm -f "$filt_bam"; step_failed=1
        fi
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$FILTER_COMPLETE"
else echo "⏩ Skipping filtering."; fi

# ─── 7. DEDUPLICATION (shared by all assays) ──────────────────────────────────
# `samtools markdup -r -s` removes duplicates and writes the statistics block
# that the metrics table parses. For pairs it uses both fragment ends; for
# single-end reads the read's own unclipped 5' position and strand.
#
# Deduplication runs BEFORE chrM removal so that the reported duplicate rate
# describes the library as sequenced, and so the mitochondrial fraction stays
# measurable from idxstats on the deduplicated BAM.
if ! marker_valid "$DUP_COMPLETE"; then
    echo "👯 Deduplicating..."
    step_failed=0
    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        filt_bam="${FILTERED_BAM_DIR}/${merged_base}.filtered.bam"
        dedup_bam="${DEDUP_BAM_DIR}/${merged_base}.dedup.bam"

        if [[ -f "$dedup_bam" ]] && samtools quickcheck "$dedup_bam" 2>/dev/null; then continue; fi
        samtools quickcheck "$filt_bam" 2>/dev/null \
            || { echo "  ❌ invalid filtered BAM: ${merged_base}"; step_failed=1; continue; }

        echo "  → ${merged_base}"
        set +e
        samtools markdup -@ "$THREADS" -r -s "$filt_bam" "$dedup_bam" \
            2> "${LOG_DIR}/${merged_base}_markdup.log"
        rc=$?
        if [[ $rc -eq 0 ]] && samtools quickcheck "$dedup_bam" 2>/dev/null; then
            samtools index -@ "$THREADS" "$dedup_bam"
            samtools idxstats "$dedup_bam" > "${LOG_DIR}/${merged_base}_idxstats.txt"
        else
            echo "  ❌ dedup failed for ${merged_base}"; rm -f "$dedup_bam"; step_failed=1
        fi
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$DUP_COMPLETE"
else echo "⏩ Skipping deduplication."; fi

# ─── 8. FINAL BAMs: chromosome restriction + blacklist ────────────────────────
# Two cleanups that the old scripts did in different places and different ways
# are done here once, identically for every assay.
#
# Chromosome restriction keeps chr1-22, chrX, chrY and therefore drops chrM,
# unplaced contigs, random scaffolds and alt haplotypes in one operation.
# Removing chrM matters most for ATAC -- Tn5 attacks the naked mitochondrial
# genome so efficiently that 20-50% of an ATAC library can be mitochondrial --
# but it is correct for every assay, and doing it identically keeps library
# sizes comparable across assays.
#
# The blacklist is subtracted at the BAM level rather than at bigWig time, so
# the final BAM, the HOMER tag directory and the bigWig all describe exactly the
# same set of reads.
if ! marker_valid "$FINALBAM_COMPLETE"; then
    echo "🚫 Restricting to chr1-22/X/Y and subtracting blacklist..."
    step_failed=0
    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        dedup_bam="${DEDUP_BAM_DIR}/${merged_base}.dedup.bam"
        chr_bam="${TMP_DIR}/${merged_base}.chrfilt.bam"
        final_bam="${FINAL_BAM_DIR}/${merged_base}.clean.final.bam"

        if [[ -f "$final_bam" ]] && samtools quickcheck "$final_bam" 2>/dev/null; then continue; fi
        samtools quickcheck "$dedup_bam" 2>/dev/null \
            || { echo "  ❌ invalid dedup BAM: ${merged_base}"; step_failed=1; continue; }

        echo "  → ${merged_base}"
        set +e
        keep_chrs=$(samtools idxstats "$dedup_bam" | cut -f1 \
                    | grep -E '^chr([0-9]{1,2}|X|Y)$' | tr '\n' ' ')
        if [[ -z "$keep_chrs" ]]; then
            echo "  ❌ no chr1-22/X/Y found in ${merged_base} -- is the index UCSC-style?"
            step_failed=1; set -e; continue
        fi
        # shellcheck disable=SC2086
        samtools view -@ "$THREADS" -b "$dedup_bam" $keep_chrs > "$chr_bam"
        bedtools intersect -v -abam "$chr_bam" -b "$BLACKLIST" \
          | samtools sort -@ "$THREADS" -T "${TMP_DIR}/${merged_base}_fs" -o "$final_bam" -
        rc=$?
        if [[ $rc -eq 0 ]] && samtools quickcheck "$final_bam" 2>/dev/null; then
            samtools index -@ "$THREADS" "$final_bam"
            rm -f "$chr_bam"
        else
            echo "  ❌ chrM/blacklist step failed for ${merged_base}"; rm -f "$final_bam"; step_failed=1
        fi
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$FINALBAM_COMPLETE"
else echo "⏩ Skipping chrM/blacklist."; fi

# ─── 9. BIGWIGS ───────────────────────────────────────────────────────────────
# (a) one raw-coverage track per library, all assays
# (b) one CPM-normalised track per sample group with biological replicates
#     merged, only for assays whose ASSAY_MERGE_REPS flag is 1 (ATAC here)
#
# Coverage is fragment-based in both cases: paired-end reads are extended to
# their mate so each fragment contributes its true measured length, single-end
# reads are extended to the per-assay fixed estimate.
if ! marker_valid "$BIGWIG_COMPLETE"; then
    echo "📈 Creating bigWigs..."
    step_failed=0

    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        akey="${GRP_TO_ASSAY[$grp]}"
        is_paired="${GRP_PAIRED[$grp]}"
        inbam="${FINAL_BAM_DIR}/${merged_base}.clean.final.bam"
        bw="${BIGWIG_DIR}/${merged_base}_raw.bw"

        [[ -f "$bw" && -f "${bw}.md5" ]] && continue
        samtools quickcheck "$inbam" 2>/dev/null \
            || { echo "  ❌ invalid final BAM: ${merged_base}"; step_failed=1; continue; }

        echo "  → raw: ${merged_base}"
        set +e
        bc=( -b "$inbam" --binSize "$BIGWIG_BIN_SIZE" -p "$THREADS" -o "$bw" )
        if [[ $is_paired -eq 1 ]]; then
            bc+=( --extendReads )
        else
            bc+=( --extendReads "${ASSAY_SE_EXTEND[$akey]}" )
        fi
        bamCoverage "${bc[@]}"
        rc=$?
        if [[ $rc -eq 0 && -f "$bw" ]]; then
            md5sum "$bw" | awk '{print $1}' > "${bw}.md5"
            cp -f "$bw" "$WEB_DIR_MASTER/" 2>/dev/null || true
        else
            echo "  ❌ bigWig failed for ${merged_base}"; step_failed=1
        fi
        set -e
    done

    for sample in "${ALL_SAMPLES[@]}"; do
        read -ra grps <<< "${SAMPLE_TO_GRPS[$sample]}"
        akey="${GRP_TO_ASSAY[${grps[0]}]}"
        [[ "${ASSAY_MERGE_REPS[$akey]}" -eq 1 ]] || continue
        [[ ${#grps[@]} -lt 2 ]] && continue
        is_paired="${GRP_PAIRED[${grps[0]}]}"
        merged_bam="${FINAL_BAM_DIR}/${sample}_merged.clean.final.bam"
        merged_bw="${BIGWIG_DIR}/${sample}_merged_CPM.bw"
        [[ -f "$merged_bw" && -f "${merged_bw}.md5" ]] && continue

        echo "  → merged CPM: ${sample} (${#grps[@]} biological replicates)"
        set +e
        parts=()
        for g in "${grps[@]}"; do
            parts+=("${FINAL_BAM_DIR}/${GRP_TO_MERGED[$g]}.clean.final.bam")
        done

        # A straight merge weights each replicate by its depth, so a much deeper
        # replicate dominates the merged track. CPM normalisation afterwards
        # fixes the overall scale, not that imbalance.
        depths=(); minD=""; maxD=""
        for p in "${parts[@]}"; do
            d=$(samtools view -@ "$THREADS" -c -F 0x4 "$p" 2>/dev/null || echo 0)
            depths+=("$d")
            [[ -z "$minD" || "$d" -lt "$minD" ]] && minD=$d
            [[ -z "$maxD" || "$d" -gt "$maxD" ]] && maxD=$d
        done
        echo "     replicate depths: ${depths[*]}"
        if [[ -n "${minD:-}" && "$minD" -gt 0 ]]; then
            ratio=$(awk -v a="$maxD" -v b="$minD" 'BEGIN{printf "%.2f", a/b}')
            echo "     max/min depth ratio: ${ratio}"
            awk -v r="$ratio" 'BEGIN{ if(r>1.5) print "     ⚠️  replicates differ by >1.5x; the deeper one dominates the merged track" }'
        fi

        if [[ ! -f "$merged_bam" ]]; then
            samtools merge -@ "$THREADS" -f "$merged_bam" "${parts[@]}"
            samtools index -@ "$THREADS" "$merged_bam"
        fi
        mbc=( -b "$merged_bam" --normalizeUsing CPM --binSize "$MERGED_BIN_SIZE"
              -p "$THREADS" -o "$merged_bw" )
        if [[ $is_paired -eq 1 ]]; then
            mbc+=( --extendReads )
        else
            mbc+=( --extendReads "${ASSAY_SE_EXTEND[$akey]}" )
        fi
        bamCoverage "${mbc[@]}"
        rc=$?
        if [[ $rc -eq 0 && -f "$merged_bw" ]]; then
            md5sum "$merged_bw" | awk '{print $1}' > "${merged_bw}.md5"
            cp -f "$merged_bw" "$WEB_DIR_MASTER/" 2>/dev/null || true
        else
            echo "  ❌ merged bigWig failed for ${sample}"; step_failed=1
        fi
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$BIGWIG_COMPLETE"
else echo "⏩ Skipping bigWigs."; fi

# ─── 10. HOMER TAG DIRECTORIES + QC PLOTS ─────────────────────────────────────
# Built for every assay, including ATAC. Three panels per library:
#   1. nucleotide frequency across HOMER's full tagFreq range, which shows the
#      Tn5 insertion preference for ATAC and CUT&Tag and adapter or bias
#      artefacts for ChIP;
#   2. the same data zoomed to +/- FREQ_ZOOM_BP around the read 5' end, where
#      the enzyme's sequence preference actually lives;
#   3. the fragment GC distribution against the hg38 background, the read-out
#      most relevant to the GC-bias question in this project.
# Panels 1 and 2 are two views of one file, so the zoom costs no extra compute.
if ! marker_valid "$HOMER_QC_COMPLETE"; then
    echo "📊 Generating HOMER tag directories and QC plots..."
    step_failed=0
    for grp in "${ALL_GROUPS[@]}"; do
        merged_base="${GRP_TO_MERGED[$grp]}"
        akey="${GRP_TO_ASSAY[$grp]}"
        inbam="${FINAL_BAM_DIR}/${merged_base}.clean.final.bam"
        tag_dir="${HOMER_QC_DIR}/${merged_base}_tagdir"
        out_pdf="${HOMER_QC_DIR}/${merged_base}_QC_Plots.pdf"

        [[ -f "$out_pdf" ]] && continue
        samtools quickcheck "$inbam" 2>/dev/null \
            || { echo "  ❌ invalid final BAM: ${merged_base}"; step_failed=1; continue; }

        echo "  → QC: ${merged_base} (${akey})"
        set +e
        if [[ ! -d "$tag_dir" ]]; then
            makeTagDirectory "$tag_dir" "$inbam" -genome "$GENOME" -checkGC \
                -fragLength 200 > "${LOG_DIR}/${merged_base}_makeTagDir.log" 2>&1
        fi

        export TAG_DIR="$tag_dir" OUT_PDF="$out_pdf" SAMPLE_NAME="${GRP_TO_PLOT[$grp]}"
        export FREQ_ZOOM_BP="$FREQ_ZOOM_BP"
        Rscript - <<'EOF'
options(warn = -1)
suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(readr); library(tidyr); library(patchwork)
})

# Font sizes are 1.5x the original QC script (base 16 -> 24, title 20 -> 30).
# Panel width is held at 9 inches each, so the canvas grows with the number of
# panels (2 panels = 18x7, 3 panels = 27x7) rather than squeezing them.
BASE_SIZE   <- 24
TITLE_SIZE  <- 30
PANEL_W_IN  <- 9
PANEL_H_IN  <- 7

zoom <- suppressWarnings(as.numeric(Sys.getenv("FREQ_ZOOM_BP")))
if (is.na(zoom) || zoom <= 0) zoom <- 50

tdir <- Sys.getenv("TAG_DIR")
if (!file.exists(file.path(tdir, "tagInfo.txt"))) quit(status = 1)

NUC_COLORS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", T = "#D62828")

# One builder used for both nucleotide panels, so the zoomed-out and zoomed-in
# views can never drift apart in styling. xlim = NULL means "show everything".
freq_panel <- function(df, ttl, xlim = NULL) {
    p <- ggplot(df, aes(x = Offset, y = Frequency, color = Nucleotide)) +
         geom_line(linewidth = 1.2) +
         scale_color_manual(values = NUC_COLORS) +
         scale_y_continuous(limits = c(0, 0.5)) +
         theme_minimal(base_size = BASE_SIZE) +
         labs(title = ttl, x = "Distance from 5' end (bp)", y = "Frequency")
    if (!is.null(xlim)) {
        # coord_cartesian zooms the view; xlim() would DROP the outside rows and
        # break the lines at the panel edge instead of clipping them cleanly.
        p <- p + geom_vline(xintercept = 0, linetype = "dashed",
                            linewidth = 0.6, color = "grey40") +
                 coord_cartesian(xlim = xlim) +
                 scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = diff(xlim) / 4))
    }
    p
}

p1 <- ggplot() + theme_void() + ggtitle("Missing tagFreq.txt")
p2 <- ggplot() + theme_void() + ggtitle("Missing tagFreq.txt")
if (file.exists(file.path(tdir, "tagFreq.txt"))) {
    f_df <- suppressMessages(read_tsv(file.path(tdir, "tagFreq.txt"), show_col_types = FALSE))[, 1:5]
    colnames(f_df) <- c("Offset", "A", "C", "G", "T")
    f_long <- pivot_longer(f_df, cols = c("A", "C", "G", "T"),
                           names_to = "Nucleotide", values_to = "Frequency")
    rng <- range(f_long$Offset, na.rm = TRUE)
    # as.integer() matters: sprintf("%+d") errors on the doubles that read_tsv
    # and as.numeric() return.
    p1 <- freq_panel(f_long,
                     sprintf("Nucleotide Frequency (%+d to %+d bp)",
                             as.integer(rng[1]), as.integer(rng[2])))
    p2 <- freq_panel(f_long,
                     sprintf("Nucleotide Frequency (%+d to %+d bp)",
                             as.integer(-zoom), as.integer(zoom)),
                     xlim = c(-zoom, zoom))
}

p3 <- ggplot() + theme_void() + ggtitle("Missing GC files")
gc_s <- file.path(tdir, "tagGCcontent.txt")
gc_g <- file.path(tdir, "genomeGCcontent.txt")
if (file.exists(gc_s) && file.exists(gc_g)) {
    s_tab <- suppressMessages(read_tsv(gc_s, show_col_types = FALSE))
    g_tab <- suppressMessages(read_tsv(gc_g, show_col_types = FALSE))
    ds <- data.frame(GC = as.numeric(s_tab[[1]]), Fraction = as.numeric(s_tab[[3]]), Type = "Sample")
    dg <- data.frame(GC = as.numeric(g_tab[[1]]), Fraction = as.numeric(g_tab[[3]]), Type = "hg38")
    p3 <- ggplot(bind_rows(ds, dg), aes(x = GC, y = Fraction, color = Type, fill = Type)) +
          geom_area(alpha = 0.2, position = "identity") +
          geom_line(linewidth = 1.5) +
          scale_color_manual(values = c(Sample = "#6a51a3", hg38 = "#969696")) +
          scale_fill_manual(values  = c(Sample = "#6a51a3", hg38 = "#969696")) +
          theme_minimal(base_size = BASE_SIZE) +
          labs(title = "GC% Distribution (200bp window)", x = "GC%", y = "Fraction")
}

panels <- list(p1, p2, p3)
final_plot <- wrap_plots(panels, nrow = 1) +
    plot_annotation(title = Sys.getenv("SAMPLE_NAME"),
                    theme = theme(plot.title = element_text(size = TITLE_SIZE,
                                                            face = "bold", hjust = 0.5)))
ggsave(Sys.getenv("OUT_PDF"), plot = final_plot,
       width = PANEL_W_IN * length(panels), height = PANEL_H_IN,
       device = "pdf", bg = "white")
EOF
        rc=$?
        [[ $rc -ne 0 ]] && { echo "  ❌ QC failed for ${merged_base}"; step_failed=1; }
        set -e
    done
    [[ $step_failed -eq 0 ]] && mark_done "$HOMER_QC_COMPLETE"
else echo "⏩ Skipping HOMER QC."; fi

# ─── 11. PROCESSING METRICS + UCSC TRACK LINES ────────────────────────────────
echo "📝 Compiling processing metrics..."

printf "IDs\tsample_name\tplot_label\tassay\tread_ends\traw_read_count\talignment_rate\taligned_read_count\tfraction_reads_pass_filtering\tpost_filter_read_count\tduplicate_rate\tduplicate_count\tfrac_mito\tfinal_mapped_fragments\tlocation_final_bam\tlocation_final_bigwig\tlocation_homer_tagdir\tbigwig_md5sum\tucsc_genome_browser_track_text\n" > "$MASTER_METRICS"

{
  echo "# UCSC genome browser track lines"
  echo "# Paste into https://genome.ucsc.edu/cgi-bin/hgCustom"
  echo ""
} > "$UCSC_TRACKS"

for grp in "${ALL_GROUPS[@]}"; do
    merged_base="${GRP_TO_MERGED[$grp]}"
    akey="${GRP_TO_ASSAY[$grp]}"
    plot_label="${GRP_TO_PLOT[$grp]}"
    is_paired="${GRP_PAIRED[$grp]}"
    ids_str=$(echo "${GRP_TO_IDS[$grp]}" | tr ' ' ','); ids_str=${ids_str%,}
    read_ends=$( [[ $is_paired -eq 1 ]] && echo "Paired" || echo "Single" )

    aligned_bam="${ALIGNED_BAM_DIR}/${merged_base}.aligned.bam"
    filt_bam="${FILTERED_BAM_DIR}/${merged_base}.filtered.bam"
    final_bam="${FINAL_BAM_DIR}/${merged_base}.clean.final.bam"
    bw="${BIGWIG_DIR}/${merged_base}_raw.bw"
    tag_dir="${HOMER_QC_DIR}/${merged_base}_tagdir"

    raw_reads="NA"; align_rate="NA"
    if [[ -f "${LOG_DIR}/${merged_base}_bowtie2.txt" ]]; then
        raw_reads=$(grep -m1 "reads; of these:" "${LOG_DIR}/${merged_base}_bowtie2.txt" | awk '{print $1}')
        align_rate=$(grep -m1 "overall alignment rate" "${LOG_DIR}/${merged_base}_bowtie2.txt" | awk '{print $1}')
        [[ -z "$raw_reads"  ]] && raw_reads="NA"
        [[ -z "$align_rate" ]] && align_rate="NA"
    fi

    aligned_reads=$(bam_mapped_reads "$aligned_bam")
    post_filt_reads=$(bam_mapped_reads "$filt_bam")
    final_reads=$(bam_mapped_reads "$final_bam")

    frac_pass="NA"
    [[ "$aligned_reads" != "NA" && "$post_filt_reads" != "NA" ]] && \
        frac_pass=$(safe_div "$post_filt_reads" "$aligned_reads")

    dup_count="NA"; dup_rate="NA"
    if [[ -f "${LOG_DIR}/${merged_base}_markdup.log" ]]; then
        dup_count=$(grep -m1 "DUPLICATE TOTAL:" "${LOG_DIR}/${merged_base}_markdup.log" | awk '{print $3}')
        read_ct=$(grep -m1 "^READ:" "${LOG_DIR}/${merged_base}_markdup.log" | awk '{print $2}')
        [[ -z "$read_ct" ]] && read_ct="$post_filt_reads"
        [[ -n "$dup_count" && "$read_ct" != "NA" ]] && dup_rate=$(safe_div "$dup_count" "$read_ct")
        [[ -z "$dup_count" ]] && dup_count="NA"
    fi

    frac_mito="NA"
    if [[ -f "${LOG_DIR}/${merged_base}_idxstats.txt" ]]; then
        frac_mito=$(awk 'BEGIN{m=0;t=0} {t+=$3} $1=="chrM"{m=$3}
                         END{ if(t>0) printf "%.4f", m/t; else printf "NA" }' \
                    "${LOG_DIR}/${merged_base}_idxstats.txt")
    fi

    # A single-end read is one fragment; a paired-end fragment is two reads.
    final_frags="NA"
    if [[ "$final_reads" != "NA" ]]; then
        if [[ $is_paired -eq 1 ]]; then final_frags=$(( final_reads / 2 )); else final_frags=$final_reads; fi
    fi

    bw_md5="NA"; [[ -f "${bw}.md5" ]] && bw_md5=$(cat "${bw}.md5")
    bw_url="${WEB_URL_MASTER}/$(basename "$bw")"
    color="${ASSAY_COLOR[$akey]}"
    track_name="${plot_label} Raw"
    ucsc_text="track type=bigWig name=\"${track_name}\" description=\"${track_name}\" bigDataUrl=${bw_url} color=${color} visibility=full yLineOnOff=on autoScale=on alwaysZero=on graphType=bar maxHeightPixels=128:75:11 windowingFunction=maximum smoothingWindow=off"
    echo "$ucsc_text" >> "$UCSC_TRACKS"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$ids_str" "$merged_base" "$plot_label" "$akey" "$read_ends" "$raw_reads" \
        "$align_rate" "$aligned_reads" "$frac_pass" "$post_filt_reads" "$dup_rate" \
        "$dup_count" "$frac_mito" "$final_frags" "$final_bam" "$bw" "$tag_dir" \
        "$bw_md5" "$ucsc_text" >> "$MASTER_METRICS"
done

# Merged biological-replicate tracks (ATAC only, by ASSAY_MERGE_REPS)
printf "sample_group\tassay\tn_bio_reps\tcontributing_libraries\tlocation_merged_bam\tlocation_merged_bigwig\tbigwig_md5sum\tucsc_genome_browser_track_text\n" > "$MERGED_METRICS"
echo "" >> "$UCSC_TRACKS"
for sample in "${ALL_SAMPLES[@]}"; do
    read -ra grps <<< "${SAMPLE_TO_GRPS[$sample]}"
    akey="${GRP_TO_ASSAY[${grps[0]}]}"
    [[ "${ASSAY_MERGE_REPS[$akey]}" -eq 1 ]] || continue
    [[ ${#grps[@]} -lt 2 ]] && continue
    merged_bam="${FINAL_BAM_DIR}/${sample}_merged.clean.final.bam"
    merged_bw="${BIGWIG_DIR}/${sample}_merged_CPM.bw"
    [[ -f "$merged_bw" ]] || continue

    base_label="${GRP_TO_PLOT[${grps[0]}]}"
    mtrack="${base_label%_*} merged (CPM)"
    murl="${WEB_URL_MASTER}/$(basename "$merged_bw")"
    mmd5="NA"; [[ -f "${merged_bw}.md5" ]] && mmd5=$(cat "${merged_bw}.md5")
    libs=""; for g in "${grps[@]}"; do libs+="${GRP_TO_MERGED[$g]},"; done; libs=${libs%,}

    mtext="track type=bigWig name=\"${mtrack}\" description=\"${sample}, ${#grps[@]} biological replicates merged, CPM-normalized\" bigDataUrl=${murl} color=${ASSAY_COLOR[$akey]} visibility=full yLineOnOff=on autoScale=on alwaysZero=on graphType=bar maxHeightPixels=128:75:11 windowingFunction=maximum smoothingWindow=off"
    echo "$mtext" >> "$UCSC_TRACKS"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$akey" "${#grps[@]}" "$libs" "$merged_bam" "$merged_bw" \
        "$mmd5" "$mtext" >> "$MERGED_METRICS"
done

echo ""
echo "✅ Processing metrics : $MASTER_METRICS"
echo "✅ Merged-track table : $MERGED_METRICS"
echo "✅ UCSC track lines   : $UCSC_TRACKS"
echo "✅ Final BAMs         : $FINAL_BAM_DIR"
echo "✅ BigWigs            : $BIGWIG_DIR"
echo "✅ HOMER QC           : $HOMER_QC_DIR"
echo "🎉 Pipeline finished."