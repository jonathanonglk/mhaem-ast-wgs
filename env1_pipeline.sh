#!/usr/bin/env bash
# =============================================================================
# env_1 — reference-based pipeline for M. haemophilum WGS
# =============================================================================
# Combines, in one conda environment (default name: env_1):
#   PART A  read trimming, QC, Snippy variant calling, snippy-core, IQ-TREE
#           (formerly mhaem_full_pipeline.sh)
#   PART B  SNP-based GWAS, targeted resistance-gene analysis, gyrA intein
#           status, Manhattan/QQ plots (formerly gwas11.sh)
#
# Environment setup: see SETUP.md (section "env_1") before running.
#
# Usage:
#   conda activate env_1          # or let the script activate it
#   bash env1_pipeline.sh [KEY=VALUE ...]
#
# Every path below can be set either as an environment variable or as a
# KEY=VALUE argument, e.g.
#   bash env1_pipeline.sh RAW_READS=/data/reads REF_SRC=/data/ref.gbff \
#                         MIC_AMI=/data/AMI.txt MIC_CIP=/data/CIP.txt \
#                         MIC_RIF=/data/RIF.txt
#
# Run "bash env1_pipeline.sh --help" for the full list.
# =============================================================================
set -euo pipefail

usage() {
    cat <<'USAGE'
env1_pipeline.sh [KEY=VALUE ...]

Settings (environment variable or KEY=VALUE argument; defaults in brackets):

  ENV_NAME     conda environment to run in                 [env_1]
  WORKDIR      working directory for all outputs           [$HOME/mhaem_wgs]
  RAW_READS    directory of raw paired FASTQ files,
               named <sample>_R1.fastq.gz / _R2.fastq.gz   [$WORKDIR/raw_reads]
  REF_SRC      annotated reference genome (GenBank .gbff)   [$HOME/GCF_000340435.2.gbff]
  MIC_AMI      amikacin MIC rank file (sample<TAB>rank)     [$HOME/AMI.txt]
  MIC_CIP      ciprofloxacin MIC rank file                  [$HOME/CIP.txt]
  MIC_RIF      rifampicin MIC rank file                     [$HOME/RIF.txt]
  OUTPUT_DEST  directory to copy final outputs to           [$WORKDIR/results]
  THREADS      CPU threads                                  [4]

Example:
  bash env1_pipeline.sh RAW_READS=~/data/fastq REF_SRC=~/data/ref.gbff \
                        MIC_AMI=~/data/AMI.txt MIC_CIP=~/data/CIP.txt \
                        MIC_RIF=~/data/RIF.txt
USAGE
}

# ---------------------------------------------------------------------------
# Command-line overrides (KEY=VALUE), applied before defaults
# ---------------------------------------------------------------------------
for arg in "$@"; do
    case "$arg" in
        -h|--help) usage; exit 0 ;;
        [A-Za-z_]*=*) declare "$arg" ;;
        *) echo "ERROR: unrecognised argument '$arg' (expected KEY=VALUE)"; exit 1 ;;
    esac
done

# ===========================================================================
# USER CONFIGURATION — edit here, or override as shown in usage above
# ===========================================================================
ENV_NAME="${ENV_NAME:-env_1}"
WORKDIR="${WORKDIR:-$HOME/mhaem_wgs}"
RAW_READS="${RAW_READS:-$WORKDIR/raw_reads}"
REF_SRC="${REF_SRC:-$HOME/GCF_000340435.2.gbff}"
MIC_AMI="${MIC_AMI:-$HOME/AMI.txt}"
MIC_CIP="${MIC_CIP:-$HOME/CIP.txt}"
MIC_RIF="${MIC_RIF:-$HOME/RIF.txt}"
OUTPUT_DEST="${OUTPUT_DEST:-$WORKDIR/results}"
THREADS="${THREADS:-4}"

# Isolates handled specially (do not edit unless the isolate list changes)

KINSHIP="$WORKDIR/phylogeny_kinship_C.tsv"              # pyseer phylogeny_distance.py --lmm --midpoint
KINSHIP_LEGACY="$WORKDIR/phylogeny_kinship_legacy.tsv"   # old 1 - d/max matrix (sensitivity analysis only)

# Clinical isolates (one per patient) — the ONLY samples analysed in the GWAS
# and in all diagnostic statistics. Identical to the CLINICAL list in
# pangenome.sh. Replicates, ATCC 29548 (CT-WNTC9548) and CT-WN17CD6059
# (contaminated; replaced by CT-WN18CA1178) are excluded so that pyseer's af
# and the 2x2 tables share the same denominator (n = 56) and same-patient
# isolates are not counted twice.
CLINICAL=(
    CT-WN15CC3812 CT-WN18BA5462 CT-WN15CB4880 CT-WN17CC6853 CT-WN16CC8097
    CT-WN15CC0814 CT-WN18CC5430 CT-WN17CA8830 CT-WN17CD0730 CT-WN18BA6316
    CT-WN15CC3097 CT-WN18CA1953 CT-WN15CB5439 CT-WN17CD1188 CT-WN18CA1163
    CT-WN15CD0483 CT-WN15CD1369 CT-WN18CA7807 CT-WN16BA7594 CT-WN16CC2935
    CT-WN16CD2223 CT-WN17CB3038 CT-WN17CA3088 CT-WN17CB5288 CT-WN18CA8070
    CT-WN18CB7320 CT-WN17BA5155 CT-WN15CC6561 CT-WN16CA5839 CT-WN17CB9744
    CT-WN16CC7074 CT-WN18CD4765 CT-WN17CD6535 CT-WN17BA6921 CT-WN18CD7756
    CT-WN15CA7460 CT-WN18BA5406 CT-WN15CC1919 CT-WN15CA4590 CT-WN17CA6664
    CT-WN18CA3189 CT-WN17CA4749 CT-WN18BA1483 CT-WN15CA8741 CT-WN17BA6280
    CT-WN18CD0191 CT-WN16CD4677 CT-WN15CC3385 CT-WN18CC4066 CT-WN15CA9423
    CT-WN15CA4114 CT-WN15CC4010 CT-WN17BA6591 CT-WN16CA1868
    CT-WN18CA7566 CT-WN18CA1178
)

EXCLUDED_FROM_CORE="CT-WN17CD6059"   # contaminated; replaced by CT-WN18CA1178
ATCC_SAMPLE="CT-WNTC9548"
REPLICATES=(
    "CT-WN18CA0775"
    "CT-WN18CA3057"
    "CT-WN18CD7757"
    "CT-WN17CC9963"
)
# ===========================================================================

# Derived paths used by both parts
REF="$WORKDIR/$(basename "$REF_SRC")"
CORE_VCF="$WORKDIR/core/core.vcf.gz"
GWAS_VCF="$WORKDIR/core/core.split.vcf.gz"
CORE_TREE="$WORKDIR/core/core.full.aln.treefile"
SNIPPY_DIR="$WORKDIR/snippy"
AMI_PHENO="$MIC_AMI"
CIP_PHENO="$MIC_CIP"
RIF_PHENO="$MIC_RIF"
KINSHIP="$WORKDIR/phylogeny_kinship_C.tsv"
KINSHIP_LEGACY="$WORKDIR/phylogeny_kinship_legacy.tsv"
CONDA_ENV="$ENV_NAME"

# Resistance rank thresholds (I grouped with R)
RIF_RESIST_RANK=5   # MIC >= 2 mg/L
AMI_RESIST_RANK=6   # MIC >= 32 mg/L
CIP_RESIST_RANK=5   # MIC >= 2 mg/L

# ---------------------------------------------------------------------------
# Activate the environment and check inputs
# ---------------------------------------------------------------------------
echo "============================================================"
echo "env_1 pipeline"
echo "  environment : $ENV_NAME"
echo "  workdir     : $WORKDIR"
echo "  raw reads   : $RAW_READS"
echo "  reference   : $REF_SRC"
echo "  MIC files   : $MIC_AMI, $MIC_CIP, $MIC_RIF"
echo "  outputs to  : $OUTPUT_DEST"
echo "============================================================"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
    if ! command -v conda >/dev/null; then
        echo "ERROR: conda not found. See SETUP.md."; exit 1
    fi
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME" || {
        echo "ERROR: conda environment '$ENV_NAME' not found. Create it as described in SETUP.md."
        exit 1
    }
fi
echo "Environment '$ENV_NAME' active."

for f in "$REF_SRC" "$MIC_AMI" "$MIC_CIP" "$MIC_RIF"; do
    [[ -f "$f" ]] || { echo "ERROR: input file not found: $f"; exit 1; }
done
[[ -d "$RAW_READS" ]] || { echo "ERROR: raw read directory not found: $RAW_READS"; exit 1; }
if ! compgen -G "$RAW_READS/*.fastq.gz" > /dev/null; then
    echo "ERROR: no *.fastq.gz files in $RAW_READS"; exit 1
fi

mkdir -p "$WORKDIR"/{raw_reads,trimmed,qc,snippy,core} "$OUTPUT_DEST"

# Raw reads are read in place; if they are not already under $WORKDIR/raw_reads
# the pipeline below uses $RAW_READS directly.
cp -n "$REF_SRC" "$REF" 2>/dev/null || cp "$REF_SRC" "$REF"

# ###########################################################################
# PART A — trimming, QC, variant calling, core alignment, phylogeny
# (formerly mhaem_full_pipeline.sh)
# ###########################################################################

# Verify Snippy installation
snippy --check >/dev/null

# fix vt normalize bug on WSL (find snippy binary dynamically regardless of mamba install location)
SNIPPY_BIN=$(command -v snippy)
# --follow-symlinks edits the real file if bin/snippy is a symlink
sed -i --follow-symlinks 's/ | vt normalize -r $refdir\/ref.fa -//' "$SNIPPY_BIN"
if grep -q "vt normalize" "$(readlink -f "$SNIPPY_BIN")"; then
    echo "  WARNING: vt normalize step still present in snippy (patch did not match)"
else
    echo "  vt normalize step removed from snippy (report in methods)"
fi

# -----------------------------------------------------------------------------
# Step 3. Create working directories
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 3: Creating working directories"
echo "============================================================"

mkdir -p "$WORKDIR"/{raw_reads,trimmed,qc,snippy,core}

# -----------------------------------------------------------------------------
# Step 4. Check inputs and copy files
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 4: Checking raw FASTQ folder and copying GenBank reference"
echo "============================================================"

# Confirm reference exists
if [[ ! -f "$REF_SRC" ]]; then
    echo "ERROR: Reference file not found:"
    echo "  $REF_SRC"
    echo ""
    echo "Please download GCF_000340435.2.gbff from NCBI and place it in:"
    echo "  E:\\M. haemophilum Project\\WGS Reference\\"
    exit 1
fi

# Copy reference
cp "$REF_SRC" "$REF"

# Raw FASTQs are read in place from $RAW_READS (not copied from E:)
if ! compgen -G "$RAW_READS/*.fastq.gz" > /dev/null; then
    echo "ERROR: No FASTQ files found in:"
    echo "  $RAW_READS"
    echo "Expected e.g. $RAW_READS/CT-WN15CA4114_R1.fastq.gz"
    exit 1
fi

echo "Reference copied; using raw FASTQs in $RAW_READS"

# -----------------------------------------------------------------------------
# Step 5. Rename *_R1_001.fastq.gz and *_R2_001.fastq.gz
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 5: Checking FASTQ file names in $RAW_READS"
echo "============================================================"

# Rename only if any Illumina-style *_R1_001 / *_R2_001 names are present
for f in "$RAW_READS/"*_R1_001.fastq.gz; do
    [[ -e "$f" ]] || continue
    mv "$f" "${f/_R1_001.fastq.gz/_R1.fastq.gz}"
done

for f in "$RAW_READS/"*_R2_001.fastq.gz; do
    [[ -e "$f" ]] || continue
    mv "$f" "${f/_R2_001.fastq.gz/_R2.fastq.gz}"
done

echo "Name check complete."

# Count FASTQ files and paired samples
FASTQ_COUNT=$(find "$RAW_READS" -maxdepth 1 -name "*.fastq.gz" | wc -l)
PAIR_COUNT=$(find "$RAW_READS" -maxdepth 1 -name "*_R1.fastq.gz" | wc -l)
echo "Found $FASTQ_COUNT FASTQ files ($PAIR_COUNT samples with an R1 file) in $RAW_READS."

if [[ "$FASTQ_COUNT" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found."
    exit 1
fi

# -----------------------------------------------------------------------------
# Step 6. Per-sample analysis
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 6: Running fastp, FastQC, and Snippy"
echo "============================================================"

for R1 in "$RAW_READS/"*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="$RAW_READS/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: Missing R2 file for $SAMPLE. Skipping."
        continue
    fi

    echo ""
    echo "========================================================"
    echo "Processing $SAMPLE"
    echo "========================================================"

    # Trim reads (skip if both trimmed files already exist)
    if [[ -f "$WORKDIR/trimmed/${SAMPLE}_R1.trim.fastq.gz" && \
          -f "$WORKDIR/trimmed/${SAMPLE}_R2.trim.fastq.gz" ]]; then
        echo "  SKIP fastp — trimmed reads already exist for $SAMPLE"
    else
        fastp \
            -i "$R1" \
            -I "$R2" \
            -o "$WORKDIR/trimmed/${SAMPLE}_R1.trim.fastq.gz" \
            -O "$WORKDIR/trimmed/${SAMPLE}_R2.trim.fastq.gz" \
            -h "$WORKDIR/qc/${SAMPLE}.fastp.html" \
            -j "$WORKDIR/qc/${SAMPLE}.fastp.json" \
            -w "$THREADS"
    fi

    # FastQC (skip if report already exists for R1)
    if [[ -f "$WORKDIR/qc/${SAMPLE}_R1.trim_fastqc.html" ]]; then
        echo "  SKIP FastQC — report already exists for $SAMPLE"
    else
        fastqc \
            -t "$THREADS" \
            -o "$WORKDIR/qc" \
            "$WORKDIR/trimmed/${SAMPLE}_R1.trim.fastq.gz" \
            "$WORKDIR/trimmed/${SAMPLE}_R2.trim.fastq.gz"
    fi

    # Mapping and SNP calling with Snippy (skip if snps.vcf already exists)
    if [[ -f "$WORKDIR/snippy/$SAMPLE/snps.vcf" ]]; then
        echo "  SKIP Snippy — output already exists for $SAMPLE"
    else
        # --force: overwrite a partial output folder left by an interrupted run
        snippy \
            --force \
            --cpus "$THREADS" \
            --outdir "$WORKDIR/snippy/$SAMPLE" \
            --ref "$REF" \
            --R1 "$WORKDIR/trimmed/${SAMPLE}_R1.trim.fastq.gz" \
            --R2 "$WORKDIR/trimmed/${SAMPLE}_R2.trim.fastq.gz"
    fi

done

# -----------------------------------------------------------------------------
# Step 7. MultiQC (excluding the contaminated isolate)
#
# The report is written to $WORKDIR/qc_61 rather than $WORKDIR/qc so that
# MultiQC never scans its own previous output, and so the 62-library report
# from earlier runs is not silently reused.
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 7: Generating MultiQC report (excluding $EXCLUDED_FROM_CORE)"
echo "============================================================"

MULTIQC_OUT="$WORKDIR/qc_61"

if [[ -f "$MULTIQC_OUT/multiqc_report.html" && "$MULTIQC_OUT/multiqc_report.html" -nt "$WORKDIR/qc" ]]; then
    echo "SKIP MultiQC — report is up to date: $MULTIQC_OUT/multiqc_report.html"
else
    multiqc \
        -f \
        -o "$MULTIQC_OUT" \
        --ignore-samples "${EXCLUDED_FROM_CORE}*" \
        "$WORKDIR/qc"
fi

# Confirm the exclusion actually took effect
if grep -q "$EXCLUDED_FROM_CORE" "$MULTIQC_OUT/multiqc_report.html"; then
    echo "  WARNING: $EXCLUDED_FROM_CORE is still present in the MultiQC report"
else
    echo "  MultiQC report excludes $EXCLUDED_FROM_CORE"
fi

# -----------------------------------------------------------------------------
# Step 7b. fastp QC summary tables
#
# Writes per-sample metrics and a summary for (i) all retained libraries and
# (ii) the clinical isolates only (replicates, the QC strain and the
# contaminated isolate excluded). The clinical summary is the one quoted in
# the manuscript.
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 7b: fastp QC summary tables"
echo "============================================================"

REPLICATES_STR=$(IFS=","; echo "${REPLICATES[*]}")

python3 - "$WORKDIR/qc" "$REPLICATES_STR" \
        "$ATCC_SAMPLE" "$EXCLUDED_FROM_CORE" \
        "$WORKDIR/qc_per_sample.tsv" "$WORKDIR/qc_summary.tsv" << 'QC_SUMMARY_EOF'
import json, glob, os, sys, statistics as st

qc_dir, rep_str, atcc, excluded, out_rows, out_summary = sys.argv[1:7]
replicates = set(r for r in rep_str.split(",") if r)

rows = []
for f in sorted(glob.glob(os.path.join(qc_dir, "*.fastp.json"))):
    s = os.path.basename(f)[:-len(".fastp.json")]
    if s == excluded:
        continue                      # contaminated isolate, not analysed
    d = json.load(open(f))
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    group = ("QC strain" if s == atcc else
             "Replicate" if s in replicates else "Clinical")
    rows.append({
        "Sample": s, "Group": group,
        "Raw_reads": b["total_reads"], "Trimmed_reads": a["total_reads"],
        "Pass_rate_pct": round(100 * a["total_reads"] / b["total_reads"], 2)
                         if b["total_reads"] else "NA",
        "Q30_before_pct": round(100 * b["q30_rate"], 2),
        "Q30_after_pct": round(100 * a["q30_rate"], 2),
        "GC_after_pct": round(100 * a["gc_content"], 2),
        "Mean_len_R1": a.get("read1_mean_length", "NA"),
        "Mean_len_R2": a.get("read2_mean_length", "NA"),
        "Duplication_pct": round(100 * d.get("duplication", {}).get("rate", float("nan")), 2),
    })

if not rows:
    sys.exit(f"ERROR: no *.fastp.json files found in {qc_dir}")

cols = list(rows[0])
with open(out_rows, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

METRICS = [("Trimmed_reads", "Trimmed reads (R1 + R2)", 1e6, "{:.2f}", " M"),
           ("Pass_rate_pct", "Pass rate", 1, "{:.1f}", " %"),
           ("Q30_before_pct", "Q30 before trimming", 1, "{:.1f}", " %"),
           ("Q30_after_pct", "Q30 after trimming", 1, "{:.1f}", " %"),
           ("GC_after_pct", "GC after trimming", 1, "{:.1f}", " %"),
           ("Duplication_pct", "Duplication", 1, "{:.1f}", " %")]

groups = [("All retained libraries", rows),
          ("Clinical isolates", [r for r in rows if r["Group"] == "Clinical"])]

with open(out_summary, "w") as fh:
    fh.write("Group\tn\tMetric\tMean\tMedian\tMin\tMax\n")
    for label, sel in groups:
        print(f"  {label} (n = {len(sel)})")
        for key, name, scale, fmt, unit in METRICS:
            v = [r[key] / scale for r in sel if isinstance(r[key], (int, float)) and r[key] == r[key]]
            if not v:
                continue
            fh.write("\t".join([label, str(len(sel)), name,
                                fmt.format(st.mean(v)), fmt.format(st.median(v)),
                                fmt.format(min(v)), fmt.format(max(v))]) + "\n")
            print(f"    {name:24s} mean {fmt.format(st.mean(v))}{unit}, "
                  f"range {fmt.format(min(v))}–{fmt.format(max(v))}{unit}")

print(f"  Per-sample table : {out_rows}")
print(f"  Summary table    : {out_summary}")
QC_SUMMARY_EOF

# -----------------------------------------------------------------------------
# Step 8. Core SNP alignment
# Excludes CT-WN17CD6059 (contaminated assembly; replaced by CT-WN18CA1178).
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 8: Running snippy-core (excluding $EXCLUDED_FROM_CORE)"
echo "============================================================"

# Build list of snippy directories, excluding the contaminated sample
SNIPPY_DIRS=()
for d in "$WORKDIR/snippy"/*/; do
    sample=$(basename "$d")
    if [[ "$sample" != "$EXCLUDED_FROM_CORE" ]]; then
        SNIPPY_DIRS+=("$d")
    fi
done

echo "  Including ${#SNIPPY_DIRS[@]} sample(s) in snippy-core (excluded: $EXCLUDED_FROM_CORE)"

if [[ -f "$WORKDIR/core/core.full.aln" ]]; then
    echo "SKIP snippy-core — core.full.aln already exists"
else
    snippy-core \
        --ref "$REF" \
        --prefix "$WORKDIR/core/core" \
        "${SNIPPY_DIRS[@]}"
fi

# bgzip + tabix index the core VCF (required by gwas11.sh)
if [[ -f "$WORKDIR/core/core.vcf.gz" ]]; then
    echo "SKIP bgzip — core.vcf.gz already exists"
else
    bgzip -c "$WORKDIR/core/core.vcf" \
        > "$WORKDIR/core/core.vcf.gz"
    tabix -p vcf "$WORKDIR/core/core.vcf.gz"
fi

# -----------------------------------------------------------------------------
# Step 9. Phylogenetic tree
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 9: Building phylogenetic tree with IQ-TREE"
echo "============================================================"

if [[ -f "$WORKDIR/core/core.full.aln.treefile" ]]; then
    echo "SKIP IQ-TREE — treefile already exists"
else
    iqtree \
        -s "$WORKDIR/core/core.full.aln" \
        -m GTR+G \
        -bb 1000 \
        -nt AUTO
fi

# -----------------------------------------------------------------------------
# Step 10. Copy final outputs to E drive
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 10: Copying final outputs to E drive"
echo "============================================================"

mkdir -p "$OUTPUT_DEST"

cp "$MULTIQC_OUT/multiqc_report.html"     "$OUTPUT_DEST/"
cp "$WORKDIR/qc_per_sample.tsv"           "$OUTPUT_DEST/"
cp "$WORKDIR/qc_summary.tsv"              "$OUTPUT_DEST/"
cp "$WORKDIR/core/core.full.aln"          "$OUTPUT_DEST/"
cp "$WORKDIR/core/core.full.aln.treefile" "$OUTPUT_DEST/"

# Optional: copy all per-sample Snippy results
# cp -r "$WORKDIR/snippy" "$OUTPUT_DEST/"

# -----------------------------------------------------------------------------
# Step 11. Extract resistance mutations
# -----------------------------------------------------------------------------
echo "============================================================"
echo "[PART A] STEP 11: Extracting resistance-associated mutations"
echo "============================================================"

#WORKDIR="$HOME/mhaem_wgs"

grep -E "B586_RS03455|B586_RS00030|B586_RS00025|B586_RS14090" \
    $(ls "$WORKDIR/snippy"/*/snps.tab | grep -v CT-WN17CD6059) \
    > "$WORKDIR/resistance_mutations_2.tsv" || true

echo "  $(wc -l < "$WORKDIR/resistance_mutations_2.tsv") variant line(s) in target genes"
cp "$WORKDIR/resistance_mutations_2.tsv" "$OUTPUT_DEST/"

# -----------------------------------------------------------------------------
# Finished
# -----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "============================================================"
echo "Results are located in:"
echo "  E:\\M. haemophilum Project\\Jon's Folder\\WGS Analysis Results"
echo "============================================================"

# ###########################################################################
# PART B — SNP-based GWAS and targeted analyses
# (formerly gwas11.sh)
# ###########################################################################

# =============================================================================
# STEP 2: Verify input files and validate phenotype sample names
# =============================================================================

echo "============================================================"
echo "[PART B] STEP 2: Verifying input files and sample name consistency"
echo "============================================================"

REQUIRED_FILES=(
    "$CORE_VCF"
    "$CORE_TREE"
    "$AMI_PHENO"
    "$CIP_PHENO"
    "$RIF_PHENO"
)

for FILE in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "$FILE" ]]; then
        echo "ERROR: Missing required file:"
        echo "  $FILE"
        exit 1
    fi
done

echo "All required input files found."
echo ""

# Extract and display VCF sample names
VCF_SAMPLES=$(bcftools query -l "$CORE_VCF")
echo "Samples detected in VCF:"
echo "$VCF_SAMPLES"
echo ""

# Warn if any VCF sample is absent from a phenotype file
# (pyseer drops missing samples silently — better to catch this early)
check_sample_names() {
    local PHENO_FILE="$1"
    local DRUG="$2"

    PHENO_SAMPLES=$(awk 'NR>1 { print $1 }' "$PHENO_FILE")

    MISSING=0
    while IFS= read -r SAMPLE; do
        # ATCC 29548 and the reference have no MIC by design — not a problem
        [[ "$SAMPLE" == "CT-WNTC9548" || "$SAMPLE" == "Reference" ]] && continue
        if ! echo "$PHENO_SAMPLES" | grep -qx "$SAMPLE"; then
            echo "  WARNING [$DRUG]: VCF sample '$SAMPLE' not found in $PHENO_FILE"
            MISSING=$((MISSING + 1))
        fi
    done <<< "$VCF_SAMPLES"

    if [[ "$MISSING" -eq 0 ]]; then
        echo "  [$DRUG] All VCF samples present in phenotype file (ATCC 29548 excluded by design). OK."
    else
        echo "  [$DRUG] $MISSING sample(s) missing — pyseer will silently drop them."
        echo "          Confirm this is intended before proceeding."
    fi
}

check_sample_names "$AMI_PHENO" "AMI"
check_sample_names "$CIP_PHENO" "CIP"
check_sample_names "$RIF_PHENO" "RIF"
echo ""

# Validate that trait column contains only numeric values
# (catches accidental '>64', '<0.12', or text left in the file)
validate_numeric() {
    local PHENO_FILE="$1"
    local DRUG="$2"

    # Strip carriage returns and normalise whitespace before checking.
    # Handles Windows-style line endings (CRLF) and mixed tab/space delimiters
    # which cause awk to mis-parse the trait column as non-numeric.
    INVALID=$(sed 's/\r//' "$PHENO_FILE" | awk 'NR>1 {
        gsub(/[ \t]+/, "\t");
        n = split($0, a, "\t");
        val = a[2];
        gsub(/^[ \t]+|[ \t]+$/, "", val);
        if (val !~ /^[0-9]+(\.[0-9]+)?$/) print NR": "$0
    }')
    if [[ -n "$INVALID" ]]; then
        echo "ERROR [$DRUG]: Non-numeric values found in trait column of $PHENO_FILE:"
        echo "$INVALID"
        echo "All trait values must be numeric ordinal ranks (e.g. 1, 2, 3...)."
        echo "Remove any '>', '<', or text values before running."
        exit 1
    else
        echo "  [$DRUG] Trait column is numeric. OK."
    fi
}

validate_numeric "$AMI_PHENO" "AMI"
validate_numeric "$CIP_PHENO" "CIP"
validate_numeric "$RIF_PHENO" "RIF"

echo ""
echo "Input validation complete."

# =============================================================================
# STEP 2a: Clinical-only phenotype files (56 isolates, one per patient)
# =============================================================================
# Writes $WORKDIR/pheno_clinical_<DRUG>.txt containing only the CLINICAL
# isolates (header "samples<TAB>trait", Windows line endings stripped).
# From here on, AMI_PHENO / CIP_PHENO / RIF_PHENO point to these files, so
# pyseer (af, LMM), the 2x2 tables (Step 7) and the candidate mutation tables
# (Step 9) all use the same 56 isolates.
# A file is only rewritten when its content changes, so unchanged inputs do
# not trigger reruns.

echo "============================================================"
echo "[PART B] STEP 2a: Writing clinical-only phenotype files (n=${#CLINICAL[@]})"
echo "============================================================"

CLINICAL_LIST="$WORKDIR/clinical_samples.txt"
printf '%s\n' "${CLINICAL[@]}" > "$CLINICAL_LIST.tmp"
if [[ -f "$CLINICAL_LIST" ]] && cmp -s "$CLINICAL_LIST.tmp" "$CLINICAL_LIST"; then
    rm -f "$CLINICAL_LIST.tmp"
else
    mv "$CLINICAL_LIST.tmp" "$CLINICAL_LIST"
fi

make_clinical_pheno() {
    local SRC="$1" DRUG="$2"
    local DEST="$WORKDIR/pheno_clinical_${DRUG}.txt"
    python3 - "$SRC" "$CLINICAL_LIST" "$DEST.tmp" "$DRUG" << 'PHENO_EOF'
import sys
src, clin_path, dest, drug = sys.argv[1:5]
clinical = [l.strip() for l in open(clin_path) if l.strip()]
clin_set = set(clinical)
kept, dropped, seen = [], [], set()
with open(src) as f:
    next(f)                                   # header
    for line in f:
        line = line.replace("\r", "").strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        s, v = parts[0], parts[1]
        if s in clin_set:
            if s in seen:
                print(f"  WARNING [{drug}]: duplicate entry for {s} — keeping first")
                continue
            seen.add(s)
            kept.append((s, v))
        else:
            dropped.append(s)
with open(dest, "w") as out:
    out.write("samples\ttrait\n")
    for s, v in kept:
        out.write(f"{s}\t{v}\n")
missing = [s for s in clinical if s not in seen]
print(f"  [{drug}] kept {len(kept)} clinical isolates; dropped {len(dropped)}: "
      f"{', '.join(dropped) if dropped else 'none'}")
if missing:
    print(f"  WARNING [{drug}]: {len(missing)} clinical isolate(s) have no MIC: {', '.join(missing)}")
PHENO_EOF
    if [[ -f "$DEST" ]] && cmp -s "$DEST.tmp" "$DEST"; then
        rm -f "$DEST.tmp"
    else
        mv "$DEST.tmp" "$DEST"
        echo "  [$DRUG] $DEST updated"
    fi
}

make_clinical_pheno "$AMI_PHENO" AMI
make_clinical_pheno "$CIP_PHENO" CIP
make_clinical_pheno "$RIF_PHENO" RIF

AMI_PHENO="$WORKDIR/pheno_clinical_AMI.txt"
CIP_PHENO="$WORKDIR/pheno_clinical_CIP.txt"
RIF_PHENO="$WORKDIR/pheno_clinical_RIF.txt"

# Clinical isolates that are missing from the core VCF (would be silently dropped)
VCF_SAMPLES_FOR_CHECK=$(bcftools query -l "$CORE_VCF")
N_MISSING_VCF=0
for S in "${CLINICAL[@]}"; do
    if ! grep -qx "$S" <<< "$VCF_SAMPLES_FOR_CHECK"; then
        echo "  WARNING: clinical isolate $S is not in $CORE_VCF"
        N_MISSING_VCF=$((N_MISSING_VCF + 1))
    fi
done
[[ "$N_MISSING_VCF" -eq 0 ]] && echo "  All ${#CLINICAL[@]} clinical isolates present in core VCF. OK."

# =============================================================================
# STEP 2b: Snippy QC table (56 clinical isolates from unique patients)
# Columns: Sample | Mean_depth | Mapping_rate_pct | Ref_coverage_pct
#          | N_SNPs | N_indels | Trimmed_reads | QC_flag
#
# Reads per-isolate snps.bam via samtools; SNP/indel counts from snps.tab.
# Pass/fail criteria:
#   FAIL if: mean depth < 10x, OR mapping rate < 80%, OR ref coverage < 85%
# =============================================================================

echo "============================================================"
echo "[PART B] STEP 2b: Snippy QC table"
echo "============================================================"

SNIPPY_QC="$WORKDIR/snippy_qc.tsv"
QC_DIR="$HOME/mhaem_wgs/qc"

if [ -f "$SNIPPY_QC" ]; then
    echo "SKIP snippy QC table (already exists: $SNIPPY_QC)"
else

# Derive clinical samples: VCF samples minus replicates and ATCC
REPLICATES_STR=$(IFS=","; echo "${REPLICATES[*]}")

python3 << STEP2B_EOF
import subprocess, os, sys, json, csv, math

SNIPPY_DIR    = "$SNIPPY_DIR"
WORKDIR       = "$WORKDIR"
CORE_VCF      = "$CORE_VCF"
QC_DIR        = "$QC_DIR"
QC_TSV        = "$SNIPPY_QC"
REPLICATES    = set(r.strip() for r in "$REPLICATES_STR".split(",") if r.strip())
ATCC_SAMPLE   = "CT-WNTC9548"

# Get all samples from VCF
try:
    vcf_samples = subprocess.check_output(
        ["bcftools", "query", "-l", CORE_VCF],
        stderr=subprocess.DEVNULL
    ).decode().strip().split("\n")
except subprocess.CalledProcessError:
    print("ERROR: could not query VCF sample names")
    sys.exit(1)

# Clinical = VCF samples minus replicates and ATCC
clinical = [s for s in vcf_samples if s not in REPLICATES and s != ATCC_SAMPLE]
print(f"  Clinical isolates: {len(clinical)}")

# ── samtools stats from snps.bam ──────────────────────────────────────────────
def bam_stats(sample):
    bam = os.path.join(SNIPPY_DIR, sample, "snps.bam")
    if not os.path.isfile(bam):
        return {"mean_depth": "NA", "mapping_rate": "NA", "ref_coverage": "NA"}

    # Mean depth
    try:
        depth_out = subprocess.check_output(
            ["samtools", "depth", "-a", bam],
            stderr=subprocess.DEVNULL
        ).decode()
        depths = [int(l.split("\t")[2]) for l in depth_out.strip().split("\n") if l]
        mean_depth = round(sum(depths) / len(depths), 2) if depths else "NA"
    except Exception:
        mean_depth = "NA"

    # Mapping rate from flagstat
    try:
        flag_out = subprocess.check_output(
            ["samtools", "flagstat", bam],
            stderr=subprocess.DEVNULL
        ).decode()
        total = mapped = None
        for line in flag_out.split("\n"):
            if "in total" in line:
                total = int(line.split()[0])
            if "primary mapped" in line and "mate mapped" not in line:
                mapped = int(line.split()[0])
        mapping_rate = round(mapped / total * 100, 2) if total and mapped else "NA"
    except Exception:
        mapping_rate = "NA"

    # Reference coverage from samtools coverage
    try:
        cov_out = subprocess.check_output(
            ["samtools", "coverage", bam],
            stderr=subprocess.DEVNULL
        ).decode()
        lines = [l for l in cov_out.strip().split("\n") if not l.startswith("#")]
        if lines:
            cols = lines[0].split("\t")
            ref_coverage = round(float(cols[5]), 2)  # coverage column
        else:
            ref_coverage = "NA"
    except Exception:
        ref_coverage = "NA"

    return {"mean_depth": mean_depth, "mapping_rate": mapping_rate, "ref_coverage": ref_coverage}

# ── SNP/indel counts from snps.tab ───────────────────────────────────────────
def snp_counts(sample):
    tab = os.path.join(SNIPPY_DIR, sample, "snps.tab")
    if not os.path.isfile(tab):
        return "NA", "NA"
    n_snp = n_indel = 0
    try:
        with open(tab) as f:
            next(f)  # skip header
            for line in f:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                vtype = parts[2].strip().lower()
                if vtype == "snp":
                    n_snp += 1
                elif vtype in ("ins", "del", "complex", "mnp"):
                    n_indel += 1
    except OSError:
        return "NA", "NA"
    return n_snp, n_indel

# ── Trimmed reads from fastp JSON ─────────────────────────────────────────────
# Pattern: {QC_DIR}/{sample}.fastp.json
def trimmed_reads(sample):
    path = os.path.join(QC_DIR, f"{sample}.fastp.json")
    if not os.path.isfile(path):
        return "NA"
    try:
        with open(path) as f:
            d = json.load(f)
        reads = (d.get("filtering_result", {}).get("passed_filter_reads") or
                 d.get("summary", {}).get("after_filtering", {}).get("total_reads"))
        return reads if reads is not None else "NA"
    except Exception:
        return "NA"

# ── Pass/fail ─────────────────────────────────────────────────────────────────
def qc_flag(depth, mapping, coverage):
    reasons = []
    try:
        if float(depth) < 10:
            reasons.append("depth<10x")
    except (ValueError, TypeError):
        pass
    try:
        if float(mapping) < 80:
            reasons.append("mapping<80%")
    except (ValueError, TypeError):
        pass
    try:
        if float(coverage) < 85:
            reasons.append("coverage<85%")
    except (ValueError, TypeError):
        pass
    return "FAIL:" + ",".join(reasons) if reasons else "PASS"

# ── Build table ───────────────────────────────────────────────────────────────
OUT_COLS = ["Sample", "Mean_depth_x", "Mapping_rate_pct", "Ref_coverage_pct",
            "N_SNPs", "N_indels", "Trimmed_reads", "QC_flag"]

rows = []
for sample in clinical:
    stats   = bam_stats(sample)
    n_snp, n_indel = snp_counts(sample)
    reads   = trimmed_reads(sample)
    flag    = qc_flag(stats["mean_depth"], stats["mapping_rate"], stats["ref_coverage"])
    rows.append({
        "Sample":           sample,
        "Mean_depth_x":     stats["mean_depth"],
        "Mapping_rate_pct": stats["mapping_rate"],
        "Ref_coverage_pct": stats["ref_coverage"],
        "N_SNPs":           n_snp,
        "N_indels":         n_indel,
        "Trimmed_reads":    reads,
        "QC_flag":          flag,
    })
    print(f"  {sample}: depth={stats['mean_depth']}x  map={stats['mapping_rate']}%  cov={stats['ref_coverage']}%  SNPs={n_snp}  {flag}")

with open(QC_TSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

n_pass = sum(1 for r in rows if r["QC_flag"] == "PASS")
n_fail = sum(1 for r in rows if r["QC_flag"].startswith("FAIL"))
print(f"\n  QC table written: {QC_TSV}")
print(f"  PASS: {n_pass}  FAIL: {n_fail}")
STEP2B_EOF

fi  # end skip

echo "Step 2b complete — Snippy QC table: $SNIPPY_QC"


# =============================================================================
# STEP 3: Phylogenetic kinship matrix (pyseer phylogeny_distance.py --lmm)
# =============================================================================
# Standard pyseer method for a phylogeny-based LMM:
#   - the IQ-TREE core SNP tree is midpoint-rooted (--midpoint)
#   - for each pair of isolates, similarity = branch length from the root to
#     their most recent common ancestor (MRCA); diagonal = root-to-tip length
#   - this is the Brownian-motion variance-covariance matrix C, which is
#     positive semi-definite by construction (a valid LMM kinship matrix)
#
# The Python script written below is pyseer's scripts/phylogeny_distance.py,
# reproduced verbatim (pyseer 1.4.x; Galardini & Lees; Apache-2.0) so the run
# does not depend on locating the script inside the conda env.
# Requires DendroPy and pandas (both installed as pyseer dependencies).
#
# Methods wording:
#   "Population structure was controlled with a linear mixed model using a
#    phylogenetic similarity matrix generated with pyseer's
#    phylogeny_distance.py (--lmm --midpoint) from the midpoint-rooted
#    core-genome SNP maximum-likelihood tree."

echo "============================================================"
echo "[PART B] STEP 3: Phylogenetic kinship matrix (pyseer phylogeny_distance.py --lmm)"
echo "============================================================"

echo "  Tree:   $CORE_TREE"
echo "  Output: $KINSHIP"

PHYLO_SCRIPT="$WORKDIR/phylogeny_distance.py"

cat > "$PHYLO_SCRIPT" << 'PHYLO_SCRIPT_EOF'
#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Extract a distance matrix from a phylogeny'''


def get_options():
    import argparse

    description = 'Extract a distance matrix from a phylogeny'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('phylogeny',
                        help='Tree file')

    parser.add_argument('--format',
                        default="newick",
                        help="Format of tree file [Default: newick]")
    parser.add_argument('--midpoint',
                        action='store_true',
                        default=False,
                        help='Midpoint root the tree before calculating '
                             'distances.')
    method_group = parser.add_mutually_exclusive_group()
    method_group.add_argument('--lmm', '--calc-C',
                              action='store_true',
                              help='Produce var-covar matrix C (as from PDDIST). '
                                   'Always uses branch lengths.')
    method_group.add_argument('--topology',
                              action='store_true',
                              default=False,
                              help='Ignore branch lengths, and only use topological '
                                   'distances')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import sys
    import re
    import pandas as pd
    import dendropy

    tree = dendropy.Tree.get(
        path=options.phylogeny,
        schema=options.format,
        preserve_underscores=True)

    if options.midpoint:
        tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)

    d = {}
    pdm = tree.phylogenetic_distance_matrix()
    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        d[taxon1.label] = d.get(taxon1.label, {})
        for taxon2 in tree.taxon_namespace:
            if taxon2.label not in d[taxon1.label].keys():
                if options.lmm:
                    mrca = pdm.mrca(taxon1, taxon2)
                    d[taxon1.label][taxon2.label] = mrca.distance_from_root()
                elif options.topology:
                    d[taxon1.label][taxon2.label] = pdm.path_edge_count(taxon1, taxon2)
                else:
                    d[taxon1.label][taxon2.label] = pdm.patristic_distance(taxon1, taxon2)

    m = pd.DataFrame(d)
    m = m.reindex(m.columns)
    m.to_csv(sys.stdout,
             sep='\t')
PHYLO_SCRIPT_EOF

if [[ -s "$KINSHIP" ]]; then
    echo "  SKIP — kinship matrix already exists: $KINSHIP"
else
    if ! python3 "$PHYLO_SCRIPT" --lmm --midpoint "$CORE_TREE" > "$KINSHIP"; then
        rm -f "$KINSHIP"
        echo "ERROR: phylogeny_distance.py failed (see traceback above)."
        echo "       Check DendroPy is installed and that $CORE_TREE is a valid Newick tree."
        exit 1
    fi
fi

if [[ ! -s "$KINSHIP" ]]; then
    echo "ERROR: Kinship matrix was not created or is empty."
    echo "       Check that DendroPy is installed: mamba install -c conda-forge dendropy"
    exit 1
fi

# Sanity checks: square, symmetric, positive semi-definite; list samples
python3 - "$KINSHIP" << 'KCHECK_EOF'
import sys
import numpy as np
import pandas as pd
K = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
print(f"  Kinship matrix: {K.shape[0]} rows x {K.shape[1]} columns")
if list(K.index) != list(K.columns):
    print("  WARNING: row and column labels differ — check tree tip labels")
v = K.values.astype(float)
if not np.allclose(v, v.T):
    print("  WARNING: matrix is not symmetric")
ev = np.linalg.eigvalsh((v + v.T) / 2)
print(f"  Smallest eigenvalue: {ev.min():.3e} "
      f"({'positive semi-definite — OK' if ev.min() > -1e-8 else 'NOT positive semi-definite — check tree'})")
print(f"  Tips: {', '.join(map(str, K.index))}")
KCHECK_EOF

echo "Kinship matrix ready."

# =============================================================================
# STEP 3b: Legacy kinship matrix (sensitivity analysis only)
# =============================================================================
# The matrix used in earlier runs: similarity = 1 - patristic distance / max.
# Not used for the main GWAS. It is rebuilt here only so Step 6c can compare
# results between the two kinship constructions.

echo "============================================================"
echo "[PART B] STEP 3b: Legacy kinship matrix (for sensitivity comparison)"
echo "============================================================"

if [[ -s "$KINSHIP_LEGACY" ]]; then
    echo "  SKIP — legacy kinship matrix already exists: $KINSHIP_LEGACY"
else
LEGACY_KINSHIP_SCRIPT="$WORKDIR/make_kinship_legacy.py"

cat > "$LEGACY_KINSHIP_SCRIPT" << 'KINSHIP_SCRIPT_EOF'
import sys
import os
import numpy as np

TREE_FILE = sys.argv[1]
KINSHIP_FILE = sys.argv[2]

tree = None
samples = []

# --- Attempt 1: ETE3 ---
try:
    from ete3 import Tree
    t = Tree(TREE_FILE)
    samples = [leaf.name for leaf in t.get_leaves()]
    n = len(samples)
    dist_matrix = np.zeros((n, n))
    leaves = t.get_leaves()
    for i, leaf1 in enumerate(leaves):
        for j, leaf2 in enumerate(leaves):
            if i < j:
                d = t.get_distance(leaf1, leaf2)
                dist_matrix[i][j] = d
                dist_matrix[j][i] = d
    print(f"  [ETE3] Tree loaded: {n} samples")
    tree = "loaded"
except ImportError:
    print("  ETE3 not available, trying Biopython...")

# --- Attempt 2: Biopython ---
if tree is None:
    try:
        from Bio import Phylo
        with open(TREE_FILE) as f:
            t = list(Phylo.parse(f, "newick"))[0]
        samples = [c.name for c in t.get_terminals()]
        n = len(samples)
        dist_matrix = np.zeros((n, n))
        for i, s1 in enumerate(samples):
            for j, s2 in enumerate(samples):
                if i < j:
                    d = t.distance(s1, s2)
                    dist_matrix[i][j] = d
                    dist_matrix[j][i] = d
        print(f"  [Biopython] Tree loaded: {n} samples")
        tree = "loaded"
    except ImportError:
        print("  Biopython not available, trying DendroPy...")

# --- Attempt 3: DendroPy ---
if tree is None:
    try:
        import dendropy
        t = dendropy.Tree.get(path=TREE_FILE, schema="newick")
        pdc = t.phylogenetic_distance_matrix()
        taxa = list(pdc.taxon_iter())
        samples = [taxon.label for taxon in taxa]
        n = len(samples)
        dist_matrix = np.zeros((n, n))
        for i, t1 in enumerate(taxa):
            for j, t2 in enumerate(taxa):
                if i < j:
                    d = pdc.distance(t1, t2)
                    dist_matrix[i][j] = d
                    dist_matrix[j][i] = d
        print(f"  [DendroPy] Tree loaded: {n} samples")
        tree = "loaded"
    except ImportError:
        pass

if tree is None:
    print("ERROR: No tree-reading library found (ete3, biopython, dendropy).")
    print("Install one with: mamba install -c conda-forge ete3")
    sys.exit(1)

# Convert distance matrix to kinship/similarity matrix
# Follows pyseer phylogeny_distance.py --lmm:
#   similarity = 1 - (distance / max_distance)
max_dist = dist_matrix.max()
if max_dist == 0:
    print("ERROR: All pairwise distances are zero — check your tree file.")
    sys.exit(1)

similarity = 1 - (dist_matrix / max_dist)

# Write pyseer-compatible tab-separated similarity matrix
with open(KINSHIP_FILE, 'w') as out:
    out.write('\t' + '\t'.join(samples) + '\n')
    for i, s in enumerate(samples):
        row = '	'.join(f"{similarity[i][j]:.6f}" for j in range(len(samples)))
        out.write(f"{s}\t{row}\n")

print(f"  Kinship matrix written: {KINSHIP_FILE}")
print(f"  Dimensions: {len(samples)} x {len(samples)}")
KINSHIP_SCRIPT_EOF

python3 "$LEGACY_KINSHIP_SCRIPT" "$CORE_TREE" "$KINSHIP_LEGACY"

if [[ ! -s "$KINSHIP_LEGACY" ]]; then
    echo "ERROR: Legacy kinship matrix was not created."
    exit 1
fi
fi  # end legacy kinship skip

# =============================================================================
# STEP 3c: Split multi-allelic sites into biallelic records
# =============================================================================
# snippy-core can report more than one ALT allele at a site (e.g. gyrA 7643
# T>G and T>C). Splitting gives one record per ALT allele, so that:
#   - each allele is tested separately by pyseer,
#   - the pyseer variant ID is always CHROM_POS_REF_ALT (parsed in Step 7),
#   - carriers of every ALT allele are genotype 1 in their own record.

echo "============================================================"
echo "[PART B] STEP 3c: Splitting multi-allelic sites (bcftools norm -m -any)"
echo "============================================================"

if [[ -s "$GWAS_VCF" && "$GWAS_VCF" -nt "$CORE_VCF" ]]; then
    echo "  SKIP — split VCF is up to date: $GWAS_VCF"
else
    bcftools norm -m -any -Oz -o "$GWAS_VCF" "$CORE_VCF"
    bcftools index -f -t "$GWAS_VCF"
fi
N_REC_IN=$(bcftools view -H "$CORE_VCF" | wc -l)
N_MULTI=$(bcftools view -H -m3 "$CORE_VCF" | wc -l)
N_REC_OUT=$(bcftools view -H "$GWAS_VCF" | wc -l)
echo "  $N_MULTI multi-allelic site(s) split: $N_REC_IN -> $N_REC_OUT records"

# =============================================================================
# STEPS 4-6: pyseer LMM GWAS — Amikacin, Ciprofloxacin, Rifampicin
# =============================================================================
# --lmm             : linear mixed model for population structure correction.
# --similarity      : phylogenetic kinship matrix C (Step 3).
# --continuous      : treat ordinal MIC ranks as a continuous phenotype.
# --min-af 0.02     : exclude variants in fewer than 2 of 56 isolates
#                     (i.e. singletons are not tested).
# Phenotypes        : clinical-only files from Step 2a (n = 56), so pyseer's
#                     af is a frequency among the 56 clinical isolates.
# --output-patterns : variant patterns for the Bonferroni threshold (Step 6b).
#
# pyseer stderr (including the h^2 estimate) is saved to pyseer_<DRUG>.log.
# Results are rerun automatically if they are older than the kinship matrix;
# when that happens, downstream hit tables and Manhattan plots are removed so
# they are regenerated from the new results.

run_pyseer() {
    local DRUG="$1" PHENO="$2" K="$3" OUT="$4" LOG="$5" PATTERNS="${6:-}"

    if [[ -s "$OUT" && "$OUT" -nt "$K" && "$OUT" -nt "$PHENO" && "$OUT" -nt "$GWAS_VCF" ]]; then
        echo "  SKIP $DRUG — $(basename "$OUT") is up to date"
        return 0
    fi
    if [[ -f "$OUT" ]]; then
        echo "  $DRUG: $(basename "$OUT") is older than the kinship matrix, phenotype file or VCF — rerunning"
    fi

    local EXTRA=()
    if [[ -n "$PATTERNS" ]]; then
        EXTRA=(--output-patterns "$PATTERNS")
        # main analysis rerun -> clear downstream outputs that depend on it
        rm -f "$WORKDIR/gwas_hits_${DRUG}.tsv" \
              "$WORKDIR/manhattan_${DRUG}.png" \
              "$WORKDIR/manhattan_combined.png" \
              "$WORKDIR/qq_${DRUG}.png" \
              "$WORKDIR/qq_plots.png" \
              "$WORKDIR/genomic_inflation.tsv"
    fi

    echo "  Running pyseer for $DRUG (kinship: $(basename "$K")) ..."
    if ! pyseer \
            --lmm \
            --phenotypes "$PHENO" \
            --vcf "$GWAS_VCF" \
            --similarity "$K" \
            --continuous \
            --min-af 0.02 \
            ${EXTRA[@]+"${EXTRA[@]}"} \
            > "$OUT" 2> "$LOG"; then
        echo "  ERROR: pyseer failed for $DRUG — last lines of $LOG:"
        tail -20 "$LOG"
        exit 1
    fi
    echo "  $DRUG: $(grep -m1 -F 'h^2' "$LOG" || echo 'h^2 not reported')"
}

echo "============================================================"
echo "STEPS 4-6: Running pyseer LMM GWAS (AMI, CIP, RIF)"
echo "============================================================"

run_pyseer AMI "$AMI_PHENO" "$KINSHIP" \
    "$WORKDIR/pyseer_results_AMI.txt" "$WORKDIR/pyseer_AMI.log" "$WORKDIR/patterns_AMI.txt"
echo "Amikacin GWAS complete."

run_pyseer CIP "$CIP_PHENO" "$KINSHIP" \
    "$WORKDIR/pyseer_results_CIP.txt" "$WORKDIR/pyseer_CIP.log" "$WORKDIR/patterns_CIP.txt"
echo "Ciprofloxacin GWAS complete."

run_pyseer RIF "$RIF_PHENO" "$KINSHIP" \
    "$WORKDIR/pyseer_results_RIF.txt" "$WORKDIR/pyseer_RIF.log" "$WORKDIR/patterns_RIF.txt"
echo "Rifampicin GWAS complete."


# =============================================================================
# STEP 6b: Bonferroni threshold from unique variant patterns
# =============================================================================
# Variants with identical presence/absence patterns across isolates are the
# same statistical test. In a clonal population, the number of unique patterns
# (not the number of variants) is the appropriate Bonferroni denominator.
# Patterns are written by pyseer --output-patterns (Steps 4-6).

echo "============================================================"
echo "[PART B] STEP 6b: Bonferroni threshold from unique variant patterns"
echo "============================================================"

for DRUG in AMI CIP RIF; do
    PAT="$WORKDIR/patterns_${DRUG}.txt"
    RES="$WORKDIR/pyseer_results_${DRUG}.txt"
    if [[ ! -s "$PAT" ]]; then
        echo "  ERROR: $PAT missing or empty."
        echo "         Delete $RES and rerun so pyseer writes the patterns file."
        exit 1
    fi
    N_PAT=$(sort -u "$PAT" | wc -l)
    N_VAR=$(( $(wc -l < "$RES") - 1 ))
    python3 -c "print(0.05/$N_PAT)" > "$WORKDIR/threshold_${DRUG}.txt"
    echo "  $DRUG: $N_VAR variants tested, $N_PAT unique patterns -> Bonferroni p < $(cat "$WORKDIR/threshold_${DRUG}.txt")"
done

# =============================================================================
# STEP 6c: Kinship sensitivity analysis (standard C matrix vs legacy matrix)
# =============================================================================
# Reruns pyseer with the legacy kinship matrix and compares, per drug:
#   h^2, Spearman correlation of -log10(p) across shared variants,
#   overlap of the top-20 variants, and counts of nominal / Bonferroni hits.
# The main analysis and all downstream steps use the standard matrix only.

echo "============================================================"
echo "[PART B] STEP 6c: Kinship sensitivity analysis (standard vs legacy)"
echo "============================================================"

for DRUG in AMI CIP RIF; do
    case "$DRUG" in
        AMI) PH="$AMI_PHENO" ;;
        CIP) PH="$CIP_PHENO" ;;
        RIF) PH="$RIF_PHENO" ;;
    esac
    run_pyseer "$DRUG" "$PH" "$KINSHIP_LEGACY" \
        "$WORKDIR/pyseer_results_${DRUG}_legacy.txt" "$WORKDIR/pyseer_${DRUG}_legacy.log"
done

python3 - "$WORKDIR" << 'KCOMP_EOF'
import os, re, sys
import numpy as np
import pandas as pd

W = sys.argv[1]
TOP_N = 20

def h2(log):
    try:
        with open(log) as f:
            for line in f:
                m = re.search(r"h\^2\s*=\s*([0-9.eE+-]+)", line)
                if m:
                    return float(m.group(1))
    except OSError:
        pass
    return float("nan")

rows = []
for d in ("AMI", "CIP", "RIF"):
    f_new = os.path.join(W, f"pyseer_results_{d}.txt")
    f_old = os.path.join(W, f"pyseer_results_{d}_legacy.txt")
    if not (os.path.isfile(f_new) and os.path.isfile(f_old)):
        print(f"  [{d}] results missing — skipped")
        continue
    new = pd.read_csv(f_new, sep="\t")[["variant", "lrt-pvalue"]]
    old = pd.read_csv(f_old, sep="\t")[["variant", "lrt-pvalue"]]
    thr_path = os.path.join(W, f"threshold_{d}.txt")
    thr = float(open(thr_path).read().strip()) if os.path.isfile(thr_path) else float("nan")

    m = new.merge(old, on="variant", suffixes=("_C", "_legacy")).dropna()
    rho = m["lrt-pvalue_C"].rank().corr(m["lrt-pvalue_legacy"].rank()) if len(m) > 2 else float("nan")
    top_new = set(new.nsmallest(TOP_N, "lrt-pvalue")["variant"])
    top_old = set(old.nsmallest(TOP_N, "lrt-pvalue")["variant"])

    rows.append({
        "drug": d,
        "h2_standard": h2(os.path.join(W, f"pyseer_{d}.log")),
        "h2_legacy": h2(os.path.join(W, f"pyseer_{d}_legacy.log")),
        "n_variants_shared": len(m),
        "spearman_rho_pvalues": round(rho, 4) if rho == rho else "NA",
        f"top{TOP_N}_overlap": len(top_new & top_old),
        "n_nominal_standard": int((new["lrt-pvalue"] < 0.05).sum()),
        "n_nominal_legacy": int((old["lrt-pvalue"] < 0.05).sum()),
        "bonferroni_threshold": thr,
        "n_bonferroni_standard": int((new["lrt-pvalue"] < thr).sum()) if thr == thr else "NA",
        "n_bonferroni_legacy": int((old["lrt-pvalue"] < thr).sum()) if thr == thr else "NA",
    })

out = os.path.join(W, "kinship_comparison.tsv")
df = pd.DataFrame(rows)
df.to_csv(out, sep="\t", index=False)
with pd.option_context("display.width", 200, "display.max_columns", None):
    print(df.to_string(index=False))
print(f"\n  Kinship comparison written: {out}")
print("  Interpretation: rho close to 1 and high top-20 overlap = results robust")
print("  to kinship construction; report the standard (C matrix) results.")
KCOMP_EOF

# =============================================================================
# STEP 7: Filter significant hits, annotate, and compute diagnostic statistics
# =============================================================================
# Filters variants with raw lrt-pvalue < 0.05, then for each hit:
#   (a) Extracts per-sample genotypes for the matching REF/ALT record from
#       the split VCF (core.split.vcf.gz) and computes a 2×2 contingency
#       table against the binary resistance phenotype (rank >= threshold).
#   (b) Annotates the variant from the per-isolate snippy snps.tab of the
#       first ALT carrier with a matching CHROM/POS/REF/ALT row.
#
# For each significant SNP:
#   SNP+  : ALT allele present (0/1 or 1/1 in VCF)
#   SNP-  : REF only (0/0)
#   R     : resistant  (ordinal rank >= threshold)
#   S     : susceptible (ordinal rank <  threshold)
#
#              R       S
#   SNP+      TP      FP
#   SNP-      FN      TN
#
# Sensitivity = TP / (TP + FN)
# Specificity = TN / (TN + FP)
# PPV         = TP / (TP + FP)
# NPV         = TN / (TN + FN)
#
# 95% confidence intervals use the Wilson score method, which performs well
# even with small counts and proportions near 0 or 1.
#
# Output columns (one combined TSV per drug):
#   [pyseer]     variant  af  filter-pvalue  lrt-pvalue  beta  beta-std-err
#                variant_h2  notes
#   [annotation] CHROM  POS  TYPE  REF  ALT  FTYPE  STRAND  NT_POS  AA_POS
#                EFFECT  LOCUS_TAG  GENE  PRODUCT
#   [direction]  direction (resistance if beta > 0, else susceptibility;
#                2x2 table flipped for susceptibility-direction variants)
#   [diagnostic] TP  FP  TN  FN  n_resistant  n_susceptible
#                sensitivity  sensitivity_ci95_lo  sensitivity_ci95_hi
#                specificity  specificity_ci95_lo  specificity_ci95_hi
#                PPV  PPV_ci95_lo  PPV_ci95_hi
#                NPV  NPV_ci95_lo  NPV_ci95_hi
#
# Resistance thresholds:
#   RIF: rank >= 5  (MIC >= 2 mg/L)  [I+R non-susceptible]
#   AMI: rank >= 6  (MIC >= 32 mg/L) [I+R non-susceptible; was rank 7]
#   CIP: rank >= 5  (MIC >= 2 mg/L)  [I+R non-susceptible; was rank 6]
#
# Replicates and ATCC 29548 are excluded from diagnostic statistics.

echo "============================================================"
echo "[PART B] STEP 7: Filtering significant hits, annotating, and computing diagnostic statistics"
echo "============================================================"

filter_and_diagnose() {
    local RESULTS_FILE="$1"
    local PHENO_FILE="$2"
    local VCF_FILE="$3"
    local RESIST_RANK="$4"
    local OUT_FILE="$5"
    local DRUG="$6"

    echo "  $DRUG: Filtering lrt-pvalue < 0.05, annotating from snps.tab, computing diagnostic stats..."
    echo "  $DRUG: Replicates and ATCC 29548 excluded from diagnostic statistics."

    # Serialise the bash REPLICATES array into a comma-separated string
    local REPLICATES_STR
    REPLICATES_STR=$(IFS=","; echo "${REPLICATES[*]}")

    python3 << PYEOF
import subprocess, sys, math, os, csv

results_file   = "$RESULTS_FILE"
pheno_file     = "$PHENO_FILE"
vcf_file       = "$VCF_FILE"
snippy_dir     = "$SNIPPY_DIR"
resist_rank    = int("$RESIST_RANK")
out_file       = "$OUT_FILE"
drug           = "$DRUG"
replicates_str = "$REPLICATES_STR"

replicates = set(r.strip() for r in replicates_str.split(",") if r.strip())
print(f"  [{drug}] Excluding {len(replicates)} replicate sample(s) from diagnostic stats.")

# -------------------------------------------------------------------------
# Wilson score 95% CI for a proportion p = k/n
# -------------------------------------------------------------------------
def wilson_ci(k, n, z=1.96):
    if n == 0:
        return float('nan'), float('nan')
    p_hat = k / n
    denom = 1 + z**2 / n
    centre = (p_hat + z**2 / (2 * n)) / denom
    margin = z * math.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2)) / denom
    return max(0.0, centre - margin), min(1.0, centre + margin)

# -------------------------------------------------------------------------
# snps.tab annotation lookup
#
# snippy writes per-isolate annotation to <snippy_dir>/<sample>/snps.tab.
# Columns (tab-separated, 0-indexed):
#   0  CHROM  1  POS  2  TYPE  3  REF  4  ALT  5  EVIDENCE
#   6  FTYPE  7  STRAND  8  NT_POS  9  AA_POS  10  EFFECT
#   11 LOCUS_TAG  12  GENE  13  PRODUCT
#
# The annotation is allele- and position-specific, so any isolate carrying
# the ALT allele will have the same annotation. We try carriers in order
# and return on the first successful match.
# -------------------------------------------------------------------------
TAB_COLS = ['CHROM', 'POS', 'TYPE', 'REF', 'ALT',
            'FTYPE', 'STRAND', 'NT_POS', 'AA_POS', 'EFFECT',
            'LOCUS_TAG', 'GENE', 'PRODUCT']

_gff_cache = {}

def load_ref_gff(carrier_samples):
    """CDS features from the snippy reference annotation (<sample>/reference/ref.gff)."""
    if 'cds' in _gff_cache:
        return _gff_cache['cds']
    cds = []
    for sample in carrier_samples:
        gff = os.path.join(snippy_dir, sample, 'reference', 'ref.gff')
        if not os.path.isfile(gff):
            continue
        with open(gff) as f:
            for line in f:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#'):
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] != 'CDS':
                    continue
                attrs = {}
                for kv in c[8].split(';'):
                    if '=' in kv:
                        k, v = kv.split('=', 1)
                        attrs[k] = v
                cds.append((c[0], int(c[3]), int(c[4]), c[6],
                            attrs.get('locus_tag', 'NA'), attrs.get('gene', 'NA'),
                            attrs.get('product', 'NA')))
        break
    _gff_cache['cds'] = cds
    return cds

def mark_intergenic(ann):
    """snippy leaves FTYPE/LOCUS_TAG empty for calls outside annotated CDS;
    label those explicitly as intergenic instead of NA."""
    if ann.get('LOCUS_TAG', 'NA') == 'NA' and ann.get('FTYPE', 'NA') == 'NA':
        ann['FTYPE'] = 'intergenic'
        ann['LOCUS_TAG'] = 'intergenic'
        ann['PRODUCT'] = 'intergenic'
        eff = ann.get('EFFECT', 'NA')
        if eff == 'NA':
            ann['EFFECT'] = 'intergenic'
        elif eff.endswith('] NA'):
            ann['EFFECT'] = eff[:-2] + 'intergenic'
    return ann

def lookup_snps_tab(chrom, pos, ref, alt, carrier_samples):
    """
    Annotate a variant, trying in order:
      (1) an exact chrom:pos:ref:alt row in a carrier's snippy snps.tab;
      (2) a multi-base snps.tab row (mnp/complex) whose span contains this
          position. pyseer tests the split single-base records, so they have
          no exact row; the parent call's gene/effect is copied, TYPE is the
          tested variant's type, and EFFECT is prefixed with the parent call;
      (3) the reference annotation (ref.gff): locus tag/gene/product of the CDS
          containing the position, or 'intergenic' if none.
    """
    pos = int(pos)
    tab_indices = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13]
    vtype = 'snp' if (len(ref) == 1 and len(alt) == 1) else (
        'ins' if len(alt) > len(ref) else ('del' if len(alt) < len(ref) else 'mnp'))
    parent = None

    for sample in carrier_samples:
        tab_path = os.path.join(snippy_dir, sample, 'snps.tab')
        if not os.path.isfile(tab_path):
            continue
        try:
            with open(tab_path) as f:
                next(f)  # skip header
                for line in f:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 14 or parts[0] != chrom:
                        continue
                    # (1) exact match
                    if parts[1] == str(pos) and parts[3] == ref and parts[4] == alt:
                        return mark_intergenic({TAB_COLS[i]: (parts[tab_indices[i]] if parts[tab_indices[i]] else 'NA')
                                                for i in range(len(TAB_COLS))})
                    # (2) remember a multi-base call spanning this position
                    if parent is None and len(parts[3]) > 1:
                        try:
                            p0 = int(parts[1])
                        except ValueError:
                            continue
                        off = pos - p0
                        if 0 <= off < len(parts[3]):
                            same_len = len(parts[3]) == len(parts[4])
                            if (not same_len) or (parts[3][off] == ref[0] and parts[4][off] == alt[0]):
                                parent = parts
        except OSError:
            continue

    if parent is not None:
        ann = {TAB_COLS[i]: (parent[tab_indices[i]] if parent[tab_indices[i]] else 'NA')
               for i in range(len(TAB_COLS))}
        tag = parent[1] + '_' + parent[3] + '_' + parent[4]
        ann.update({'CHROM': chrom, 'POS': str(pos), 'REF': ref, 'ALT': alt, 'TYPE': vtype,
                    'EFFECT': '[from multi-base call ' + tag + '] ' + ann.get('EFFECT', 'NA')})
        return mark_intergenic(ann)

    # (3) reference annotation fallback
    ann = {c: 'NA' for c in TAB_COLS}
    ann.update({'CHROM': chrom, 'POS': str(pos), 'REF': ref, 'ALT': alt, 'TYPE': vtype})
    hits = [c for c in load_ref_gff(carrier_samples) if c[0] == chrom and c[1] <= pos <= c[2]]
    if hits:
        _, start, end, strand, lt, gene, product = hits[0]
        ann.update({'FTYPE': 'CDS', 'STRAND': strand, 'LOCUS_TAG': lt, 'GENE': gene,
                    'PRODUCT': product,
                    'EFFECT': '[from reference annotation; no snps.tab call]'})
    else:
        ann.update({'FTYPE': 'intergenic', 'LOCUS_TAG': 'intergenic', 'PRODUCT': 'intergenic',
                    'EFFECT': 'intergenic'})
    return ann

# -------------------------------------------------------------------------
# Parse pyseer output: filter to lrt-pvalue < 0.05, keep all columns
# -------------------------------------------------------------------------
PYSEER_COLS = ['variant', 'af', 'filter-pvalue', 'lrt-pvalue',
               'beta', 'beta-std-err', 'variant_h2', 'notes']

header_in = None
pcol      = None
sig_rows  = []

with open(results_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        line = line.rstrip('\n')

        if header_in is None:
            header_in = line.split('\t')
            try:
                pcol = header_in.index('lrt-pvalue')
            except ValueError:
                print(f"ERROR [{drug}]: 'lrt-pvalue' column not found in {results_file}")
                sys.exit(1)
            continue

        parts = line.split('\t')
        if len(parts) <= pcol:
            continue

        try:
            pval = float(parts[pcol])
        except ValueError:
            continue

        if pval >= 0.05:
            continue

        row = {}
        for col in PYSEER_COLS:
            if col in header_in:
                idx = header_in.index(col)
                row[col] = parts[idx] if idx < len(parts) else 'NA'
            else:
                row[col] = 'NA'

        sig_rows.append(row)

sig_rows.sort(key=lambda r: float(r['lrt-pvalue']) if r['lrt-pvalue'] != 'NA' else 1)

n_hits = len(sig_rows)
print(f"  [{drug}] {n_hits} significant hit(s) (lrt-pvalue < 0.05).")

if n_hits == 0:
    print(f"  [{drug}] No output file written.")
    sys.exit(0)

# -------------------------------------------------------------------------
# Load phenotype: sample -> ordinal rank (replicates and ATCC 29548 excluded)
# -------------------------------------------------------------------------
pheno = {}
with open(pheno_file) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            sample = parts[0].strip()
            if sample not in replicates and sample != "CT-WNTC9548":
                pheno[sample] = float(parts[1])

# -------------------------------------------------------------------------
# Column definitions
# -------------------------------------------------------------------------
# TAB_COLS defined above (matches snps.tab column order exactly)
DIAG_COLS = [
    'direction',
    'TP', 'FP', 'TN', 'FN', 'n_resistant', 'n_susceptible',
    'sensitivity', 'sensitivity_ci95_lo', 'sensitivity_ci95_hi',
    'specificity', 'specificity_ci95_lo', 'specificity_ci95_hi',
    'PPV',         'PPV_ci95_lo',         'PPV_ci95_hi',
    'NPV',         'NPV_ci95_lo',         'NPV_ci95_hi',
]

def fmt(v):
    if isinstance(v, float):
        return 'NA' if math.isnan(v) else f'{v:.4f}'
    return str(v)

results = []

for row in sig_rows:
    variant_id = row['variant']

    # Parse CHROM_POS_REF_ALT from the pyseer variant ID
    try:
        chrom, pos, ref, alt = variant_id.rsplit('_', 3)
    except ValueError:
        print(f"  WARNING: Cannot parse variant ID '{variant_id}' — skipping.")
        ann  = {c: 'NA' for c in TAB_COLS}
        diag = {c: 'NA' for c in DIAG_COLS}
        results.append({**row, **ann, **diag})
        continue

    region = f"{chrom}:{pos}-{pos}"

    # -----------------------------------------------------------------------
    # (a) Extract per-sample genotypes from core.vcf.gz
    #     We do this first so we have the carrier list for snps.tab lookup.
    # -----------------------------------------------------------------------
    try:
        gt_output = subprocess.check_output(
            ["bcftools", "query", "-r", region, "-f",
             "%REF\t%ALT[\t%SAMPLE=%GT]\n", vcf_file],
            stderr=subprocess.DEVNULL
        ).decode().strip()
    except subprocess.CalledProcessError:
        print(f"  WARNING: bcftools query failed for {variant_id} — skipping.")
        ann  = {c: 'NA' for c in TAB_COLS}
        ann.update({'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt})
        diag = {c: 'NA' for c in DIAG_COLS}
        results.append({**row, **ann, **diag})
        continue

    if not gt_output:
        print(f"  WARNING: No genotype data for {variant_id} — skipping.")
        ann  = {c: 'NA' for c in TAB_COLS}
        ann.update({'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt})
        diag = {c: 'NA' for c in DIAG_COLS}
        results.append({**row, **ann, **diag})
        continue

    # Build sample -> has_alt from the record whose REF/ALT match this
    # variant (after Step 3c a position can have one record per ALT allele);
    # collect carriers for the snps.tab lookup
    has_alt  = {}
    carriers = []
    matched  = False
    for rec in gt_output.split('\n'):
        fields = rec.split('\t')
        if len(fields) < 3 or fields[0] != ref or fields[1] != alt:
            continue
        matched = True
        for sg in fields[2:]:
            sample, _, gt = sg.rpartition('=')
            is_alt = '1' in gt.replace('|', '/').split('/')
            has_alt[sample] = is_alt
            if is_alt:
                carriers.append(sample)
        break

    if not matched:
        print(f"  WARNING: No VCF record matching {variant_id} — skipping.")
        ann  = {c: 'NA' for c in TAB_COLS}
        ann.update({'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt})
        diag = {c: 'NA' for c in DIAG_COLS}
        results.append({**row, **ann, **diag})
        continue

    # -----------------------------------------------------------------------
    # (b) Annotation lookup from per-isolate snps.tab
    #     Uses the first carrier whose snps.tab contains a matching row.
    # -----------------------------------------------------------------------
    ann_dict = lookup_snps_tab(chrom, pos, ref, alt, carriers)

    eff = ann_dict.get('EFFECT', 'NA')
    if eff.startswith('[from reference annotation') or eff == 'intergenic':
        print(f"  NOTE: {variant_id} not in snps.tab — annotated from ref.gff "
              f"({ann_dict.get('LOCUS_TAG', 'NA')}).")

    # -----------------------------------------------------------------------
    # (c) 2x2 contingency table
    #     Replicates already absent from pheno dict.
    # -----------------------------------------------------------------------
    TP = FP = TN = FN = 0
    for sample, rank in pheno.items():
        if sample not in has_alt:
            continue
        resistant = (rank >= resist_rank)
        snp_pos   = has_alt[sample]
        if snp_pos and resistant:
            TP += 1
        elif snp_pos and not resistant:
            FP += 1
        elif not snp_pos and resistant:
            FN += 1
        else:
            TN += 1

    total = TP + FP + TN + FN
    if total == 0:
        print(f"  WARNING: No overlapping samples for {variant_id} — skipping.")
        diag = {c: 'NA' for c in DIAG_COLS}
        results.append({**row, **ann_dict, **diag})
        continue

    # -----------------------------------------------------------------------
    # (d) Determine association direction from beta sign.
    #     The reference genome (ATCC 29548 / CT-WNTC9548) is resistant to
    #     CIP, so variants are called relative to a resistant background.
    #     Positive beta  → ALT allele associated with HIGHER MIC → resistance marker
    #     Negative beta  → ALT allele associated with LOWER  MIC → susceptibility marker
    #     For susceptibility-direction variants, flip the 2×2 table so that
    #     TP/FP/TN/FN and all derived stats (sens/spec/PPV/NPV) reflect the
    #     ability of the variant to PREDICT SUSCEPTIBILITY, not resistance.
    #       Flip: TP_s=FP_r, FP_s=TP_r, TN_s=FN_r, FN_s=TN_r
    # -----------------------------------------------------------------------
    try:
        beta_val = float(row.get('beta', 'nan'))
    except (ValueError, TypeError):
        beta_val = float('nan')

    if not math.isnan(beta_val) and beta_val < 0:
        direction = "susceptibility"
        # Susceptibility prediction: ALT present + susceptible = TP
        TP_d, FP_d, TN_d, FN_d = FP, TP, FN, TN
    else:
        direction = "resistance"
        TP_d, FP_d, TN_d, FN_d = TP, FP, TN, FN

    # n_resistant / n_susceptible are phenotype counts — unchanged by direction
    n_resistant   = TP + FN   # total resistant samples in pheno
    n_susceptible = TN + FP   # total susceptible samples in pheno

    # Point estimates (from direction-adjusted contingency)
    sens = TP_d / (TP_d + FN_d) if (TP_d + FN_d) > 0 else float('nan')
    spec = TN_d / (TN_d + FP_d) if (TN_d + FP_d) > 0 else float('nan')
    ppv  = TP_d / (TP_d + FP_d) if (TP_d + FP_d) > 0 else float('nan')
    npv  = TN_d / (TN_d + FN_d) if (TN_d + FN_d) > 0 else float('nan')

    # Wilson score 95% CIs
    sens_lo, sens_hi = wilson_ci(TP_d, TP_d + FN_d)
    spec_lo, spec_hi = wilson_ci(TN_d, TN_d + FP_d)
    ppv_lo,  ppv_hi  = wilson_ci(TP_d, TP_d + FP_d)
    npv_lo,  npv_hi  = wilson_ci(TN_d, TN_d + FN_d)

    diag = {
        'direction':            direction,
        'TP':                   TP_d,
        'FP':                   FP_d,
        'TN':                   TN_d,
        'FN':                   FN_d,
        'n_resistant':          n_resistant,
        'n_susceptible':        n_susceptible,
        'sensitivity':          sens,
        'sensitivity_ci95_lo':  sens_lo,
        'sensitivity_ci95_hi':  sens_hi,
        'specificity':          spec,
        'specificity_ci95_lo':  spec_lo,
        'specificity_ci95_hi':  spec_hi,
        'PPV':                  ppv,
        'PPV_ci95_lo':          ppv_lo,
        'PPV_ci95_hi':          ppv_hi,
        'NPV':                  npv,
        'NPV_ci95_lo':          npv_lo,
        'NPV_ci95_hi':          npv_hi,
    }
    results.append({**row, **ann_dict, **diag})

# -------------------------------------------------------------------------
# Write combined TSV
# Column order: pyseer | snps.tab annotation | diagnostic stats
# -------------------------------------------------------------------------
out_cols = PYSEER_COLS + TAB_COLS + DIAG_COLS

with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_ALL)
    writer.writerow(out_cols)
    for r in results:
        writer.writerow([fmt(r.get(c, 'NA')) for c in out_cols])

print(f"  [{drug}] Combined hits + annotation + diagnostic stats written to {out_file}")
print(f"  [{drug}] Best lrt-pvalue: {sig_rows[0]['lrt-pvalue']}")
PYEOF

    echo ""
}

if [[ -f "$WORKDIR/gwas_hits_AMI.tsv" ]]; then
    echo "  SKIP AMI — gwas_hits_AMI.tsv already exists"
else
    filter_and_diagnose \
        "$WORKDIR/pyseer_results_AMI.txt" \
        "$AMI_PHENO" \
        "$GWAS_VCF" \
        "$AMI_RESIST_RANK" \
        "$WORKDIR/gwas_hits_AMI.tsv" \
        "AMI"
fi

if [[ -f "$WORKDIR/gwas_hits_CIP.tsv" ]]; then
    echo "  SKIP CIP — gwas_hits_CIP.tsv already exists"
else
    filter_and_diagnose \
        "$WORKDIR/pyseer_results_CIP.txt" \
        "$CIP_PHENO" \
        "$GWAS_VCF" \
        "$CIP_RESIST_RANK" \
        "$WORKDIR/gwas_hits_CIP.tsv" \
        "CIP"
fi

if [[ -f "$WORKDIR/gwas_hits_RIF.tsv" ]]; then
    echo "  SKIP RIF — gwas_hits_RIF.tsv already exists"
else
    filter_and_diagnose \
        "$WORKDIR/pyseer_results_RIF.txt" \
        "$RIF_PHENO" \
        "$GWAS_VCF" \
        "$RIF_RESIST_RANK" \
        "$WORKDIR/gwas_hits_RIF.tsv" \
        "RIF"
fi

# =============================================================================
# STEP 8: Generate Manhattan plots
# =============================================================================
# Produces one per-drug plot and one combined three-panel figure.
#
# Layout per panel:
#   x-axis : genomic position (parsed from pyseer variant ID CHROM_POS_REF_ALT)
#   y-axis : -log10(lrt-pvalue)
#   colours: alternating shades per chromosome/contig for visual separation
#   lines  : suggestive threshold (blue dashed, p = 0.05)
#             Bonferroni threshold (red dashed, p < 0.05 / n unique patterns,
#             read from threshold_<DRUG>.txt written in Step 6b)
#   labels : top hits annotated with variant ID (non-overlapping, above points)
#
# Variants where lrt-pvalue == 0 or is NA are excluded (cannot log-transform).
# Requires matplotlib and numpy (both standard in bioinformatics conda envs).

echo "============================================================"
echo "[PART B] STEP 8: Generating Manhattan plots"
echo "============================================================"

MANHATTAN_SCRIPT="$WORKDIR/make_manhattan.py"

cat > "$MANHATTAN_SCRIPT" << 'MANHATTAN_SCRIPT_EOF'
"""
Manhattan plot generator for pyseer LMM GWAS output.

Usage:
    python3 make_manhattan.py \
        <workdir> \
        <ami_results> <cip_results> <rif_results>

Outputs (written to workdir):
    manhattan_AMI.png
    manhattan_CIP.png
    manhattan_RIF.png
    manhattan_combined.png
"""

import sys
import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend — no display required
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

WORKDIR      = sys.argv[1]
RESULTS_AMI  = sys.argv[2]
RESULTS_CIP  = sys.argv[3]
RESULTS_RIF  = sys.argv[4]

# Colour palettes per drug (alternating light/dark pairs per contig)
PALETTES = {
    'AMI': ['#1a6b8a', '#4db3d4'],   # teal tones
    'CIP': ['#7b2d8b', '#b36ec7'],   # purple tones
    'RIF': ['#8b3a1a', '#d4804d'],   # amber tones
}

# Number of top hits to label on each plot
N_LABEL = 0

# Minimum -log10(p) for a point to be eligible for labelling
MIN_LOG10P_LABEL = 1.3  # corresponds to p < 0.05

# Figure dimensions
FIG_SINGLE_W = 14
FIG_SINGLE_H = 5
FIG_COMBINED_W = 14
FIG_COMBINED_H = 13

DPI = 180

# Bonferroni thresholds from unique variant patterns (written by Step 6b).
# Falls back to 0.05 / n_variants if a threshold file is missing.
BONF = {}
for _d in ('AMI', 'CIP', 'RIF'):
    _p = os.path.join(WORKDIR, f"threshold_{_d}.txt")
    if os.path.isfile(_p):
        with open(_p) as _f:
            BONF[_d] = float(_f.read().strip())
        print(f"  [{_d}] Bonferroni threshold (unique patterns): {BONF[_d]:.2e}")
    else:
        print(f"  [{_d}] WARNING: {_p} not found — using 0.05 / n_variants")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_pyseer_results(results_file):
    """
    Parse a pyseer output file.

    Returns a list of dicts with keys:
        variant_id, chrom, pos, neg_log10p
    Variants with missing, zero, or non-finite p-values are skipped.
    """
    rows = []
    header = None
    pcol = None

    with open(results_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            if header is None:
                header = line.split('\t')
                try:
                    pcol = header.index('lrt-pvalue')
                except ValueError:
                    print(f"  ERROR: 'lrt-pvalue' column not found in {results_file}")
                    sys.exit(1)
                continue

            parts = line.split('\t')
            if len(parts) <= pcol:
                continue

            variant_id = parts[0]

            try:
                pval = float(parts[pcol])
            except ValueError:
                continue

            if pval <= 0 or not math.isfinite(pval):
                continue

            # Parse CHROM_POS_REF_ALT  (rsplit from right, 3 splits)
            try:
                chrom, pos_str, ref, alt = variant_id.rsplit('_', 3)
                pos = int(pos_str)
            except (ValueError, TypeError):
                # If variant ID format is unexpected, use index as position
                chrom = 'unknown'
                pos   = len(rows)

            rows.append({
                'variant_id':  variant_id,
                'chrom':       chrom,
                'pos':         pos,
                'pval':        pval,
                'neg_log10p':  -math.log10(pval),
            })

    return rows


def assign_plot_positions(rows):
    """
    Convert genomic positions to monotonically increasing plot x-coordinates,
    adding a fixed gap between chromosomes/contigs.

    Returns:
        x_coords   : list of float plot-x values (same order as rows)
        chrom_mids : dict {chrom: mid_x} for x-axis tick labels
        chrom_list : ordered list of unique chroms
    """
    # Collect unique chroms in order of first appearance
    seen = []
    for r in rows:
        if r['chrom'] not in seen:
            seen.append(r['chrom'])
    chrom_list = seen

    # Determine max position per chrom
    chrom_max = {}
    for r in rows:
        c = r['chrom']
        chrom_max[c] = max(chrom_max.get(c, 0), r['pos'])

    # Gap between chroms = 5% of total genome span
    total_span = sum(chrom_max.values())
    gap = max(total_span * 0.05, 1)

    # Build cumulative offset per chrom
    offsets = {}
    cursor = 0
    for c in chrom_list:
        offsets[c] = cursor
        cursor += chrom_max[c] + gap

    # Assign x-coordinates
    x_coords = [offsets[r['chrom']] + r['pos'] for r in rows]

    # Chromosome midpoints for tick labels
    chrom_mids = {}
    for c in chrom_list:
        # Find min and max x for this chrom
        xs = [offsets[c] + r['pos'] for r in rows if r['chrom'] == c]
        chrom_mids[c] = (min(xs) + max(xs)) / 2

    return x_coords, chrom_mids, chrom_list


def draw_manhattan(ax, rows, x_coords, chrom_mids, chrom_list, palette,
                   drug, n_variants):
    """
    Draw a single Manhattan panel onto ax.
    """
    if not rows:
        ax.text(0.5, 0.5, f'No data for {drug}',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
        ax.set_title(drug, fontsize=13, fontweight='bold')
        return

    neg_log10p = np.array([r['neg_log10p'] for r in rows])
    x          = np.array(x_coords)

    # Assign colour by chrom index (alternating)
    chrom_idx  = {c: i for i, c in enumerate(chrom_list)}
    colors     = [palette[chrom_idx[r['chrom']] % len(palette)] for r in rows]

    # Scatter plot — all variants
    ax.scatter(x, neg_log10p, c=colors, s=8, alpha=0.75,
               linewidths=0, rasterized=True)

    # Significance thresholds — Bonferroni from unique patterns (Step 6b)
    suggestive_thresh  = -math.log10(0.05)
    bonf_pval          = BONF.get(drug, 0.05 / n_variants)
    bonferroni_thresh  = -math.log10(bonf_pval)

    ax.axhline(suggestive_thresh, color='steelblue', linestyle='--',
               linewidth=1.0, alpha=0.8, label='Suggestive (p < 0.05)')

    ax.axhline(bonferroni_thresh, color='crimson', linestyle='--',
               linewidth=1.0, alpha=0.9,
               label=f'Bonferroni (p < {bonf_pval:.2e})'.replace('e-06', ' × 10⁻⁶').replace('e-07', ' × 10⁻⁷').replace('e-05', ' × 10⁻⁵'))

    # Label top N hits
    # Sort by neg_log10p descending; label those above MIN_LOG10P_LABEL
    labellable = sorted(
        [(i, r) for i, r in enumerate(rows) if r['neg_log10p'] >= MIN_LOG10P_LABEL],
        key=lambda t: t[1]['neg_log10p'],
        reverse=True
    )[:N_LABEL]

    label_positions = []
    for idx, r in labellable:
        xi = x[idx]
        yi = r['neg_log10p']
        # Simple upward offset; enough for non-overlapping in most cases
        y_text = yi + 0.15 * (neg_log10p.max() - neg_log10p.min() + 1)
        ax.annotate(
            r['variant_id'],
            xy=(xi, yi),
            xytext=(xi, y_text),
            fontsize=5.5,
            color='#222222',
            ha='center',
            va='bottom',
            arrowprops=dict(arrowstyle='-', color='#888888', lw=0.6),
            clip_on=True,
        )

    # X-axis: chromosome tick labels
    # Only label if there are multiple chroms; otherwise show genomic position
    if len(chrom_list) > 1:
        ax.set_xticks([chrom_mids[c] for c in chrom_list])
        ax.set_xticklabels(
            [c.replace('Chromosome_', 'Chr').replace('chromosome_', 'Chr')
             for c in chrom_list],
            fontsize=7, rotation=45, ha='right'
        )
        ax.set_xlabel('Chromosome / Contig', fontsize=9)
    else:
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f'{int(v):,}')
        )
        ax.set_xlabel('Genomic Position (bp)', fontsize=9)

    # Tidy y-axis
    y_max = max(neg_log10p.max() * 1.15, bonferroni_thresh * 1.1)
    ax.set_ylim(bottom=0, top=y_max)
    ax.set_ylabel('−log₁₀(p)', fontsize=9)

    # Drug label as panel title
    drug_labels = {'AMI': 'Amikacin (AMI)', 'CIP': 'Ciprofloxacin (CIP)',
                   'RIF': 'Rifampicin (RIF)'}
    ax.set_title(drug_labels.get(drug, drug), fontsize=11, fontweight='bold', pad=6)

    # Minimal spine style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', labelsize=8)

    # Legend (placed upper right, small font)
    legend_elements = [
        Line2D([0], [0], color='steelblue', linestyle='--',
               lw=1.2, label='Suggestive  p < 0.05'),
        Line2D([0], [0], color='crimson', linestyle='--',
               lw=1.2, label=f'Bonferroni  p < {bonf_pval:.2e}'.replace('e-06', ' × 10⁻⁶').replace('e-07', ' × 10⁻⁷').replace('e-05', ' × 10⁻⁵')),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='upper right',
              framealpha=0.7, edgecolor='#cccccc')


def plot_single(drug, rows, palette, out_path):
    """Save a single-drug Manhattan plot."""
    if not rows:
        print(f"  [{drug}] No plottable variants — skipping individual plot.")
        return

    x_coords, chrom_mids, chrom_list = assign_plot_positions(rows)
    n_variants = len(rows)

    fig, ax = plt.subplots(figsize=(FIG_SINGLE_W, FIG_SINGLE_H))
    fig.patch.set_facecolor('white')

    draw_manhattan(ax, rows, x_coords, chrom_mids, chrom_list,
                   palette, drug, n_variants)

    plt.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  [{drug}] Individual Manhattan plot → {out_path}")


def plot_combined(drug_data, out_path):
    """Save a three-panel combined Manhattan plot (AMI / CIP / RIF)."""
    drugs = ['AMI', 'CIP', 'RIF']

    fig, axes = plt.subplots(3, 1, figsize=(FIG_COMBINED_W, FIG_COMBINED_H))
    fig.patch.set_facecolor('white')

    for ax, drug in zip(axes, drugs):
        rows    = drug_data[drug]['rows']
        palette = PALETTES[drug]

        if not rows:
            ax.text(0.5, 0.5, f'No data for {drug}',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        x_coords, chrom_mids, chrom_list = assign_plot_positions(rows)
        draw_manhattan(ax, rows, x_coords, chrom_mids, chrom_list,
                       palette, drug, len(rows))

    fig.suptitle(
        'M. haemophilum GWAS — Manhattan Plots\n'
        'Amikacin · Ciprofloxacin · Rifampicin',
        fontsize=13, fontweight='bold', y=1.01
    )

    plt.tight_layout(h_pad=2.5)
    fig.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Combined Manhattan plot → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

drug_files = {
    'AMI': RESULTS_AMI,
    'CIP': RESULTS_CIP,
    'RIF': RESULTS_RIF,
}

drug_data = {}

for drug, results_file in drug_files.items():
    print(f"  [{drug}] Parsing {results_file} ...")
    if not os.path.isfile(results_file):
        print(f"  [{drug}] File not found — skipping.")
        drug_data[drug] = {'rows': []}
        continue

    rows = parse_pyseer_results(results_file)
    print(f"  [{drug}] {len(rows):,} plottable variants loaded.")
    drug_data[drug] = {'rows': rows}

    # Per-drug plot
    out_single = os.path.join(WORKDIR, f"manhattan_{drug}.png")
    plot_single(drug, rows, PALETTES[drug], out_single)

# Combined three-panel plot
out_combined = os.path.join(WORKDIR, "manhattan_combined.png")
plot_combined(drug_data, out_combined)

print("  Manhattan plots complete.")
MANHATTAN_SCRIPT_EOF

if [[ -f "$WORKDIR/manhattan_combined.png" ]]; then
    echo "SKIP — Manhattan plots already exist"
else
    python3 "$MANHATTAN_SCRIPT" \
        "$WORKDIR" \
        "$WORKDIR/pyseer_results_AMI.txt" \
        "$WORKDIR/pyseer_results_CIP.txt" \
        "$WORKDIR/pyseer_results_RIF.txt"
fi

echo "Manhattan plots written."

# =============================================================================
# STEP 8b: QQ plots and genomic inflation factor (lambda)
# =============================================================================
# For each drug, observed -log10(lrt-pvalue) is plotted against the value
# expected under the null (uniform p-values; expected p_i = (i - 0.5) / n).
#   lambda = median(chi2_1 quantile of observed p) / 0.4549 (median of chi2_1)
# Dashed diagonal: expected; red dashed line: Bonferroni threshold from
# threshold_<DRUG>.txt (Step 6b). Horizontal plateaus are expected in a clonal
# bacterial population: many variants share one presence/absence pattern and
# therefore receive identical p-values, which also inflates lambda.
# Variants with lrt-pvalue NA or 0 are excluded.
# Outputs: qq_<DRUG>.png, qq_plots.png (three panels), genomic_inflation.tsv

echo "============================================================"
echo "[PART B] STEP 8b: QQ plots and genomic inflation factor"
echo "============================================================"

QQ_SCRIPT="$WORKDIR/make_qq.py"

cat > "$QQ_SCRIPT" << 'QQ_SCRIPT_EOF'
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

W = sys.argv[1]
DRUGS = [("AMI", "Amikacin (AMI)", "#1f5fa8"),
         ("CIP", "Ciprofloxacin (CIP)", "#d9533e"),
         ("RIF", "Rifampicin (RIF)", "#3f9e2a")]
MEDIAN_CHI2_1 = chi2.ppf(0.5, 1)   # 0.4549

def load(drug):
    f = os.path.join(W, f"pyseer_results_{drug}.txt")
    p = pd.read_csv(f, sep="\t")["lrt-pvalue"]
    p = pd.to_numeric(p, errors="coerce").dropna()
    p = p[(p > 0) & (p <= 1)].values
    thr_f = os.path.join(W, f"threshold_{drug}.txt")
    thr = float(open(thr_f).read().strip()) if os.path.isfile(thr_f) else float("nan")
    return p, thr

def draw(ax, drug, name, colour, p, thr):
    n = len(p)
    obs = -np.log10(np.sort(p))
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    lam = np.median(chi2.isf(p, 1)) / MEDIAN_CHI2_1
    ax.scatter(exp, obs, s=9, c=colour, alpha=0.6, linewidths=0)
    top = exp.max() * 1.1
    ax.plot([0, top], [0, top], "k--", lw=1.1, label="Expected")
    ymax = obs.max()
    if thr == thr:
        ax.axhline(-np.log10(thr), color="red", ls="--", lw=1,
                   label=f"Bonferroni (p < {thr:.2e})")
        ymax = max(ymax, -np.log10(thr))
    ax.set_xlim(0, top)
    ax.set_ylim(0, ymax * 1.12)
    ax.set_title(f"{name}\nλ = {lam:.3f}", fontsize=12)
    ax.set_xlabel("Expected −log₁₀(p)", fontsize=11)
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8.5, loc="lower right")
    return n, lam

rows, data = [], {}
for drug, name, colour in DRUGS:
    f = os.path.join(W, f"pyseer_results_{drug}.txt")
    if not os.path.isfile(f):
        print(f"  [{drug}] {f} missing — skipped")
        continue
    data[drug] = load(drug)

if not data:
    sys.exit("  ERROR: no pyseer results found for QQ plots")

# per-drug figures
for drug, name, colour in DRUGS:
    if drug not in data:
        continue
    p, thr = data[drug]
    fig, ax = plt.subplots(figsize=(5.5, 5.2), dpi=150)
    n, lam = draw(ax, drug, name, colour, p, thr)
    fig.tight_layout()
    fig.savefig(os.path.join(W, f"qq_{drug}.png"))
    plt.close(fig)
    rows.append({"drug": drug, "n_tested": n, "lambda": round(lam, 3),
                 "min_p": float(p.min()), "bonferroni_threshold": thr})
    print(f"  {drug}: {n} variants, lambda = {lam:.3f} -> qq_{drug}.png")

# combined figure
present = [d for d in DRUGS if d[0] in data]
fig, axs = plt.subplots(1, len(present), figsize=(5.2 * len(present), 5.35), dpi=150, squeeze=False)
for ax, (drug, name, colour) in zip(axs[0], present):
    draw(ax, drug, name, colour, *data[drug])
fig.suptitle("Reference-based SNP GWAS: Quantile-Quantile Plots\n"
             "M. haemophilum (n=56 clinical isolates)", fontweight="bold", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(W, "qq_plots.png"))
plt.close(fig)

pd.DataFrame(rows).to_csv(os.path.join(W, "genomic_inflation.tsv"), sep="\t", index=False)
print("  QQ plots complete: qq_plots.png, genomic_inflation.tsv")
QQ_SCRIPT_EOF

QQ_STALE=0
for DRUG in AMI CIP RIF; do
    if [[ ! -s "$WORKDIR/qq_${DRUG}.png" \
          || "$WORKDIR/pyseer_results_${DRUG}.txt" -nt "$WORKDIR/qq_${DRUG}.png" \
          || "$WORKDIR/threshold_${DRUG}.txt" -nt "$WORKDIR/qq_${DRUG}.png" ]]; then
        QQ_STALE=1
    fi
done
if [[ "$QQ_STALE" -eq 0 && -s "$WORKDIR/qq_plots.png" && -s "$WORKDIR/genomic_inflation.tsv" ]]; then
    echo "SKIP — QQ plots are up to date"
    column -t -s $'\t' "$WORKDIR/genomic_inflation.tsv" 2>/dev/null || cat "$WORKDIR/genomic_inflation.tsv"
else
    python3 "$QQ_SCRIPT" "$WORKDIR"
fi

echo "QQ plots written."

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "============================================================"
echo "GWAS STEPS 1-8b COMPLETE (Steps 9-10 follow)"
echo "============================================================"
echo ""
echo "Phylogenetic kinship matrix (pyseer phylogeny_distance.py --lmm --midpoint):"
echo "  $KINSHIP"
echo ""
echo "Kinship sensitivity analysis (standard vs legacy matrix):"
echo "  $WORKDIR/kinship_comparison.tsv"
echo ""
echo "Bonferroni thresholds (0.05 / n unique patterns):"
echo "  $WORKDIR/threshold_AMI.txt  ($(cat "$WORKDIR/threshold_AMI.txt" 2>/dev/null || echo NA))"
echo "  $WORKDIR/threshold_CIP.txt  ($(cat "$WORKDIR/threshold_CIP.txt" 2>/dev/null || echo NA))"
echo "  $WORKDIR/threshold_RIF.txt  ($(cat "$WORKDIR/threshold_RIF.txt" 2>/dev/null || echo NA))"
echo ""
echo "Full GWAS results:"
echo "  $WORKDIR/pyseer_results_AMI.txt"
echo "  $WORKDIR/pyseer_results_CIP.txt"
echo "  $WORKDIR/pyseer_results_RIF.txt"
echo ""
echo "Significant hits + annotation + diagnostic statistics (combined, one file per drug):"
echo "  $WORKDIR/gwas_hits_AMI.tsv  [resistance threshold: rank >= $AMI_RESIST_RANK, MIC >= 32 mg/L (I=R)]"
echo "  $WORKDIR/gwas_hits_CIP.tsv  [resistance threshold: rank >= $CIP_RESIST_RANK, MIC >= 2 mg/L (I=R)]"
echo "  $WORKDIR/gwas_hits_RIF.tsv  [resistance threshold: rank >= $RIF_RESIST_RANK, MIC >= 2 mg/L]"
echo ""
echo "  Columns:"
echo "    [pyseer]     variant | af | filter-pvalue | lrt-pvalue | beta | beta-std-err |"
echo "                 variant_h2 | notes"
echo "    [annotation] CHROM | POS | TYPE | REF | ALT | FTYPE | STRAND | NT_POS | AA_POS |"
echo "                 EFFECT | LOCUS_TAG | GENE | PRODUCT | direction"
echo "    [diagnostic] TP | FP | TN | FN | n_resistant | n_susceptible |"
echo "                 sensitivity | sensitivity_ci95_lo | sensitivity_ci95_hi |"
echo "                 specificity | specificity_ci95_lo | specificity_ci95_hi |"
echo "                 PPV | PPV_ci95_lo | PPV_ci95_hi |"
echo "                 NPV | NPV_ci95_lo | NPV_ci95_hi"
echo "  Annotation source: per-isolate snps.tab files from snippy ($SNIPPY_DIR/<sample>/snps.tab)."
echo "  For each hit, the first ALT-carrier sample with a matching snps.tab row is used."
echo "  CIs: Wilson score 95% confidence intervals."
echo ""
echo "  Note: Replicates and ATCC 29548 excluded from diagnostic stats:"
for R in "${REPLICATES[@]}"; do echo "    $R"; done
echo ""
echo "Manhattan plots:"
echo "  $WORKDIR/manhattan_AMI.png"
echo "  $WORKDIR/manhattan_CIP.png"
echo "  $WORKDIR/manhattan_RIF.png"
echo "  $WORKDIR/manhattan_combined.png"
echo ""
echo "QQ plots and genomic inflation (lambda):"
echo "  $WORKDIR/qq_AMI.png"
echo "  $WORKDIR/qq_CIP.png"
echo "  $WORKDIR/qq_RIF.png"
echo "  $WORKDIR/qq_plots.png"
echo "  $WORKDIR/genomic_inflation.tsv"
echo ""
echo "IMPORTANT — post-run checks:"
echo "  1. Confirm sample name matching from Step 2 output."
echo "     Tree tip labels, VCF sample names, and phenotype file sample"
echo "     names must all match exactly."
echo "  2. Check LMM heritability (h²) per drug in pyseer_<DRUG>.log."
echo "     Low h² may indicate limited genetic basis or insufficient power."
echo "  3. With 56 isolates, power is limited. Treat significant hits as"
echo "     candidates requiring experimental validation."
echo "  4. Diagnostic stats are based on SNP presence/absence vs binary"
echo "     resistance classification. Interpret PPV/NPV in the context of"
echo "     your cohort's resistance prevalence (which affects both metrics)."
echo "  5. Manhattan plots show all variants from the full pyseer output."
echo "     Blue dashed line = suggestive threshold (p < 0.05)."
echo "     Red dashed line  = Bonferroni threshold (p < 0.05 / n unique patterns; see threshold_<DRUG>.txt)."
echo ""
echo "============================================================"
# =============================================================================
# STEP 9: Candidate resistance mutation summary tables
# =============================================================================
# Reads resistance_mutations_2.tsv (from mhaem_full_pipeline.sh Step 11) and
# produces one summary TSV per drug (RIF / AMI / CIP).
#
# Gene → drug groupings:
#   RIF : B586_RS03455  (rpoB)
#   AMI : B586_RS14090  (16S rRNA)
#   CIP : B586_RS00025  (gyrB) + B586_RS00030  (gyrA)
#
# Output columns per row (one row = one unique mutation):
#   Mutation              CHROM_POS_REF_ALT  (pyseer-style variant ID)
#   Gene / Product / Effect
#   N_with_mutation       clinical isolates (n=56 set) carrying ALT
#   Note                  flags: present in all isolates (uninformative),
#                         pooled from multiple representations, or split
#                         from a multi-base call
#   Source_calls          original snippy representations + carrier counts
#   N_missing / Missing_isolates
#                         isolates with N or - in core.full.aln at this
#                         position (no/low coverage) — excluded from the 2x2
#   Coverage_check        core.full.aln (SNPs) or snps.tab only (indels)
#   Conflicts             isolates where snps.tab and core.full.aln disagree
#
# Variant normalisation (replaces vt normalize, removed from snippy):
#   multi-base substitutions (e.g. ACAC>GCAT) are split into single-base SNPs,
#   indels are trimmed of shared flanking bases, and carriers are pooled per
#   normalised variant so each isolate is counted once per variant.
#   MIC50                 MIC (mg/L) inhibiting >=50% of carriers: the
#                         ceil(n/2)-th value of the carriers' sorted MICs.
#                         Always a tested concentration (no averaging or
#                         rounding); censored values keep "<=" / ">"
#   MIC_range             lowest - highest MIC among carriers (e.g. "1 - >8")
#   Sensitivity / Specificity / PPV / NPV  with Wilson 95% CIs
#
# Resistance thresholds:
#   RIF: rank >= 5  (MIC >= 2 mg/L)  [I+R non-susceptible]
#   AMI: rank >= 6  (MIC >= 32 mg/L) [I+R non-susceptible; was rank 7]
#   CIP: rank >= 5  (MIC >= 2 mg/L)  [I+R non-susceptible; was rank 6]
#
# Replicates excluded (n=4, as defined in REPLICATES array above).
# ATCC 29548 excluded (reference, no MIC data).
#
# Output format: tab-separated with all fields double-quoted (QUOTE_ALL).
# This prevents Excel from auto-converting fields like "4/436" (AA_POS)
# into dates when the file is opened directly.
# =============================================================================

echo "============================================================"
echo "[PART B] STEP 9: Candidate resistance mutation summary tables"
echo "============================================================"

RES_MUT_TSV="$WORKDIR/resistance_mutations_2.tsv"

if [[ ! -f "$RES_MUT_TSV" ]]; then
    echo "ERROR: $RES_MUT_TSV not found."
    echo "       Run mhaem_full_pipeline.sh Step 11 first."
    exit 1
fi

# Serialise replicates array into comma-separated string for Python
REPLICATES_STR=$(IFS=","; echo "${REPLICATES[*]}")

# Skip if all three output TSVs exist and are newer than their inputs
# (resistance_mutations_2.tsv and the clinical-only phenotype files)
STEP9_UP_TO_DATE=1
for D in RIF AMI CIP; do
    OUT9="$WORKDIR/resistance_mutations_${D}.tsv"
    if [[ ! -f "$OUT9" || "$RES_MUT_TSV" -nt "$OUT9" || \
          "$WORKDIR/pheno_clinical_${D}.txt" -nt "$OUT9" || \
          "$WORKDIR/core/core.full.aln" -nt "$OUT9" || "$0" -nt "$OUT9" ]]; then
        STEP9_UP_TO_DATE=0
    fi
done
if [[ "$STEP9_UP_TO_DATE" -eq 1 ]]; then
    echo "  SKIP — resistance mutation TSVs already exist"
else

# NOTE: heredoc is UNQUOTED so bash expands $WORKDIR, $REPLICATES_STR, etc.
python3 << STEP9_EOF

import os, sys, math, csv
from collections import defaultdict

# ---------------------------------------------------------------------------
# Paths and parameters — injected by bash (unquoted heredoc)
# ---------------------------------------------------------------------------
WORKDIR         = "$WORKDIR"
RES_MUT_TSV     = os.path.join(WORKDIR, "resistance_mutations_2.tsv")
AMI_PHENO       = "$AMI_PHENO"
CIP_PHENO       = "$CIP_PHENO"
RIF_PHENO       = "$RIF_PHENO"
REPLICATES_STR  = "$REPLICATES_STR"

replicates = set(r.strip() for r in REPLICATES_STR.split(",") if r.strip())
print(f"  Excluded from carrier counts + diagnostic stats: {sorted(replicates)} + CT-WNTC9548")

RIF_RESIST_RANK = $RIF_RESIST_RANK
AMI_RESIST_RANK = $AMI_RESIST_RANK
CIP_RESIST_RANK = $CIP_RESIST_RANK

# ---------------------------------------------------------------------------
# Locus-tag → drug grouping
# ---------------------------------------------------------------------------
LOCUS_DRUG = {
    "B586_RS03455": "RIF",
    "B586_RS14090": "AMI",
    "B586_RS00025": "CIP",
    "B586_RS00030": "CIP",
}

DRUG_PHENO = {
    "RIF": RIF_PHENO,
    "AMI": AMI_PHENO,
    "CIP": CIP_PHENO,
}

DRUG_RESIST_RANK = {
    "RIF": RIF_RESIST_RANK,
    "AMI": AMI_RESIST_RANK,
    "CIP": CIP_RESIST_RANK,
}

# MIC rank -> reported MIC (mg/L). Censored end ranks keep their qualifier
# ("<=" written in ASCII so Excel opens the TSV without encoding problems).
RIF_RANK_TO_MIC = {1: "<=0.12", 2: "0.25", 3: "0.5", 4: "1", 5: "2", 6: "4", 7: "8", 8: ">8"}
AMI_RANK_TO_MIC = {1: "<=1", 2: "2", 3: "4", 4: "8", 5: "16", 6: "32", 7: "64", 8: ">64"}
CIP_RANK_TO_MIC = {1: "<=0.12", 2: "0.25", 3: "0.5", 4: "1", 5: "2", 6: "4", 7: "8", 8: "16", 9: ">16"}

DRUG_RANK_TO_MIC = {"RIF": RIF_RANK_TO_MIC, "AMI": AMI_RANK_TO_MIC, "CIP": CIP_RANK_TO_MIC}

# ---------------------------------------------------------------------------
# Wilson score 95% CI
# ---------------------------------------------------------------------------
def wilson_ci(k, n, z=1.96):
    if n == 0:
        return float('nan'), float('nan')
    p_hat = k / n
    denom = 1 + z**2 / n
    centre = (p_hat + z**2 / (2 * n)) / denom
    margin = z * math.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2)) / denom
    return max(0.0, centre - margin), min(1.0, centre + margin)

def fmt(v, decimals=4):
    if isinstance(v, float):
        return "NA" if math.isnan(v) else f"{v:.{decimals}f}"
    if v is None:
        return "NA"
    return str(v)

# ---------------------------------------------------------------------------
# Load phenotype files: sample -> ordinal rank
# Excludes replicates AND CT-WNTC9548
# ---------------------------------------------------------------------------
def load_pheno(path):
    pheno = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sample = parts[0].strip()
                if sample not in replicates and sample != "CT-WNTC9548":
                    try:
                        pheno[sample] = float(parts[1])
                    except ValueError:
                        pass
    return pheno

pheno = {
    "RIF": load_pheno(RIF_PHENO),
    "AMI": load_pheno(AMI_PHENO),
    "CIP": load_pheno(CIP_PHENO),
}

print(f"  Phenotype sample counts after exclusions:")
for drug in ("RIF", "AMI", "CIP"):
    print(f"    {drug}: {len(pheno[drug])} samples")

# ---------------------------------------------------------------------------
# Parse resistance_mutations_2.tsv
# ---------------------------------------------------------------------------
mutation_carriers = defaultdict(lambda: defaultdict(set))
mutation_meta     = {}
mutation_sources  = defaultdict(lambda: defaultdict(set))   # mut_id -> original call -> samples

# ---------------------------------------------------------------------------
# Variant normalisation (replaces the vt normalize step removed from snippy)
#   - Equal-length REF/ALT (snp / mnp / substitution-type complex): split into
#     one single-base SNP per differing position. E.g. 7439 CGCGA>TGCGA
#     becomes 7439 C>T; 5623 ACAC>GCAT becomes 5623 A>G + 5626 C>T.
#   - Unequal-length REF/ALT (ins / del / complex indel): trim shared bases
#     from the right, then from the left (keeping at least one base), so the
#     same indel gets one representation across isolates.
# Carriers are then pooled per normalised variant, so an isolate is counted
# once per variant regardless of how snippy reported it.
# ---------------------------------------------------------------------------
def normalise(pos, ref, alt):
    pos = int(pos)
    if len(ref) == len(alt):
        return [(pos + i, r, a) for i, (r, a) in enumerate(zip(ref, alt)) if r != a]
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return [(pos, ref, alt)]

EXPECTED_LOCI = set(LOCUS_DRUG.keys())
n_lines = n_skipped = 0

with open(RES_MUT_TSV) as f:
    for raw_line in f:
        line = raw_line.rstrip('\n')
        if not line:
            continue
        if ':' not in line:
            n_skipped += 1
            continue

        file_path, tsv_row = line.split(':', 1)

        parts_path = file_path.replace('\\\\', '/').split('/')
        try:
            snippy_idx = next(i for i, p in enumerate(parts_path) if p == 'snippy')
            sample = parts_path[snippy_idx + 1]
        except (StopIteration, IndexError):
            print(f"  WARNING: Cannot parse sample from path: {file_path}")
            n_skipped += 1
            continue

        cols = tsv_row.split('\t')
        if len(cols) < 14:
            n_skipped += 1
            continue

        chrom     = cols[0].strip()
        pos       = cols[1].strip()
        ref       = cols[3].strip()
        alt       = cols[4].strip()
        effect    = cols[10].strip()
        locus_tag = cols[11].strip()
        gene      = cols[12].strip()
        product   = cols[13].strip()

        if locus_tag not in EXPECTED_LOCI:
            n_skipped += 1
            continue

        drug      = LOCUS_DRUG[locus_tag]
        orig_id   = f"{pos}_{ref}_{alt}"
        atoms     = normalise(pos, ref, alt)
        is_direct = (len(atoms) == 1 and atoms[0] == (int(pos), ref, alt))

        for a_pos, a_ref, a_alt in atoms:
            mut_id = f"{chrom}_{a_pos}_{a_ref}_{a_alt}"
            meta = mutation_meta.get(mut_id)
            # Prefer the annotation from a call reported directly as this
            # variant; a split multi-base call carries the effect of the whole
            # block, so it is used only if no direct call exists (flagged).
            if meta is None or (is_direct and not meta["DIRECT"]):
                mutation_meta[mut_id] = {
                    "CHROM": chrom, "POS": str(a_pos), "REF": a_ref, "ALT": a_alt,
                    "LOCUS_TAG": locus_tag, "GENE": gene, "PRODUCT": product,
                    "EFFECT": effect if is_direct else f"[from multi-base call {orig_id}] {effect}",
                    "DRUG": drug, "DIRECT": is_direct,
                }
            if sample not in replicates and sample != "CT-WNTC9548":
                mutation_carriers[drug][mut_id].add(sample)
                mutation_sources[mut_id][orig_id].add(sample)

        n_lines += 1

print(f"  Parsed {n_lines} lines; {n_skipped} skipped.")
for drug in ("RIF", "AMI", "CIP"):
    n_pooled = sum(1 for m in mutation_carriers[drug] if len(mutation_sources[m]) > 1)
    print(f"    {drug}: {len(mutation_carriers[drug])} unique normalised variant(s); "
          f"{n_pooled} pooled from more than one snippy representation")

# ---------------------------------------------------------------------------
# Coverage-aware genotypes from snippy-core core.full.aln
#
# snps.tab only lists variants that were called, so an isolate with no row
# could be wild type OR uncovered. core.full.aln has every isolate's base at
# every reference position: the base itself, N (low coverage) or - (no
# coverage). For each single-base SNP the isolate is classified as:
#   carrier  : base == ALT (or called in snps.tab)
#   missing  : N / - (excluded from that variant's 2x2 table)
#   wild type: any other base
# Indels are not in the alignment and keep the snps.tab-only classification.
#
# Safety checks (fall back to snps.tab-only if either fails):
#   - the Reference row must carry each variant's REF base (coordinates OK)
#   - snps.tab carriers must show the ALT base in the alignment (>= 95%)
# ---------------------------------------------------------------------------
CORE_ALN = os.path.join(WORKDIR, "core", "core.full.aln")
ALN_OK   = False
aln_bases, ref_bases = {}, {}

all_clinical = set()
for d in pheno:
    all_clinical |= set(pheno[d])

def is_single(m):
    return len(m["REF"]) == 1 and len(m["ALT"]) == 1

snp_positions = sorted({int(m["POS"]) for m in mutation_meta.values() if is_single(m)})

if os.path.isfile(CORE_ALN) and snp_positions:
    def _store(name, chunks):
        if name is None or not (name == "Reference" or name in all_clinical):
            return
        seq = "".join(chunks)
        bases = {p: (seq[p - 1].upper() if p <= len(seq) else "-") for p in snp_positions}
        if name == "Reference":
            ref_bases.update(bases)
        else:
            aln_bases[name] = bases

    name, chunks = None, []
    with open(CORE_ALN) as f:
        for line in f:
            if line.startswith(">"):
                _store(name, chunks)
                name, chunks = line[1:].strip().split()[0], []
            elif name == "Reference" or name in all_clinical:
                chunks.append(line.strip())
    _store(name, chunks)

    ref_mism = [m for m in mutation_meta.values()
                if is_single(m) and ref_bases.get(int(m["POS"])) != m["REF"].upper()]
    n_obs = n_agree = 0
    for d in mutation_carriers:
        for mid, samples in mutation_carriers[d].items():
            m = mutation_meta[mid]
            if not is_single(m):
                continue
            for s in samples:
                if s in aln_bases:
                    n_obs += 1
                    n_agree += aln_bases[s].get(int(m["POS"])) == m["ALT"].upper()
    concord = n_agree / n_obs if n_obs else 0.0

    if not ref_bases:
        print("  WARNING: no 'Reference' row in core.full.aln — coverage NOT checked (snps.tab only).")
    elif ref_mism:
        print(f"  WARNING: core.full.aln Reference row disagrees with REF at {len(ref_mism)} position(s) "
              "— coverage NOT checked (snps.tab only).")
    elif concord < 0.95:
        print(f"  WARNING: only {concord:.1%} of snps.tab carriers show ALT in core.full.aln "
              "— coverage NOT checked (snps.tab only).")
    else:
        ALN_OK = True
        missing_iso = [s for s in all_clinical if s not in aln_bases]
        print(f"  core.full.aln: {len(aln_bases)} clinical isolates, {len(snp_positions)} SNP position(s); "
              f"snps.tab/alignment concordance {concord:.1%}")
        if missing_iso:
            print(f"  WARNING: not in core.full.aln (treated as missing): {', '.join(sorted(missing_iso))}")
else:
    print("  WARNING: core.full.aln not found — coverage NOT checked (snps.tab only).")

# ---------------------------------------------------------------------------
# Output columns
# ---------------------------------------------------------------------------
OUT_COLS = [
    "Mutation", "Gene", "Product", "Effect",
    "N_with_mutation", "MIC50", "MIC_range", "Note", "Source_calls",
    "N_missing", "Missing_isolates", "Coverage_check", "Conflicts",
    # Resistance-direction stats (ALT allele predicts resistance)
    "TP", "FP", "TN", "FN", "n_resistant", "n_susceptible",
    "Sensitivity", "Sensitivity_CI95_lo", "Sensitivity_CI95_hi",
    "Specificity", "Specificity_CI95_lo", "Specificity_CI95_hi",
    "PPV", "PPV_CI95_lo", "PPV_CI95_hi",
    "NPV", "NPV_CI95_lo", "NPV_CI95_hi",
    # Susceptibility-direction stats (ALT allele predicts susceptibility)
    # Relevant when the reference genome is resistant (e.g. CIP gyrA/gyrB variants).
    # 2x2 table is flipped: TP_S=FP, FP_S=TP, TN_S=FN, FN_S=TN
    "TP_S", "FP_S", "TN_S", "FN_S",
    "Sensitivity_S", "Sensitivity_S_CI95_lo", "Sensitivity_S_CI95_hi",
    "Specificity_S", "Specificity_S_CI95_lo", "Specificity_S_CI95_hi",
    "PPV_S", "PPV_S_CI95_lo", "PPV_S_CI95_hi",
    "NPV_S", "NPV_S_CI95_lo", "NPV_S_CI95_hi",
]

# ---------------------------------------------------------------------------
# Compute stats and write one TSV per drug
# All fields are double-quoted (csv.QUOTE_ALL) to prevent Excel from
# auto-converting values like "4/436" in AA_POS into dates.
# ---------------------------------------------------------------------------
for drug in ("RIF", "AMI", "CIP"):
    out_path       = os.path.join(WORKDIR, f"resistance_mutations_{drug}.tsv")
    drug_pheno     = pheno[drug]
    resist_rank    = DRUG_RESIST_RANK[drug]
    rank_to_mic    = DRUG_RANK_TO_MIC[drug]
    drug_mutations = mutation_carriers[drug]

    if not drug_mutations:
        print(f"  [{drug}] No mutations found — skipping output.")
        continue

    rows = []

    for mut_id in sorted(drug_mutations.keys(),
                         key=lambda m: int(mutation_meta[m]["POS"])):

        meta     = mutation_meta[mut_id]
        tab_carriers = drug_mutations[mut_id] & set(drug_pheno)   # called in snps.tab
        use_aln  = ALN_OK and is_single(meta)
        pos_i    = int(meta["POS"])
        alt_b    = meta["ALT"].upper()

        # Per-isolate genotype: carrier / wt / missing
        gt, conflicts = {}, []
        for s in drug_pheno:
            in_tab = s in tab_carriers
            if not use_aln:
                gt[s] = "carrier" if in_tab else "wt"
                continue
            b = aln_bases.get(s, {}).get(pos_i)
            if b is None:
                gt[s] = "carrier" if in_tab else "missing"
            elif b == alt_b or in_tab:
                gt[s] = "carrier"
            elif b in ("N", "-", "X", "?"):
                gt[s] = "missing"
            else:
                gt[s] = "wt"
            if b is not None and in_tab != (b == alt_b):
                conflicts.append(f"{s}(tab={'ALT' if in_tab else 'none'},aln={b})")

        carriers  = {s for s, g in gt.items() if g == "carrier"}
        missing   = sorted(s for s, g in gt.items() if g == "missing")
        evaluable = [s for s, g in gt.items() if g != "missing"]
        n_carriers = len(carriers)
        if n_carriers == 0:
            continue
        coverage_check = ("core.full.aln" if use_aln else
                          "snps.tab only (indel)" if ALN_OK else
                          "snps.tab only (alignment unavailable)")

        # How snippy originally reported this variant, with carrier counts
        src = mutation_sources[mut_id]
        source_calls = "; ".join(
            f"{k} (n={len(src[k] & carriers)})"
            for k in sorted(src, key=lambda k: -len(src[k] & carriers))
            if src[k] & carriers)

        notes = []
        if n_carriers == len(evaluable):
            notes.append("present in all evaluable clinical isolates (differs from reference only; uninformative)")
        if len(src) > 1:
            notes.append("pooled from multiple snippy representations")
        if not meta["DIRECT"]:
            notes.append("split from multi-base call")
        if missing:
            notes.append(f"{len(missing)} isolate(s) without coverage excluded")
        if conflicts:
            notes.append("snps.tab/alignment disagreement — check")
        note = "; ".join(notes)

        # MIC50 = MIC inhibiting >=50% of carriers: the ceil(n/2)-th value of
        # the sorted carrier MICs (a tested concentration; no averaging/rounding)
        carrier_ranks = sorted(int(round(drug_pheno[s])) for s in carriers if s in drug_pheno)
        if carrier_ranks:
            n = len(carrier_ranks)
            mic50 = rank_to_mic.get(carrier_ranks[math.ceil(n / 2) - 1], "NA")
            lo, hi = rank_to_mic.get(carrier_ranks[0], "NA"), rank_to_mic.get(carrier_ranks[-1], "NA")
            mic_range = lo if lo == hi else f"{lo} - {hi}"
        else:
            mic50, mic_range = "NA", "NA"

        TP = FP = TN = FN = 0
        for sample in evaluable:                      # isolates without coverage excluded
            resistant = (drug_pheno[sample] >= resist_rank)
            snp_pos   = (sample in carriers)
            if   snp_pos and     resistant: TP += 1
            elif snp_pos and not resistant: FP += 1
            elif not snp_pos and resistant: FN += 1
            else:                           TN += 1

        total = TP + FP + TN + FN
        if total == 0:
            print(f"  [{drug}] WARNING: No overlapping samples for {mut_id} — skipping.")
            continue

        sens = TP / (TP + FN) if (TP + FN) > 0 else float('nan')
        spec = TN / (TN + FP) if (TN + FP) > 0 else float('nan')
        ppv  = TP / (TP + FP) if (TP + FP) > 0 else float('nan')
        npv  = TN / (TN + FN) if (TN + FN) > 0 else float('nan')

        sens_lo, sens_hi = wilson_ci(TP, TP + FN)
        spec_lo, spec_hi = wilson_ci(TN, TN + FP)
        ppv_lo,  ppv_hi  = wilson_ci(TP, TP + FP)
        npv_lo,  npv_hi  = wilson_ci(TN, TN + FN)

        # Susceptibility-direction stats (flipped 2x2 table)
        # ALT allele in susceptible isolate = TP_S; ALT in resistant = FP_S
        TP_S, FP_S, TN_S, FN_S = FP, TP, FN, TN

        sens_S = TP_S / (TP_S + FN_S) if (TP_S + FN_S) > 0 else float('nan')
        spec_S = TN_S / (TN_S + FP_S) if (TN_S + FP_S) > 0 else float('nan')
        ppv_S  = TP_S / (TP_S + FP_S) if (TP_S + FP_S) > 0 else float('nan')
        npv_S  = TN_S / (TN_S + FN_S) if (TN_S + FN_S) > 0 else float('nan')

        sens_S_lo, sens_S_hi = wilson_ci(TP_S, TP_S + FN_S)
        spec_S_lo, spec_S_hi = wilson_ci(TN_S, TN_S + FP_S)
        ppv_S_lo,  ppv_S_hi  = wilson_ci(TP_S, TP_S + FP_S)
        npv_S_lo,  npv_S_hi  = wilson_ci(TN_S, TN_S + FN_S)

        rows.append({
            "Mutation":               mut_id,
            "Gene":                   meta["GENE"],
            "Product":                meta["PRODUCT"],
            "Effect":                 meta["EFFECT"],
            "N_with_mutation":        n_carriers,
            "Note":                   note,
            "Source_calls":           source_calls,
            "N_missing":              len(missing),
            "Missing_isolates":       ", ".join(missing),
            "Coverage_check":         coverage_check,
            "Conflicts":              "; ".join(conflicts),
            "MIC50":                  mic50,
            "MIC_range":              mic_range,
            # Resistance direction
            "TP":                     TP,
            "FP":                     FP,
            "TN":                     TN,
            "FN":                     FN,
            "n_resistant":            TP + FN,
            "n_susceptible":          TN + FP,
            "Sensitivity":            fmt(sens),
            "Sensitivity_CI95_lo":    fmt(sens_lo),
            "Sensitivity_CI95_hi":    fmt(sens_hi),
            "Specificity":            fmt(spec),
            "Specificity_CI95_lo":    fmt(spec_lo),
            "Specificity_CI95_hi":    fmt(spec_hi),
            "PPV":                    fmt(ppv),
            "PPV_CI95_lo":            fmt(ppv_lo),
            "PPV_CI95_hi":            fmt(ppv_hi),
            "NPV":                    fmt(npv),
            "NPV_CI95_lo":            fmt(npv_lo),
            "NPV_CI95_hi":            fmt(npv_hi),
            # Susceptibility direction
            "TP_S":                   TP_S,
            "FP_S":                   FP_S,
            "TN_S":                   TN_S,
            "FN_S":                   FN_S,
            "Sensitivity_S":          fmt(sens_S),
            "Sensitivity_S_CI95_lo":  fmt(sens_S_lo),
            "Sensitivity_S_CI95_hi":  fmt(sens_S_hi),
            "Specificity_S":          fmt(spec_S),
            "Specificity_S_CI95_lo":  fmt(spec_S_lo),
            "Specificity_S_CI95_hi":  fmt(spec_S_hi),
            "PPV_S":                  fmt(ppv_S),
            "PPV_S_CI95_lo":          fmt(ppv_S_lo),
            "PPV_S_CI95_hi":          fmt(ppv_S_hi),
            "NPV_S":                  fmt(npv_S),
            "NPV_S_CI95_lo":          fmt(npv_S_lo),
            "NPV_S_CI95_hi":          fmt(npv_S_hi),
        })

        print(f"  [{drug}] {mut_id}  n={n_carriers}  missing={len(missing)}  "
              f"TP={TP} FP={FP} TN={TN} FN={FN}  "
              f"sens={fmt(sens)} spec={fmt(spec)} PPV={fmt(ppv)} NPV={fmt(npv)}")

    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_ALL)
        writer.writerow(OUT_COLS)
        for r in rows:
            writer.writerow([fmt(r.get(c, 'NA')) if isinstance(r.get(c), float)
                             else str(r.get(c, 'NA')) for c in OUT_COLS])

    print(f"  [{drug}] {len(rows)} mutation(s) written -> {out_path}")

print("")
print("Step 9 complete.")
STEP9_EOF

fi  # end Step 9 skip

echo ""
echo "Resistance mutation summary tables:"
for D in RIF AMI CIP; do
    if [[ -f "$WORKDIR/resistance_mutations_${D}.tsv" ]]; then
        echo "  $WORKDIR/resistance_mutations_${D}.tsv"
    else
        echo "  [$D] no variants in target gene(s) — no table written"
    fi
done
echo ""
echo "Columns: Mutation | Gene | Product | Effect | N_with_mutation | MIC50 | MIC_range | Note | Source_calls |"
echo "         N_missing | Missing_isolates | Coverage_check | Conflicts |"
echo "         [Resistance direction] TP | FP | TN | FN | n_resistant | n_susceptible |"
echo "         Sensitivity (95% CI) | Specificity (95% CI) | PPV (95% CI) | NPV (95% CI) |"
echo "         [Susceptibility direction] TP_S | FP_S | TN_S | FN_S |"
echo "         Sensitivity_S (95% CI) | Specificity_S (95% CI) | PPV_S (95% CI) | NPV_S (95% CI)"
echo "  Susceptibility-direction stats use a flipped 2x2 table (TP_S=FP, FP_S=TP, TN_S=FN, FN_S=TN)."
echo "  Relevant when the reference genome is resistant (e.g. CIP gyrA/gyrB variants)."
echo "CIs: Wilson score 95% (equivalent to MedCalc Diagnostic Test Calculator)."
echo "Denominators: 56 clinical isolates (4 replicates + ATCC 29548 excluded)."
echo "Format: tab-separated, all fields double-quoted to prevent Excel date autoformat."

# =============================================================================
# STEP 10: gyrA intein status and lineage MIC comparison (exploratory)
# =============================================================================
# The reference gyrA (ATCC 29548) carries a ~420-aa intein. Isolates lacking it
# have no read coverage across the intein, so it marks the two main lineages.
#
# Intein status is called from two independent sources and cross-checked:
#   (a) Prokka gyrA length in each de novo assembly (pangenome.sh output):
#         > 3000 bp = intein present; 2300-2800 bp = absent;
#         anything else or a split gene = unresolved
#   (b) core.full.aln: % of the reference gyrA CDS that is uncovered (- or N):
#         < 10% = present; > 20% = absent; otherwise unresolved
#   Final status = the shared call when (a) and (b) agree, else "discordant".
#
# The intein boundaries are also estimated from the alignment: reference
# positions within gyrA that are uncovered in >= 90% of intein-negative
# clinical isolates.
#
# Outputs:
#   intein_status.tsv        one row per isolate (clinical, replicate, ATCC):
#                            gyrA length, % uncovered, calls, final status,
#                            and MICs for clinical isolates
#   intein_lineage_MIC.tsv   clinical isolates with a concordant call:
#                            n, MIC50, MIC90, range and % non-susceptible by
#                            intein status, with Mann-Whitney (MIC
#                            distribution) and Fisher exact (non-susceptible
#                            proportion) p-values
#
# INTERPRETATION: exploratory. Isolates within a lineage are not independent,
# so these p-values overstate the evidence (the LMM GWAS, which accounts for
# relatedness, found no significant variant). Treat as a lineage-level
# observation, not a resistance mechanism.
# =============================================================================

echo ""
echo "============================================================"
echo "[PART B] STEP 10: gyrA intein status and lineage MIC comparison"
echo "============================================================"

PROKKA_DIR="$HOME/pangenome/prokka"
REF_GBFF="$WORKDIR/GCF_000340435.2.gbff"
INTEIN_STATUS="$WORKDIR/intein_status.tsv"
INTEIN_MIC="$WORKDIR/intein_lineage_MIC.tsv"

STEP10_UP_TO_DATE=1
for F in "$INTEIN_STATUS" "$INTEIN_MIC"; do
    if [[ ! -f "$F" || "$0" -nt "$F" || "$WORKDIR/core/core.full.aln" -nt "$F" || \
          "$AMI_PHENO" -nt "$F" || "$CIP_PHENO" -nt "$F" || "$RIF_PHENO" -nt "$F" ]]; then
        STEP10_UP_TO_DATE=0
    fi
done

if [[ "$STEP10_UP_TO_DATE" -eq 1 ]]; then
    echo "  SKIP — intein outputs are up to date"
else
REPLICATES_STR=$(IFS=","; echo "${REPLICATES[*]}")
python3 - "$WORKDIR" "$PROKKA_DIR" "$REF_GBFF" "$CLINICAL_LIST" "$REPLICATES_STR" \
          "$AMI_PHENO" "$CIP_PHENO" "$RIF_PHENO" "$INTEIN_STATUS" "$INTEIN_MIC" \
          "$AMI_RESIST_RANK" "$CIP_RESIST_RANK" "$RIF_RESIST_RANK" << 'STEP10_EOF'
import os, re, sys, csv, math, glob

(workdir, prokka_dir, gbff, clin_list, rep_str, ami_p, cip_p, rif_p,
 out_status, out_mic, ami_r, cip_r, rif_r) = sys.argv[1:14]

ATCC = "CT-WNTC9548"
clinical   = [l.strip() for l in open(clin_list) if l.strip()]
replicates = [r for r in rep_str.split(",") if r]
samples    = clinical + replicates + [ATCC]
group      = {**{s: "clinical" for s in clinical},
              **{s: "replicate" for s in replicates}, ATCC: "ATCC 29548"}

LAB = {"RIF": {1:"<=0.12",2:"0.25",3:"0.5",4:"1",5:"2",6:"4",7:"8",8:">8"},
       "AMI": {1:"<=1",2:"2",3:"4",4:"8",5:"16",6:"32",7:"64",8:">64"},
       "CIP": {1:"<=0.12",2:"0.25",3:"0.5",4:"1",5:"2",6:"4",7:"8",8:"16",9:">16"}}
RES  = {"AMI": int(ami_r), "CIP": int(cip_r), "RIF": int(rif_r)}
PHEN = {"AMI": ami_p, "CIP": cip_p, "RIF": rif_p}

pheno = {}
for d, path in PHEN.items():
    pheno[d] = {}
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 2:
                pheno[d][p[0]] = int(round(float(p[1])))

# ---- (a) Prokka gyrA length -------------------------------------------------
def prokka_call(s):
    tsv = os.path.join(prokka_dir, s, f"{s}.tsv")
    if not os.path.isfile(tsv):
        return "NA", "unresolved"
    lens = []
    with open(tsv) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) > 3 and c[3] == "gyrA":
                try:
                    lens.append(int(c[2]))
                except ValueError:
                    pass
    if len(lens) != 1:
        return "/".join(map(str, lens)) or "NA", "unresolved"
    L = lens[0]
    if L > 3000:
        return str(L), "present"
    if 2300 <= L <= 2800:
        return str(L), "absent"
    return str(L), "unresolved"

# ---- reference gyrA CDS coordinates from the GenBank file --------------------
gyra = None
if os.path.isfile(gbff):
    loc = None
    with open(gbff) as f:
        for line in f:
            m = re.match(r"\s{5}CDS\s+(\S+)", line)
            if m:
                loc = m.group(1)
            elif loc and 'intein-containing DNA gyrase subunit A' in line:
                nums = list(map(int, re.findall(r"\d+", loc)))
                gyra = (min(nums), max(nums))
                break
            elif line.startswith("ORIGIN"):
                loc = None
if gyra:
    print(f"  Reference gyrA CDS: {gyra[0]}-{gyra[1]} ({gyra[1]-gyra[0]+1} bp)")
else:
    print("  WARNING: intein-containing gyrA not found in the GenBank file — alignment call skipped")

# ---- (b) alignment coverage across gyrA -------------------------------------
aln_seg = {}
aln_path = os.path.join(workdir, "core", "core.full.aln")
if gyra and os.path.isfile(aln_path):
    a, b = gyra
    name, chunks = None, []
    def store():
        if name in group:
            aln_seg[name] = "".join(chunks)[a - 1:b].upper()
    with open(aln_path) as f:
        for line in f:
            if line.startswith(">"):
                store()
                name, chunks = line[1:].strip().split()[0], []
            elif name in group:
                chunks.append(line.strip())
    store()

def aln_call(s):
    seg = aln_seg.get(s)
    if not seg:
        return "NA", "unresolved"
    pct = 100.0 * sum(ch in "-N" for ch in seg) / len(seg)
    if pct < 10:
        return f"{pct:.1f}", "present"
    if pct > 20:
        return f"{pct:.1f}", "absent"
    return f"{pct:.1f}", "unresolved"

# ---- per-isolate table -------------------------------------------------------
rows, status = [], {}
for s in samples:
    plen, pc = prokka_call(s)
    apct, ac = aln_call(s)
    if pc == ac and pc != "unresolved":
        final = pc
    elif "unresolved" in (pc, ac) and {pc, ac} - {"unresolved"}:
        final = ({pc, ac} - {"unresolved"}).pop() + " (single source)"
    else:
        final = "discordant" if pc != ac else "unresolved"
    status[s] = final
    row = {"Sample": s, "Group": group[s], "gyrA_bp_prokka": plen,
           "gyrA_pct_uncovered_aln": apct, "Intein_prokka": pc,
           "Intein_alignment": ac, "Intein_status": final}
    for d in ("CIP", "RIF", "AMI"):
        r = pheno[d].get(s)
        row[f"{d}_MIC"] = LAB[d][r] if (r is not None and group[s] == "clinical") else ""
    rows.append(row)

cols = ["Sample", "Group", "gyrA_bp_prokka", "gyrA_pct_uncovered_aln",
        "Intein_prokka", "Intein_alignment", "Intein_status",
        "CIP_MIC", "RIF_MIC", "AMI_MIC"]
with open(out_status, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", quoting=csv.QUOTE_ALL)
    w.writeheader()
    w.writerows(rows)

cl = [s for s in clinical]
n_pos = sum(status[s] == "present" for s in cl)
n_neg = sum(status[s] == "absent" for s in cl)
other = [f"{s}={status[s]}" for s in cl if status[s] not in ("present", "absent")]
print(f"  Clinical isolates: intein present {n_pos}, absent {n_neg}"
      + (f"; not concordant: {', '.join(other)}" if other else "; all calls concordant"))
print(f"  ATCC 29548: {status[ATCC]}")
print(f"  Replicates: " + ", ".join(f"{r}={status[r]}" for r in replicates))

# ---- intein boundaries estimated from the alignment ---------------------------
neg = [aln_seg[s] for s in cl if status[s] == "absent" and s in aln_seg]
if gyra and neg:
    L = len(neg[0])
    unc = [sum(seq[i] in "-N" for seq in neg) / len(neg) >= 0.9 for i in range(L)]
    runs, start = [], None
    for i, u in enumerate(unc + [False]):
        if u and start is None:
            start = i
        elif not u and start is not None:
            runs.append((start, i - 1)); start = None
    if runs:
        s0, e0 = max(runs, key=lambda r: r[1] - r[0])
        a = gyra[0]
        print(f"  Intein region (uncovered in >=90% of intein-negative isolates): "
              f"{a + s0}-{a + e0} ({e0 - s0 + 1} bp; gyrA c.{s0 + 1}-c.{e0 + 1})")

# ---- lineage MIC comparison (clinical, concordant calls only) ---------------
try:
    from scipy.stats import mannwhitneyu, fisher_exact
except ImportError:
    mannwhitneyu = fisher_exact = None
    print("  WARNING: scipy not available — p-values not computed")

def q(r, p):
    r = sorted(r)
    return r[math.ceil(p * len(r)) - 1]

out = []
for d in ("CIP", "RIF", "AMI"):
    g = {k: [pheno[d][s] for s in cl if status[s] == k and s in pheno[d]]
         for k in ("present", "absent")}
    if not g["present"] or not g["absent"]:
        continue
    ns = {k: sum(r >= RES[d] for r in g[k]) for k in g}
    mw = fp = float("nan")
    if mannwhitneyu:
        mw = mannwhitneyu(g["present"], g["absent"], alternative="two-sided").pvalue
        fp = fisher_exact([[ns["present"], len(g["present"]) - ns["present"]],
                           [ns["absent"],  len(g["absent"])  - ns["absent"]]])[1]
    for k, name in (("present", "intein+"), ("absent", "intein-")):
        r = g[k]
        out.append({
            "Drug": d, "Group": name, "n": len(r),
            "MIC50": LAB[d][q(r, .5)], "MIC90": LAB[d][q(r, .9)],
            "MIC_range": f"{LAB[d][min(r)]} - {LAB[d][max(r)]}",
            "NS_threshold": f">={LAB[d][RES[d]]}",
            "n_NS": ns[k], "pct_NS": f"{100 * ns[k] / len(r):.0f}",
            "MannWhitney_p": f"{mw:.3g}" if k == "present" else "",
            "Fisher_NS_p":   f"{fp:.3g}" if k == "present" else "",
        })

mcols = ["Drug", "Group", "n", "MIC50", "MIC90", "MIC_range", "NS_threshold",
         "n_NS", "pct_NS", "MannWhitney_p", "Fisher_NS_p"]
with open(out_mic, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=mcols, delimiter="\t", quoting=csv.QUOTE_ALL)
    w.writeheader()
    w.writerows(out)

print("")
print("  " + "\t".join(mcols))
for r in out:
    print("  " + "\t".join(str(r[c]) for c in mcols))
print(f"\n  Per-isolate status : {out_status}")
print(f"  Lineage MIC table  : {out_mic}")
print("  Exploratory: isolates within a lineage are not independent, so these")
print("  p-values overstate the evidence; report as a lineage-level observation.")
STEP10_EOF
fi

echo ""
echo "============================================================"
echo "GWAS PIPELINE COMPLETED SUCCESSFULLY (Steps 1-10)"
echo "============================================================"
