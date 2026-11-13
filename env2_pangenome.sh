#!/usr/bin/env bash
# =============================================================================
# env_2 — pan-genome analysis, Figure 4 and verification checks
# =============================================================================
# Combines, in one conda environment (default name: env_2):
#   PART A  SPAdes, QUAST, Prokka, Panaroo, IQ-TREE, pyseer pan-genome GWAS,
#           Manhattan/QQ plots            (formerly pangenome.sh)
#   PART B  Figure 4: MIC distributions by gyrA intein status
#                                         (formerly make_figure4.sh)
#   PART C  verification checks behind the figures quoted in the manuscript
#                                         (formerly adhoc_checks.sh)
#
# env_1 must have been run first: this script uses its trimmed reads, Snippy
# output, core alignment and intein status table.
#
# Environment setup: see SETUP.md (section "env_2") before running.
#
# Usage:
#   bash env2_pangenome.sh [KEY=VALUE ...]
#   bash env2_pangenome.sh --help
# =============================================================================
set -euo pipefail

usage() {
    cat <<'USAGE'
env2_pangenome.sh [KEY=VALUE ...]

Settings (environment variable or KEY=VALUE argument; defaults in brackets):

  ENV_NAME        conda environment to run in                 [env_2]
  MHAEM_WORKDIR   env_1 working directory (its outputs are
                  the inputs here)                            [$HOME/mhaem_wgs]
  OUTDIR          output directory for this script            [$HOME/pangenome]
  REF_GBFF        annotated reference genome (GenBank .gbff)  [$MHAEM_WORKDIR/GCF_000340435.2.gbff]
  MIC_AMI         amikacin MIC rank file (sample<TAB>rank)    [$HOME/AMI.txt]
  MIC_CIP         ciprofloxacin MIC rank file                 [$HOME/CIP.txt]
  MIC_RIF         rifampicin MIC rank file                    [$HOME/RIF.txt]
  THREADS_PER_JOB CPU threads per SPAdes/Prokka job           [4]
  PARALLEL_JOBS   jobs run simultaneously                     [1]
  SPADES_MEM      RAM (GB) per SPAdes job                     [16]
  ADHOC_SECTIONS  which PART C checks to run                  [A B C D E F G H I J]

Example:
  bash env2_pangenome.sh MHAEM_WORKDIR=~/mhaem_wgs OUTDIR=~/pangenome \
                         MIC_AMI=~/data/AMI.txt MIC_CIP=~/data/CIP.txt \
                         MIC_RIF=~/data/RIF.txt THREADS_PER_JOB=8
USAGE
}

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
ENV_NAME="${ENV_NAME:-env_2}"
MHAEM_WORKDIR="${MHAEM_WORKDIR:-$HOME/mhaem_wgs}"
OUTDIR="${OUTDIR:-$HOME/pangenome}"
REF_GBFF="${REF_GBFF:-$MHAEM_WORKDIR/GCF_000340435.2.gbff}"
MIC_AMI="${MIC_AMI:-$HOME/AMI.txt}"
MIC_CIP="${MIC_CIP:-$HOME/CIP.txt}"
MIC_RIF="${MIC_RIF:-$HOME/RIF.txt}"

THREADS_PER_JOB="${THREADS_PER_JOB:-4}"
PARALLEL_JOBS="${PARALLEL_JOBS:-1}"
SPADES_MEM="${SPADES_MEM:-16}"
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-500}"
CORE_THRESHOLD="${CORE_THRESHOLD:-0.98}"
ADHOC_SECTIONS="${ADHOC_SECTIONS:-A B C D E F G H I J}"
# ===========================================================================

# Derived paths
TRIMMED="$MHAEM_WORKDIR/trimmed"
MHAEM_QC="$MHAEM_WORKDIR/qc"
SPADES_OUT=$OUTDIR/spades
QUAST_OUT=$OUTDIR/quast
PROKKA_OUT=$OUTDIR/prokka
PANAROO_OUT=$OUTDIR/panaroo
IQTREE_OUT=$OUTDIR/iqtree
PYSEER_OUT=$OUTDIR/pyseer
ANTIBIOTICS=(AMI CIP RIF)
INTEIN_STATUS="$MHAEM_WORKDIR/intein_status.tsv"
FIGURE4_PNG="${FIGURE4_PNG:-$OUTDIR/figure4_intein_mic.png}"

echo "============================================================"
echo "env_2 pan-genome analysis"
echo "  environment    : $ENV_NAME"
echo "  env_1 workdir  : $MHAEM_WORKDIR"
echo "  output dir     : $OUTDIR"
echo "  trimmed reads  : $TRIMMED"
echo "  MIC files      : $MIC_AMI, $MIC_CIP, $MIC_RIF"
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

for f in "$MIC_AMI" "$MIC_CIP" "$MIC_RIF"; do
    [[ -f "$f" ]] || { echo "ERROR: MIC file not found: $f"; exit 1; }
done
[[ -d "$TRIMMED" ]] || {
    echo "ERROR: trimmed reads not found at $TRIMMED — run env1_pipeline.sh first"; exit 1; }

# =============================================================================
# SAMPLE GROUPS — do not edit unless isolate list changes
# =============================================================================

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

REPLICATES=(
    CT-WN18CA0775 CT-WN18CA3057 CT-WN18CD7757 CT-WN17CC9963
)

ATCC="CT-WNTC9548"

# Derived groups
ALL_SAMPLES=("${CLINICAL[@]}" "${REPLICATES[@]}" "$ATCC")
PANAROO_SAMPLES=("${CLINICAL[@]}" "$ATCC")   # 56 clinical + ATCC 29548
# pyseer uses CLINICAL only — ATCC excluded because it has no phenotype entry;
# pyseer automatically uses the intersection of phenotype file and genotype matrix

# =============================================================================
# SETUP
# =============================================================================

PANGENOME_OUT=$OUTDIR

mkdir -p "$SPADES_OUT" "$QUAST_OUT" "$PROKKA_OUT" "$PANAROO_OUT" \
         "$IQTREE_OUT" "$PYSEER_OUT"

# Export all variables needed by Python heredocs
export SPADES_OUT QUAST_OUT PROKKA_OUT PANAROO_OUT IQTREE_OUT PYSEER_OUT
export PANGENOME_OUT MIC_AMI MIC_CIP MIC_RIF THREADS_PER_JOB PARALLEL_JOBS

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "Sample counts:"
log "  Clinical isolates : ${#CLINICAL[@]}"
log "  Replicates        : ${#REPLICATES[@]}"
log "  ATCC 29548        : 1 ($ATCC)"
log "  Total for SPAdes  : ${#ALL_SAMPLES[@]}"
log "  Panaroo input     : ${#PANAROO_SAMPLES[@]} (clinical + ATCC)"

# =============================================================================
# STEP 1: De novo assembly with SPAdes (ALL samples)
# =============================================================================
log "=== [PART A] STEP 1: SPAdes assembly (all ${#ALL_SAMPLES[@]} samples) ==="

run_spades() {
    local sample=$1
    local r1=$TRIMMED/${sample}_R1.trim.fastq.gz
    local r2=$TRIMMED/${sample}_R2.trim.fastq.gz
    local outdir=$SPADES_OUT/$sample

    mkdir -p "$outdir"

    if [ -f "$outdir/contigs.fasta" ]; then
        echo "[SPAdes] SKIP $sample (already done)"
        return 0
    fi

    echo "[SPAdes] Assembling $sample ..."
    spades.py \
        -1 "$r1" \
        -2 "$r2" \
        -o "$outdir" \
        --careful \
        --cov-cutoff auto \
        -t "$THREADS_PER_JOB" \
        -m "$SPADES_MEM" \
        > "$outdir/spades.log" 2>&1 \
        && echo "[SPAdes] DONE  $sample" \
        || echo "[SPAdes] FAIL  $sample — check $outdir/spades.log"
}

export -f run_spades
export TRIMMED SPADES_OUT THREADS_PER_JOB SPADES_MEM

job_count=0
for sample in "${ALL_SAMPLES[@]}"; do
    run_spades "$sample" &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$PARALLEL_JOBS" ]; then
        wait; job_count=0
    fi
done
wait

log "Step 1 complete"

# =============================================================================
# STEP 1b: Filter short contigs (ALL samples)
# =============================================================================
log "=== [PART A] STEP 1b: Filtering contigs < ${MIN_CONTIG_LEN}bp ==="

filter_contigs() {
    local sample=$1
    local raw=$SPADES_OUT/$sample/contigs.fasta
    local filtered=$SPADES_OUT/$sample/contigs_filtered.fasta

    if [ ! -f "$raw" ]; then
        echo "[Filter] No assembly for $sample, skipping"
        return 0
    fi

    if [ -f "$filtered" ]; then
        echo "[Filter] SKIP $sample"
        return 0
    fi

    awk -v min="$MIN_CONTIG_LEN" '
        /^>/ { if (seq!="" && length(seq)>=min) { print header; print seq } header=$0; seq="" }
        !/^>/ { seq=seq$0 }
        END   { if (seq!="" && length(seq)>=min) { print header; print seq } }
    ' "$raw" > "$filtered"

    contig_count=$(grep -c "^>" "$filtered" || true)
    echo "[Filter] $sample: $contig_count contigs kept"
}

export -f filter_contigs
export SPADES_OUT MIN_CONTIG_LEN

job_count=0
for sample in "${ALL_SAMPLES[@]}"; do
    filter_contigs "$sample" &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$PARALLEL_JOBS" ]; then
        wait; job_count=0
    fi
done
wait

log "Contig filtering complete"

# =============================================================================
# STEP 2: Assembly QC with QUAST (ALL samples)
# =============================================================================
log "=== [PART A] STEP 2: QUAST assembly QC ==="

# Collect all filtered assemblies for a single multi-sample QUAST run
ASSEMBLY_FILES=()
ASSEMBLY_LABELS=()
for sample in "${ALL_SAMPLES[@]}"; do
    f=$SPADES_OUT/$sample/contigs_filtered.fasta
    if [ -f "$f" ]; then
        ASSEMBLY_FILES+=("$f")
        ASSEMBLY_LABELS+=("$sample")
    fi
done

# Build --labels string
LABELS_STR=$(IFS=','; echo "${ASSEMBLY_LABELS[*]}")

if [ -f "$QUAST_OUT/report.tsv" ]; then
    log "SKIP QUAST (report.tsv already exists)"
else
    log "Running QUAST on ${#ASSEMBLY_FILES[@]} assemblies ..."
    quast.py \
        "${ASSEMBLY_FILES[@]}" \
        --labels "$LABELS_STR" \
        -o "$QUAST_OUT" \
        --threads $((THREADS_PER_JOB * PARALLEL_JOBS)) \
        > "$QUAST_OUT/quast.log" 2>&1
    log "QUAST complete — summary: $QUAST_OUT/report.tsv"
fi
log "  Key metrics to review: N50, # contigs, total length, largest contig"
log "  Flag any sample with N50 < 20 kb or > 500 contigs for exclusion"

# =============================================================================
# STEP 3: Gene annotation with Prokka (ALL samples)
# =============================================================================
log "=== [PART A] STEP 3: Prokka annotation (all ${#ALL_SAMPLES[@]} samples) ==="

run_prokka() {
    local sample=$1
    local assembly=$SPADES_OUT/$sample/contigs_filtered.fasta
    local outdir=$PROKKA_OUT/$sample

    mkdir -p "$outdir"

    if [ ! -f "$assembly" ]; then
        echo "[Prokka] No assembly for $sample, skipping"
        return 0
    fi

    if [ -f "$outdir/${sample}.gff" ]; then
        echo "[Prokka] SKIP $sample (already done)"
        return 0
    fi

    echo "[Prokka] Annotating $sample ..."
    prokka \
        --outdir "$outdir" \
        --prefix "$sample" \
        --genus Mycobacterium \
        --species haemophilum \
        --kingdom Bacteria \
        --cpus "$THREADS_PER_JOB" \
        --force \
        "$assembly" \
        > "$outdir/prokka.log" 2>&1 \
        && echo "[Prokka] DONE  $sample" \
        || echo "[Prokka] FAIL  $sample — check $outdir/prokka.log"
}

export -f run_prokka
export PROKKA_OUT SPADES_OUT THREADS_PER_JOB

job_count=0
for sample in "${ALL_SAMPLES[@]}"; do
    run_prokka "$sample" &
    job_count=$((job_count + 1))
    if [ "$job_count" -ge "$PARALLEL_JOBS" ]; then
        wait; job_count=0
    fi
done
wait

log "Step 3 complete"

# =============================================================================
# STEP 3b: Pan-genome QC table (one row per clinical isolate)
# Columns: Sample | N50 | Contigs_ge500bp | Total_size_Mb | Largest_contig_bp
#          | GC_pct | Prokka_genes | Trimmed_reads | QC_flag
#
# Pass/fail criteria:
#   FAIL if any of:
#     - Total assembled size > 6.0 Mb (possible contamination)
#     - N50 < 20,000 bp (fragmented)
#     - > 500 contigs (>= 500 bp)
#     - GC content outside 60–74%
# =============================================================================
log "=== [PART A] STEP 3b: Pan-genome QC table ==="

QC_TSV=$OUTDIR/pangenome_qc.tsv
QUAST_TRANSPOSED=$QUAST_OUT/transposed_report.tsv

if [ -f "$QC_TSV" ]; then
    log "SKIP QC table (already exists: $QC_TSV)"
else

# Export CLINICAL array as newline-separated string for Python
CLINICAL_STR=$(printf '%s\n' "${CLINICAL[@]}")

python3 - << STEP3B_EOF
import os, sys, math, csv

SPADES_OUT  = "$SPADES_OUT"
PROKKA_OUT  = "$PROKKA_OUT"
QUAST_TSV   = "$QUAST_TRANSPOSED"
QC_TSV      = "$QC_TSV"
QC_DIR      = "$MHAEM_QC"
CLINICAL_STR = """$CLINICAL_STR"""

clinical = [s.strip() for s in CLINICAL_STR.strip().split("\n") if s.strip()]
print(f"  Building QC table for {len(clinical)} clinical isolates...")

# ---------------------------------------------------------------------------
# Parse QUAST transposed_report.tsv
# ---------------------------------------------------------------------------
quast = {}   # sample -> {N50, contigs, total_length, largest_contig}
if os.path.isfile(QUAST_TSV):
    with open(QUAST_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            s = row.get("Assembly", "").strip()
            if not s:
                continue
            quast[s] = {
                "N50":            row.get("N50", "NA").strip(),
                "contigs":        row.get("# contigs (>= 500 bp)", row.get("# contigs", "NA")).strip(),
                "total_length":   row.get("Total length (>= 500 bp)", row.get("Total length", "NA")).strip(),
                "largest_contig": row.get("Largest contig", "NA").strip(),
            }
    print(f"  QUAST: loaded {len(quast)} entries")
else:
    print(f"  WARNING: QUAST transposed report not found at {QUAST_TSV}")

# ---------------------------------------------------------------------------
# GC content from filtered FASTA
# ---------------------------------------------------------------------------
def calc_gc(fasta_path):
    gc = total = 0
    seq = []
    try:
        with open(fasta_path) as f:
            for line in f:
                if line.startswith(">"):
                    continue
                seq.append(line.strip().upper())
        s = "".join(seq)
        gc    = s.count("G") + s.count("C")
        total = len(s) - s.count("N")
        return round(gc / total * 100, 2) if total > 0 else float("nan")
    except OSError:
        return float("nan")

# ---------------------------------------------------------------------------
# Prokka gene count (CDS lines in .gff)
# ---------------------------------------------------------------------------
def prokka_genes(sample):
    gff = os.path.join(PROKKA_OUT, sample, f"{sample}.gff")
    if not os.path.isfile(gff):
        return "NA"
    count = 0
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if "\tCDS\t" in line:
                count += 1
    return count

# ---------------------------------------------------------------------------
# Trimmed read count from fastp JSON
# Pattern: {QC_DIR}/{sample}.fastp.json
# ---------------------------------------------------------------------------
def trimmed_reads(sample):
    path = os.path.join(QC_DIR, f"{sample}.fastp.json")
    if not os.path.isfile(path):
        return "NA"
    try:
        import json
        with open(path) as f:
            d = json.load(f)
        reads = (d.get("filtering_result", {}).get("passed_filter_reads") or
                 d.get("summary", {}).get("after_filtering", {}).get("total_reads"))
        return reads if reads is not None else "NA"
    except Exception:
        return "NA"

# ---------------------------------------------------------------------------
# Pass/fail logic
# ---------------------------------------------------------------------------
def qc_flag(total_mb, n50, contigs, gc):
    reasons = []
    try:
        if float(total_mb) > 6.0:
            reasons.append("size>6Mb")
    except (ValueError, TypeError):
        pass
    try:
        if float(n50) < 20000:
            reasons.append("N50<20kb")
    except (ValueError, TypeError):
        pass
    try:
        if int(contigs) > 500:
            reasons.append("contigs>500")
    except (ValueError, TypeError):
        pass
    try:
        if not (60.0 <= float(gc) <= 74.0):
            reasons.append(f"GC={gc}%")
    except (ValueError, TypeError):
        pass
    return "FAIL:" + ",".join(reasons) if reasons else "PASS"

# ---------------------------------------------------------------------------
# Build rows
# ---------------------------------------------------------------------------
OUT_COLS = [
    "Sample", "N50_bp", "Contigs_ge500bp", "Total_size_Mb",
    "Largest_contig_bp", "GC_pct", "Prokka_genes", "Trimmed_reads", "QC_flag"
]

rows = []
for sample in clinical:
    q = quast.get(sample, {})
    fasta = os.path.join(SPADES_OUT, sample, "contigs_filtered.fasta")

    n50     = q.get("N50", "NA")
    contigs = q.get("contigs", "NA")
    largest = q.get("largest_contig", "NA")

    # Total size in Mb
    try:
        total_mb = round(int(q.get("total_length", 0)) / 1e6, 3)
    except (ValueError, TypeError):
        total_mb = "NA"

    gc      = calc_gc(fasta)
    genes   = prokka_genes(sample)
    reads   = trimmed_reads(sample)
    flag    = qc_flag(total_mb, n50, contigs, gc)

    rows.append({
        "Sample":           sample,
        "N50_bp":           n50,
        "Contigs_ge500bp":  contigs,
        "Total_size_Mb":    total_mb,
        "Largest_contig_bp": largest,
        "GC_pct":           gc if not (isinstance(gc, float) and math.isnan(gc)) else "NA",
        "Prokka_genes":     genes,
        "Trimmed_reads":    reads,
        "QC_flag":          flag,
    })

with open(QC_TSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

n_pass = sum(1 for r in rows if r["QC_flag"] == "PASS")
n_fail = sum(1 for r in rows if r["QC_flag"].startswith("FAIL"))
print(f"  QC table written: {QC_TSV}")
print(f"  PASS: {n_pass}  FAIL: {n_fail}")
for r in rows:
    if r["QC_flag"].startswith("FAIL"):
        print(f"    FAIL — {r['Sample']}: {r['QC_flag']}")
STEP3B_EOF

fi  # end skip check

log "Step 3b complete — QC table: $QC_TSV"


# =============================================================================
# STEP 4: Pan-genome with Panaroo (56 clinical + ATCC only)
# =============================================================================
log "=== [PART A] STEP 4: Panaroo pan-genome (${#PANAROO_SAMPLES[@]} samples) ==="
log "  Excluded from Panaroo: replicates (${REPLICATES[*]})"

GFF_FILES=()
MISSING_GFF=()
for sample in "${PANAROO_SAMPLES[@]}"; do
    gff=$PROKKA_OUT/$sample/${sample}.gff
    if [ -f "$gff" ]; then
        GFF_FILES+=("$gff")
    else
        MISSING_GFF+=("$sample")
    fi
done

if [ ${#MISSING_GFF[@]} -gt 0 ]; then
    log "WARNING: Missing GFF for: ${MISSING_GFF[*]}"
fi

if [ -f "$PANAROO_OUT/gene_presence_absence.Rtab" ]; then
    log "SKIP Panaroo (gene_presence_absence.Rtab already exists)"
else
log "Running Panaroo with ${#GFF_FILES[@]} isolates ..."
panaroo \
    -i "${GFF_FILES[@]}" \
    -o "$PANAROO_OUT" \
    --clean-mode moderate \
    -a core \
    --core_threshold "$CORE_THRESHOLD" \
    -t $((THREADS_PER_JOB * PARALLEL_JOBS)) \
    > "$PANAROO_OUT/panaroo.log" 2>&1
fi

log "Step 4 complete"
log "  Gene presence/absence : $PANAROO_OUT/gene_presence_absence.Rtab"
log "  Core genome alignment : $PANAROO_OUT/core_gene_alignment.aln"

# =============================================================================
# STEP 5: Phylogenetic tree with IQ-TREE (on Panaroo core genome alignment)
# =============================================================================
log "=== [PART A] STEP 5: IQ-TREE phylogenetic tree ==="

CORE_ALN=$PANAROO_OUT/core_gene_alignment.aln

if [ ! -f "$CORE_ALN" ]; then
    log "ERROR: core_gene_alignment.aln not found — ensure Panaroo ran with -a core"
    exit 1
fi

TREE=$IQTREE_OUT/core_genome.treefile

if [ -f "$TREE" ]; then
    log "SKIP IQ-TREE (treefile already exists: $TREE)"
else
    log "Running IQ-TREE on core genome alignment ..."
    IQTREE_BIN=$(command -v iqtree2 || command -v iqtree || echo "")
    if [ -z "$IQTREE_BIN" ]; then
        log "ERROR: iqtree2/iqtree not found — run: conda install -c bioconda iqtree"
        exit 1
    fi

    "$IQTREE_BIN" \
        -s "$CORE_ALN" \
        -m GTR+G \
        -nt $((THREADS_PER_JOB * PARALLEL_JOBS)) \
        -pre "$IQTREE_OUT/core_genome" \
        -bb 1000 \
        > "$IQTREE_OUT/iqtree.log" 2>&1

    if [ ! -f "$TREE" ]; then
        log "ERROR: IQ-TREE did not produce a treefile at $TREE"
        exit 1
    fi
    log "IQ-TREE complete — tree: $TREE"
fi
log "  ATCC 29548 ($ATCC) is included as a phylogenetic outgroup/anchor"

# =============================================================================
# STEP 6: Kinship matrix from pan-genome presence/absence (similarity_pyseer)
# =============================================================================
log "=== [PART A] STEP 6: Kinship matrix from pan-genome presence/absence ==="

KINSHIP=$PYSEER_OUT/kinship_matrix.tsv
SAMPLES_LIST=$PYSEER_OUT/samples_list.txt
PRESENCE_ABSENCE=$PANAROO_OUT/gene_presence_absence.Rtab

# Write clinical-only sample list (ATCC excluded — no phenotype)
printf '%s\n' "${CLINICAL[@]}" > "$SAMPLES_LIST"
log "  Sample list written: ${#CLINICAL[@]} clinical isolates"

if [ -f "$KINSHIP" ]; then
    log "SKIP kinship matrix (already exists: $KINSHIP)"
else
    similarity_pyseer \
        --pres "$PRESENCE_ABSENCE" \
        "$SAMPLES_LIST" \
        > "$KINSHIP"
    log "Kinship matrix written to $KINSHIP"
fi
log "  Matrix based on pan-genome gene presence/absence (56 clinical isolates)"

# =============================================================================
# STEP 7: Phenotype file preparation (ranks already pre-computed)
# Input files have a header row and include replicates — strip both.
# =============================================================================
log "=== [PART A] STEP 7: Phenotype file preparation (56 clinical isolates) ==="

python3 - <<'PYEOF'
import sys, os

CLINICAL = {
    "CT-WN15CC3812","CT-WN18BA5462","CT-WN15CB4880","CT-WN17CC6853","CT-WN16CC8097",
    "CT-WN15CC0814","CT-WN18CC5430","CT-WN17CA8830","CT-WN17CD0730","CT-WN18BA6316",
    "CT-WN15CC3097","CT-WN18CA1953","CT-WN15CB5439","CT-WN17CD1188","CT-WN18CA1163",
    "CT-WN15CD0483","CT-WN15CD1369","CT-WN18CA7807","CT-WN16BA7594","CT-WN16CC2935",
    "CT-WN16CD2223","CT-WN17CB3038","CT-WN17CA3088","CT-WN17CB5288","CT-WN18CA8070",
    "CT-WN18CB7320","CT-WN17BA5155","CT-WN15CC6561","CT-WN16CA5839","CT-WN17CB9744",
    "CT-WN16CC7074","CT-WN18CD4765","CT-WN17CD6535","CT-WN17BA6921","CT-WN18CD7756",
    "CT-WN15CA7460","CT-WN18BA5406","CT-WN15CC1919","CT-WN15CA4590","CT-WN17CA6664",
    "CT-WN18CA3189","CT-WN17CA4749","CT-WN18BA1483","CT-WN15CA8741","CT-WN17BA6280",
    "CT-WN18CD0191","CT-WN16CD4677","CT-WN15CC3385","CT-WN18CC4066","CT-WN15CA9423",
    "CT-WN15CA4114","CT-WN15CC4010","CT-WN17BA6591","CT-WN16CA1868",
    "CT-WN18CA7566","CT-WN18CA1178"
}

outdir = os.environ['PYSEER_OUT']

jobs = [
    (os.environ['MIC_RIF'], os.path.join(outdir, 'phenotypes_RIF.tsv'), 'RIF'),
    (os.environ['MIC_CIP'], os.path.join(outdir, 'phenotypes_CIP.tsv'), 'CIP'),
    (os.environ['MIC_AMI'], os.path.join(outdir, 'phenotypes_AMI.tsv'), 'AMI'),
]

for infile, outfile, name in jobs:
    if not os.path.isfile(infile):
        print(f"ERROR: phenotype file not found for {name}: {infile}", file=sys.stderr)
        sys.exit(1)
    if os.path.isfile(outfile):
        print(f"  SKIP {name} (phenotype file already exists)")
        continue
    kept = 0
    skipped_header = 0
    skipped_replicate = 0
    with open(infile) as f, open(outfile, 'w') as out:
        out.write("sample\ttrait\n")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                print(f"  WARNING: skipping malformed line: {line}", file=sys.stderr)
                continue
            sample, rank_str = parts[0], parts[1].strip()
            # Skip header row
            if sample.lower() in ('samples', 'sample', 'id'):
                skipped_header += 1
                continue
            # Skip non-clinical samples (replicates, ATCC)
            if sample not in CLINICAL:
                skipped_replicate += 1
                continue
            # Validate rank is an integer
            try:
                int(rank_str)
            except ValueError:
                print(f"  WARNING: non-integer rank '{rank_str}' for {sample} — skipping", file=sys.stderr)
                continue
            out.write(f"{sample}\t{rank_str}\n")
            kept += 1
    print(f"  {name}: {kept} clinical isolates written, {skipped_replicate} non-clinical skipped, {skipped_header} header row stripped")
    if kept != 56:
        print(f"  WARNING: expected 56 clinical isolates for {name}, got {kept}", file=sys.stderr)

print("Phenotype file preparation complete.")
PYEOF

log "Step 7 complete"

# =============================================================================
# STEP 8: Pyseer GWAS (56 clinical isolates only)
# =============================================================================
# --min-af 0.02 : gene clusters present in fewer than 2 of 56 isolates are not
#                 tested (matches the SNP GWAS in gwas11.sh).
log "=== [PART A] STEP 8: Pyseer GWAS (56 clinical, ATCC excluded via phenotype file) ==="

PRESENCE_ABSENCE=$PANAROO_OUT/gene_presence_absence.Rtab
PRESENCE_ABSENCE_CLINICAL=$PYSEER_OUT/gene_presence_absence_clinical.Rtab

# Create clinical-only Rtab (exclude CT-WNTC9548) so pyseer af is based on 56 clinical isolates only
export PRESENCE_ABSENCE PRESENCE_ABSENCE_CLINICAL
if [ -f "$PRESENCE_ABSENCE_CLINICAL" ]; then
    log "SKIP clinical Rtab (already exists)"
else
    python3 - <<'RTAB_EOF'
import os, sys
infile  = os.environ["PRESENCE_ABSENCE"]
outfile = os.environ["PRESENCE_ABSENCE_CLINICAL"]
exclude = {"CT-WNTC9548"}
with open(infile) as fin, open(outfile, "w") as fout:
    header = fin.readline().rstrip("\n").split("\t")
    keep_idx = [i for i, s in enumerate(header) if s not in exclude]
    fout.write("\t".join(header[i] for i in keep_idx) + "\n")
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        fout.write("\t".join(parts[i] for i in keep_idx) + "\n")
print(f"  Clinical Rtab written: {outfile}")
RTAB_EOF
fi

for antibiotic in "${ANTIBIOTICS[@]}"; do
    pheno=$PYSEER_OUT/phenotypes_${antibiotic}.tsv
    results=$PYSEER_OUT/${antibiotic}_results.txt
    threshold_file=$PYSEER_OUT/${antibiotic}_threshold.txt
    patterns_file=$PYSEER_OUT/${antibiotic}_patterns.txt

    if [ ! -f "$pheno" ]; then
        log "ERROR: Ranked phenotype file not found: $pheno — check Step 7"
        continue
    fi

    if [ -f "$results" ]; then
        log "  SKIP $antibiotic Pyseer (results already exist)"
        continue
    fi

    log "Running Pyseer for $antibiotic ..."
    pyseer \
        --lmm \
        --phenotypes "$pheno" \
        --pres "$PRESENCE_ABSENCE_CLINICAL" \
        --similarity "$KINSHIP" \
        --min-af 0.02 \
        --cpu $((THREADS_PER_JOB * PARALLEL_JOBS)) \
        --output-patterns "$patterns_file" \
        > "$results" \
        2> "$PYSEER_OUT/${antibiotic}_pyseer.log"

    log "  $antibiotic: $(( $(wc -l < "$results") - 1 )) gene clusters tested"

    # Bonferroni threshold = 0.05 / number of unique presence/absence patterns
    if [ ! -s "$patterns_file" ]; then
        log "ERROR: $patterns_file missing or empty — pyseer did not write patterns"
        exit 1
    fi
    N_PAT=$(sort -u "$patterns_file" | wc -l)
    python3 -c "print(0.05/$N_PAT)" > "$threshold_file"
    log "  $N_PAT unique patterns -> Bonferroni p < $(cat "$threshold_file")"

    log "  Top 20 hits for $antibiotic:"
    head -n 1 "$results"
    awk 'NR > 1' "$results" | sort -g -k4,4 | head -n 20 || true
done

# =============================================================================
# STEP 9: GWAS summary — diagnostic statistics (Python, embedded)
# Requires: pandas openpyxl scipy
#   conda install -c conda-forge pandas openpyxl scipy
# =============================================================================
log "=== [PART A] STEP 9: GWAS summary + diagnostic statistics ==="

# Check required Python packages
python3 - <<'PYEOF'
import importlib, sys
missing = [pkg for pkg in ("pandas","openpyxl","scipy") if importlib.util.find_spec(pkg) is None]
if missing:
    print(f"ERROR: missing Python packages: {', '.join(missing)}", file=sys.stderr)
    print("Run: conda install -c conda-forge pandas openpyxl scipy", file=sys.stderr)
    sys.exit(1)
PYEOF

GWAS_XLSX=$PANGENOME_OUT/gwas_summary.xlsx
if [ -f "$GWAS_XLSX" ]; then
    log "SKIP GWAS summary (gwas_summary.xlsx already exists)"
else

python3 - <<'PYEOF'
"""
GWAS post-processing for M. haemophilum pan-genome GWAS.
Reads Pyseer results for AMI, CIP, RIF; keeps nominal hits (lrt-pvalue < 0.05)
and flags those passing Bonferroni (0.05 / n unique patterns) in 'bonferroni_sig';
computes sensitivity, specificity, PPV, NPV (Wilson 95% CI) per gene cluster;
direction follows beta sign (positive = resistance, negative = susceptibility).
Writes gwas_summary.xlsx — one sheet per antibiotic.
"""

import os, math, sys
import pandas as pd

PYSEER_OUT  = os.environ["PYSEER_OUT"]
PANAROO_OUT = os.environ["PANAROO_OUT"]
OUTPUT_XLSX = os.path.join(os.environ["PANGENOME_OUT"], "gwas_summary.xlsx")
ANTIBIOTICS = ["AMI", "CIP", "RIF"]

# Resistance threshold: ordinal rank at/above which isolate is non-susceptible
# (intermediate MIC is grouped with resistant)
RESISTANCE_THRESHOLD = {"RIF": 5, "CIP": 5, "AMI": 6}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def wilson_ci(k, n, z=1.96):
    """Wilson score 95% CI for proportion k/n. Returns (proportion, lo, hi)."""
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    p = k / n
    centre = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    margin  = (z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / (1 + z**2 / n)
    return (round(p, 4),
            round(max(0.0, centre - margin), 4),
            round(min(1.0, centre + margin), 4))


def diagnostic_stats(gene_vec, resistant_vec, direction):
    """
    2×2 table + diagnostic statistics.

    direction == 'resistance'  : gene present predicts RESISTANT
        TP = present & resistant   FP = present & susceptible
        TN = absent  & susceptible FN = absent  & resistant

    direction == 'susceptible' : gene present predicts SUSCEPTIBLE (beta < 0)
        Flip: susceptible is the positive condition
        TP = present & susceptible FP = present & resistant
        TN = absent  & resistant   FN = absent  & susceptible
    """
    if direction == "resistance":
        tp = sum(g == 1 and r == 1 for g, r in zip(gene_vec, resistant_vec))
        fp = sum(g == 1 and r == 0 for g, r in zip(gene_vec, resistant_vec))
        tn = sum(g == 0 and r == 0 for g, r in zip(gene_vec, resistant_vec))
        fn = sum(g == 0 and r == 1 for g, r in zip(gene_vec, resistant_vec))
    else:
        tp = sum(g == 1 and r == 0 for g, r in zip(gene_vec, resistant_vec))
        fp = sum(g == 1 and r == 1 for g, r in zip(gene_vec, resistant_vec))
        tn = sum(g == 0 and r == 1 for g, r in zip(gene_vec, resistant_vec))
        fn = sum(g == 0 and r == 0 for g, r in zip(gene_vec, resistant_vec))

    n_resistant   = sum(resistant_vec)
    n_susceptible = len(resistant_vec) - n_resistant

    sens, sens_lo, sens_hi = wilson_ci(tp, tp + fn)
    spec, spec_lo, spec_hi = wilson_ci(tn, tn + fp)
    ppv,  ppv_lo,  ppv_hi  = wilson_ci(tp, tp + fp)
    npv,  npv_lo,  npv_hi  = wilson_ci(tn, tn + fn)

    return dict(
        TP=tp, FP=fp, TN=tn, FN=fn,
        n_resistant=n_resistant, n_susceptible=n_susceptible,
        sensitivity=sens,    sensitivity_ci95_lo=sens_lo,  sensitivity_ci95_hi=sens_hi,
        specificity=spec,    specificity_ci95_lo=spec_lo,  specificity_ci95_hi=spec_hi,
        PPV=ppv,             PPV_ci95_lo=ppv_lo,            PPV_ci95_hi=ppv_hi,
        NPV=npv,             NPV_ci95_lo=npv_lo,            NPV_ci95_hi=npv_hi,
    )


COLUMNS = [
    "variant", "af", "filter-pvalue", "lrt-pvalue", "beta", "beta-std-err",
    "variant_h2", "notes", "bonferroni_threshold", "bonferroni_sig",
    "CHROM", "POS", "TYPE", "REF", "ALT", "FTYPE", "STRAND",
    "NT_POS", "AA_POS", "EFFECT", "LOCUS_TAG", "GENE", "PRODUCT",
    "direction",
    "TP", "FP", "TN", "FN",
    "n_resistant", "n_susceptible",
    "sensitivity", "sensitivity_ci95_lo", "sensitivity_ci95_hi",
    "specificity", "specificity_ci95_lo", "specificity_ci95_hi",
    "PPV", "PPV_ci95_lo", "PPV_ci95_hi",
    "NPV", "NPV_ci95_lo", "NPV_ci95_hi",
]

# ---------------------------------------------------------------------------
# Load gene presence/absence matrix (Panaroo Rtab)
# ---------------------------------------------------------------------------
rtab_path = os.path.join(PANAROO_OUT, "gene_presence_absence.Rtab")
print(f"Loading gene presence/absence matrix: {rtab_path}")
if not os.path.isfile(rtab_path):
    print(f"ERROR: Rtab matrix not found — ensure Step 4 (Panaroo) completed", file=sys.stderr)
    sys.exit(1)

rtab = pd.read_csv(rtab_path, sep="\t", index_col=0)
print(f"  {rtab.shape[0]} gene clusters × {rtab.shape[1]} samples")

# ---------------------------------------------------------------------------
# Process each antibiotic
# ---------------------------------------------------------------------------
writer = pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl")

for ab in ANTIBIOTICS:
    print(f"\n{'='*60}\nProcessing {ab} ...")

    results_path   = os.path.join(PYSEER_OUT, f"{ab}_results.txt")
    threshold_path = os.path.join(PYSEER_OUT, f"{ab}_threshold.txt")
    pheno_path     = os.path.join(PYSEER_OUT, f"phenotypes_{ab}.tsv")

    # --- Bonferroni threshold ---
    threshold = None
    if os.path.isfile(threshold_path):
        with open(threshold_path) as f:
            for line in f:
                line = line.strip()
                if line.lower().startswith("threshold"):
                    try:
                        threshold = float(line.split()[-1])
                        break
                    except ValueError:
                        pass
                else:
                    try:
                        threshold = float(line.split()[-1])
                        break
                    except ValueError:
                        continue
    if threshold is None:
        # Fall back to Bonferroni based on number of variants in results file
        try:
            n_tests = sum(1 for _ in open(results_path)) - 1
            threshold = 0.05 / n_tests
            print(f"  Bonferroni threshold calculated from results: {threshold:.2e} (0.05/{n_tests})")
        except Exception:
            threshold = 0.05
            print(f"  WARNING: could not parse threshold, using 0.05")
    print(f"  Bonferroni threshold: {threshold:.2e}")

    # --- Pyseer results ---
    if not os.path.isfile(results_path):
        print(f"  ERROR: results file not found: {results_path}")
        no_hit_df = pd.DataFrame(columns=COLUMNS)
        no_hit_df.loc[0, "variant"] = "no hits"
        no_hit_df.to_excel(writer, sheet_name=ab, index=False)
        continue

    res = pd.read_csv(results_path, sep="\t")
    res.columns = [c.strip() for c in res.columns]
    if "beta-std-err" not in res.columns and "beta_std-err" in res.columns:
        res.rename(columns={"beta_std-err": "beta-std-err"}, inplace=True)

    hits = res[res["lrt-pvalue"] < 0.05].copy()  # exploratory: all nominal hits
    hits["bonferroni_threshold"] = threshold
    hits["bonferroni_sig"] = hits["lrt-pvalue"] < threshold
    print(f"  Nominal hits (p < 0.05): {len(hits)}; "
          f"Bonferroni-significant (p < {threshold:.2e}): {int(hits['bonferroni_sig'].sum())}")

    if hits.empty:
        print(f"  No nominal associations (p < 0.05) for {ab}")
        no_hit_df = pd.DataFrame(columns=COLUMNS)
        no_hit_df.loc[0, "variant"] = "no hits"
        no_hit_df.to_excel(writer, sheet_name=ab, index=False)
        continue

    # --- phenotype (ordinal ranks → resistant/susceptible binary) ---
    pheno = pd.read_csv(pheno_path, sep="\t", names=["sample", "rank"], skiprows=1)
    res_threshold = RESISTANCE_THRESHOLD[ab]
    pheno["resistant"] = (pheno["rank"] >= res_threshold).astype(int)
    pheno_dict = dict(zip(pheno["sample"], pheno["resistant"]))
    print(f"  Phenotype: {pheno['resistant'].sum()} resistant, "
          f"{(pheno['resistant'] == 0).sum()} susceptible "
          f"(non-susceptible threshold: rank >= {res_threshold})")

    # --- build output rows ---
    rows = []
    for _, hit in hits.iterrows():
        gene = hit["variant"]
        beta = hit.get("beta", float("nan"))

        row = {col: hit.get(col, float("nan")) for col in COLUMNS if col in hit.index}

        # annotation columns absent in pan-genome output
        for col in ["CHROM","POS","TYPE","REF","ALT","FTYPE","STRAND",
                    "NT_POS","AA_POS","EFFECT","LOCUS_TAG","GENE","PRODUCT"]:
            if col not in row or (not isinstance(row.get(col), str) and math.isnan(float(row.get(col, float("nan"))))):
                row[col] = "NA"

        direction = "resistance" if (not math.isnan(float(beta)) and float(beta) > 0) else "susceptible"
        row["direction"] = direction

        if gene in rtab.index:
            gene_row = rtab.loc[gene]
            shared = [s for s in pheno_dict.keys() if s in gene_row.index]
            gene_vec  = [int(gene_row[s]) for s in shared]
            resist_vec = [pheno_dict[s] for s in shared]
            dstats = diagnostic_stats(gene_vec, resist_vec, direction)
        else:
            print(f"  WARNING: gene '{gene}' not found in Rtab matrix")
            dstats = {k: float("nan") for k in [
                "TP","FP","TN","FN","n_resistant","n_susceptible",
                "sensitivity","sensitivity_ci95_lo","sensitivity_ci95_hi",
                "specificity","specificity_ci95_lo","specificity_ci95_hi",
                "PPV","PPV_ci95_lo","PPV_ci95_hi",
                "NPV","NPV_ci95_lo","NPV_ci95_hi"]}

        row.update(dstats)
        rows.append(row)

    out_df = pd.DataFrame(rows, columns=COLUMNS)
    out_df.to_excel(writer, sheet_name=ab, index=False)
    print(f"  Written {len(out_df)} rows to sheet '{ab}'")

writer.close()
print(f"\nDone. GWAS summary written to: {OUTPUT_XLSX}")
PYEOF

fi  # end skip check for GWAS summary

log "Step 9 complete — gwas_summary.xlsx written to $PANGENOME_OUT/gwas_summary.xlsx"

# =============================================================================
# STEP 10: Manhattan plots and QQ plots (embedded Python)
# Requires: matplotlib numpy
#   conda install -c conda-forge matplotlib numpy
# =============================================================================
log "=== [PART A] STEP 10: Manhattan and QQ plots ==="

MANHATTAN_PNG=$PANGENOME_OUT/manhattan_plots.png
QQ_PNG=$PANGENOME_OUT/qq_plots.png

if [ -f "$MANHATTAN_PNG" ] && [ -f "$QQ_PNG" ]; then
    log "SKIP plots (manhattan_plots.png and qq_plots.png already exist)"
else

python3 - <<'PYEOF'
import os, sys, math
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
except ImportError:
    print("ERROR: matplotlib not found. Run: conda install -c conda-forge matplotlib", file=sys.stderr)
    sys.exit(1)

PYSEER_OUT  = os.environ["PYSEER_OUT"]
PANGENOME_OUT = os.environ["PANGENOME_OUT"]
ANTIBIOTICS = ["AMI", "CIP", "RIF"]
AB_COLORS   = {"AMI": "#2166ac", "CIP": "#d6604d", "RIF": "#4dac26"}
AB_LABELS   = {"AMI": "Amikacin (AMI)", "CIP": "Ciprofloxacin (CIP)", "RIF": "Rifampicin (RIF)"}

# -----------------------------------------------------------------------
# Load results
# -----------------------------------------------------------------------
data = {}
for ab in ANTIBIOTICS:
    fpath = os.path.join(PYSEER_OUT, f"{ab}_results.txt")
    if not os.path.isfile(fpath):
        print(f"WARNING: {fpath} not found — skipping {ab}")
        continue
    import csv
    pvals = []
    with open(fpath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                pvals.append(float(row["lrt-pvalue"]))
            except (KeyError, ValueError):
                continue
    if pvals:
        data[ab] = sorted(pvals)
    print(f"  {ab}: {len(pvals)} p-values loaded")

if not data:
    print("ERROR: No pyseer results found — run Steps 1-8 first", file=sys.stderr)
    sys.exit(1)

# Per-drug Bonferroni threshold = 0.05 / n unique patterns (written in Step 8).
# Falls back to 0.05 / n tested variants if the threshold file is missing.
bonf_thr = {}
for ab, v in data.items():
    tpath = os.path.join(PYSEER_OUT, f"{ab}_threshold.txt")
    if os.path.isfile(tpath):
        with open(tpath) as tf:
            bonf_thr[ab] = float(tf.read().strip())
        print(f"  Bonferroni threshold [{ab}]: {bonf_thr[ab]:.2e}  (0.05 / unique patterns)")
    else:
        bonf_thr[ab] = 0.05 / len(v)
        print(f"  WARNING: {tpath} missing — Bonferroni [{ab}]: {bonf_thr[ab]:.2e} (0.05 / {len(v)} variants)")

# -----------------------------------------------------------------------
# Load gene names for labelling top hits on Manhattan
# -----------------------------------------------------------------------
gene_data = {}  # ab -> list of (pval, gene_name) sorted by pval ascending
for ab in ANTIBIOTICS:
    fpath = os.path.join(PYSEER_OUT, f"{ab}_results.txt")
    if not os.path.isfile(fpath):
        continue
    import csv
    rows = []
    with open(fpath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                p = float(row["lrt-pvalue"])
                g = row.get("variant", "")
                rows.append((p, g))
            except (KeyError, ValueError):
                continue
    gene_data[ab] = sorted(rows, key=lambda x: x[0])

# -----------------------------------------------------------------------
# FIGURE 1: Manhattan plots (3 stacked panels)
# -----------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=False)
fig.suptitle("Pan-Genome GWAS: Manhattan Plots\n"
             "M. haemophilum (n=56 clinical isolates)",
             fontsize=13, fontweight='bold', y=0.98)

for ax, ab in zip(axes, ANTIBIOTICS):
    if ab not in gene_data:
        ax.set_title(f"{AB_LABELS[ab]} — no data")
        continue

    rows = gene_data[ab]
    pvals_sorted = [r[0] for r in rows]
    genes_sorted = [r[1] for r in rows]

    # Clamp very small p-values for display
    logp = [-math.log10(max(p, 1e-300)) for p in pvals_sorted]
    x    = list(range(len(logp)))

    # Alternating colors for visual separation
    colors = [AB_COLORS[ab] if i % 2 == 0 else "#999999" for i in range(len(logp))]

    ax.scatter(x, logp, c=colors, s=8, alpha=0.7, linewidths=0, rasterized=True)

    # Reference lines — suggestive first so it appears on top in the legend
    ab_bonf   = bonf_thr.get(ab, 0.05 / len(logp))
    bonf_line = -math.log10(ab_bonf)
    nom_line  = -math.log10(0.05)
    ax.axhline(nom_line,  color='orange', linestyle=':', linewidth=1.0,
               label="Suggestive p < 0.05")
    ax.axhline(bonf_line, color='red',   linestyle='--', linewidth=1.2,
               label=f"Bonferroni p < {ab_bonf:.2e}".replace("e-05", " × 10⁻⁵").replace("e-04", " × 10⁻⁴"))

    ax.set_ylabel("−log₁₀(p)", fontsize=9)
    ax.set_title(f"{AB_LABELS[ab]}", fontsize=10)
    ax.set_xlim(-10, len(logp) + 10)
    ax.set_ylim(0, max(max(logp) * 1.15, bonf_line * 1.2))
    ax.legend(fontsize=7, loc='upper right')
    ax.tick_params(axis='x', labelbottom=False)  # hide x tick labels (cluster indices)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[-1].set_xlabel("Gene cluster (sorted by association p-value)", fontsize=9)
axes[-1].tick_params(axis='x', labelbottom=True)

plt.tight_layout(rect=[0, 0, 1, 0.96])
out_manhattan = os.path.join(PANGENOME_OUT, "manhattan_plots.png")
plt.savefig(out_manhattan, dpi=180, bbox_inches='tight')
plt.close()
print(f"  Manhattan plot saved: {out_manhattan}")

# Individual Manhattan plots per drug
for ab in ANTIBIOTICS:
    if ab not in gene_data:
        continue
    rows_ab = gene_data[ab]
    logp_ab = [-math.log10(max(r[0], 1e-300)) for r in rows_ab]
    colors_ab = [AB_COLORS[ab] if i % 2 == 0 else "#999999" for i in range(len(logp_ab))]
    ab_bonf_i = bonf_thr.get(ab, 0.05 / len(logp_ab))
    bonf_line = -math.log10(ab_bonf_i)
    nom_line  = -math.log10(0.05)

    fig_i, ax_i = plt.subplots(figsize=(12, 4))
    ax_i.scatter(range(len(logp_ab)), logp_ab, c=colors_ab, s=8, alpha=0.7,
                 linewidths=0, rasterized=True)
    ax_i.axhline(nom_line,  color='orange', linestyle=':', linewidth=1.0,
                 label="Suggestive p < 0.05")
    ax_i.axhline(bonf_line, color='red', linestyle='--', linewidth=1.2,
                 label=f"Bonferroni p < {ab_bonf_i:.2e}".replace("e-05", " × 10⁻⁵").replace("e-04", " × 10⁻⁴"))
    ax_i.set_xlabel("Gene cluster (sorted by association p-value)", fontsize=9)
    ax_i.set_ylabel("−log₁₀(p)", fontsize=9)
    ax_i.set_title(f"{AB_LABELS[ab]}", fontsize=11, fontweight='bold')
    ax_i.set_xlim(-10, len(logp_ab) + 10)
    ax_i.set_ylim(0, max(max(logp_ab) * 1.15, bonf_line * 1.2))
    ax_i.legend(fontsize=8, loc='upper right')
    ax_i.spines['top'].set_visible(False)
    ax_i.spines['right'].set_visible(False)
    plt.tight_layout()
    out_i = os.path.join(PANGENOME_OUT, f"manhattan_{ab}.png")
    fig_i.savefig(out_i, dpi=180, bbox_inches='tight')
    plt.close(fig_i)
    print(f"  Individual Manhattan plot saved: {out_i}")

# -----------------------------------------------------------------------
# FIGURE 2: QQ plots (3 panels)
# -----------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 3, figsize=(13, 4.5))
fig2.suptitle("Pan-Genome GWAS: Quantile-Quantile Plots\n"
              "M. haemophilum (n=56 clinical isolates)",
              fontsize=13, fontweight='bold')

for ax, ab in zip(axes2, ANTIBIOTICS):
    if ab not in data:
        ax.set_title(f"{AB_LABELS[ab]} — no data")
        continue

    obs  = sorted(data[ab])           # ascending p-values
    n    = len(obs)
    exp  = [(i + 0.5) / n for i in range(n)]   # uniform expected quantiles

    obs_logp = [-math.log10(max(p, 1e-300)) for p in obs]
    exp_logp = [-math.log10(max(p, 1e-300)) for p in exp]

    # Genomic inflation factor lambda (using median)
    median_obs = np.median(obs)
    median_exp = np.median(exp)
    lam = math.log(median_obs) / math.log(median_exp) if median_exp != 1.0 else 1.0
    # More standard: lambda = chi2_obs_median / chi2_expected_median
    # Use scipy if available, else approximate
    try:
        from scipy.stats import chi2
        chi2_obs = chi2.isf(median_obs, df=1)
        chi2_exp = chi2.isf(0.5, df=1)
        lam = chi2_obs / chi2_exp
    except ImportError:
        pass

    ax.scatter(exp_logp, obs_logp, s=10, alpha=0.6,
               color=AB_COLORS[ab], linewidths=0, rasterized=True)

    # Diagonal reference line
    max_val = max(max(exp_logp), max(obs_logp)) * 1.05
    ax.plot([0, max_val], [0, max_val], 'k--', linewidth=1.0, label='Expected')

    # Bonferroni reference (per-drug threshold)
    ab_bonf_qq = bonf_thr.get(ab, 0.05)
    bonf_line_qq = -math.log10(ab_bonf_qq)
    ax.axhline(bonf_line_qq, color='red', linestyle='--',
               linewidth=0.9, alpha=0.7,
               label=f"Bonferroni (p < {ab_bonf_qq:.2e})".replace("e-05", " × 10⁻⁵").replace("e-04", " × 10⁻⁴"))

    ax.set_xlabel("Expected −log₁₀(p)", fontsize=9)
    ax.set_ylabel("Observed −log₁₀(p)", fontsize=9)
    ax.set_title(f"{AB_LABELS[ab]}\nλ = {lam:.3f}", fontsize=10)
    ax.set_xlim(0, max(exp_logp) * 1.1)
    # y-axis always includes the Bonferroni line so it is visible even when
    # no cluster approaches it (matches the SNP QQ plots from gwas11.sh Step 8b)
    ax.set_ylim(0, max(max(obs_logp), bonf_line_qq) * 1.12)
    ax.legend(fontsize=7, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
out_qq = os.path.join(PANGENOME_OUT, "qq_plots.png")
plt.savefig(out_qq, dpi=180, bbox_inches='tight')
plt.close()
print(f"  QQ plot saved: {out_qq}")
print("\nStep 10 complete.")
PYEOF

fi  # end skip check for plots

log "Step 10 complete — plots written to $PANGENOME_OUT/"

# =============================================================================
# DONE
# =============================================================================
log "=== Pipeline complete ==="
echo ""
echo "Output summary:"
echo "  Assemblies (all 61)    : $SPADES_OUT/<sample>/contigs_filtered.fasta"
echo "  QUAST QC               : $QUAST_OUT/report.tsv"
echo "  Annotations (all 61)   : $PROKKA_OUT/<sample>/<sample>.gff"
echo "  Pan-genome (57)        : $PANAROO_OUT/gene_presence_absence.Rtab"
echo "  Core alignment (57)    : $PANAROO_OUT/core_gene_alignment.aln"
echo "  IQ-TREE tree (57)      : $IQTREE_OUT/core_genome.treefile"
echo "  Kinship matrix (56)    : $PYSEER_OUT/kinship_matrix.tsv"
echo "  GWAS results (56)      : $PYSEER_OUT/<antibiotic>_results.txt"
echo "  GWAS patterns          : $PYSEER_OUT/<antibiotic>_patterns.txt"
echo "  GWAS thresholds        : $PYSEER_OUT/<antibiotic>_threshold.txt (0.05 / unique patterns)"
echo "  GWAS summary (xlsx)    : $PANGENOME_OUT/gwas_summary.xlsx"
echo "  Manhattan plots        : $PANGENOME_OUT/manhattan_plots.png"
echo "  QQ plots               : $PANGENOME_OUT/qq_plots.png"
echo ""
echo "Replicate QC (separate step — not in this pipeline):"
echo "  Compare replicate pairs using Mash distance or SNP distance"
echo "  from the core genome alignment to confirm within-patient concordance."

# ###########################################################################
# PART B — Figure 4: MIC distributions by gyrA intein status
# (formerly make_figure4.sh)
# ###########################################################################

INPUT="$INTEIN_STATUS"
OUTPUT="$FIGURE4_PNG"

if [[ ! -s "$INPUT" ]]; then
    echo "SKIP Figure 4 — $INPUT not found (run env1_pipeline.sh PART B Step 10 first)"
else
    echo "============================================================"
    echo "[PART B] Figure 4: MIC distributions by gyrA intein status"
    echo "  input : $INPUT"
    echo "  output: $OUTPUT"
    echo "============================================================"
python3 - "$INPUT" "$OUTPUT" << 'FIG4_EOF'
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = sys.argv[1]
out = sys.argv[2]

df = pd.read_csv(src, sep="\t", dtype=str)
df = df[df["Group"] == "clinical"]                       # 56 primary isolates only
pos = df[df["Intein_status"] == "present"]
neg = df[df["Intein_status"] == "absent"]

# MIC scale of the SLOMYCO panel (Supplementary Table 1) and CLSI M24S
# non-susceptibility breakpoint (I + R) for each drug
PANELS = [
    ("CIP", "Ciprofloxacin (CIP)", ["≤0.12", "0.25", "0.5", "1", "2", "4", "8", "16", ">16"], "2"),
    ("AMI", "Amikacin (AMI)",      ["≤1", "2", "4", "8", "16", "32", "64", ">64"],            "32"),
    ("RIF", "Rifampicin (RIF)",    ["≤0.12", "0.25", "0.5", "1", "2", "4", "8", ">8"],        "2"),
]
C_POS, C_NEG = "#2a78d6", "#eb6834"

def norm(v):
    """Normalise MIC strings so '4', '4.0' and '<=1' match the panel labels."""
    v = str(v).strip().replace("<=", "≤").replace(">=", "≥")
    try:
        f = float(v); return f"{f:g}"
    except ValueError:
        return v

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), dpi=300, sharey=True)
for ax, (drug, title, scale, bp) in zip(axes, PANELS):
    cp = pos[f"{drug}_MIC"].map(norm).value_counts()
    cn = neg[f"{drug}_MIC"].map(norm).value_counts()
    unknown = (set(cp.index) | set(cn.index)) - set(scale)
    if unknown:
        sys.exit(f"{drug}: MIC values not on the panel scale: {unknown}")
    x = range(len(scale)); w = 0.4
    for off, counts, n, col, lab in [(-w / 2, cp, len(pos), C_POS, f"Intein-positive (n={len(pos)})"),
                                     (+w / 2, cn, len(neg), C_NEG, f"Intein-negative (n={len(neg)})")]:
        k = [int(counts.get(m, 0)) for m in scale]
        pct = [100 * v / n for v in k]
        bars = ax.bar([i + off for i in x], pct, width=w * 0.95, color=col, label=lab,
                      edgecolor="white", linewidth=0.8)
        for b, v in zip(bars, k):
            if v:
                ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5, str(v),
                        ha="center", va="bottom", fontsize=8, color="#333333")
    bpx = scale.index(bp) - 0.5                        # line between S and NS categories
    ax.axvline(bpx, ls="--", color="#444444", lw=1)
    ax.text(bpx + 0.08, 97, "Non-susceptible →", fontsize=8, color="#444444", va="top")
    ax.set_xticks(list(x)); ax.set_xticklabels(scale, fontsize=9)
    ax.set_xlabel("MIC (µg/mL)", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", color="#e5e5e5", lw=0.6); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("Isolates in group (%)", fontsize=10)
# legend inside the first panel, upper left, kept compact so it clears the
# breakpoint line; the clearance is checked below
leg = axes[0].legend(frameon=False, fontsize=9, loc="upper left",
                     bbox_to_anchor=(0.01, 0.97),
                     handlelength=1.1, handletextpad=0.5,
                     labelspacing=0.35, borderpad=0.2, borderaxespad=0.0)
fig.tight_layout()
fig.canvas.draw()
# verify the legend does not reach the leftmost breakpoint line (CIP)
_bp_x = PANELS[0][2].index(PANELS[0][3]) - 0.5
_line_px = axes[0].transData.transform((_bp_x, 0))[0]
_leg_px = leg.get_window_extent().x1
if _leg_px >= _line_px:
    print(f"  WARNING: legend reaches the CIP breakpoint line "
          f"({_leg_px:.0f}px vs {_line_px:.0f}px) - widen figsize or reduce fontsize")
else:
    print(f"  legend clears the breakpoint line by {_line_px - _leg_px:.0f} px")
fig.savefig(out)
print(f"Saved {out}: intein-positive n={len(pos)}, intein-negative n={len(neg)}")
FIG4_EOF
    echo "[PART B] Figure 4 complete: $OUTPUT"
fi

# ###########################################################################
# PART C — verification checks
# (formerly adhoc_checks.sh)
#
# These are read-only checks; a failing check reports and the script carries
# on, so -e is relaxed for this part.
# ###########################################################################
set +e

WORKDIR="$MHAEM_WORKDIR"
PANGENOME="$OUTDIR"
OUT="$WORKDIR/adhoc"
mkdir -p "$OUT"

QC="$WORKDIR/qc"
SNIPPY="$WORKDIR/snippy"
GBFF="$REF_GBFF"
CORE_ALN="$WORKDIR/core/core.full.aln"

EXCLUDED=""                              # optional: a library to exclude
ATCC="CT-WNTC9548"
REPLICATES="CT-WN18CA0775,CT-WN18CA3057,CT-WN18CD7757,CT-WN17CC9963"

# rrs coordinates in NZ_CP011883.2 (B586_RS14090), confirmed in section H
RRS_START=3008819
RRS_END=3010350

# Clusters cross-tabulated against intein status in section J
CLUSTERS_J="${CLUSTERS_J:-group_727 group_1448 group_2225 group_3128 group_3133}"

# gyrA coordinates (B586_RS00030) and the intein within it
GYRA_START=7266
GYRA_END=11045

WANTED="$ADHOC_SECTIONS"
run_section() { [[ " $WANTED " == *" $1 "* ]]; }
hdr() { echo; echo "============================================================"; echo "$1"; echo "============================================================"; }

# -----------------------------------------------------------------------------
if run_section A; then
hdr "A. Software versions"
echo "-- recorded in output files (the versions that produced the results) --"
B1=$(ls -d "$SNIPPY"/*/ 2>/dev/null | head -1)
[[ -n "${B1:-}" ]] && samtools view -H "$B1/snps.bam" 2>/dev/null | grep '^@PG' | grep -o 'PN:bwa.*VN:[^[:space:]]*' | head -1
[[ -n "${B1:-}" ]] && grep -m1 -i '^##source=freebayes' "$B1/snps.raw.vcf" 2>/dev/null
[[ -n "${B1:-}" ]] && grep -m1 -i '^##SnpEffVersion' "$B1/snps.vcf" 2>/dev/null
[[ -f "$CORE_ALN.log" ]] && head -1 "$CORE_ALN.log"
[[ -f "$QC/multiqc_report.html" ]] && grep -m1 -o 'MultiQC v[0-9.]*' "$QC/multiqc_report.html"
F1=$(ls "$QC"/*_fastqc.zip 2>/dev/null | head -1)
[[ -n "${F1:-}" ]] && python3 -c "
import zipfile,sys
z=zipfile.ZipFile(sys.argv[1]); print(z.read(z.namelist()[0].split('/')[0]+'/fastqc_data.txt').decode().splitlines()[0])" "$F1"
echo "-- currently installed (verify these match the run) --"
for t in "fastp --version" "samtools --version" "bcftools --version" "pyseer --version"; do
    $t 2>&1 | head -1
done
python3 -c "import Bio; print('biopython', Bio.__version__)" 2>/dev/null || \
  conda run -n "$ENV_NAME" python -c "import Bio; print('biopython', Bio.__version__, '(pangenome env)')" 2>/dev/null
fi

# -----------------------------------------------------------------------------
if run_section B; then
hdr "B. Read QC summary from fastp JSON (clinical isolates)"
python3 - "$QC" "$REPLICATES" "$ATCC" "$EXCLUDED" "$OUT/read_qc_per_sample.tsv" <<'PY'
import json, glob, os, sys, statistics as st
qc, reps, atcc, excl, out = sys.argv[1:6]
reps = set(r for r in reps.split(",") if r)
rows = []
for f in sorted(glob.glob(os.path.join(qc, "*.fastp.json"))):
    s = os.path.basename(f)[:-len(".fastp.json")]
    if s == excl: continue
    d = json.load(open(f))
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    rows.append(dict(Sample=s,
        Group="QC strain" if s == atcc else "Replicate" if s in reps else "Clinical",
        Raw_reads=b["total_reads"], Trimmed_reads=a["total_reads"],
        Pass_rate_pct=round(100*a["total_reads"]/b["total_reads"], 2),
        Q30_before_pct=round(100*b["q30_rate"], 2), Q30_after_pct=round(100*a["q30_rate"], 2),
        GC_after_pct=round(100*a["gc_content"], 2),
        Duplication_pct=round(100*d.get("duplication", {}).get("rate", float("nan")), 2),
        Raw_len_R1=b.get("read1_mean_length"), Trim_len_R1=a.get("read1_mean_length")))
if not rows: sys.exit(f"no *.fastp.json in {qc}")
cols = list(rows[0])
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows: fh.write("\t".join(str(r[c]) for c in cols) + "\n")
clin = [r for r in rows if r["Group"] == "Clinical"]
print(f"  clinical isolates: {len(clin)}  (all retained libraries: {len(rows)})")
for k, lab, sc, fmt in [("Trimmed_reads","Trimmed reads (M)",1e6,"{:.2f}"),
                        ("Pass_rate_pct","Pass rate (%)",1,"{:.1f}"),
                        ("Q30_after_pct","Q30 after (%)",1,"{:.1f}"),
                        ("GC_after_pct","GC after (%)",1,"{:.1f}"),
                        ("Duplication_pct","Duplication (%)",1,"{:.1f}")]:
    v = [r[k]/sc for r in clin]
    print(f"    {lab:20s} mean {fmt.format(st.mean(v))}, range {fmt.format(min(v))}-{fmt.format(max(v))}")
print(f"  per-sample table: {out}")
PY
fi

# -----------------------------------------------------------------------------
if run_section C; then
hdr "C. FastQC modules and per-base-sequence-content failures"
python3 - "$QC" "$OUT/fastqc_modules.tsv" <<'PY'
import zipfile, glob, os, sys, collections
qc, out = sys.argv[1:3]
mods, pos, nfail, nfiles = collections.defaultdict(collections.Counter), collections.Counter(), 0, 0
for f in sorted(glob.glob(os.path.join(qc, "*_fastqc.zip"))):
    z = zipfile.ZipFile(f); root = z.namelist()[0].split("/")[0]
    txt = z.read(root + "/fastqc_data.txt").decode(); nfiles += 1
    for line in z.read(root + "/summary.txt").decode().splitlines():
        status, module, _ = line.split("\t"); mods[module][status] += 1
    blk = txt.split(">>Per base sequence content")[1].split(">>END_MODULE")[0]
    if not blk.lstrip().startswith("fail"): continue
    nfail += 1
    rows = [r for r in (l.split("\t") for l in blk.splitlines()
            if l and not l.startswith("#") and "\t" in l) if len(r) == 5]
    for r in rows:
        if abs(float(r[1])-float(r[4])) > 20 or abs(float(r[2])-float(r[3])) > 20:
            pos[r[0]] += 1
print(f"  files assessed: {nfiles}")
with open(out, "w") as fh:
    fh.write("Module\tPASS\tWARN\tFAIL\n")
    for m, c in mods.items():
        fh.write(f"{m}\t{c['PASS']}\t{c['WARN']}\t{c['FAIL']}\n")
        print(f"    {m:40s} PASS {c['PASS']:4d}  WARN {c['WARN']:4d}  FAIL {c['FAIL']:4d}")
print(f"  per base sequence content: {nfail} file(s) FAIL; positions breaching the 20% limit:")
for p, n in pos.most_common(): print(f"    base {p}: {n} file(s)")
print(f"  module table: {out}")
PY
fi

# -----------------------------------------------------------------------------
if run_section D; then
hdr "D. Sequencing batches (raw read length)"
python3 - "$QC" "$EXCLUDED" <<'PY'
import json, glob, os, sys, collections
qc, excl = sys.argv[1:3]
g = collections.Counter()
for f in sorted(glob.glob(os.path.join(qc, "*.fastp.json"))):
    if os.path.basename(f).startswith(excl): continue
    b = json.load(open(f))["summary"]["before_filtering"]
    g[(b.get("read1_mean_length"), b.get("read2_mean_length"))] += 1
for k, v in sorted(g.items()):
    print(f"    raw mean length R1 x R2 = {k[0]} x {k[1]}: {v} libraries")
print("  note: mean length can fall below the nominal cycle count if reads were")
print("        trimmed before delivery; check the maximum length in a raw FASTQ.")
PY
fi

# -----------------------------------------------------------------------------
if run_section E; then
hdr "E. Alignment depth, mapping rate and reference coverage"
QCTAB="$WORKDIR/snippy_qc.tsv"
if [[ -f "$QCTAB" ]]; then
python3 - "$QCTAB" <<'PY'
import csv, sys, statistics as st
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
def col(name):
    return [float(r[name]) for r in rows if r.get(name) not in (None, "", "NA")]
print(f"  rows in snippy_qc.tsv: {len(rows)}")
for name, lab in [("Mean_depth","Mean depth (x)"), ("Mapping_rate_pct","Mapping rate (%)"),
                  ("Ref_coverage_pct","Reference coverage (%)")]:
    v = col(name)
    if v: print(f"    {lab:24s} mean {st.mean(v):.2f}  median {st.median(v):.2f}  "
                f"range {min(v):.1f}-{max(v):.1f}")
PY
else
    echo "  $QCTAB not found — run gwas11.sh Step 2b first"
fi
if [[ -f "$SNIPPY/$ATCC/snps.bam" ]]; then
    printf "    %-24s " "ATCC 29548 mean depth"
    samtools depth -a "$SNIPPY/$ATCC/snps.bam" | awk '{t+=$3} END{printf "%.1fx\n", t/NR}'
fi
fi

# -----------------------------------------------------------------------------
if run_section F; then
hdr "F. rrs (B586_RS14090): depth and unfiltered variant calls"
FIRST_BAM=$(ls -d "$SNIPPY"/*/ 2>/dev/null | head -1)
if [[ -z "${FIRST_BAM:-}" || ! -f "$FIRST_BAM/snps.bam" ]]; then
  echo "  no Snippy BAM files under $SNIPPY — skipping"
else
CHR=$(samtools view -H "$FIRST_BAM/snps.bam" | grep -m1 '^@SQ' | sed 's/.*SN:\([^\t]*\).*/\1/')
echo "  region: $CHR:$RRS_START-$RRS_END"
echo "  16S rRNA copies in the reference: $(grep -c '/product="16S ribosomal RNA"' "$GBFF")"
{ echo -e "Sample\tMin_depth\tMean_depth"
  for d in "$SNIPPY"/*/; do
      s=$(basename "$d"); [[ "$s" == "$EXCLUDED" ]] && continue
      samtools depth -a -r "$CHR:$RRS_START-$RRS_END" "$d/snps.bam" 2>/dev/null |
        awk -v s="$s" '{t+=$3; if(NR==1||$3<m)m=$3} END{if(NR)printf "%s\t%d\t%.0f\n", s, m, t/NR}'
  done; } > "$OUT/rrs_depth.tsv"
echo "  lowest minimum depths:"
tail -n +2 "$OUT/rrs_depth.tsv" | sort -k2,2n | head -5 | sed 's/^/    /'
echo "  unfiltered (snps.raw.vcf) calls within rrs:"
n=0
for d in "$SNIPPY"/*/; do
    s=$(basename "$d"); f="$d/snps.raw.vcf"; [[ -f "$f" ]] || continue
    while read -r line; do
        [[ -z "$line" ]] && continue
        echo "    $s  $line"; n=$((n+1))
    done < <(awk -v c="$CHR" -v a="$RRS_START" -v b="$RRS_END" 'BEGIN{OFS="\t"}
             !/^#/ && $1==c && $2>=a && $2<=b {
               match($8,/DP=[0-9]+/); dp=substr($8,RSTART+3,RLENGTH-3)
               match($8,/AO=[0-9,]+/); ao=substr($8,RSTART+3,RLENGTH-3)
               print $2,$4">"$5,"QUAL="int($6),"DP="dp,"AO="ao}' "$f")
done
echo "    (calls listed above include $EXCLUDED if present; study isolates should have none)"
echo "  depth table: $OUT/rrs_depth.tsv"
fi
fi

# -----------------------------------------------------------------------------
if run_section G; then
hdr "G. gyrA intein boundaries"
echo "-- G1. protein alignment of an intein-positive vs intein-negative assembly --"
POS_FAA="$PANGENOME/prokka/CT-WN15CA4114/CT-WN15CA4114.faa"
NEG_FAA="$PANGENOME/prokka/CT-WN15CC3812/CT-WN15CC3812.faa"
if [[ -f "$POS_FAA" && -f "$NEG_FAA" ]]; then
cat > "$OUT/intein_align.py" <<'PY'
import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
def gyrA(f):
    for r in SeqIO.parse(f, "fasta"):
        if "gyrase subunit A" in r.description or " gyrA" in r.description:
            return r.id, str(r.seq).rstrip("*")
    sys.exit(f"no gyrA in {f}")
pid, pos = gyrA(sys.argv[1]); nid, neg = gyrA(sys.argv[2])
al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                     substitution_matrix=substitution_matrices.load("BLOSUM62"))
a = al.align(pos, neg)[0]
blocks = a.aligned[0]
gaps = [(blocks[i][1], blocks[i+1][0]) for i in range(len(blocks)-1)]
b1, b2 = max(gaps, key=lambda g: g[1]-g[0])
gyra_start = int(sys.argv[3])
print(f"    intein+ {pid}: {len(pos)} aa | intein- {nid}: {len(neg)} aa")
print(f"    intein = residues {b1+1}-{b2} ({b2-b1} aa); preceded by {pos[b1-1]}{b1}")
print(f"    reference nt {gyra_start + b1*3}-{gyra_start + b2*3 - 1}")
print(f"    junctions: start {pos[b1]}  end {pos[b2-2:b2]}  +1 {pos[b2]}")
PY
  (python3 "$OUT/intein_align.py" "$POS_FAA" "$NEG_FAA" "$GYRA_START" 2>/dev/null ||
   conda run -n "$ENV_NAME" python "$OUT/intein_align.py" "$POS_FAA" "$NEG_FAA" "$GYRA_START")
else
  echo "    Prokka .faa files not found under $PANGENOME/prokka"
fi
echo "-- G2. interval uncovered in >=90% of intein-negative isolates --"
python3 - "$CORE_ALN" "$WORKDIR/intein_status.tsv" "$GYRA_START" "$GYRA_END" <<'PY'
import csv, os, sys
aln, status, gs, ge = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
if not (os.path.isfile(aln) and os.path.isfile(status)):
    sys.exit("    core.full.aln or intein_status.tsv not found (run gwas11.sh Step 10)")
seqs, name = {}, None
for line in open(aln):
    if line.startswith(">"): name = line[1:].strip(); seqs[name] = []
    else: seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
rows = list(csv.DictReader(open(status), delimiter="\t"))
cands = [c for c in rows[0] if any(str(r[c]).strip().lower() in ("absent", "present") for r in rows)]
col = cands[-1]
neg = [r["Sample"] for r in rows if str(r[col]).strip().lower() == "absent" and r["Sample"] in seqs]
print(f"    column '{col}': {len(neg)} intein-negative isolates")
unc = [p for p in range(gs-1, ge) if sum(seqs[s][p] in "-N" for s in neg) >= 0.9*len(neg)]
if unc: print(f"    uncovered in >=90%: {unc[0]+1}-{unc[-1]+1} ({len(unc)} bp)")
PY
fi

# -----------------------------------------------------------------------------
if run_section H; then
hdr "H. Target locus tags vs the reference annotation"
if [[ ! -f "$GBFF" ]]; then echo "  reference GenBank file not found: $GBFF"; else
python3 - "$GBFF" <<'PY'
import re, sys
tags = {"B586_RS03455": "rpoB (RIF)", "B586_RS00030": "gyrA (CIP)",
        "B586_RS00025": "gyrB (CIP)", "B586_RS14090": "rrs / 16S rRNA (AMI)"}
feats, cur, key = [], None, None
for line in open(sys.argv[1]):
    if line.startswith("ORIGIN"): cur = None; continue
    m = re.match(r"^ {5}(\S+)\s+(\S+)", line)
    if m:
        cur = {"type": m.group(1), "loc": m.group(2)}; feats.append(cur); key = None; continue
    if cur is None: continue
    q = re.match(r'^ {21}/(\w+)=?"?(.*?)"?$', line.rstrip())
    if q: key = q.group(1); cur[key] = q.group(2)
    elif key and line.startswith(" " * 21): cur[key] += " " + line.strip().strip('"')
for t, exp in tags.items():
    hits = [f for f in feats if f.get("locus_tag") == t and f["type"] != "gene"] or \
           [f for f in feats if f.get("locus_tag") == t]
    for f in hits or [None]:
        print(f"    {t:14s} {exp:22s} " +
              (f"{f['type']:6s} {f.get('gene','-'):8s} {f['loc']:28s} {f.get('product','-')}"
               if f else "NOT FOUND"))
PY
fi
fi

# -----------------------------------------------------------------------------
if run_section I; then
hdr "I. Intein status by collection year: exact test"
python3 - "$WORKDIR/intein_status.tsv" <<'PY'
import csv, os, sys, collections
from math import comb
f = sys.argv[1]
if not os.path.isfile(f): sys.exit("    intein_status.tsv not found (run gwas11.sh Step 10)")
rows = [r for r in csv.DictReader(open(f), delimiter="\t")
        if "clinical" in str(r.get("Group", "")).lower() and "replicate" not in str(r.get("Group","")).lower()]
cands = [c for c in rows[0] if any(str(r[c]).strip().lower() in ("absent","present") for r in rows)]
col = cands[-1]
years = sorted({r["Sample"][5:7] for r in rows})
n = [sum(1 for r in rows if r["Sample"][5:7] == y) for y in years]
k = [sum(1 for r in rows if r["Sample"][5:7] == y and
         str(r[col]).strip().lower() == "present") for y in years]
N, K = sum(n), sum(k)
print(f"    intein-positive by year: " +
      "; ".join(f"20{y}: {a}/{b}" for y, a, b in zip(years, k, n)))
def prob(x):
    r = 1
    for xi, ni in zip(x, n): r *= comb(ni, xi)
    return r / comb(N, K)
obs, tot = prob(k), 0.0
import itertools
for a in itertools.product(*[range(ni+1) for ni in n[:-1]]):
    last = K - sum(a)
    if 0 <= last <= n[-1]:
        p = prob(list(a) + [last])
        if p <= obs * (1 + 1e-9): tot += p
print(f"    Fisher-Freeman-Halton exact p = {tot:.3f}")
try:
    from scipy.stats import chi2_contingency
    import numpy as np
    c, p, d, e = chi2_contingency(np.array([k, [a-b for a, b in zip(n, k)]]))
    print(f"    chi-square p = {p:.3f} (expected counts {e[0].round(2).tolist()}; "
          f"exact test preferred when any expected count < 5)")
except Exception:
    pass
PY
fi

# -----------------------------------------------------------------------------
if run_section J; then
hdr "J. Pan-genome clusters vs gyrA intein status"
RTAB="$PANGENOME/panaroo/gene_presence_absence.Rtab"
python3 - "$RTAB" "$WORKDIR/intein_status.tsv" "$CLUSTERS_J" \
         "$OUT/cluster_intein_crosstab.tsv" <<'PY'
import csv, os, sys
rtab, status, wanted, out = sys.argv[1:5]
if not os.path.isfile(rtab):
    sys.exit(f"    {rtab} not found (run the pan-genome steps first)")
if not os.path.isfile(status):
    sys.exit(f"    {status} not found (run gwas11.sh / PART B Step 10 first)")

rows = list(csv.DictReader(open(status), delimiter="\t"))
cands = [c for c in rows[0] if any(str(r[c]).strip().lower() in ("absent", "present") for r in rows)]
col = cands[-1]
clin = [r for r in rows if str(r.get("Group", "")).strip().lower().startswith("clinical")]
pos = {r["Sample"] for r in clin if str(r[col]).strip().lower() == "present"}
neg = {r["Sample"] for r in clin if str(r[col]).strip().lower() == "absent"}
print(f"    clinical isolates: {len(clin)}  intein-positive: {len(pos)}  intein-negative: {len(neg)}")
if not pos or not neg:
    sys.exit("    could not determine intein groups — check intein_status.tsv")

with open(rtab) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    samples = header[1:]
    want = [w for w in wanted.split() if w]
    found = {}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[0] in want:
            found[p[0]] = {s for s, v in zip(samples, p[1:]) if v.strip() not in ("", "0")}

try:
    from scipy.stats import fisher_exact
except Exception:
    fisher_exact = None

with open(out, "w") as fh:
    fh.write("Cluster\tIntein_pos_carrying\tIntein_pos_total\t"
             "Intein_neg_carrying\tIntein_neg_total\tFisher_p\n")
    for g in want:
        if g not in found:
            print(f"    {g:14s} NOT FOUND in {os.path.basename(rtab)}")
            continue
        a, b = len(found[g] & pos), len(found[g] & neg)
        p = ""
        if fisher_exact:
            p = f"{fisher_exact([[a, len(pos)-a], [b, len(neg)-b]])[1]:.3g}"
        fh.write(f"{g}\t{a}\t{len(pos)}\t{b}\t{len(neg)}\t{p}\n")
        print(f"    {g:14s} intein-positive {a:2d}/{len(pos)}   "
              f"intein-negative {b:2d}/{len(neg)}" + (f"   Fisher p = {p}" if p else ""))
print(f"    table: {out}")
PY
fi

echo
echo "============================================================"
echo "Done. Output files in $OUT"
echo "============================================================"
