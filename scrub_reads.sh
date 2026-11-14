#!/usr/bin/env bash
# =============================================================================
# scrub_reads.sh — run NCBI sra-human-scrubber over every raw FASTQ, write the
# scrubbed copies to a parallel directory, and tabulate how many reads were
# removed from each file.
#
# Setup (once):
#   conda create -y -n scrub -c conda-forge -c bioconda sra-human-scrubber
#   conda activate scrub
#   init_db.sh                     # downloads the human k-mer database (~3 GB)
#
# Run:
#   bash scrub_reads.sh
#   bash scrub_reads.sh IN=~/mhaem_wgs/raw_reads OUT=~/mhaem_wgs/raw_reads_scrubbed
#   bash scrub_reads.sh LIMIT=2    # do the first 2 files only, to check it works
#
# Originals are never modified. Already-scrubbed outputs are skipped, so the
# script can be interrupted and re-run.
#
# Output:
#   $OUT/<sample>_R<n>.fastq.gz    scrubbed, gzip-compressed
#   $OUT/scrub_report.tsv          file, reads_before, reads_after,
#                                  reads_removed, percent_removed
# =============================================================================
set -uo pipefail
for arg in "$@"; do declare "$arg"; done

IN="${IN:-$HOME/mhaem_wgs/raw_reads}"
OUT="${OUT:-$HOME/mhaem_wgs/raw_reads_scrubbed}"
EXCLUDE="${EXCLUDE:-CT-WN17CD6059}"   # sample prefix(es) to skip, space separated
TMP="${TMP:-${TMPDIR:-/tmp}}"
LIMIT="${LIMIT:-0}"                   # 0 = all files
REPORT="$OUT/scrub_report.tsv"

command -v scrub.sh >/dev/null || {
    echo "ERROR: scrub.sh not on PATH. Activate the environment that has"
    echo "       sra-human-scrubber installed, e.g.  conda activate scrub"
    exit 1
}
[ -d "$IN" ] || { echo "ERROR: input directory not found: $IN"; exit 1; }
mkdir -p "$OUT"

# ---------------------------------------------------------------- file list
mapfile -t FILES < <(find "$IN" -maxdepth 1 -name '*.fastq.gz' -printf '%f\n' | sort)
if [ ${#FILES[@]} -eq 0 ]; then echo "ERROR: no *.fastq.gz in $IN"; exit 1; fi

KEEP=()
for f in "${FILES[@]}"; do
    skip=0
    for e in $EXCLUDE; do [[ "$f" == "$e"* ]] && skip=1; done
    [ $skip -eq 0 ] && KEEP+=("$f")
done
echo "Found ${#FILES[@]} FASTQ files; ${#KEEP[@]} to process "
echo "(excluded prefix(es): $EXCLUDE)"
[ "$LIMIT" -gt 0 ] && { KEEP=("${KEEP[@]:0:$LIMIT}"); echo "LIMIT=$LIMIT — processing first ${#KEEP[@]} only"; }

[ -s "$REPORT" ] || printf 'file\treads_before\treads_after\treads_removed\tpercent_removed\n' > "$REPORT"

# ---------------------------------------------------------------- main loop
n=0
for f in "${KEEP[@]}"; do
    n=$((n+1))
    dest="$OUT/$f"
    if [ -s "$dest" ] && grep -q "^${f}	" "$REPORT"; then
        echo "[$n/${#KEEP[@]}] skip (done): $f"
        continue
    fi

    echo "[$n/${#KEEP[@]}] scrubbing: $f"
    tin="$TMP/scrub_in_$$.fastq"
    tout="$TMP/scrub_out_$$.fastq"
    trap 'rm -f "$tin" "$tout"' EXIT

    if ! zcat "$IN/$f" > "$tin"; then
        echo "  ERROR: could not decompress $f — skipped"; rm -f "$tin"; continue
    fi
    before=$(( $(wc -l < "$tin") / 4 ))

    # scrub.sh accepts -i/-o in current versions; fall back to stdin/stdout
    if ! scrub.sh -i "$tin" -o "$tout" -r >/dev/null 2>"$TMP/scrub_err_$$"; then
        if ! scrub.sh -r < "$tin" > "$tout" 2>"$TMP/scrub_err_$$"; then
            echo "  ERROR: scrub.sh failed on $f:"; sed 's/^/    /' "$TMP/scrub_err_$$"
            rm -f "$tin" "$tout" "$TMP/scrub_err_$$"; continue
        fi
    fi
    rm -f "$TMP/scrub_err_$$"

    if [ ! -s "$tout" ]; then
        echo "  ERROR: scrubber produced an empty file for $f — skipped"
        rm -f "$tin" "$tout"; continue
    fi
    after=$(( $(wc -l < "$tout") / 4 ))
    removed=$(( before - after ))
    pct=$(awk -v r="$removed" -v b="$before" 'BEGIN{ printf (b>0 ? "%.6f" : "NA"), (b>0 ? 100*r/b : 0) }')

    gzip -c "$tout" > "$dest.part" && mv "$dest.part" "$dest"
    printf '%s\t%d\t%d\t%d\t%s\n' "$f" "$before" "$after" "$removed" "$pct" >> "$REPORT"
    echo "  $before -> $after reads  (removed $removed, ${pct}%)"
    rm -f "$tin" "$tout"
done

# ---------------------------------------------------------------- summary
echo
echo "Report: $REPORT"
awk -F'\t' 'NR>1 {b+=$2; a+=$3; r+=$4; n++}
END {
  if (n==0) { print "no files processed"; exit }
  printf "Files processed      : %d\n", n
  printf "Total reads before   : %d\n", b
  printf "Total reads after    : %d\n", a
  printf "Total reads removed  : %d (%.6f%% of all reads)\n", r, (b>0?100*r/b:0)
}' "$REPORT"
awk -F'\t' 'NR>1 && $4>0 {print "  " $1 ": " $4 " reads (" $5 "%)"}' "$REPORT" \
  | head -20 | { grep . && echo "  (files with any removal, first 20)" || echo "  no reads removed from any file"; }
