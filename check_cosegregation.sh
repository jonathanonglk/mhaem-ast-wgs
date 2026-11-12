#!/usr/bin/env bash
# =============================================================================
# check_cosegregation.sh — test the Discussion claim that the top CIP SNP
# signals "share their pattern with many co-segregating variants elsewhere in
# the genome".
#
# Run:   bash check_cosegregation.sh
#        bash check_cosegregation.sh RESULTS=/home/jonat/mhaem_wgs/pyseer_results_CIP.txt
#        bash check_cosegregation.sh DRUG=AMI TOPN=3
#
# Method: variants that share a presence/absence pattern necessarily receive
# identical LRT p-values and identical beta from the LMM. So grouping the
# pyseer results by (lrt-pvalue, beta) recovers the pattern groups without
# needing the genotype matrix. For the group(s) containing the lowest p-value,
# the script reports how many variants share that pattern and how they are
# distributed along the chromosome.
#
# Reports:
#   1. Size of the pattern group containing the top hit, and its coordinate
#      span, spread across genome deciles, and nearest/furthest members
#   2. The same for the largest pattern groups overall
#   3. Genome-wide summary: how many variants, how many distinct patterns,
#      median/max group size  (the "24,294 variants = 826 patterns" fact)
# =============================================================================
set -uo pipefail
for arg in "$@"; do declare "$arg"; done

WORKDIR="${WORKDIR:-$HOME/mhaem_wgs}"
DRUG="${DRUG:-CIP}"
RESULTS="${RESULTS:-$WORKDIR/pyseer_results_${DRUG}.txt}"
TOPN="${TOPN:-5}"          # how many top-ranked pattern groups to detail
GENE_LO="${GENE_LO:-0}"    # optional: highlight members inside a locus
GENE_HI="${GENE_HI:-0}"

python3 - "$RESULTS" "$DRUG" "$TOPN" "$GENE_LO" "$GENE_HI" <<'PYEOF'
import sys, os, re, gzip
from collections import defaultdict

path, drug, topn, glo, ghi = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
if not os.path.isfile(path):
    sys.exit(f"results file not found: {path}")

op = gzip.open if path.endswith(".gz") else open
with op(path, "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    try:
        ip = hdr.index("lrt-pvalue")
    except ValueError:
        sys.exit("no lrt-pvalue column; is this a pyseer results table?")
    ib = hdr.index("beta") if "beta" in hdr else None
    rows = []
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) <= ip: continue
        try: p = float(f[ip])
        except ValueError: continue
        b = None
        if ib is not None and len(f) > ib:
            try: b = float(f[ib])
            except ValueError: b = None
        rows.append((f[0], p, b))

if not rows:
    sys.exit("no variants parsed")

def coord(name):
    """Pull the genomic position out of a pyseer variant id.
    Handles NZ_CP011883_9984_A_G, CHROM_POS_REF_ALT, and POS_REF_ALT."""
    m = re.findall(r'_(\d+)_[ACGTNacgtn*\-]+_[ACGTNacgtn*\-]+$', name)
    if m: return int(m[0])
    m = re.match(r'^(\d+)_', name)
    if m: return int(m.group(1))
    m = re.findall(r'(\d{3,})', name)
    return int(m[-1]) if m else None

# group by (p, beta) == same presence/absence pattern
groups = defaultdict(list)
for name, p, b in rows:
    groups[(round(p, 12), None if b is None else round(b, 10))].append(name)

sizes = sorted(len(v) for v in groups.values())
gmax = max(coord(n) or 0 for n, _, _ in rows)

print(f"\n=== {drug}: {os.path.basename(path)} ===")
print(f"Variants tested          : {len(rows):,}")
print(f"Distinct (p, beta) groups: {len(groups):,}   <- pattern-equivalence classes")
print(f"Group size median / max  : {sizes[len(sizes)//2]} / {sizes[-1]}")
print(f"Largest observed coord   : {gmax:,} bp")

def describe(key, label):
    members = groups[key]
    cs = sorted(c for c in (coord(n) for n in members) if c is not None)
    p, b = key
    print(f"\n-- {label}")
    print(f"   LRT p = {p:.3e}" + (f", beta = {b:+.4g}" if b is not None else ""))
    print(f"   variants sharing this pattern: {len(members):,}")
    if not cs:
        print("   (could not parse coordinates from variant ids)"); return
    print(f"   coordinate range: {cs[0]:,} – {cs[-1]:,} bp  (span {cs[-1]-cs[0]:,} bp)")
    if gmax:
        dec = [0]*10
        for c in cs:
            dec[min(9, int(10 * c / (gmax + 1)))] += 1
        print("   spread across genome deciles: " + " ".join(f"{d:>4d}" for d in dec))
        occupied = sum(1 for d in dec if d)
        print(f"   deciles occupied: {occupied}/10 "
              + ("-> scattered genome-wide" if occupied >= 5
                 else "-> clustered in part of the genome"))
    gaps = [b_ - a_ for a_, b_ in zip(cs, cs[1:])]
    if gaps:
        gaps.sort()
        print(f"   median gap between consecutive members: {gaps[len(gaps)//2]:,} bp")
    if glo and ghi:
        inside = [c for c in cs if glo <= c <= ghi]
        print(f"   members inside {glo:,}-{ghi:,}: {len(inside)} of {len(cs)}"
              f"  -> {len(cs)-len(inside)} lie OUTSIDE that locus")

# the top hit(s)
best = min(groups, key=lambda k: k[0])
describe(best, f"Pattern group containing the LOWEST p-value")

print("\n-- Next lowest-p pattern groups")
for key in sorted(groups, key=lambda k: k[0])[1:topn]:
    members = groups[key]
    cs = sorted(c for c in (coord(n) for n in members) if c is not None)
    dec = [0]*10
    for c in cs: dec[min(9, int(10 * c / (gmax + 1)))] += 1
    print(f"   p = {key[0]:.3e}  n = {len(members):>5,}  "
          f"span = {(cs[-1]-cs[0]) if cs else 0:>9,} bp  deciles = {sum(1 for d in dec if d)}/10")

print("\n-- Largest pattern groups overall")
for key in sorted(groups, key=lambda k: -len(groups[k]))[:topn]:
    members = groups[key]
    cs = sorted(c for c in (coord(n) for n in members) if c is not None)
    dec = [0]*10
    for c in cs: dec[min(9, int(10 * c / (gmax + 1)))] += 1
    print(f"   p = {key[0]:.3e}  n = {len(members):>5,}  "
          f"span = {(cs[-1]-cs[0]) if cs else 0:>9,} bp  deciles = {sum(1 for d in dec if d)}/10")

print("\nInterpretation: if the top-hit group contains many variants spread over "
      "most genome deciles, the signal marks a clonal lineage rather than a locus, "
      "because in a clonal organism co-inherited variants are scattered genome-wide "
      "rather than physically linked.\n")
PYEOF
