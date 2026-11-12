#!/usr/bin/env bash
# =============================================================================
# compare_abc_snp_mic.sh — compare CIP MICs between isolates carrying the
# reference and alternative alleles at the B586_RS18800 ABC transporter
# permease variants (the block-1 CIP GWAS signal).
#
# Run:   bash compare_abc_snp_mic.sh
#        bash compare_abc_snp_mic.sh VCF=~/mhaem_wgs/core/core.vcf MIC=~/mhaem_wgs/CIP.txt
#        bash compare_abc_snp_mic.sh POS=4021393,4021694,4021771 DRUG=CIP
#
# Inputs
#   VCF  : the multi-sample VCF used for the GWAS (snippy-core output)
#   MIC  : whitespace/tab-separated, one line per isolate: <sample> <MIC>
#          MIC may be written as 4, >8, <=0.25, 0.5/9.5 etc.
#   POS  : comma-separated reference positions to genotype
#
# Output
#   - genotype x MIC cross-tabulation for each position and for the haplotype
#   - MIC50 / range / full distribution per allele group
#   - Mann-Whitney U on ordinal MIC rank (two-sided)
#   - Fisher exact on susceptible vs non-susceptible (CIP MIC >= 2 = NS)
#   - per-isolate table written to abc_snp_mic.tsv
# =============================================================================
set -uo pipefail
for arg in "$@"; do declare "$arg"; done
if [ -z "${MIC:-}" ]; then
  echo "MIC=<file> is required (ordinal MIC rank file for the drug)." >&2
  echo "e.g. bash $0 MIC=~/data/CIP.txt VCF=~/mhaem_wgs/core/core.vcf" >&2
  exit 2
fi

WORKDIR="${WORKDIR:-$HOME/mhaem_wgs}"
DRUG="${DRUG:-CIP}"
VCF="${VCF:-$WORKDIR/core/core.vcf}"
MIC="${MIC:-}"                 # REQUIRED: the ordinal MIC rank file for $DRUG
POS="${POS:-4021393,4021694,4021771}"
SNIPPY="${SNIPPY:-$WORKDIR/snippy}"
# Phenotype files hold ORDINAL RANKS, not ug/mL. Non-susceptible cut on the rank
# scale: CIP rank >= 5 (MIC >= 2 ug/mL), AMI rank >= 6 (MIC >= 32 ug/mL).
case "$DRUG" in
  CIP) DEF_CUT=5 ;;
  AMI) DEF_CUT=6 ;;
  *)   DEF_CUT=5 ;;
esac
NS_CUT="${NS_CUT:-$DEF_CUT}"   # non-susceptible if rank >= this
EXCLUDE_SAMPLES="${EXCLUDE_SAMPLES:-CT-WN17CD6059 CT-WN18CD7757 CT-WN17CC9963 CT-WN18CA0775 CT-WN18CA3057 CT-WNTC9548}"

python3 - "$VCF" "$MIC" "$POS" "$SNIPPY" "$NS_CUT" "$DRUG" "$EXCLUDE_SAMPLES" <<'PYEOF'
import sys, os, re, gzip, glob, itertools
from collections import defaultdict

vcf, micf, posarg, snippy, ns_cut, drug, excl = sys.argv[1:8]
POS = [int(p) for p in posarg.split(",") if p.strip()]
NS  = float(ns_cut)
EXCL = set(excl.split())

def op(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

# ---------------------------------------------------------------- MIC values
def parse_mic(s):
    """'>8' -> 8 with flag '>' ; '<=0.25' -> 0.25 with '<=' ; '0.5/9.5' -> 0.5"""
    s = str(s).strip().replace("≤", "<=").replace("≥", ">=")
    m = re.match(r'^\s*(<=|>=|<|>)?\s*([0-9.]+)', s.split("/")[0])
    if not m: return None, None
    return float(m.group(2)), (m.group(1) or "")

if not os.path.isfile(micf):
    sys.exit(f"MIC file not found: {micf}  (pass MIC=<file>)")
mic, flag = {}, {}
with open(micf) as fh:
    for line in fh:
        f = line.split()
        if len(f) < 2: continue
        if f[0].lower() in ("sample", "samples", "strain", "isolate"): continue
        v, fl = parse_mic(f[1])
        if v is None: continue
        mic[f[0]] = v; flag[f[0]] = fl
print(f"MIC values read: {len(mic)} isolates from {os.path.basename(micf)}")

# ---------------------------------------------------------------- genotypes
geno = defaultdict(dict)   # geno[sample][pos] = 'REF' | 'ALT' | 'NA'
samples = []

if os.path.isfile(vcf):
    with op(vcf) as fh:
        for line in fh:
            if line.startswith("##"): continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]; continue
            if not samples: continue
            try: p = int(f[1])
            except ValueError: continue
            if p not in POS: continue
            for s, cell in zip(samples, f[9:]):
                gt = cell.split(":")[0].replace("|", "/")
                a = gt.split("/")[0]
                geno[s][p] = "NA" if a in (".", "") else ("REF" if a == "0" else "ALT")
    print(f"Genotypes from VCF: {len(samples)} samples, "
          f"{sum(1 for s in samples if any(p in geno[s] for p in POS))} with calls at target positions")

# fall back to / fill gaps from snippy per-sample snps.tab
if not samples or any(p not in geno[s] for s in samples for p in POS):
    tabs = glob.glob(os.path.join(snippy, "*", "snps.tab"))
    if tabs:
        print(f"Filling from {len(tabs)} snippy snps.tab files "
              "(absence of a call at a position = REF)")
        for t in tabs:
            s = os.path.basename(os.path.dirname(t))
            if s not in samples: samples.append(s)
            found = set()
            with open(t) as fh:
                for line in fh:
                    f = line.split("\t")
                    if len(f) < 2: continue
                    try: p = int(f[1])
                    except ValueError: continue
                    if p in POS: found.add(p)
            for p in POS:
                geno[s].setdefault(p, "ALT" if p in found else "REF")

samples = [s for s in samples if s not in EXCL]
samples = [s for s in samples if s in mic]
if not samples:
    sys.exit("No samples with both genotype and MIC — check sample naming in the MIC file.")
print(f"Isolates analysed: {len(samples)} (excluded: {', '.join(sorted(EXCL))})\n")

# ---------------------------------------------------------------- stats
def mic50(vals):
    v = sorted(vals)
    return v[(len(v) + 1)//2 - 1] if v else float('nan')

def fmt(v, fl=""): return f"{fl}{v:g}"

def mannwhitney(a, b):
    """two-sided Mann-Whitney U with tie correction, normal approximation"""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0: return float('nan'), float('nan')
    allv = sorted(a + b)
    ranks, i = {}, 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j+1] == allv[i]: j += 1
        r = (i + j) / 2 + 1
        ranks[allv[i]] = r
        i = j + 1
    R1 = sum(ranks[x] for x in a)
    U1 = R1 - n1 * (n1 + 1) / 2
    U = min(U1, n1 * n2 - U1)
    mu = n1 * n2 / 2
    cnt = defaultdict(int)
    for x in allv: cnt[x] += 1
    N = n1 + n2
    tie = sum(c**3 - c for c in cnt.values())
    sd = ((n1 * n2 / (N * (N - 1))) * ((N**3 - N - tie) / 12)) ** 0.5
    if sd == 0: return U, float('nan')
    import math
    z = (abs(U - mu) - 0.5) / sd
    return U, math.erfc(z / 2**0.5)

def fisher(a, b, c, d):
    """two-sided Fisher exact on [[a,b],[c,d]]"""
    from math import comb
    n = a + b + c + d
    r1, c1 = a + b, a + c
    def pr(x): return comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
    p0 = pr(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p0 * (1 + 1e-9))

out = open("abc_snp_mic.tsv", "w")
out.write("sample\t" + "\t".join(f"pos_{p}" for p in POS) + "\thaplotype\tMIC\tcategory\n")

for p in POS + ["HAPLOTYPE"]:
    if p == "HAPLOTYPE":
        grp = defaultdict(list)
        for s in samples:
            g = {geno[s].get(q, "NA") for q in POS}
            lab = "ALT(all)" if g == {"ALT"} else ("REF(all)" if g == {"REF"} else "mixed/NA")
            grp[lab].append(s)
        title = f"All {len(POS)} positions combined"
    else:
        grp = defaultdict(list)
        for s in samples: grp[geno[s].get(p, "NA")].append(s)
        title = f"Position {p:,}"

    print(f"== {title} ==")
    rows = []
    for lab in sorted(grp, key=lambda k: (k.startswith("mixed"), k)):
        ss = grp[lab]
        vs = [mic[s] for s in ss]
        ns = sum(1 for v in vs if v >= NS)
        dist = defaultdict(int)
        for s in ss: dist[fmt(mic[s], flag[s])] += 1
        print(f"  {lab:10s} n={len(ss):>3}  MIC50={mic50(vs):>6g}  "
              f"range={min(vs):g}-{max(vs):g}  "
              f"non-susceptible={ns}/{len(ss)} ({100*ns/len(ss):.0f}%)")
        print(f"             MICs: " + ", ".join(f"{k}x{v}" for k, v in
              sorted(dist.items(), key=lambda kv: parse_mic(kv[0])[0])))
        rows.append((lab, vs, ns))

    ref = next((r for r in rows if r[0].startswith("REF")), None)
    alt = next((r for r in rows if r[0].startswith("ALT")), None)
    if ref and alt:
        U, pmw = mannwhitney(ref[1], alt[1])
        pf = fisher(ref[2], len(ref[1]) - ref[2], alt[2], len(alt[1]) - alt[2])
        print(f"  Mann-Whitney U = {U:g}, two-sided p = {pmw:.4g}   "
              f"(REF vs ALT MIC distribution)")
        print(f"  Fisher exact   p = {pf:.4g}   "
              f"(non-susceptible: REF {ref[2]}/{len(ref[1])} vs ALT {alt[2]}/{len(alt[1])})")
    print()

for s in sorted(samples):
    g = [geno[s].get(p, "NA") for p in POS]
    hap = "ALT(all)" if set(g) == {"ALT"} else ("REF(all)" if set(g) == {"REF"} else "mixed/NA")
    cat = "non-susceptible" if mic[s] >= NS else "susceptible"
    out.write(f"{s}\t" + "\t".join(g) + f"\t{hap}\t{fmt(mic[s], flag[s])}\t{cat}\n")
out.close()
print(f"Per-isolate table written to {os.path.abspath('abc_snp_mic.tsv')}")
print(f"\nNote: {drug} non-susceptible defined as MIC >= {NS:g} ug/mL (CLSI M24S Table 5).")
print("MIC values at panel boundaries are treated at their nominal value; "
      "Mann-Whitney is on ranks, so boundary handling does not affect the test.")
PYEOF
