#!/usr/bin/env bash
# =============================================================================
# test_kinship_absorption.sh — distinguish "the kinship correction absorbed a
# real lineage signal" from "the association was never strong" for the gene
# clusters that track the gyrA intein lineage.
#
# Run:  bash test_kinship_absorption.sh
#       bash test_kinship_absorption.sh PANDIR=/home/jonat/pangenome/pyseer
#
# For every cluster whose presence/absence closely matches intein status it
# reports:
#   carriers                 - n isolates carrying it (and whether it passed
#                              pyseer's 2..55 frequency filter, i.e. testable)
#   in results file?         - TESTED vs NOT TESTED  (settles the first question)
#   filter-pvalue            - pyseer's association WITHOUT kinship correction
#   lrt-pvalue               - pyseer's association WITH the LMM kinship term
#   attenuation              - lrt / filter, how much the correction cost
#   raw Mann-Whitney p       - independent uncorrected test on MIC rank
#
# Verdict logic
#   raw/filter p small AND lrt p large      -> correction absorbed the signal
#   raw/filter p large  AND lrt p large     -> association was never strong
#   cluster absent from results file        -> never tested (filter excluded it)
# =============================================================================
set -uo pipefail
for arg in "$@"; do declare "$arg"; done

WORKDIR="${WORKDIR:-$HOME/mhaem_wgs}"       # env_1 WORKDIR
PANDIR="${PANDIR:-$HOME/pangenome/pyseer}"  # env_2 OUTDIR/pyseer
RTAB="${RTAB:-$PANDIR/gene_presence_absence_clinical.Rtab}"
INTEIN="${INTEIN:-$WORKDIR/intein_status.tsv}"
DRUGS="${DRUGS:-AMI CIP RIF}"
MIN_MATCH="${MIN_MATCH:-0.90}"   # min concordance with intein status to be listed
MAXN="${MAXN:-12}"               # max clusters to report

python3 - "$PANDIR" "$RTAB" "$INTEIN" "$DRUGS" "$MIN_MATCH" "$MAXN" <<'PYEOF'
import sys, os, csv, math
from collections import defaultdict

pandir, rtab, inteinf, drugs, minmatch, maxn = sys.argv[1:7]
drugs = drugs.split(); minmatch = float(minmatch); maxn = int(maxn)

# ---------------------------------------------------------------- intein status
if not os.path.isfile(inteinf): sys.exit(f"intein file not found: {inteinf}")
intein = {}
with open(inteinf) as fh:
    for r in csv.DictReader(fh, delimiter="\t", quotechar='"'):
        s = r.get("Sample") or r.get("sample")
        v = (r.get("Intein_status") or "").strip().lower()
        if s and v in ("present", "absent"):
            intein[s] = 1 if v == "present" else 0

# ---------------------------------------------------------------- Rtab
if not os.path.isfile(rtab): sys.exit(f"Rtab not found: {rtab}")
with open(rtab) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    samples = hdr[1:]
    pres = {}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 2: continue
        pres[f[0]] = {s: int(v) for s, v in zip(samples, f[1:]) if v in ("0", "1")}

shared = [s for s in samples if s in intein]
pos = [s for s in shared if intein[s] == 1]
neg = [s for s in shared if intein[s] == 0]
print(f"Rtab: {len(pres):,} clusters x {len(samples)} samples; "
      f"{len(shared)} with intein status ({len(pos)} positive, {len(neg)} negative)\n")

# clusters tracking the intein lineage
cands = []
for g, d in pres.items():
    inp = sum(d.get(s, 0) for s in pos)
    inn = sum(d.get(s, 0) for s in neg)
    match = (inp + (len(neg) - inn)) / len(shared)          # concordance with intein+
    if match >= minmatch and inp >= 1:
        cands.append((match, g, inp, inn, inp + inn))
cands.sort(key=lambda x: (-x[0], x[1]))
if not cands:
    sys.exit(f"No clusters matched intein status at >= {minmatch:.0%} concordance.")

# ---------------------------------------------------------------- pyseer results
res = {}
for drug in drugs:
    path = None
    for cand in (f"{drug}_results.txt", f"pyseer_results_{drug}.txt", f"{drug}.assoc.txt"):
        if os.path.isfile(os.path.join(pandir, cand)):
            path = os.path.join(pandir, cand); break
    if not path:
        print(f"!! no results file for {drug} in {pandir}"); continue
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        idx = {n: k for k, n in enumerate(h)}
        tbl = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f: continue
            def num(col):
                try: return float(f[idx[col]])
                except Exception: return None
            tbl[f[0]] = dict(filter_p=num("filter-pvalue"), lrt_p=num("lrt-pvalue"),
                             beta=num("beta"), af=num("af"))
    res[drug] = tbl
    print(f"{drug}: {len(tbl):,} clusters in {os.path.basename(path)}")

# ---------------------------------------------------------------- phenotypes
pheno = {}
for drug in drugs:
    for cand in (f"phenotypes_{drug}.tsv", f"{drug}.txt", f"pheno_{drug}.tsv"):
        p = os.path.join(pandir, cand)
        if not os.path.isfile(p):
            p = os.path.join(os.path.dirname(pandir), cand)
        if os.path.isfile(p):
            d = {}
            with open(p) as fh:
                for line in fh:
                    f = line.split()
                    if len(f) < 2 or f[0].lower() in ("sample", "samples"): continue
                    try: d[f[0]] = float(f[1])
                    except ValueError: pass
            if d: pheno[drug] = d; break

def mannwhitney(a, b):
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0: return float('nan')
    allv = sorted(a + b); ranks = {}; i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j+1] == allv[i]: j += 1
        ranks[allv[i]] = (i + j) / 2 + 1; i = j + 1
    R1 = sum(ranks[x] for x in a)
    U1 = R1 - n1*(n1+1)/2; U = min(U1, n1*n2 - U1)
    mu = n1*n2/2
    cnt = defaultdict(int)
    for x in allv: cnt[x] += 1
    N = n1 + n2
    tie = sum(c**3 - c for c in cnt.values())
    sd = ((n1*n2/(N*(N-1))) * ((N**3 - N - tie)/12)) ** 0.5
    if sd == 0: return float('nan')
    z = (abs(U - mu) - 0.5)/sd
    return math.erfc(z / 2**0.5)

# ---------------------------------------------------------------- report
print(f"\nClusters tracking intein status at >= {minmatch:.0%} concordance "
      f"(showing up to {maxn}):\n")
for match, g, inp, inn, n in cands[:maxn]:
    testable = 2 <= n <= len(samples) - 1
    print(f"== {g}   carriers {n}/{len(shared)}  "
          f"(intein+ {inp}/{len(pos)}, intein- {inn}/{len(neg)})  "
          f"concordance {match:.0%}")
    print(f"   pyseer frequency filter (2..{len(samples)-1} carriers): "
          + ("PASSES -> should have been tested" if testable else "EXCLUDED -> not testable"))
    for drug in drugs:
        row = res.get(drug, {}).get(g)
        if row is None:
            print(f"   {drug}: NOT TESTED (absent from results file)")
            continue
        fp, lp = row['filter_p'], row['lrt_p']
        att = (lp / fp) if (fp and lp and fp > 0) else None
        ph = pheno.get(drug, {})
        carr = [ph[s] for s in shared if ph.get(s) is not None and pres[g].get(s) == 1]
        non  = [ph[s] for s in shared if ph.get(s) is not None and pres[g].get(s) == 0]
        raw = mannwhitney(carr, non) if carr and non else float('nan')
        verdict = ""
        if fp is not None and lp is not None:
            if fp < 0.01 and lp > 0.05:   verdict = "  <- strong before correction, null after: ABSORBED"
            elif fp > 0.05 and lp > 0.05: verdict = "  <- weak before and after: NEVER STRONG"
            elif lp <= 0.05:              verdict = "  <- survives correction"
        print(f"   {drug}: TESTED   filter-p = {fp:<10.3g} lrt-p = {lp:<10.3g} "
              f"attenuation = {('%.0fx' % att) if att else 'NA'}   raw MW p = {raw:.3g}{verdict}")
    print()

print("filter-pvalue is pyseer's association test WITHOUT the kinship term;")
print("lrt-pvalue is the same cluster under the LMM WITH kinship. A large gap")
print("means the correction removed the signal; similar values mean it was never there.")
PYEOF
