#!/usr/bin/env bash
# =============================================================================
# check_gwas_table5.sh — re-derive every number in Table 5 from the pyseer
# output files, independently of the scripts that produced them.
#
# Run:   bash check_gwas_table5.sh
#        bash check_gwas_table5.sh WORKDIR=~/mhaem_wgs PANOUT=~/pangenome
#
# Recomputes, per drug and per approach:
#   Features tested      = number of result rows
#   Unique patterns      = from the pyseer --output-patterns file (count_patterns.py)
#   Bonferroni threshold = 0.05 / unique patterns
#   Lowest p             = min(lrt-pvalue)
#   Nominal hits         = n(lrt-pvalue < 0.05)
#   Significant hits     = n(lrt-pvalue < Bonferroni threshold)
#   lambda               = median(chi2.isf(lrt-pvalue, df=1)) / 0.455
#   h^2                  = parsed from the pyseer run log (NOT recomputable from
#                          the results table — FaST-LMM reports it at fit time)
#
# NOTE on h^2: pyseer prints "h^2 = ..." to stderr while setting up the LMM.
# If you did not tee the stderr to a log file, h^2 cannot be recovered from the
# results table; rerun pyseer with `2>&1 | tee <drug>.log` to capture it.
# =============================================================================
set -uo pipefail

for arg in "$@"; do declare "$arg"; done

WORKDIR="${WORKDIR:-$HOME/mhaem_wgs}"
PANOUT="${PANOUT:-$HOME/pangenome}"           # env2_pangenome.sh OUTDIR
GWAS="${GWAS:-$WORKDIR}"                      # env_1 writes pyseer_results_<DRUG>.txt here
PANGWAS="${PANGWAS:-$PANOUT/pyseer}"          # env_2 writes <DRUG>_results.txt here
DRUGS="${DRUGS:-AMI CIP RIF}"
# filename substrings that are never a results table (set EXCLUDE=... to override;
# "legacy" is excluded by default so the standard run is checked)
EXCLUDE="${EXCLUDE:-pattern|hits|legacy|threshold|phenotype|kinship|samples|stderr|\.log$|\.err$|\.Rtab$}"

python3 - "$GWAS" "$PANGWAS" "$DRUGS" "$EXCLUDE" <<'PYEOF'
import sys, os, glob, math, gzip, re

GWAS, PANGWAS, DRUGS = sys.argv[1], sys.argv[2], sys.argv[3].split()
EXCLUDE = sys.argv[4] if len(sys.argv) > 4 else r"pattern|hits|legacy|threshold|phenotype|kinship|samples|stderr|\.log$|\.err$|\.Rtab$"

try:
    from scipy.stats import chi2
    def isf(p): return chi2.isf(p, 1)
except ImportError:
    # exact inverse of the 1-df chi-square survival function without scipy:
    # P(X>x) = erfc(sqrt(x/2))  =>  x = 2 * erfcinv(p)^2
    def _erfcinv(y):
        # Newton refinement on an initial rational approximation
        if y <= 0 or y >= 2: return float('inf') if y <= 0 else float('-inf')
        x = 0.0
        lo, hi = -10.0, 10.0
        for _ in range(200):                      # bisection: robust, fast enough
            x = (lo + hi) / 2
            if math.erfc(x) > y: lo = x
            else: hi = x
        return x
    def isf(p): return 2.0 * _erfcinv(p) ** 2

def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def find(d, drug, kinds):
    """Locate a pyseer output / patterns / log file for this drug."""
    hits = []
    for pat in kinds:
        hits += glob.glob(os.path.join(d, pat.format(drug=drug)))
        hits += glob.glob(os.path.join(d, pat.format(drug=drug.lower())))
    return sorted(set(hits))

def read_pvals(path):
    """Return (n_rows, [lrt-pvalues]) from a pyseer results table."""
    pv, n = [], 0
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            col = header.index("lrt-pvalue")
        except ValueError:
            cands = [i for i, h in enumerate(header) if "pvalue" in h and "filter" not in h]
            if not cands:
                return None, None
            col = cands[0]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= col: continue
            n += 1
            try: v = float(f[col])
            except ValueError: continue
            if v > 0: pv.append(v)
    return n, pv

def count_patterns(path):
    """Unique presence/absence patterns = unique lines in the patterns file."""
    if not path or not os.path.isfile(path): return None
    seen = set()
    with opener(path) as fh:
        for line in fh:
            line = line.strip()
            if line: seen.add(line)
    return len(seen)

def grab_h2(paths):
    """Parse 'h^2 = 0.38' (or 'heritability') out of a pyseer log."""
    rx = re.compile(r"h\^?2\s*[=:]\s*([0-9.]+)|heritability\s*[=:]\s*([0-9.]+)", re.I)
    for p in paths:
        try: txt = open(p, errors="ignore").read()
        except OSError: continue
        m = rx.search(txt)
        if m: return float(m.group(1) or m.group(2))
    return None

def lam(pv):
    if not pv: return float('nan')
    s = sorted(isf(p) for p in pv)
    n = len(s)
    med = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    return med / 0.455

rows = []
for label, d in (("SNP-based (Snippy)", GWAS), ("Pan-genome (Panaroo)", PANGWAS)):
    if not os.path.isdir(d):
        print(f"!! directory not found: {d}")
        continue
    for drug in DRUGS:
        res = find(d, drug, ["*{drug}*.assoc*", "*{drug}*.txt", "*{drug}*.tsv", "{drug}*"])
        res = [f for f in res if not re.search(EXCLUDE, os.path.basename(f), re.I)]
        # prefer a file that actually looks like the full pyseer results table
        pref = [f for f in res if "result" in os.path.basename(f).lower()]
        if pref: res = pref
        if not res:
            print(f"!! no results file for {label} / {drug} in {d}")
            continue
        path = res[0]
        if len(res) > 1:
            print(f"   note: {len(res)} candidate files for {label}/{drug}; using {os.path.basename(path)}")
        n, pv = read_pvals(path)
        if n is None:
            print(f"!! no lrt-pvalue column in {path}")
            continue
        pat_f = find(d, drug, ["*{drug}*pattern*", "*pattern*{drug}*"])
        npat  = count_patterns(pat_f[0] if pat_f else None)
        bonf  = 0.05 / npat if npat else None
        h2    = grab_h2([f for f in find(d, drug,
                 ["*{drug}*.log", "*{drug}*.err", "*{drug}*stderr*",
                  "*log*{drug}*", "*{drug}*.out"])
                 if "legacy" not in os.path.basename(f).lower()])
        rows.append(dict(
            approach=label, drug=drug, file=os.path.basename(path),
            n=n, npat=npat, bonf=bonf,
            lowest=min(pv) if pv else float('nan'),
            nominal=sum(1 for p in pv if p < 0.05),
            sig=sum(1 for p in pv if bonf and p < bonf) if bonf else None,
            lam=lam(pv), h2=h2,
        ))

if not rows:
    sys.exit("No pyseer outputs found — pass GWAS=<dir> PANGWAS=<dir> explicitly.")

hdr = f"{'Approach':22s} {'Drug':4s} {'Features':>8s} {'Patterns':>8s} {'Bonferroni':>11s} " \
      f"{'Lowest p':>10s} {'Nominal':>7s} {'Sig':>4s} {'h^2':>5s} {'lambda':>7s}"
print("\n" + hdr); print("-" * len(hdr))
for r in rows:
    s_bonf = "%.2e" % r['bonf'] if r['bonf'] else "NA"
    s_h2   = "%.2f" % r['h2']   if r['h2'] is not None else "NA"
    s_sig  = str(r['sig'])      if r['sig'] is not None else "NA"
    s_pat  = str(r['npat'])     if r['npat'] else "NA"
    print(f"{r['approach']:22s} {r['drug']:4s} {r['n']:8d} {s_pat:>8s} "
          f"{s_bonf:>11s} {r['lowest']:10.2e} {r['nominal']:7d} "
          f"{s_sig:>4s} {s_h2:>5s} {r['lam']:7.3f}")

print("\nCompare against Table 5:")
exp = {
 ("SNP-based (Snippy)","AMI"):   (24071, 821, 6.09e-5, 1.10e-3, 166, 0, 0.38, 1.82),
 ("SNP-based (Snippy)","CIP"):   (24294, 826, 6.05e-5, 8.24e-3, 254, 0, 0.44, 3.14),
 ("SNP-based (Snippy)","RIF"):   (24294, 826, 6.05e-5, 1.39e-3,  74, 0, 0.50, 1.28),
 ("Pan-genome (Panaroo)","AMI"): ( 1539, 328, 1.52e-4, 4.57e-3,   7, 0, 0.35, 2.91),
 ("Pan-genome (Panaroo)","CIP"): ( 1562, 334, 1.50e-4, 7.29e-3,  21, 0, 0.61, 2.28),
 ("Pan-genome (Panaroo)","RIF"): ( 1562, 334, 1.50e-4, 1.14e-2,  22, 0, 0.65, 1.83),
}
ok = True
for r in rows:
    e = exp.get((r['approach'], r['drug']))
    if not e: continue
    en, epat, ebonf, elow, enom, esig, eh2, elam = e
    bad = []
    if r['n']    != en:   bad.append(f"features {r['n']} vs {en}")
    if r['npat'] and r['npat'] != epat: bad.append(f"patterns {r['npat']} vs {epat}")
    if abs(r['lam'] - elam) > 0.005:    bad.append(f"lambda {r['lam']:.3f} vs {elam}")
    if r['nominal'] != enom:            bad.append(f"nominal {r['nominal']} vs {enom}")
    if abs(math.log10(r['lowest']) - math.log10(elow)) > 0.01:
        bad.append(f"lowest p {r['lowest']:.2e} vs {elow:.2e}")
    if r['h2'] is not None and abs(r['h2'] - eh2) > 0.005:
        bad.append(f"h^2 {r['h2']:.2f} vs {eh2}")
    if bad:
        ok = False
        print(f"  MISMATCH {r['approach']}/{r['drug']}: " + "; ".join(bad))
    else:
        note = "" if r['h2'] is not None else "  (h^2 not in log — not checked)"
        print(f"  OK       {r['approach']}/{r['drug']}{note}")
print("\nAll checked values match Table 5." if ok else "\nSee mismatches above.")
PYEOF
