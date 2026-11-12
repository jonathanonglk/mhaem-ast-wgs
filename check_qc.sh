#!/usr/bin/env bash
# =============================================================================
# check_qc.sh — check fastp, FastQC and MultiQC outputs for the M. haemophilum
# WGS run (mhaem_full_pipeline.sh, Steps 6–7)
#
# Run:   bash ~/check_qc.sh            (any env with python3; stdlib only)
#
# Reads: ~/mhaem_wgs/raw_reads/<sample>_R1/_R2.fastq.gz   (sample list)
#        ~/mhaem_wgs/qc/<sample>.fastp.json / .html
#        ~/mhaem_wgs/qc/<sample>_R{1,2}.trim_fastqc.zip / .html
#        ~/mhaem_wgs/qc/multiqc_report.html + multiqc_data/
#
# Writes: ~/mhaem_wgs/qc/qc_check_per_sample.tsv
#         ~/mhaem_wgs/qc/qc_check_fastqc_modules.tsv
#
# Checks:
#   1. Every sample has fastp json/html and FastQC R1/R2 zip/html
#   2. fastp per-sample metrics + summary (all, and 56 clinical isolates only)
#   3. fastp R1 == R2 read counts after trimming (pairs intact)
#   4. FastQC module PASS/WARN/FAIL counts; FastQC "Total Sequences" vs fastp
#   5. MultiQC report exists, is newer than all inputs, and includes every sample
# =============================================================================
set -euo pipefail
WORKDIR="${WORKDIR:-$HOME/mhaem_wgs}"

python3 - "$WORKDIR" <<'PYEOF'
import sys, os, json, glob, zipfile, statistics as st

W = sys.argv[1]
QC, RAW = os.path.join(W, "qc"), os.path.join(W, "raw_reads")

# Replicates and QC strain — excluded from the "clinical" summary (as in the manuscript)
NON_CLINICAL = {"CT-WN18CD7757", "CT-WN17CC9963", "CT-WN18CA0775", "CT-WN18CA3057", "CT-WNTC9548"}

samples = sorted(os.path.basename(f)[:-len("_R1.fastq.gz")]
                 for f in glob.glob(os.path.join(RAW, "*_R1.fastq.gz")))
if not samples:
    sys.exit(f"No *_R1.fastq.gz found in {RAW}")
print(f"Samples in raw_reads: {len(samples)}  "
      f"(clinical: {len([s for s in samples if s not in NON_CLINICAL])})\n")

problems = []

# ---------------------------------------------------------------- 1. files present
print("== 1. Output files present ==")
expected = {
    "fastp json": "{s}.fastp.json", "fastp html": "{s}.fastp.html",
    "FastQC R1 zip": "{s}_R1.trim_fastqc.zip", "FastQC R2 zip": "{s}_R2.trim_fastqc.zip",
    "FastQC R1 html": "{s}_R1.trim_fastqc.html", "FastQC R2 html": "{s}_R2.trim_fastqc.html",
}
for label, pat in expected.items():
    missing = [s for s in samples if not os.path.isfile(os.path.join(QC, pat.format(s=s)))]
    print(f"  {label:15s}: {len(samples)-len(missing)}/{len(samples)}"
          + (f"  MISSING: {', '.join(missing)}" if missing else "  OK"))
    if missing: problems.append(f"{label} missing for {len(missing)} sample(s)")

# ---------------------------------------------------------------- 2–3. fastp
print("\n== 2. fastp metrics ==")
rows = []
for s in samples:
    f = os.path.join(QC, f"{s}.fastp.json")
    if not os.path.isfile(f):
        continue
    d = json.load(open(f))
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    fr = d.get("filtering_result", {})
    r1_after = d.get("read1_after_filtering", {}).get("total_reads")
    r2_after = d.get("read2_after_filtering", {}).get("total_reads")
    rows.append(dict(
        sample=s, clinical="yes" if s not in NON_CLINICAL else "no",
        raw_reads=b["total_reads"], trimmed_reads=a["total_reads"],
        pass_pct=100 * a["total_reads"] / b["total_reads"] if b["total_reads"] else float("nan"),
        q30_before_pct=100 * b["q30_rate"], q30_after_pct=100 * a["q30_rate"],
        gc_after_pct=100 * a["gc_content"],
        mean_len_r1_after=a.get("read1_mean_length"), mean_len_r2_after=a.get("read2_mean_length"),
        duplication_pct=100 * d.get("duplication", {}).get("rate", float("nan")),
        adapter_trimmed_reads=d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
        low_quality=fr.get("low_quality_reads"), too_many_N=fr.get("too_many_N_reads"),
        too_short=fr.get("too_short_reads"),
        insert_peak=d.get("insert_size", {}).get("peak"),
        r1_after=r1_after, r2_after=r2_after,
    ))

def summ(vals, fmt):
    vals = [v for v in vals if v is not None and v == v]
    if not vals: return "NA"
    return (f"mean {fmt.format(st.mean(vals))}, median {fmt.format(st.median(vals))}, "
            f"range {fmt.format(min(vals))}–{fmt.format(max(vals))}")

for label, sel in [("All samples", rows), ("Clinical isolates", [r for r in rows if r["clinical"] == "yes"])]:
    print(f"  {label} (n = {len(sel)})")
    print(f"    Raw reads (R1+R2)     : {summ([r['raw_reads']/1e6 for r in sel], '{:.2f}')} M")
    print(f"    Trimmed reads (R1+R2) : {summ([r['trimmed_reads']/1e6 for r in sel], '{:.2f}')} M")
    print(f"    Pass rate             : {summ([r['pass_pct'] for r in sel], '{:.1f}')} %")
    print(f"    Q30 after trimming    : {summ([r['q30_after_pct'] for r in sel], '{:.1f}')} %")
    print(f"    GC after trimming     : {summ([r['gc_after_pct'] for r in sel], '{:.1f}')} %")
    print(f"    Duplication           : {summ([r['duplication_pct'] for r in sel], '{:.1f}')} %")
    print(f"    Mean R1 length after  : {summ([r['mean_len_r1_after'] for r in sel], '{:.0f}')} bp")

# Flags (adjust thresholds as needed)
print("\n  Flags (pass rate <90%, Q30 <80%, trimmed reads <1 M, GC outside 62–65%):")
nflag = 0
for r in rows:
    why = []
    if r["pass_pct"] < 90: why.append(f"pass {r['pass_pct']:.1f}%")
    if r["q30_after_pct"] < 80: why.append(f"Q30 {r['q30_after_pct']:.1f}%")
    if r["trimmed_reads"] < 1e6: why.append(f"{r['trimmed_reads']:,} reads")
    if not 62 <= r["gc_after_pct"] <= 65: why.append(f"GC {r['gc_after_pct']:.1f}%")
    if why:
        nflag += 1; print(f"    {r['sample']}: {'; '.join(why)}")
if not nflag: print("    none")

print("\n== 3. Read pairs intact after fastp (R1 == R2) ==")
bad = [r["sample"] for r in rows if r["r1_after"] is not None and r["r1_after"] != r["r2_after"]]
odd = [r["sample"] for r in rows if r["r1_after"] is not None and r["r1_after"] + r["r2_after"] != r["trimmed_reads"]]
print("  R1 == R2 for all samples: " + ("OK" if not bad else "NO — " + ", ".join(bad)))
print("  R1 + R2 == fastp total_reads: " + ("OK" if not odd else "NO — " + ", ".join(odd)))
if bad or odd: problems.append("fastp R1/R2 counts inconsistent")

# ---------------------------------------------------------------- 4. FastQC
print("\n== 4. FastQC (trimmed reads) ==")
modules, fq_total = {}, {}
for z in sorted(glob.glob(os.path.join(QC, "*_fastqc.zip"))):
    name = os.path.basename(z)[:-len("_fastqc.zip")]           # e.g. CT-WN15CA4114_R1.trim
    try:
        with zipfile.ZipFile(z) as zf:
            root = zf.namelist()[0].split("/")[0]
            for line in zf.read(f"{root}/summary.txt").decode().splitlines():
                status, module, _ = line.split("\t")
                modules.setdefault(module, {}).setdefault(status, []).append(name)
            for line in zf.read(f"{root}/fastqc_data.txt").decode().splitlines():
                if line.startswith("Total Sequences"):
                    fq_total[name] = int(line.split("\t")[1]); break
    except Exception as e:
        problems.append(f"cannot read {os.path.basename(z)}: {e}")
print(f"  FastQC reports parsed: {len(fq_total)}")
mod_rows = []
print(f"  {'Module':42s} PASS  WARN  FAIL")
for m, sd in modules.items():
    p, w, f = (len(sd.get(k, [])) for k in ("PASS", "WARN", "FAIL"))
    mod_rows.append((m, p, w, f, ", ".join(sd.get("FAIL", []))))
    print(f"  {m:42s} {p:4d}  {w:4d}  {f:4d}")
print("  (Expected for bacterial WGS: 'Per base sequence content', 'Sequence Duplication Levels'\n"
      "   and 'Per sequence GC content' often WARN/FAIL without indicating a problem; a high-GC\n"
      "   genome such as M. haemophilum (~64%) commonly trips the GC module.)")

# FastQC total vs fastp per read
mism = []
for r in rows:
    for rd, key in (("R1", "r1_after"), ("R2", "r2_after")):
        fq = fq_total.get(f"{r['sample']}_{rd}.trim")
        if fq is not None and r[key] is not None and fq != r[key]:
            mism.append(f"{r['sample']}_{rd} (FastQC {fq:,} vs fastp {r[key]:,})")
print("  FastQC Total Sequences == fastp per-read counts: " + ("OK" if not mism else "NO — " + "; ".join(mism)))
if mism: problems.append("FastQC and fastp read counts differ")

# ---------------------------------------------------------------- 5. MultiQC
print("\n== 5. MultiQC ==")
rep = os.path.join(QC, "multiqc_report.html")
if not os.path.isfile(rep):
    print("  multiqc_report.html MISSING"); problems.append("MultiQC report missing")
else:
    inputs = glob.glob(os.path.join(QC, "*.fastp.json")) + glob.glob(os.path.join(QC, "*_fastqc.zip"))
    newest = max(os.path.getmtime(f) for f in inputs) if inputs else 0
    stale = os.path.getmtime(rep) < newest
    print("  Report present; " + ("STALE — older than some fastp/FastQC outputs, rerun multiqc"
                                  if stale else "newer than all inputs: OK"))
    if stale: problems.append("MultiQC report older than inputs")
    ddir = os.path.join(QC, "multiqc_data")
    for fname, label in (("multiqc_fastp.txt", "fastp"), ("multiqc_fastqc.txt", "FastQC")):
        f = os.path.join(ddir, fname)
        if not os.path.isfile(f):
            print(f"  {fname}: not found"); continue
        names = [l.split("\t")[0] for l in open(f).read().splitlines()[1:] if l.strip()]
        absent = [s for s in samples if not any(n.startswith(s) for n in names)]
        print(f"  {label:6s} entries in MultiQC: {len(names)}"
              + ("  — all samples included" if not absent else f"  — MISSING: {', '.join(absent)}"))
        if absent: problems.append(f"MultiQC {label} missing {len(absent)} sample(s)")

# ---------------------------------------------------------------- outputs
cols = list(rows[0].keys()) if rows else []
with open(os.path.join(QC, "qc_check_per_sample.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(f"{r[c]:.2f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
with open(os.path.join(QC, "qc_check_fastqc_modules.tsv"), "w") as fh:
    fh.write("module\tPASS\tWARN\tFAIL\tFAIL_files\n")
    for m in mod_rows: fh.write("\t".join(map(str, m)) + "\n")

print("\n== Summary ==")
print("  No problems found." if not problems else "  " + "\n  ".join(problems))
print(f"  Per-sample table : {os.path.join(QC, 'qc_check_per_sample.tsv')}")
print(f"  FastQC modules   : {os.path.join(QC, 'qc_check_fastqc_modules.tsv')}")
PYEOF
