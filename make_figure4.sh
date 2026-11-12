#!/usr/bin/env bash
# =============================================================================
# make_figure4.sh — Figure 4: MIC distributions of CIP, AMI and RIF by gyrA
# intein status (56 clinical M. haemophilum isolates from unique patients)
# =============================================================================
# Input : ~/mhaem_wgs/intein_status.tsv  (written by gwas11.sh Step 10)
# Output: ~/mhaem_wgs/figure4_intein_mic.png  (3900 x 1260 px, 300 dpi)
#
# Bars  : % of isolates within each group (intein-positive n=14, blue;
#         intein-negative n=42, orange) at each MIC; numbers above bars are counts
# Line  : CLSI M24S non-susceptibility breakpoint (CIP >=2, AMI >=32, RIF >=2 ug/mL)
#
# Requires python3 with pandas + matplotlib (e.g. the mhaem_wgs environment).
#
# Run with defaults:
#   bash make_figure4.sh
# or with your own paths:
#   bash make_figure4.sh /path/to/intein_status.tsv /path/to/figure4.png
# =============================================================================
set -euo pipefail

INPUT="${1:-$HOME/mhaem_wgs/intein_status.tsv}"
OUTPUT="${2:-$HOME/mhaem_wgs/figure4_intein_mic.png}"

if [[ ! -s "$INPUT" ]]; then
    echo "ERROR: $INPUT not found or empty — run gwas11.sh (Step 10) first."
    exit 1
fi

echo "============================================================"
echo "Figure 4: MIC distributions by gyrA intein status"
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

echo "Figure 4 complete: $OUTPUT"
