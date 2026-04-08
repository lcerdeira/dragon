#!/usr/bin/env python3
"""Generate Figure 1: Dragon Architecture Diagram (updated for Phases 1-4)."""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

fig, ax = plt.subplots(1, 1, figsize=(14, 9))
ax.set_xlim(0, 14)
ax.set_ylim(0, 9.2)
ax.axis("off")

# ============================================================================
# Colour scheme
# ============================================================================
COL_INPUT = "#E8F4FD"
COL_INDEX = "#D5E8D4"
COL_QUERY = "#FFF2CC"
COL_CORE = "#DAE8FC"
COL_OUTPUT = "#F8CECC"
COL_HEADER = "#4472C4"
COL_ARROW = "#333333"
COL_BORDER_IDX = "#82B366"
COL_BORDER_QRY = "#D6B656"
COL_BORDER_IN = "#6C8EBF"
COL_BORDER_OUT = "#B85450"
COL_DATA_FLOW = "#6C8EBF"
COL_NEW = "#E6B8AF"          # Phase 1-4 additions
COL_BORDER_NEW = "#CC4125"
COL_ML = "#D5A6BD"           # ML components
COL_BORDER_ML = "#A64D79"
COL_SURV = "#B4A7D6"         # Surveillance
COL_BORDER_SURV = "#674EA7"

def draw_box(x, y, w, h, text, facecolor, edgecolor="#666666", fontsize=8,
             bold=False, text_color="black", lw=1.5, alpha=1.0):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.08",
        facecolor=facecolor, edgecolor=edgecolor,
        linewidth=lw, alpha=alpha,
    )
    ax.add_patch(box)
    weight = "bold" if bold else "normal"
    ax.text(
        x + w / 2, y + h / 2, text,
        ha="center", va="center", fontsize=fontsize,
        fontweight=weight, color=text_color,
        wrap=True,
    )

def draw_arrow(x1, y1, x2, y2, color=COL_ARROW, style="-|>", lw=1.5,
               connectionstyle="arc3,rad=0"):
    arrow = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle=style,
        color=color,
        linewidth=lw,
        connectionstyle=connectionstyle,
        mutation_scale=15,
    )
    ax.add_patch(arrow)

# ============================================================================
# SECTION HEADERS
# ============================================================================
draw_box(0.2, 8.4, 5.8, 0.5, "INDEX CONSTRUCTION (offline, one-time)",
         COL_HEADER, COL_HEADER, fontsize=10, bold=True, text_color="white")

draw_box(7.0, 8.4, 6.8, 0.5, "QUERY PIPELINE (online, <4 GB RAM)",
         COL_HEADER, COL_HEADER, fontsize=10, bold=True, text_color="white")

# ============================================================================
# LEFT PANEL: INDEX CONSTRUCTION
# ============================================================================

# Input: Genome FASTA files
draw_box(1.3, 7.3, 2.7, 0.65, "Genome FASTA files\n(N genomes)",
         COL_INPUT, COL_BORDER_IN, fontsize=8.5, bold=True)

# GGCAT
draw_box(0.3, 6.1, 2.3, 0.7, "GGCAT\n(ccdBG builder)",
         COL_INDEX, COL_BORDER_IDX, fontsize=8.5, bold=True)

# ccdBG output
draw_box(3.0, 6.1, 2.8, 0.7, "Coloured compacted\nde Bruijn graph",
         COL_CORE, COL_BORDER_IN, fontsize=8)

# Arrow: input -> GGCAT
draw_arrow(2.65, 7.3, 1.45, 6.85, lw=2)
# Arrow: GGCAT -> ccdBG
draw_arrow(2.6, 6.45, 3.0, 6.45, lw=2)

# Four index components
# FM-index
draw_box(0.2, 4.5, 1.7, 0.95, "Run-length\nFM-index\n(r-index)",
         COL_INDEX, COL_BORDER_IDX, fontsize=7.5, bold=True)

# Colour index
draw_box(2.1, 4.5, 1.7, 0.95, "Colour Index\n(Roaring\nBitmaps)",
         COL_INDEX, COL_BORDER_IDX, fontsize=7.5, bold=True)

# Genome paths
draw_box(4.0, 4.5, 1.7, 0.95, "Genome Path\nIndex\n(delta-coded)",
         COL_INDEX, COL_BORDER_IDX, fontsize=7.5, bold=True)

# NEW: Specificity index (Phase 1C)
draw_box(0.2, 3.2, 2.4, 0.85, "Specificity Index\n(private unitigs\nper genome)",
         COL_NEW, COL_BORDER_NEW, fontsize=7.5, bold=True)

# Arrows from ccdBG down to four components
draw_arrow(1.8, 6.1, 1.05, 5.5, lw=1.5)
draw_arrow(3.5, 6.1, 2.95, 5.5, lw=1.5)
draw_arrow(5.0, 6.1, 4.85, 5.5, lw=1.5)

# Arrow: Colour Index -> Specificity
draw_arrow(2.95, 4.5, 1.4, 4.1, lw=1.2, color=COL_BORDER_NEW)

# Elias-Fano sub-component
draw_box(3.8, 3.2, 1.9, 0.55, "Elias-Fano\nposition index",
         "#E1D5E7", "#9673A6", fontsize=7)
draw_arrow(1.05, 4.5, 4.75, 3.8, lw=1.0, color="#9673A6",
           connectionstyle="arc3,rad=-0.1")

# On-disk storage
draw_box(0.5, 1.8, 5.3, 0.8, "On-disk index (~100 GB for 2M genomes)\nmemory-mapped  |  laptop/workstation profiles",
         "#F5F5F5", "#999999", fontsize=8)

# NEW: Overlay box (Phase 3B)
draw_box(0.5, 0.7, 5.3, 0.7, "Incremental overlays (dragon update)\noverlay_001/ ... overlay_N/ → dragon compact",
         COL_NEW, COL_BORDER_NEW, fontsize=7.5)

# Arrows: components -> disk
draw_arrow(1.05, 3.2, 2.0, 2.65, lw=1.0, color="#999999")
draw_arrow(2.95, 4.5, 3.0, 2.65, lw=1.0, color="#999999")
draw_arrow(4.85, 4.5, 4.5, 2.65, lw=1.0, color="#999999")

# ============================================================================
# RIGHT PANEL: QUERY PIPELINE
# ============================================================================

# Input query
draw_box(8.8, 7.3, 2.8, 0.65, "Query FASTA\n(gene / plasmid / read)",
         COL_INPUT, COL_BORDER_IN, fontsize=8.5, bold=True)

# Stage 1
draw_box(7.5, 6.2, 4.5, 0.65, "Stage 1: FM-index backward search\n(variable-length seed matching)",
         COL_QUERY, COL_BORDER_QRY, fontsize=8, bold=True)

# Stage 2 (updated: Phase 1A)
draw_box(7.5, 5.15, 4.5, 0.65, "Stage 2: Info-weighted candidate filtering\n(IC = log2(N/|color|) per unitig vote)",
         COL_NEW, COL_BORDER_NEW, fontsize=7.5, bold=True)

# Stage 3 (updated: Phase 1B + 1C)
draw_box(7.5, 3.9, 4.5, 0.85, "Stage 3: ML-scored colinear chaining\n(10-feature logistic regression scorer)\n+ specificity re-ranking (private unitigs)",
         COL_ML, COL_BORDER_ML, fontsize=7.5, bold=True)

# Stage 4
draw_box(7.5, 2.9, 4.5, 0.65, "Stage 4: Banded alignment + filtering\n(identity, coverage, score-ratio, dedup)",
         COL_QUERY, COL_BORDER_QRY, fontsize=8, bold=True)

# Output block (Phase 2 + 4)
draw_box(7.3, 0.7, 2.8, 1.7, "PAF\nBLAST6\nGFA + Bandage\ncolours",
         COL_OUTPUT, COL_BORDER_OUT, fontsize=8, bold=True)

# Surveillance output (Phase 2)
draw_box(10.5, 0.7, 3.0, 1.7, "dragon summarize\n\nSpecies prevalence\nANI distributions\nGenome-level detail\n(TSV / JSON)",
         COL_SURV, COL_BORDER_SURV, fontsize=7.5, bold=True)

# Vertical arrows between stages (right panel)
draw_arrow(10.2, 7.3, 9.75, 6.9, lw=2)
draw_arrow(9.75, 6.2, 9.75, 5.85, lw=2)
draw_arrow(9.75, 5.15, 9.75, 4.8, lw=2)
draw_arrow(9.75, 3.9, 9.75, 3.6, lw=2)
draw_arrow(9.75, 2.9, 9.0, 2.45, lw=2)
draw_arrow(9.75, 2.9, 12.0, 2.45, lw=2)

# ML training loop arrow (Phase 1B)
ax.annotate(
    "", xy=(12.3, 4.75), xytext=(12.3, 3.55),
    arrowprops=dict(arrowstyle="-|>", color=COL_BORDER_ML, lw=1.5,
                    connectionstyle="arc3,rad=-0.5"),
)
ax.text(12.7, 4.15, "train\nweights", fontsize=6.5, color=COL_BORDER_ML,
        ha="center", fontstyle="italic",
        bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                  edgecolor="none", alpha=0.9))
draw_box(12.1, 4.8, 1.6, 0.55, "seed_scorer.py\n(logistic reg.)",
         COL_ML, COL_BORDER_ML, fontsize=6.5)
draw_arrow(12.05, 4.35, 12.05, 4.8, lw=1.2, color=COL_BORDER_ML)

# ============================================================================
# CROSS-PANEL DATA FLOW ARROWS
# ============================================================================
# FM-index -> Stage 1
draw_arrow(1.9, 4.98, 7.45, 6.52, lw=1.2, color=COL_DATA_FLOW,
           style="-|>", connectionstyle="arc3,rad=-0.12")

# Colour Index -> Stage 2
draw_arrow(3.8, 4.98, 7.45, 5.47, lw=1.2, color=COL_DATA_FLOW,
           style="-|>", connectionstyle="arc3,rad=-0.05")

# Genome Path Index -> Stage 3
draw_arrow(5.7, 4.98, 7.45, 4.32, lw=1.2, color=COL_DATA_FLOW,
           style="-|>", connectionstyle="arc3,rad=0.1")

# Specificity Index -> Stage 3
draw_arrow(2.6, 3.62, 7.45, 4.1, lw=1.2, color=COL_BORDER_NEW,
           style="-|>", connectionstyle="arc3,rad=0.2")

# Labels for data flow arrows
ax.text(4.0, 6.3, "seeds", fontsize=7, color=COL_DATA_FLOW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                  edgecolor="none", alpha=0.85))
ax.text(5.2, 5.35, "colours", fontsize=7, color=COL_DATA_FLOW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                  edgecolor="none", alpha=0.85))
ax.text(6.3, 4.35, "paths", fontsize=7, color=COL_DATA_FLOW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                  edgecolor="none", alpha=0.85))
ax.text(4.3, 3.65, "specificity", fontsize=6.5, color=COL_BORDER_NEW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                  edgecolor="none", alpha=0.85))

# ============================================================================
# Legend for new components
# ============================================================================
leg_y = 0.15
ax.text(0.5, leg_y, "Legend:", fontsize=7, fontweight="bold")
for x, color, border, label in [
    (1.6, COL_INDEX, COL_BORDER_IDX, "Original"),
    (3.3, COL_NEW, COL_BORDER_NEW, "Phase 1/3 (accuracy + updates)"),
    (6.5, COL_ML, COL_BORDER_ML, "Phase 1B (ML scoring)"),
    (9.5, COL_SURV, COL_BORDER_SURV, "Phase 2 (surveillance)"),
]:
    box = FancyBboxPatch((x, leg_y - 0.08), 0.3, 0.25,
                         boxstyle="round,pad=0.04",
                         facecolor=color, edgecolor=border, lw=1.2)
    ax.add_patch(box)
    ax.text(x + 0.4, leg_y + 0.04, label, fontsize=6.5, va="center")

# ============================================================================
# Save
# ============================================================================
output_dir = "figures"
fig.savefig(f"{output_dir}/fig1_architecture.pdf", bbox_inches="tight")
fig.savefig(f"{output_dir}/fig1_architecture.png", bbox_inches="tight")
plt.close(fig)
print(f"Figure 1 saved to {output_dir}/fig1_architecture.{{pdf,png}}")
