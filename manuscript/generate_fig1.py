#!/usr/bin/env python3
"""Generate Figure 1: Dragon Architecture Diagram."""

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

fig, ax = plt.subplots(1, 1, figsize=(13, 7.5))
ax.set_xlim(0, 13)
ax.set_ylim(0, 7.8)
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
draw_box(0.2, 7.1, 5.8, 0.5, "INDEX CONSTRUCTION (offline, one-time)",
         COL_HEADER, COL_HEADER, fontsize=10, bold=True, text_color="white")

draw_box(7.0, 7.1, 5.8, 0.5, "QUERY PIPELINE (online, <4 GB RAM)",
         COL_HEADER, COL_HEADER, fontsize=10, bold=True, text_color="white")

# ============================================================================
# LEFT PANEL: INDEX CONSTRUCTION
# ============================================================================

# Input: Genome FASTA files
draw_box(1.3, 6.1, 2.7, 0.65, "Genome FASTA files\n(N genomes)",
         COL_INPUT, COL_BORDER_IN, fontsize=8.5, bold=True)

# GGCAT
draw_box(0.3, 4.9, 2.3, 0.7, "GGCAT\n(ccdBG builder)",
         COL_INDEX, COL_BORDER_IDX, fontsize=8.5, bold=True)

# ccdBG output
draw_box(3.0, 4.9, 2.8, 0.7, "Coloured compacted\nde Bruijn graph",
         COL_CORE, COL_BORDER_IN, fontsize=8)

# Arrow: input -> GGCAT
draw_arrow(2.65, 6.1, 1.45, 5.65, lw=2)
# Arrow: GGCAT -> ccdBG
draw_arrow(2.6, 5.25, 3.0, 5.25, lw=2)

# Three index components — spread wider to avoid crowding
# FM-index
draw_box(0.2, 3.3, 2.0, 0.95, "Run-length\nFM-index\n(r-index)",
         COL_INDEX, COL_BORDER_IDX, fontsize=8, bold=True)

# Colour index
draw_box(2.5, 3.3, 2.0, 0.95, "Colour Index\n(Roaring\nBitmaps)",
         COL_INDEX, COL_BORDER_IDX, fontsize=8, bold=True)

# Genome paths
draw_box(4.8, 3.3, 1.8, 0.95, "Genome Path\nIndex\n(delta-coded)",
         COL_INDEX, COL_BORDER_IDX, fontsize=8, bold=True)

# Arrows from ccdBG down to three components
draw_arrow(1.8, 4.9, 1.2, 4.3, lw=1.5)
draw_arrow(4.0, 4.9, 3.5, 4.3, lw=1.5)
draw_arrow(5.4, 4.9, 5.7, 4.3, lw=1.5)

# Elias-Fano sub-component
draw_box(0.2, 2.2, 2.0, 0.6, "Elias-Fano\nposition index",
         "#E1D5E7", "#9673A6", fontsize=7.5)

draw_arrow(1.2, 3.3, 1.2, 2.85, lw=1.2)

# On-disk storage
draw_box(1.0, 0.8, 4.8, 0.8, "On-disk index (~100 GB for 2M genomes)\nmemory-mapped at query time",
         "#F5F5F5", "#999999", fontsize=8)

# Arrows: components -> disk
draw_arrow(1.2, 2.2, 2.2, 1.65, lw=1.0, color="#999999")
draw_arrow(3.5, 3.3, 3.4, 1.65, lw=1.0, color="#999999")
draw_arrow(5.7, 3.3, 4.8, 1.65, lw=1.0, color="#999999")

# ============================================================================
# RIGHT PANEL: QUERY PIPELINE
# ============================================================================

# Input query
draw_box(8.5, 6.1, 2.5, 0.65, "Query FASTA\n(gene / plasmid / read)",
         COL_INPUT, COL_BORDER_IN, fontsize=8.5, bold=True)

# Stage 1
draw_box(7.5, 5.0, 4.0, 0.7, "Stage 1: FM-index backward search\n(variable-length seed matching)",
         COL_QUERY, COL_BORDER_QRY, fontsize=8, bold=True)

# Stage 2
draw_box(7.5, 3.9, 4.0, 0.7, "Stage 2: Candidate genome filtering\n(colour-based voting, Roaring Bitmaps)",
         COL_QUERY, COL_BORDER_QRY, fontsize=8, bold=True)

# Stage 3
draw_box(7.5, 2.8, 4.0, 0.7, "Stage 3: Graph-aware colinear chaining\n(Fenwick tree DP, O(h log h))",
         COL_QUERY, COL_BORDER_QRY, fontsize=8, bold=True)

# Stage 4
draw_box(7.5, 1.7, 4.0, 0.7, "Stage 4: Banded wavefront alignment\n(WFA along genome paths)",
         COL_QUERY, COL_BORDER_QRY, fontsize=8, bold=True)

# Output
draw_box(8.5, 0.6, 2.5, 0.65, "Alignment\nresults",
         COL_OUTPUT, COL_BORDER_OUT, fontsize=8.5, bold=True)

# Vertical arrows between stages (right panel)
draw_arrow(9.75, 6.1, 9.5, 5.75, lw=2)
draw_arrow(9.5, 5.0, 9.5, 4.65, lw=2)
draw_arrow(9.5, 3.9, 9.5, 3.55, lw=2)
draw_arrow(9.5, 2.8, 9.5, 2.45, lw=2)
draw_arrow(9.5, 1.7, 9.75, 1.3, lw=2)

# ============================================================================
# CROSS-PANEL DATA FLOW ARROWS (curved to avoid overlaps)
# ============================================================================
# FM-index -> Stage 1 (seeds) — curve upward
draw_arrow(2.2, 3.78, 7.45, 5.35, lw=1.2, color=COL_DATA_FLOW,
           style="-|>", connectionstyle="arc3,rad=-0.15")

# Colour Index -> Stage 2 (colours) — straight
draw_arrow(4.5, 3.78, 7.45, 4.25, lw=1.2, color=COL_DATA_FLOW,
           style="-|>", connectionstyle="arc3,rad=0.0")

# Genome Path Index -> Stage 3 (paths) — curve downward
draw_arrow(6.6, 3.78, 7.45, 3.15, lw=1.2, color=COL_DATA_FLOW,
           style="-|>", connectionstyle="arc3,rad=0.15")

# Labels for data flow arrows (positioned to not overlap)
ax.text(4.3, 5.05, "seeds", fontsize=7, color=COL_DATA_FLOW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
ax.text(5.6, 4.2, "colours", fontsize=7, color=COL_DATA_FLOW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
ax.text(6.7, 3.25, "paths", fontsize=7, color=COL_DATA_FLOW,
        fontstyle="italic", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))

# ============================================================================
# Save
# ============================================================================
output_dir = "/Users/lshlt19/GitHub/Dragon/manuscript/figures"
fig.savefig(f"{output_dir}/fig1_architecture.pdf", bbox_inches="tight")
fig.savefig(f"{output_dir}/fig1_architecture.png", bbox_inches="tight")
plt.close(fig)
print(f"Figure 1 saved to {output_dir}/fig1_architecture.{{pdf,png}}")
