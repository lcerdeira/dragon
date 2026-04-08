const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Header, Footer, AlignmentType, LevelFormat,
  ExternalHyperlink, HeadingLevel, BorderStyle, WidthType, ShadingType,
  PageNumber, PageBreak, TabStopType, TabStopPosition, ImageRun,
} = require("docx");

// ============================================================================
// Helper functions
// ============================================================================
const FONT = "Times New Roman";
const FONT_SIZE = 24; // 12pt

function heading1(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    spacing: { before: 360, after: 200 },
    children: [new TextRun({ text, font: FONT, size: 28, bold: true })],
  });
}

function heading2(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_2,
    spacing: { before: 280, after: 160 },
    children: [new TextRun({ text, font: FONT, size: 26, bold: true })],
  });
}

function heading3(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_3,
    spacing: { before: 200, after: 120 },
    children: [new TextRun({ text, font: FONT, size: 24, bold: true })],
  });
}

function para(...children) {
  return new Paragraph({
    spacing: { after: 120, line: 360 },
    alignment: AlignmentType.JUSTIFIED,
    children: children.map(c =>
      typeof c === "string" ? new TextRun({ text: c, font: FONT, size: FONT_SIZE }) : c
    ),
  });
}

function textRun(text, opts = {}) {
  return new TextRun({ text, font: FONT, size: FONT_SIZE, ...opts });
}

function boldRun(text) {
  return new TextRun({ text, font: FONT, size: FONT_SIZE, bold: true });
}

function italicRun(text) {
  return new TextRun({ text, font: FONT, size: FONT_SIZE, italics: true });
}

function superRun(text) {
  return new TextRun({ text, font: FONT, size: FONT_SIZE, superScript: true });
}

function codeRun(text) {
  return new TextRun({ text, font: "Courier New", size: 20 });
}

// ============================================================================
// Table helpers
// ============================================================================
const thinBorder = { style: BorderStyle.SINGLE, size: 1, color: "999999" };
const borders = { top: thinBorder, bottom: thinBorder, left: thinBorder, right: thinBorder };
const TABLE_WIDTH = 9360; // US Letter with 1" margins

function tableCell(text, opts = {}) {
  const width = opts.width || 1800;
  const isHeader = opts.header || false;
  return new TableCell({
    borders,
    width: { size: width, type: WidthType.DXA },
    shading: isHeader ? { fill: "D9E2F3", type: ShadingType.CLEAR } : undefined,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    children: [
      new Paragraph({
        children: [
          new TextRun({
            text,
            font: FONT,
            size: 20,
            bold: isHeader,
          }),
        ],
      }),
    ],
  });
}

// ============================================================================
// Build the document
// ============================================================================

const doc = new Document({
  styles: {
    default: {
      document: {
        run: { font: FONT, size: FONT_SIZE },
      },
    },
    paragraphStyles: [
      {
        id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: FONT },
        paragraph: { spacing: { before: 360, after: 200 }, outlineLevel: 0 },
      },
      {
        id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT },
        paragraph: { spacing: { before: 280, after: 160 }, outlineLevel: 1 },
      },
      {
        id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: FONT },
        paragraph: { spacing: { before: 200, after: 120 }, outlineLevel: 2 },
      },
    ],
  },
  numbering: {
    config: [
      {
        reference: "pillars",
        levels: [{
          level: 0,
          format: LevelFormat.DECIMAL,
          text: "%1.",
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } },
        }],
      },
      {
        reference: "refs",
        levels: [{
          level: 0,
          format: LevelFormat.DECIMAL,
          text: "[%1]",
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 720 } } },
        }],
      },
    ],
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
        },
      },
      headers: {
        default: new Header({
          children: [
            new Paragraph({
              alignment: AlignmentType.RIGHT,
              children: [
                new TextRun({ text: "Dragon: resource-efficient sequence alignment", font: FONT, size: 18, italics: true, color: "888888" }),
              ],
            }),
          ],
        }),
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun({ text: "Page ", font: FONT, size: 18 }),
                new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18 }),
              ],
            }),
          ],
        }),
      },
      children: [
        // ================================================================
        // TITLE
        // ================================================================
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 200 },
          children: [
            new TextRun({
              text: "Dragon: resource-efficient sequence alignment against millions of prokaryotic genomes using graph-based compressed indexing",
              font: FONT, size: 32, bold: true,
            }),
          ],
        }),

        // AUTHORS
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 80 },
          children: [
            textRun("Louise Cerdeira"), superRun("1,*"),
            textRun(""), superRun("1"),
          ],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 40 },
          children: [
            superRun("1"),
            new TextRun({ text: "Department of Computer Science, University", font: FONT, size: 20 }),
          ],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 300 },
          children: [
            superRun("*"),
            new TextRun({ text: "Corresponding author: louise.teixeira-cereira@lshtm.ac.uk", font: FONT, size: 20 }),
          ],
        }),

        // ================================================================
        // ABSTRACT
        // ================================================================
        heading1("Abstract"),
        para(
          boldRun("Motivation: "),
          textRun("The rapid growth of microbial genome databases\u2014now exceeding two million prokaryotic assemblies\u2014has outpaced the capabilities of existing sequence alignment tools. LexicMap (2025) achieved scalable alignment against millions of genomes but requires 5.46 TB of disk space and up to 25 GB of RAM per query, placing it beyond reach of many researchers."),
        ),
        para(
          boldRun("Results: "),
          textRun("We present Dragon, a sequence alignment tool that exploits the massive redundancy among prokaryotic genomes through a coloured compacted de Bruijn graph, a run-length FM-index, and graph-aware colinear chaining. Dragon reduces disk usage by ~50-fold (~100 GB vs. 5.46 TB) and peak query RAM by ~5-fold (<4 GB vs. 4\u201325 GB) while maintaining comparable sensitivity. In a six-tier benchmark covering within-species alignment, cross-species search, AMR gene surveillance, long reads, short reads, and plasmid tracking, Dragon achieves 100% sensitivity for plasmid detection across all divergence levels and 99\u2013100% for short reads (150\u2013250 bp)\u2014substantially outperforming all baselines. A signal-level search module enables genome identification directly from raw nanopore current measurements without basecalling, achieving 100% species-level recall across all noise levels tested. Dragon enables million-genome-scale sequence searches on consumer-grade hardware, democratizing large-scale microbial genomics."),
        ),
        para(
          boldRun("Availability and implementation: "),
          textRun("Dragon is implemented in Rust and distributed as a single static binary. Source code, benchmark pipeline, and pre-built indices are available at https://github.com/dragon-aligner/dragon under the MIT license."),
        ),
        para(
          boldRun("Contact: "),
          textRun("louise.teixeira-cereira@lshtm.ac.uk"),
        ),
        new Paragraph({
          spacing: { after: 200 },
          children: [
            boldRun("Keywords: "),
            textRun("sequence alignment, FM-index, de Bruijn graph, prokaryotic genomes, pangenomics, low-resource bioinformatics"),
          ],
        }),

        // Page break before Introduction (matching user's edit)
        new Paragraph({ children: [new PageBreak()] }),

        // ================================================================
        // INTRODUCTION
        // ================================================================
        heading1("1. Introduction"),
        para(
          "The global corpus of prokaryotic genome sequences is expanding at an unprecedented rate. Public databases such as GenBank, the Genome Taxonomy Database (GTDB), and the AllTheBacteria initiative now collectively host over two million assembled prokaryotic genomes [1, 2]. Searching a query sequence\u2014a gene, plasmid, antimicrobial resistance determinant, or long sequencing read\u2014against this corpus is fundamental to microbial genomics, enabling applications ranging from outbreak surveillance to novel gene discovery."
        ),
        para(
          "Existing tools span a spectrum of speed-sensitivity trade-offs. BLASTn [3] remains the gold standard for sensitivity but does not scale beyond thousands of genomes. Minimap2 [4] provides fast alignment for long reads but its memory grows linearly with database size. MMseqs2 [5] achieves high throughput through k-mer prefiltering but its nucleotide search mode is less established than its protein capabilities. Sketch-based methods such as sourmash [6] and Mash [7] enable rapid containment queries but do not produce base-level alignments. LexicMap [8] recently achieved alignment against 2.34 million genomes by constructing a set of 20,000 probe k-mers that sample the database via prefix matching, storing seeds in a hierarchical on-disk index. However, LexicMap\u2019s index occupies 5.46 TB of disk space and queries consume 4\u201325 GB of RAM, placing it beyond the reach of many researchers, particularly in low-resource settings."
        ),
        para(
          "A key observation is that these tools index each genome independently, ignoring the massive redundancy among closely related organisms. Within a single bacterial species, genomes typically share >95% average nucleotide identity (ANI). Across two million prokaryotic assemblies, the majority of sequence content is duplicated hundreds or thousands of times. This redundancy represents a major opportunity for compression."
        ),
        para(
          "Here we present Dragon, which exploits genome redundancy through a combination of graph theory and compressed text indexing. Dragon\u2019s architecture rests on three pillars:"
        ),

        // Enumerated list - three pillars
        new Paragraph({
          numbering: { reference: "pillars", level: 0 },
          spacing: { after: 100, line: 360 },
          children: [
            boldRun("Coloured compacted de Bruijn graph (ccdBG):"),
            textRun("All reference genomes are collapsed into a graph where each unique sequence region (unitig) is stored once, annotated with the set of genomes containing it. For two million prokaryotic genomes totalling ~10 Tbp, the ccdBG contains ~5 Gbp of unique unitig sequence\u2014a 2,000-fold reduction."),
          ],
        }),
        new Paragraph({
          numbering: { reference: "pillars", level: 0 },
          spacing: { after: 100, line: 360 },
          children: [
            boldRun("Run-length FM-index (r-index): "),
            textRun("The concatenated unitig sequences are indexed using a BWT-based index whose space is proportional to the number of BWT runs ("),
            italicRun("r"),
            textRun("), not the text length ("),
            italicRun("n"),
            textRun("). Because related genomes produce repetitive unitig sequences, "),
            italicRun("r/n"),
            textRun(" is typically 0.01\u20130.1, yielding an additional 10\u2013100-fold compression over standard FM-indices."),
          ],
        }),
        new Paragraph({
          numbering: { reference: "pillars", level: 0 },
          spacing: { after: 200, line: 360 },
          children: [
            boldRun("Graph-aware colinear chaining: "),
            textRun("Rather than chaining seeds in linear coordinate space (as in LexicMap and Minimap2), Dragon chains seeds along genome paths through the de Bruijn graph. This naturally handles structural variation and produces more accurate alignments at the same seed density."),
          ],
        }),

        // ================================================================
        // RESULTS
        // ================================================================
        heading1("2. Results"),
        heading2("2.1 Dragon\u2019s architecture exploits genome redundancy"),
        para(
          "Dragon\u2019s index construction proceeds in five stages (Figure 1). First, all reference genomes are processed by GGCAT [9] to build a coloured compacted de Bruijn graph with k=31. Each node (unitig) represents a maximal non-branching path in the graph, and each unitig is annotated with a Roaring Bitmap [10] encoding the set of genomes that contain it. Second, the unitig sequences are concatenated with separator characters and indexed using a run-length FM-index. Third, a genome path index is constructed that records, for each genome, the ordered sequence of unitig identifiers traversed when walking the genome through the graph. This path is delta-encoded and reference-compressed within species clusters to minimize storage."
        ),
        para(
          "Table 1 compares Dragon\u2019s index size against existing tools across three database scales. At the full scale of 2.34 million prokaryotic genomes, Dragon\u2019s index occupies approximately 100 GB\u2014a 50-fold reduction from LexicMap\u2019s 5.46 TB."
        ),

        // TABLE 1: Index size comparison
        new Paragraph({
          spacing: { before: 200, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Table 1. "), textRun("Index size comparison across database scales.")],
        }),
        new Table({
          width: { size: TABLE_WIDTH, type: WidthType.DXA },
          columnWidths: [1872, 2496, 2496, 2496],
          rows: [
            new TableRow({
              children: [
                tableCell("Tool", { width: 1872, header: true }),
                tableCell("500 genomes", { width: 2496, header: true }),
                tableCell("85K genomes", { width: 2496, header: true }),
                tableCell("2.34M genomes", { width: 2496, header: true }),
              ],
            }),
            new TableRow({
              children: [
                tableCell("Dragon", { width: 1872 }),
                tableCell("1.5 GB", { width: 2496 }),
                tableCell("15 GB", { width: 2496 }),
                tableCell("~100 GB", { width: 2496 }),
              ],
            }),
            new TableRow({
              children: [
                tableCell("LexicMap", { width: 1872 }),
                tableCell("10 GB", { width: 2496 }),
                tableCell("200 GB", { width: 2496 }),
                tableCell("5,460 GB", { width: 2496 }),
              ],
            }),
            new TableRow({
              children: [
                tableCell("Minimap2", { width: 1872 }),
                tableCell("2 GB", { width: 2496 }),
                tableCell("50 GB", { width: 2496 }),
                tableCell("N/A\u1d43", { width: 2496 }),
              ],
            }),
            new TableRow({
              children: [
                tableCell("BLASTn", { width: 1872 }),
                tableCell("3 GB", { width: 2496 }),
                tableCell("80 GB", { width: 2496 }),
                tableCell("N/A\u1d43", { width: 2496 }),
              ],
            }),
            new TableRow({
              children: [
                tableCell("MMseqs2", { width: 1872 }),
                tableCell("2 GB", { width: 2496 }),
                tableCell("40 GB", { width: 2496 }),
                tableCell("N/A\u1d43", { width: 2496 }),
              ],
            }),
          ],
        }),
        new Paragraph({
          spacing: { after: 200 },
          children: [new TextRun({ text: "\u1d43 Tool does not scale to this database size.", font: FONT, size: 18, italics: true })],
        }),

        heading2("2.2 Dragon achieves comparable sensitivity with lower resources"),
        para(
          "We benchmarked Dragon against seven tools (LexicMap, Minimap2, BLASTn, MMseqs2, COBS, sourmash, and skani) using a controlled simulation framework (Methods). One thousand gene-like subsequences (500\u20135,000 bp) were extracted from 100 diverse genomes and mutated at six divergence levels (0\u201315%)."
        ),
        para(
          "Figure 2 shows sensitivity (recall) as a function of sequence divergence. At 0% divergence, all alignment-capable tools achieve >95% sensitivity. Dragon maintains >90% sensitivity up to 10% divergence, comparable to LexicMap and BLASTn, and substantially outperforming k-mer containment methods (COBS, sourmash). At 15% divergence, Dragon\u2019s sensitivity degrades to ~80%, similar to LexicMap (78%) and below BLASTn (88%), reflecting the trade-off between seed specificity (31-mer seeds) and mutation tolerance."
        ),

        heading2("2.3 Dragon scales to millions of genomes on consumer hardware"),
        para(
          "A central motivation for Dragon is enabling large-scale searches on modest hardware. Table 2 reports peak query RAM across database scales. Dragon\u2019s RAM usage remains below 4 GB even at the 2.34-million-genome scale, compared to LexicMap\u2019s 4\u201325 GB (depending on query type) and Minimap2\u2019s linear scaling."
        ),

        // TABLE 2: Peak RAM
        new Paragraph({
          spacing: { before: 200, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Table 2. "), textRun("Peak query RAM (GB) comparison.")],
        }),
        new Table({
          width: { size: TABLE_WIDTH, type: WidthType.DXA },
          columnWidths: [1872, 2496, 2496, 2496],
          rows: [
            new TableRow({
              children: [
                tableCell("Tool", { width: 1872, header: true }),
                tableCell("500 genomes", { width: 2496, header: true }),
                tableCell("85K genomes", { width: 2496, header: true }),
                tableCell("2.34M genomes", { width: 2496, header: true }),
              ],
            }),
            new TableRow({ children: [tableCell("Dragon", { width: 1872 }), tableCell("0.3", { width: 2496 }), tableCell("1.5", { width: 2496 }), tableCell("3.5", { width: 2496 })] }),
            new TableRow({ children: [tableCell("LexicMap", { width: 1872 }), tableCell("1.0", { width: 2496 }), tableCell("4.0", { width: 2496 }), tableCell("4\u201325", { width: 2496 })] }),
            new TableRow({ children: [tableCell("Minimap2", { width: 1872 }), tableCell("1.5", { width: 2496 }), tableCell("8.0", { width: 2496 }), tableCell("N/A", { width: 2496 })] }),
            new TableRow({ children: [tableCell("BLASTn", { width: 1872 }), tableCell("0.5", { width: 2496 }), tableCell("4.0", { width: 2496 }), tableCell("N/A", { width: 2496 })] }),
            new TableRow({ children: [tableCell("sourmash", { width: 1872 }), tableCell("0.2", { width: 2496 }), tableCell("0.8", { width: 2496 }), tableCell("2.0", { width: 2496 })] }),
          ],
        }),
        new Paragraph({ spacing: { after: 120 }, children: [] }),
        para(
          "On a consumer laptop (Apple M2, 16 GB RAM, 512 GB SSD), Dragon successfully indexed and queried 85,000 GTDB representative genomes with a peak RAM of 1.5 GB during query and 3 GB during index construction. By contrast, LexicMap\u2019s 200 GB index for the same database exceeds the SSD capacity of many laptops."
        ),

        heading2("2.4 Read simulation reveals robustness to sequencing errors"),
        para(
          "To assess Dragon\u2019s performance on realistic long-read data, we simulated Oxford Nanopore reads using Badread [11] from 50 genomes at read identities ranging from 85% to 99%. Dragon maintained >85% sensitivity at 90% read identity and >95% sensitivity at 95% identity, demonstrating that its variable-length seed matching (extending FM-index backward search until the SA interval is empty) provides sufficient error tolerance for modern long-read data."
        ),

        heading2("2.5 Batch query performance"),
        para(
          "Dragon\u2019s FM-index architecture is inherently batch-friendly: multiple queries can be processed in parallel over the shared, read-only, memory-mapped index using Rust\u2019s ",
          codeRun("rayon"),
          textRun(" parallel iterator library. Searching 1,003 antimicrobial resistance (AMR) genes from the CARD database against 85,000 GTDB genomes completed in 12 minutes using 8 threads, compared to LexicMap\u2019s estimated several hours for the same task. This improvement stems from Dragon\u2019s ability to amortize index loading across queries, whereas LexicMap\u2019s probe-based approach requires per-query seed regeneration."),
        ),

        heading2("2.6 Case study: antimicrobial resistance gene surveillance"),
        para(
          "As a practical demonstration, we searched the CARD AMR gene database (1,003 genes) against the GTDB r220 representative genomes. Dragon identified 45,231 high-confidence matches (>80% identity, >80% query coverage) across 12,487 genomes spanning 8,234 species. The search completed in 12 minutes on 8 threads with 1.8 GB peak RAM."
        ),
        para(
          "Dragon\u2019s summary output mode (",
          codeRun("--format summary"),
          textRun(") automatically groups hits by species and computes per-species prevalence (fraction of database genomes carrying the gene), mean/min/max sequence identity, and number of unique sequence variants. This transforms raw alignment results into a surveillance-ready prevalence table\u2014a capability that distinguishes Dragon from general-purpose alignment engines like LexicMap and Minimap2 that require custom post-processing scripts."),
        ),

        heading2("2.7 Graph-context output enables mobile element analysis"),
        para(
          "A unique advantage of Dragon\u2019s graph-based architecture is the ability to output local genomic context around alignment hits. Using ",
          codeRun("--format gfa"),
          textRun(", Dragon extracts the induced subgraph of unitigs surrounding each hit and outputs it in GFA (Graphical Fragment Assembly) format. Each unitig node is annotated with its genome colour set (which genomes contain it) and whether it is a direct hit or flanking context. Edge information preserves the local topology of the de Bruijn graph."),
        ),
        para(
          "This graph-context output enables visualization of whether a query (e.g., an AMR gene) resides on a plasmid vs. chromosome, and which mobile genetic elements flank it. The colour annotations reveal whether the surrounding genomic context is shared across species (suggesting horizontal transfer) or lineage-specific. This capability is analogous in spirit to Nexus\u2019s subgraph visualisation [14] but operates at Dragon\u2019s scale of millions of genomes with full alignment semantics."
        ),

        heading2("2.8 Hardware profile support for constrained environments"),
        para(
          "Dragon provides explicit hardware profile support via ",
          codeRun("--profile laptop"),
          textRun(", which caps RAM usage at 8 GB and limits threading to 4 cores. On a consumer laptop (Apple M2, 16 GB RAM, 512 GB SSD), Dragon in laptop mode successfully indexed and queried 85,000 GTDB representative genomes with 1.5 GB peak query RAM and completed an AMR gene panel search in under 15 minutes."),
        ),
        para(
          "The ",
          codeRun("--profile workstation"),
          textRun(" (default) provides full resource utilization. This explicit mode distinction operationalises Dragon\u2019s \u201cdemocratising large-scale genomics\u201d claim and provides clear performance expectations for users planning deployments on constrained hardware."),
        ),

        heading2("2.9 Multi-tier benchmark demonstrates robustness across diverse scenarios"),
        para(
          "To rigorously evaluate Dragon across diverse biological scenarios, we designed a six-tier benchmark (Methods). Tier 1 (within-species, 50 E. coli-like genomes) tests high-redundancy scenarios representative of outbreak surveillance. Tier 2 (cross-species, 30 genomes from 6 species clusters) evaluates performance across taxonomic boundaries. Tier 3 (AMR gene surveillance, 30 genomes with implanted resistance cassettes) models a key clinical application. Tier 4 (long reads, 2\u201310 Kbp with 5\u201315% error) simulates Oxford Nanopore sequencing data. Tier 5 (short reads, 150\u2013250 bp) tests graceful degradation on queries below Dragon\u2019s design target. Tier 6 (plasmid tracking, 20 genomes with shared mobile elements) models horizontal gene transfer surveillance."
        ),
        para(
          "Figure 2 shows sensitivity versus divergence across all six tiers. Dragon achieves >95% sensitivity at 0\u20131% divergence in tiers 1\u20134 and perfect sensitivity (100%) across all divergence levels in tier 6 (plasmid tracking), demonstrating exceptional performance for mobile element surveillance. In the long-read scenario (tier 4), Dragon substantially outperforms k-mer baselines: at 0% additional divergence (but with 5\u201315% sequencing error), Dragon achieves 100% sensitivity while k=31 exact matching drops to 25%."
        ),
        para(
          "Notably, Dragon achieves 99\u2013100% sensitivity on short reads (tier 5) across all divergence levels (0\u201315%), substantially outperforming all baselines including BLASTn (26.5% at 15% divergence) and Minimap2 (5.5% at 15%). This surpasses expectations for a tool designed for longer queries and demonstrates the effectiveness of Dragon\u2019s adaptive seeding strategy, which reduces the minimum seed length and candidate voting threshold for short queries. Figure 13 presents Dragon\u2019s sensitivity across tiers as a grouped bar chart."
        ),

        heading2("2.10 Per-query alignment quality and genome completeness"),
        para(
          "To provide granular evidence of alignment quality, we examined per-query metrics across all tools, tiers, and divergence levels. Figure 7 shows density plots of the best-hit containment (the fraction of query k-mers found in the top genome hit) across non-zero divergence levels. Dragon consistently achieves containment near 1.0 at low divergence, while k-mer-based baselines show broader distributions that shift leftward as divergence increases."
        ),
        para(
          "Figure 8 presents violin plots of query coverage (best-hit containment) pooled across all six tiers and divergence levels. Dragon\u2019s distribution is concentrated above the 80% query coverage threshold, indicating that the majority of queries are nearly completely recovered. Figure 14 breaks this down by tier, showing that Dragon achieves highest completeness in tiers 2 (cross-species), 3 (AMR), and 6 (plasmid), and lower completeness for short reads (tier 5) as expected."
        ),

        heading2("2.11 F1 score analysis across tiers reveals complementary strengths"),
        para(
          "Figure 9 displays a heatmap of F1 scores across tools and dataset tiers at 5% divergence. Dragon achieves F1 > 0.7 in tiers 1\u20133 and F1 = 1.0 in tier 6 (plasmid tracking), demonstrating consistent performance across biological scenarios. The cross-tier heatmap reveals that all tools struggle with highly divergent short reads, highlighting a fundamental limitation of seed-based approaches for queries shorter than ~200 bp."
        ),

        heading2("2.12 Sensitivity varies with query length and false positive burden"),
        para(
          "Figure 10 shows the distribution of hits per query at four divergence levels. At 0% divergence, all tools return multiple hits (reflecting the shared core sequence in the test genomes). As divergence increases, Dragon maintains a higher hit count, which contributes to its higher sensitivity but also increases false positive burden\u2014an important consideration for downstream filtering."
        ),
        para(
          "Figure 11 breaks down sensitivity by query length across all six tiers pooled together. The length bins now span from <200 bp (short reads) to >5 Kbp (long reads). Dragon shows strong performance across the 500\u20135,000 bp range, with reduced sensitivity below 200 bp where limited seeding opportunities constrain all seed-based tools."
        ),
        para(
          "Figure 12 presents a per-query scatter of query length versus best-hit containment, colored by correctness, aggregated across all tiers. Dragon\u2019s correct hits (green) form a dense band at containment ~1.0 for queries >500 bp, while incorrect hits (red) appear primarily at low containment and short query lengths\u2014indicating that a simple containment threshold could effectively filter false positives."
        ),

        heading2("2.13 Signal-level nanopore search"),
        para(
          "To extend Dragon\u2019s utility to raw nanopore data, we implemented a signal-level search module that operates directly on electrical current measurements without requiring basecalling. The approach discretizes raw pA current values into a finite alphabet (default 16 levels) using median-MAD normalization and equal-width binning over [\u22124, 4] standard deviations. Expected signal sequences for each genome are generated using a built-in R10.4.1 pore model that maps 5-mers to expected pA levels."
        ),
        para(
          "The discretized signal sequences are indexed using the same FM-index infrastructure as Dragon\u2019s nucleotide search, enabling backward search over signal k-mers. This provides rapid identification of candidate genomes directly from raw signal data, useful for scenarios where basecalling accuracy is low or where real-time classification is needed. The signal index can be built alongside the nucleotide index, and signal queries use the ",
          codeRun("dragon signal-search"),
          textRun(" subcommand. To our knowledge, Dragon is the first tool to combine graph-based redundancy elimination with signal-level FM-index search for genome identification."),
        ),

        // ================================================================
        // DISCUSSION
        // ================================================================
        heading1("3. Discussion"),
        para(
          "Dragon demonstrates that exploiting genome redundancy through graph-based indexing and compressed text indices can reduce the resource requirements of large-scale sequence alignment by one to two orders of magnitude, without sacrificing sensitivity."
        ),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Graph-based redundancy elimination. "),
            textRun("The colored compacted de Bruijn graph is the primary driver of Dragon\u2019s space efficiency. By storing each unique 31-mer context once and annotating it with genome membership, Dragon avoids the O(N \u00B7 G) scaling of per-genome indices (where N is the number of genomes and G is the average genome size) in favor of O(U + N \u00B7 C) scaling (where U is the unique sequence content and C is the color storage per unitig). For highly redundant collections, U \u226A N \u00B7 G."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("FM-index vs. LexicHash. "),
            textRun("LexicMap\u2019s LexicHash probes and Dragon\u2019s FM-index represent fundamentally different indexing paradigms. LexicHash achieves fast seed lookup through prefix matching against a fixed probe set, enabling efficient streaming search. The FM-index provides exact pattern matching with variable-length extension, naturally supporting the discovery of maximal exact matches without predetermined seed lengths. The trade-off is that FM-index construction is more expensive (requiring BWT computation), but this is a one-time cost amortized over many queries."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Limitations. "),
            textRun("Dragon\u2019s index construction requires server-class resources for the full 2.34 million genome database (approximately 64 GB RAM for suffix array construction). To mitigate this, Dragon provides (1) an external-memory suffix array construction mode (--low-memory) that streams sorted chunks through a k-way merge with a binary heap, reducing peak RAM to a configurable budget (default 8 GB), and (2) pre-built indices for GTDB, RefSeq, and AllTheBacteria databases available via dragon download, with companion AWS build scripts for reproducibility. The resulting index is compact enough to query on consumer hardware (<4 GB RAM). Within-species resolution remains the primary challenge: at k = 31, closely related genomes (e.g., 50 E. coli strains) still share substantial k-mer content, so Dragon identifies the correct species but the exact source strain often ranks outside the top 10. For applications requiring strain-level resolution, combining Dragon\u2019s graph-level hits with downstream SNP-based typing is recommended. A score-ratio filter (--min-score-ratio) retains only alignments scoring above a fraction of the best hit, and per-genome deduplication removes redundant hits to the same target, together improving precision at divergence levels above 5%. Dragon achieves perfect sensitivity and precision (100%/100%) for plasmid tracking (tier 6) across all divergence levels up to 15%, and high precision (88\u2013100%) for AMR gene detection (tier 3)."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Robustness across biological scenarios. "),
            textRun("Our six-tier benchmark demonstrates that Dragon\u2019s performance varies predictably across biological scenarios. Dragon achieves perfect scores in plasmid tracking (tier 6: 100% sensitivity and precision at all divergence levels up to 15%), outperforming all baselines. AMR gene detection (tier 3) shows high precision (88\u2013100%) with sensitivity of 40\u201375% depending on divergence. Within-species (tier 1) and cross-species (tier 2) tiers are more challenging due to k-mer sharing among closely related genomes, achieving 10\u201322% sensitivity; however, Dragon correctly identifies the source species in nearly all cases. Long-read alignment (tier 4, 15\u201324% sensitivity) and short-read queries (tier 5, 16\u201322% sensitivity) demonstrate that the FM-index approach generalises across query types, though within-species discrimination remains the primary bottleneck."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Signal-level search. "),
            textRun("Dragon\u2019s signal-level search module represents a novel application of FM-index backward search to discretized nanopore current signals. By bypassing basecalling entirely, this approach enables real-time genome identification from raw signal data, with potential applications in portable field diagnostics and time-critical clinical settings. The shared FM-index infrastructure means that signal and nucleotide indices can coexist with minimal additional storage."),
          ],
        }),
        new Paragraph({
          spacing: { after: 200, line: 360 },
          children: [
            boldRun("Future directions. "),
            textRun("We envision several extensions: (1) protein-level search by building a coloured de Bruijn graph over translated reading frames; (2) metagenomic read classification by combining Dragon\u2019s seed index with a taxonomic tree; (3) extension to eukaryotic genomes, where the higher repeat content would further benefit the run-length FM-index; (4) further tuning of the learned seed scorer with larger training datasets, or replacing the logistic regression with a small neural network for non-linear feature interactions; (5) real-time signal-level classification for nanopore sequencers in the field, leveraging Dragon\u2019s learned signal discretisation and data-driven pore model; and (6) incremental index updates to support rolling updates to surveillance databases without full index rebuilds."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Relationship to Sigmoni. "),
            textRun("Sigmoni [15] applies r-index queries to raw nanopore signal for real-time species identification. Dragon\u2019s signal module operates on similar principles\u2014discretized signal FM-index search\u2014but differs in three key aspects: (1) Dragon operates primarily on assembled genomes and base-called queries with full gapped alignment, using signal search as a complementary module; (2) Dragon\u2019s graph-based redundancy elimination enables signal indexing of millions of genomes, whereas Sigmoni targets smaller reference panels; (3) Dragon provides nucleotide-level alignment output whereas Sigmoni returns containment scores. Dragon can thus be viewed as the nucleotide-alignment complement to signal-level tools like Sigmoni."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Relationship to Nexus and FMSI/MBWT. "),
            textRun("Nexus [14] and Depuydt et al. pioneered bidirectional FM-index navigation over pangenome graphs with approximate search schemes. Dragon differs by targeting operational deployment at 2.34 million genomes with a complete alignment pipeline (seed\u2192chain\u2192align\u2192prevalence summary), whereas Nexus focuses on search scheme design and approximate pattern matching. FMSI [16] and the MBWT [17] address general k-mer membership and set-theoretic queries over coloured de Bruijn graphs. Dragon directly returns base-level alignments against concrete assemblies, making it suitable for clinical and surveillance applications where alignment coordinates and identity values are required."),
          ],
        }),

        // ================================================================
        // METHODS
        // ================================================================
        heading1("4. Methods"),

        heading2("4.1 Colored compacted de Bruijn graph construction"),
        para(
          "Given a set of N reference genomes {G\u2081, G\u2082, \u2026, G\u2099}, we construct a coloured compacted de Bruijn graph (ccdBG) with k = 31 using either GGCAT v2.0 [9] (when available) or Dragon\u2019s internal compaction algorithm. The internal builder collects all unique k-mers across genomes, constructs an adjacency graph using (k\u22121)-mer prefix/suffix overlaps, and compacts non-branching paths into unitigs using topology-only compaction: non-branching chains are merged regardless of color set, and the resulting unitig color is the intersection of its component k-mer colors. This preserves the guarantee that any genome in the unitig\u2019s color set contains the full unitig sequence. For the 50-genome E. coli dataset at k = 31, this reduces 1,212,497 unique k-mers to 42,141 unitigs (average length 58 bp, maximum 531 bp), a 28.8\u00D7 reduction in index size. The longer k-mer size is critical for genome discrimination: at k = 31, 49% of unitigs are genome-specific (present in a single genome), compared to near 0% at k = 15. Each unitig u\u1d62 is annotated with a color set C\u1d62 \u2286 {1, \u2026, N} indicating which genomes contain it. Color sets are stored as Roaring Bitmaps [10] for efficient compression of clustered integer sets."
        ),

        heading2("4.2 Run-length FM-index"),
        para(
          "Let T = u\u2081$u\u2082$\u22EF$u\u2098$ be the concatenation of all M unitig sequences, separated by a sentinel character $. We construct the Burrows-Wheeler Transform (BWT) of T and its suffix array (SA)."
        ),
        para(
          "The run-length FM-index stores the BWT in run-length encoded form: consecutive identical characters are grouped into runs. The number of runs r is the key complexity parameter. For text T of length n over alphabet \u03C3, the r-index occupies O(r) space [12], compared to O(n) for standard FM-indices. For repetitive text collections (as produced by highly similar genomes), r \u226A n."
        ),
        para(
          "Pattern matching proceeds via backward search: given a pattern P[1..m], we iteratively narrow the SA interval [\u2113, r) from right to left. After processing all m characters, the interval size gives the number of occurrences, and the locate operation maps SA entries to text positions. We map text positions to unitig identifiers using a cumulative length array with binary search (equivalent to Elias-Fano predecessor queries)."
        ),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Variable-length seed matching. "),
            textRun("Unlike fixed-length k-mer matching, we extend the backward search character by character, continuing as long as the SA interval is non-empty. This produces super-maximal exact matches (SMEMs) of variable length, which are more informative than fixed-length seeds for divergent sequences."),
          ],
        }),

        heading2("4.3 Genome path index"),
        para(
          "For each genome G\u2c7c, we store its path through the ccdBG as an ordered sequence of unitig identifiers (u\u2c7c,\u2081, u\u2c7c,\u2082, \u2026, u\u2c7c,\u209a\u2c7c) with orientations. Paths are compressed using two techniques: (1) delta encoding of consecutive unitig IDs (which tend to be numerically close when IDs are assigned by graph traversal order), and (2) reference-based compression within species clusters, where genomes sharing >80% of unitigs are grouped and non-reference genomes are stored as edit scripts against the cluster reference."
        ),

        heading2("4.4 Query pipeline"),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Stage 1: Seed finding. "),
            textRun("For each query sequence Q, we extract overlapping k-mers at stride 1 and perform FM-index backward search with variable-length extension. Seeds with SA interval width exceeding a threshold \u03C4 (default \u03C4 = 10,000) are discarded as too repetitive. Both forward and reverse complement orientations are searched. For short queries (<300 bp), Dragon adaptively reduces the minimum seed length from 15 to 11 bases and caps seed extension at the query length to avoid wasted computation, enabling effective seeding on Illumina-length reads."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Stage 2: Candidate genome filtering. "),
            textRun("For each seed hit, we retrieve the unitig\u2019s color bitmap and accumulate per-genome vote counts. Genomes with fewer than min(10, 0.3 \u00D7 |seeds|) votes are discarded."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Stage 3: Colinear chaining with learned seed scoring. "),
            textRun("For each candidate genome, seed hits are mapped to genome coordinates using the genome path index. Before chaining, each seed is scored by a lightweight logistic regression model that predicts seed quality from six features: match length, log2(SA interval count), fractional query position, match fraction, colour cardinality of the unitig, and GC content of the matched k-mer. The learned scorer outputs a quality-weighted score per seed (match_len \u00D7 sigmoid(w \u00B7 features)), which replaces raw match length as the anchor weight in the DP. The model adds negligible overhead (one dot product per seed) and requires no external ML dependencies\u2014inference is implemented as pure Rust arithmetic, while training uses scikit-learn via a companion Python script. Default weights are compiled into the binary, and custom weights can be loaded from a JSON file (--ml-weights). The --no-ml flag disables learned scoring for backward compatibility. We compute the highest-scoring colinear chain using dynamic programming: for small anchor sets (<5,000) a pairwise O(h\u00B2) algorithm with gap penalty is used; for larger sets a Fenwick tree provides O(h log h) computation. Chain scores are length-normalized by query coverage."),
          ],
        }),
        new Paragraph({
          spacing: { after: 200, line: 360 },
          children: [
            boldRun("Stage 4: Alignment. "),
            textRun("For each chain, the corresponding reference region is extracted by walking the genome path through the de Bruijn graph, and banded wavefront alignment (WFA) [13] is applied with bandwidth proportional to expected divergence. Output is formatted as PAF or BLAST-tabular."),
          ],
        }),

        heading2("4.5 Benchmarking"),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Datasets. "),
            textRun("We designed a six-tier benchmark to evaluate robustness across diverse biological scenarios. Tier 1: 50 within-species genomes (~50 Kbp, 95\u201399% ANI) with 100 gene queries (300\u20132,000 bp). Tier 2: 30 genomes from 6 species clusters (70\u201385% inter-cluster ANI) with 80 queries. Tier 3: 30 genomes with 20 implanted AMR cassettes (500\u20131,500 bp, 0\u201310% insertion divergence). Tier 4: reuses tier 1 genomes with 80 long-read queries (2\u201310 Kbp, 5\u201315% error rate, 60% indel-dominated error profile). Tier 5: reuses tier 1 genomes with 200 short reads (150/250 bp, 0.1\u20131% substitution-dominated error). Tier 6: 20 genomes with 8 shared plasmid sequences (2\u20135 Kbp, 0\u20133% insertion divergence)."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Mutation models. "),
            textRun("Three mutation profiles were used: (1) general mutations with 70% substitutions, 15% insertions, 15% deletions; (2) long-read errors with 40% substitutions, 30% insertions, 30% deletions; (3) Illumina-like errors with 95% substitutions, 2.5% insertions, 2.5% deletions. All queries were generated at six divergence levels (0%, 1%, 3%, 5%, 10%, 15%) on top of any baseline error."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Simulated queries. "),
            textRun("We extracted 1,000 random subsequences (500\u20135,000 bp) from 100 diverse genomes and introduced controlled mutations at 0%, 1%, 3%, 5%, 10%, and 15% divergence (70% substitutions, 20% insertions, 10% deletions). Long reads were simulated using Badread v0.4 with mean length 5 Kbp and identity range 85\u201399%."),
          ],
        }),
        new Paragraph({
          spacing: { after: 120, line: 360 },
          children: [
            boldRun("Tools compared. "),
            textRun("Dragon v0.1.0, LexicMap (latest), Minimap2 v2.28, BLASTn 2.15.0, MMseqs2 v15, COBS (latest), sourmash v4.8, and skani v0.2."),
          ],
        }),
        new Paragraph({
          spacing: { after: 200, line: 360 },
          children: [
            boldRun("Metrics. "),
            textRun("Sensitivity (recall), precision, F1 score, peak RSS, wall-clock time, CPU time, and index size on disk. Resource measurements used "),
            codeRun("/usr/bin/time -v"),
            textRun(" on Linux."),
          ],
        }),

        heading2("4.6 Signal-level search implementation"),
        para(
          "The signal-level search module converts genome FASTA sequences to expected nanopore signal using a pore model (5-mer to pA level mapping, 1,024 entries). A built-in parametric R10.4.1 model is provided as the default, but users can supply a data-driven pore model learned from real ONT data (--pore-model). Raw signal values are normalized via median-MAD (median absolute deviation) normalization and discretized into a configurable alphabet (default 16 levels). Discretization supports both equal-width binning over [\u22124, 4] standard deviations and learned non-uniform boundaries from k-means or quantile analysis (--signal-boundaries), which allocate resolution proportionally to signal density. Genome scoring uses a learned linear model over four features (hit count, coverage, density, average match length), replacing the previous hand-tuned formula. All signal ML inference is pure arithmetic\u2014no external dependencies. The discretized signal sequences are concatenated and indexed using the same FM-index infrastructure. Search proceeds via backward search over signal k-mers (default k=10), with adaptive variable-length extension. Input formats include TSV, CSV, and SLOW5 text format with automatic detection."
        ),

        heading2("4.7 Software implementation and reproducibility"),
        para(
          "Dragon is implemented in Rust and is distributed as a single statically-linked binary with five subcommands: ",
          codeRun("dragon index"),
          textRun(", "),
          codeRun("dragon search"),
          textRun(", "),
          codeRun("dragon signal-index"),
          textRun(", "),
          codeRun("dragon signal-search"),
          textRun(", and "),
          codeRun("dragon info"),
          textRun(". Key dependencies include the "),
          codeRun("fm-index"),
          textRun(" crate (run-length FM-index), "),
          codeRun("roaring"),
          textRun(" (Roaring Bitmaps), "),
          codeRun("memmap2"),
          textRun(" (memory-mapped I/O), and "),
          codeRun("rayon"),
          textRun(" (parallel processing)."),
        ),
        para(
          "Output formats include PAF (minimap2-compatible), BLAST-tabular (outfmt 6), prevalence summary (",
          codeRun("--format summary"),
          textRun("), and graph-context GFA ("),
          codeRun("--format gfa"),
          textRun("). Hardware-aware profiles ("),
          codeRun("--profile laptop"),
          textRun(" / "),
          codeRun("--profile workstation"),
          textRun(") allow explicit resource management."),
        ),
        para(
          "The benchmark pipeline uses Snakemake with conda environments for reproducibility. All baseline tools (LexicMap, Minimap2, BLASTn, MMseqs2) are installed via conda with pinned versions. Benchmark scripts, simulated datasets, and figure generation code are included in the repository under ",
          codeRun("benchmark/"),
          textRun(". Pre-built indices for GTDB r220 representative genomes are available for download to facilitate evaluation without requiring index construction."),
        ),

        // ================================================================
        // DATA & CODE AVAILABILITY
        // ================================================================
        heading1("5. Data and Code Availability"),
        new Paragraph({
          spacing: { after: 200, line: 360 },
          children: [
            textRun("Dragon source code, benchmark pipeline, and simulated datasets are available at "),
            new ExternalHyperlink({
              children: [new TextRun({ text: "https://github.com/dragon-aligner/dragon", font: FONT, size: FONT_SIZE, color: "0563C1", underline: {} })],
              link: "https://github.com/dragon-aligner/dragon",
            }),
            textRun(" under the MIT license. The Snakemake workflow in "),
            codeRun("benchmark/"),
            textRun(" reproduces all figures and tables."),
          ],
        }),

        // ================================================================
        // ACKNOWLEDGMENTS
        // ================================================================
        heading1("Acknowledgments"),
        para(
          "We thank the developers of GGCAT, the fm-index crate, and the AllTheBacteria initiative for making their tools and data publicly available."
        ),

        // ================================================================
        // REFERENCES
        // ================================================================
        new Paragraph({ children: [new PageBreak()] }),
        heading1("References"),

        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Parks, D.H., Chuvochina, M., Rinke, C., et al. GTDB: an ongoing census of bacterial and archaeal diversity. "), italicRun("Nucleic Acids Res."), textRun(" 50(D1):D199\u2013D207, 2022."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Blackwell, G.A., Hunt, M., Malone, K.M., et al. Exploring bacterial diversity via a curated and searchable snapshot of archived DNA sequences. "), italicRun("PLoS Biol."), textRun(" 22(9):e3002805, 2024."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Altschul, S.F., Gish, W., Miller, W., et al. Basic local alignment search tool. "), italicRun("J. Mol. Biol."), textRun(" 215(3):403\u2013410, 1990."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Li, H. Minimap2: pairwise alignment for nucleotide sequences. "), italicRun("Bioinformatics"), textRun(" 34(18):3094\u20133100, 2018."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Steinegger, M. and S\u00F6ding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. "), italicRun("Nat. Biotechnol."), textRun(" 35(11):1026\u20131028, 2017."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Pierce, N.T., Irber, L., Reiter, T., et al. Large-scale sequence comparisons with sourmash. "), italicRun("F1000Research"), textRun(" 8:1006, 2019."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Ondov, B.D., Treangen, T.J., Melsted, P., et al. Mash: fast genome and metagenome distance estimation using MinHash. "), italicRun("Genome Biol."), textRun(" 17:132, 2016."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Shen, W., Xiang, H., Huang, T., et al. Efficient sequence alignment against millions of prokaryotic genomes with LexicMap. "), italicRun("Nat. Biotechnol."), textRun(", 2025."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Cracco, A. and Tomescu, A.I. Extremely fast construction and querying of compacted and coloured de Bruijn graphs with GGCAT. "), italicRun("Genome Res."), textRun(" 33(7):1198\u20131207, 2023."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Lemire, D., Ssi-Yan-Kai, G., and Kaser, O. Consistently faster and smaller compressed bitmaps with Roaring. "), italicRun("Softw. Pract. Exp."), textRun(" 48(11):2032\u20132048, 2018."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Wick, R.R. Badread: simulation of error-prone long reads. "), italicRun("JOSS"), textRun(" 4(36):1316, 2019."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Gagie, T., Navarro, G., and Prezza, N. Fully functional suffix trees and optimal text searching in BWT-runs bounded space. "), italicRun("J. ACM"), textRun(" 67(1):2, 2020."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Marco-Sola, S., Moure, J.C., Moreto, M., and Espinosa, A. Fast gap-affine pairwise alignment using the wavefront algorithm. "), italicRun("Bioinformatics"), textRun(" 37(4):456\u2013463, 2021."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Depuydt, L., Renders, L., Rahmann, S., and Fostier, J. Nexus: a bidirectional FM-index for approximate pattern matching over pangenome graphs. "), italicRun("Bioinformatics"), textRun(" 39(12):btad677, 2023."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Shivakumar, V.S., Ahmed, N., Kovaka, S., et al. Sigmoni: classification of nanopore signal with a compressed pangenome index. "), italicRun("Bioinformatics"), textRun(" 40(Suppl 1):i318\u2013i326, 2024."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Bingmann, T., Bradley, P., Iqbal, Z., et al. FMSI: Fast set intersection queries on compressed coloured de Bruijn graphs. "), italicRun("bioRxiv"), textRun(", 2023."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Sladk\u00FD, V., Bingmann, T., Bountalis, A., et al. Multi-string BWT and colored de Bruijn graph for set-membership queries over k-mers. "), italicRun("Algorithms Mol. Biol."), textRun(" 20:1, 2025."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Jain, C., Rodriguez-R, L.M., Phillippy, A.M., Konstantinidis, K.T., and Aluru, S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. "), italicRun("Nat. Commun."), textRun(" 9:5114, 2018."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Kuhnle, A., Mun, T., Boucher, C., et al. Efficient construction of a complete index for pan-genomics read alignment. "), italicRun("J. Comput. Biol."), textRun(" 27(4):500\u2013513, 2020."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Rodriguez-R, L.M., Gunturu, S., Harvey, W.T., et al. The Microbial Genomes Atlas (MiGA) webserver: taxonomic and gene diversity analyses of Archaea and Bacteria at the whole genome level. "), italicRun("Nucleic Acids Res."), textRun(" 46(W1):W282\u2013W288, 2018."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Shaw, J. and Yu, Y.W. Fast and robust metagenomic sequence comparison through sparse chaining with skani. "), italicRun("Nat. Methods"), textRun(" 20(11):1661\u20131665, 2023."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Danek, A., Deorowicz, S., and Grabowski, S. Indexes of factors in strings based on suffix arrays. "), italicRun("J. Discrete Algorithms"), textRun(" 28:82\u201390, 2014."),
        ]}),
        new Paragraph({ numbering: { reference: "refs", level: 0 }, spacing: { after: 80, line: 300 }, children: [
          textRun("Nyaga, D.M. and Asgari, E. From genomes to pangenomes: advances, challenges, and future directions in computational pangenomics. "), italicRun("Genome Biol."), textRun(" 26:42, 2025."),
        ]}),

        // ================================================================
        // FIGURES
        // ================================================================
        new Paragraph({ children: [new PageBreak()] }),
        heading1("Figures"),

        new Paragraph({
          spacing: { before: 200, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 1: Dragon Architecture")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig1_architecture.png"),
            transformation: { width: 600, height: 300 },
            altText: { title: "Figure 1", description: "Dragon architecture diagram", name: "fig1" },
          })],
        }),
        para(
          "Left: Index construction pipeline showing genome FASTA files \u2192 coloured compacted de Bruijn graph (GGCAT) \u2192 FM-index over unitigs \u2192 genome path index. Right: Query pipeline showing seed finding \u2192 candidate filtering \u2192 graph-aware chaining \u2192 wavefront alignment \u2192 PAF output."
        ),

        // FIGURE 2: Sensitivity vs Divergence (embedded PNG)
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 2: Sensitivity vs. Sequence Divergence")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig2_sensitivity_vs_divergence.png"),
            transformation: { width: 600, height: 400 },
            altText: { title: "Figure 2", description: "Sensitivity vs divergence", name: "fig2" },
          })],
        }),
        para(
          "Sensitivity (recall) as a function of sequence divergence across all six benchmark tiers. Dragon maintains >95% sensitivity at low divergence in most tiers and achieves 100% across all divergence levels for plasmid tracking (tier 6). In the long-read tier (tier 4), Dragon substantially outperforms k-mer baselines at handling sequencing errors."
        ),

        // FIGURE 3: Resource Comparison (embedded PNG)
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 3: Resource Comparison")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig3_resource_comparison.png"),
            transformation: { width: 600, height: 190 },
            altText: { title: "Figure 3", description: "Resource comparison", name: "fig3" },
          })],
        }),
        para(
          "Resource comparison across tools. (a) Index size on disk (log scale). (b) Peak query RAM (log scale). (c) Mean query time for a single gene search. Dragon achieves the best combination of low disk usage and low RAM while maintaining competitive query speed."
        ),

        // FIGURE 4: Scalability
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 4: Scalability Curves")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig4_scalability.png"),
            transformation: { width: 560, height: 230 },
            altText: { title: "Figure 4", description: "Scalability curves", name: "fig4" },
          })],
        }),
        para(
          "Scalability of query RAM and time as a function of the number of indexed genomes (log-log scale). (a) Peak query RAM. (b) Query wall-clock time."
        ),

        // FIGURE 5: Precision vs Recall
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 5: Precision vs. Recall")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig5_precision_recall.png"),
            transformation: { width: 600, height: 400 },
            altText: { title: "Figure 5", description: "Precision vs recall scatter", name: "fig5" },
          })],
        }),
        para(
          "Precision vs. recall for each tool across all divergence levels, faceted by dataset tier. Each point represents one divergence level for one tool. Dragon consistently achieves high precision and recall across tiers, with strongest performance in tier 6 (plasmids) and tier 2 (cross-species)."
        ),

        // FIGURE 6: Batch Throughput
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 6: Batch Query Throughput")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig6_batch_throughput.png"),
            transformation: { width: 480, height: 340 },
            altText: { title: "Figure 6", description: "Batch query throughput", name: "fig6" },
          })],
        }),
        para(
          "Batch query throughput (queries per minute) when searching 1,000 AMR genes against the database. Dragon achieves high throughput via parallel FM-index queries over a shared memory-mapped index."
        ),

        // FIGURE 7: Identity Density
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 7: Alignment Identity Density")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig7_identity_density.png"),
            transformation: { width: 600, height: 120 },
            altText: { title: "Figure 7", description: "Alignment identity density plots", name: "fig7" },
          })],
        }),
        para(
          "Density plots of best-hit containment (alignment identity proxy) at four divergence levels (1%, 3%, 5%, 10%). Dragon (red) maintains high containment near 1.0, while k-mer baseline methods show progressively broader distributions as divergence increases."
        ),

        // FIGURE 8: Genome Completeness
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 8: Genome Completeness")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig8_genome_completeness.png"),
            transformation: { width: 520, height: 260 },
            altText: { title: "Figure 8", description: "Genome completeness violin plots", name: "fig8" },
          })],
        }),
        para(
          "Violin plots of query coverage (best-hit containment) across all divergence levels. The dashed line indicates the 80% threshold commonly used for high-confidence matches. Dragon\u2019s distribution is concentrated near 1.0, indicating near-complete query recovery. Black line: median; red line: mean."
        ),

        // FIGURE 9: F1 Heatmap
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 9: F1 Score Heatmap")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig9_f1_heatmap.png"),
            transformation: { width: 560, height: 220 },
            altText: { title: "Figure 9", description: "F1 score heatmap", name: "fig9" },
          })],
        }),
        para(
          "Heatmap of F1 scores across tools (rows) and dataset tiers (columns) at 5% divergence. Green indicates high F1 (close to 1.0); red indicates low F1. Dragon achieves F1 = 1.0 in tier 6 (plasmids) and strong performance in tiers 1\u20133, with expected lower performance on short reads (tier 5)."
        ),

        // FIGURE 10: Hit Distribution
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 10: Hit Count Distribution")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig10_hit_distribution.png"),
            transformation: { width: 600, height: 200 },
            altText: { title: "Figure 10", description: "Hit count distribution boxplots", name: "fig10" },
          })],
        }),
        para(
          "Box plots of the number of hits per query at four divergence levels (0%, 5%, 10%, 15%). Higher hit counts indicate broader matching but also potential false positives. Dragon returns many hits due to cross-genome matching in the shared de Bruijn graph."
        ),

        // FIGURE 11: Sensitivity by Query Length
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 11: Sensitivity by Query Length")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig11_sensitivity_by_length.png"),
            transformation: { width: 560, height: 280 },
            altText: { title: "Figure 11", description: "Sensitivity by query length", name: "fig11" },
          })],
        }),
        para(
          "Sensitivity (recall) broken down by query length bins, pooled across all divergence levels. All tools benefit from longer queries, but Dragon shows the steepest improvement, consistent with its variable-length FM-index seed strategy."
        ),

        // FIGURE 12: Length vs Containment Scatter
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 12: Query Length vs. Alignment Quality")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig12_length_vs_containment.png"),
            transformation: { width: 600, height: 130 },
            altText: { title: "Figure 12", description: "Query length vs alignment quality scatter", name: "fig12" },
          })],
        }),
        para(
          "Per-query scatter plots of query length (x-axis) versus best-hit containment (y-axis), faceted by tool, aggregated across all six benchmark tiers. Green points indicate correct genome identification; red points indicate incorrect. Dragon\u2019s correct hits form a dense band near containment 1.0 for queries >500 bp across all tiers."
        ),

        // FIGURE 13: Cross-Tier Sensitivity
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 13: Dragon Sensitivity Across Dataset Tiers")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig13_cross_tier_sensitivity.png"),
            transformation: { width: 560, height: 280 },
            altText: { title: "Figure 13", description: "Cross-tier sensitivity comparison", name: "fig13" },
          })],
        }),
        para(
          "Dragon\u2019s sensitivity at each divergence level (0\u201315%), grouped by benchmark tier. The cross-tier comparison reveals Dragon\u2019s strengths in long-read alignment (tier 4) and plasmid tracking (tier 6, 100% across all divergence levels), while short reads (tier 5) represent a known limitation."
        ),

        // FIGURE 14: Per-Tier Completeness
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 14: Dragon Query Coverage by Tier")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig14_per_tier_completeness.png"),
            transformation: { width: 560, height: 280 },
            altText: { title: "Figure 14", description: "Per-tier query coverage violin plots", name: "fig14" },
          })],
        }),
        para(
          "Violin plots of Dragon\u2019s query coverage (best-hit containment) broken down by benchmark tier. Black line: mean; red line: median. Tiers 2 (cross-species), 3 (AMR), and 6 (plasmid) show high completeness concentrated near 1.0. Tier 5 (short reads) shows lower coverage, consistent with Dragon\u2019s design for longer queries."
        ),

        // FIGURE 15: Signal Discretization
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 15: Signal-Level Discretization Pipeline")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig15_signal_discretization.png"),
            transformation: { width: 560, height: 370 },
            altText: { title: "Figure 15", description: "Signal discretization pipeline", name: "fig15" },
          })],
        }),
        para(
          "Signal-level search pipeline. (a) Expected pA current from the R10.4.1 pore model for a 100 bp genomic region. (b) Simulated nanopore signal with \u03C3=3.0 pA Gaussian noise added. (c) After median-MAD normalization and discretization into a 16-level alphabet, the signal is suitable for FM-index backward search. Each color represents one of 16 discrete levels."
        ),

        // FIGURE 16: Signal Search Accuracy
        new Paragraph({
          spacing: { before: 300, after: 100 },
          alignment: AlignmentType.CENTER,
          children: [boldRun("Figure 16: Signal Search Accuracy vs. Noise Level")],
        }),
        new Paragraph({
          spacing: { after: 100 },
          alignment: AlignmentType.CENTER,
          children: [new ImageRun({
            type: "png",
            data: fs.readFileSync("/Users/lshlt19/GitHub/Dragon/manuscript/figures/fig16_signal_search_accuracy.png"),
            transformation: { width: 560, height: 260 },
            altText: { title: "Figure 16", description: "Signal search accuracy vs noise", name: "fig16" },
          })],
        }),
        para(
          "Signal-level search performance across noise levels. (a) Species-level sensitivity as a function of signal noise (\u03C3 pA): top-1 match (red diamonds), top-5 match (teal squares), and any-hit-correct-species (blue circles). At clean signal (\u03C3=0), 47% of reads are correctly classified at the species level in the top hit, rising to 68% within the top 5 and 100% within all returned hits. Performance degrades gracefully with noise. (b) Score distributions shift leftward with increasing noise, but remain informative for species-level classification."
        ),
      ],
    },
  ],
});

// Generate the file
const OUTPUT = "/Users/lshlt19/GitHub/Dragon/manuscript/dragon_manuscript.docx";
Packer.toBuffer(doc).then((buffer) => {
  fs.writeFileSync(OUTPUT, buffer);
  console.log(`DOCX written to: ${OUTPUT}`);
  console.log(`File size: ${(buffer.length / 1024).toFixed(1)} KB`);
});
