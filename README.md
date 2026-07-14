# Bulk RNA-seq Analysis of COVID-19 Leukocytes
### Transcriptomic Profiling of COVID-19 vs Non-COVID Hospitalized Patients

[![Pipeline](https://img.shields.io/badge/Pipeline-FastQC%20→%20HISAT2%20→%20DESeq2-blue)]()
[![Genome](https://img.shields.io/badge/Genome-GRCh38%20(Ensembl%20109)-green)]()
[![Dataset](https://img.shields.io/badge/Dataset-GSE157103-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103)
[![Language](https://img.shields.io/badge/Language-Bash%20%7C%20R-red)]()

---

## Overview

This project implements a complete **bulk RNA-sequencing (RNA-seq) pipeline** from raw sequencing reads to differential gene expression results and biological interpretation. The analysis was carried out entirely independently on Ubuntu/WSL using standard bioinformatics tools.

**Biological Question:** Which genes are significantly more or less expressed in the blood of COVID-19 positive patients compared to non-COVID-19 hospitalized patients?

**Key Finding:** A robust interferon-stimulated gene (ISG) signature was identified in COVID-19 patient leukocytes, with IFI27 showing the highest upregulation (log2FC = 11.15, padj = 1.84×10⁻¹⁹). Haemoglobin genes HBA1 and HBA2 were the most significantly downregulated genes (padj = 1.96×10⁻⁴² and 9.57×10⁻²⁴ respectively), consistent with published reports of COVID-19's impact on haematological function.

---

## Dataset

| Parameter | Details |
|---|---|
| **GEO Accession** | [GSE157103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103) |
| **SRA BioProject** | PRJNA660067 |
| **Reference Paper** | Overmyer et al. (2021) *Cell Systems* 12(1):23-40 |
| **Tissue** | Leukocytes (peripheral blood white blood cells) |
| **Sequencing Platform** | Illumina NovaSeq 6000 |
| **Library Layout** | Paired-end (2 × 102 bp) |
| **Samples Used** | 6 (3 COVID-19 positive, 3 non-COVID hospitalized) |
| **Reference Genome** | GRCh38, Ensembl release 109 |

**Samples:**

| SRR Accession | Condition | Library Size |
|---|---|---|
| SRR12544419 | COVID-19 Positive | ~40M reads |
| SRR12544420 | COVID-19 Positive | ~39M reads |
| SRR12544421 | COVID-19 Positive | ~23M reads |
| SRR12544527 | Non-COVID Hospitalized | ~39M reads |
| SRR12544528 | Non-COVID Hospitalized | ~38M reads |
| SRR12544529 | Non-COVID Hospitalized | ~43M reads |

> **Note:** The non-COVID group are other hospitalized patients — not healthy controls. This is a COVID-19 vs other-illness comparison, making detected signals specific to SARS-CoV-2 infection rather than generalized illness.

---

## Pipeline Overview

```
Raw FASTQ files (SRA)
        │
        ▼
┌─────────────────┐
│  Quality Control │  FastQC
│  (raw reads)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    Trimming     │  Trimmomatic PE
│  (adapters +   │  TruSeq3-PE adapters
│  low quality)   │  SLIDINGWINDOW:4:15
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Quality Control │  FastQC
│  (trimmed reads) │  (post-trim verification)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Alignment     │  HISAT2 (splice-aware)
│  GRCh38 genome  │  → SAMtools sort + index
└────────┬────────┘
         │  BAM files
         ▼
┌─────────────────┐
│ Quantification  │  featureCounts
│  (reads/gene)   │  Ensembl 109 GTF
└────────┬────────┘
         │  Count matrix
         ▼
┌─────────────────┐
│  Differential   │  DESeq2 (R)
│  Expression     │  Negative binomial model
│  Analysis       │  BH FDR correction
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Visualizations  │  Volcano plot
│                 │  Heatmap
│                 │  PCA plot
└─────────────────┘
```

---

## Repository Structure

```
covid19-rnaseq-pipeline/
│
├── README.md                          # This file
│
├── scripts/
│   ├── pipeline.sh                    # Full bash pipeline (FastQC → featureCounts)
│   └── deseq2_analysis.R              # DESeq2 differential expression + visualizations
│
├── results/
│   ├── fastqc_reports/
│   │   ├── raw/                       # Pre-trimming FastQC HTML reports
│   │   └── trimmed/                   # Post-trimming FastQC HTML reports
│   ├── counts/
│   │   ├── raw_counts.txt             # featureCounts output (count matrix)
│   │   ├── raw_counts.txt.summary     # featureCounts assignment summary
│   │   ├── deseq2_results.csv         # Full DESeq2 results (all genes)
│   │   ├── volcano_plot.png           # Volcano plot
│   │   ├── heatmap.png                # Top 30 DE genes heatmap
│   │   └── pca_plot.png               # Sample clustering PCA
│   └── alignment_stats/               # HISAT2 alignment rate per sample
│
└── metadata/
    └── SraRunTable.csv                # Sample metadata from NCBI SRA
```

> **Note on raw data:** Raw FASTQ files (~15-25 GB total) and intermediate BAM files are not included in this repository due to size. Download instructions are provided below.

---

## Requirements

### System
- Ubuntu 20.04 LTS or later (or WSL2 on Windows)
- Minimum 50 GB free disk space
- Minimum 16 GB RAM (recommended for HISAT2 index building)

### Software (install via conda)

```bash
conda create -n rnaseq_pipeline python=3.9 -y
conda activate rnaseq_pipeline

conda install -c bioconda -c conda-forge \
  fastqc trimmomatic hisat2 samtools subread sra-tools -y
```

### R packages

```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages(c("ggplot2", "pheatmap", "ggrepel", "RColorBrewer"))
```

---

## Reproducing This Analysis

### Step 1 — Download raw data from SRA

```bash
# Create project structure
mkdir -p ~/rnaseq_project/{raw_data,fastqc_reports/{raw,trimmed},trimmed,aligned,counts,reference,scripts}
cd ~/rnaseq_project

# Download samples
prefetch SRR12544419 SRR12544420 SRR12544421 SRR12544527 SRR12544528 SRR12544529

# Convert to FASTQ (paired-end)
while read accession; do
    fastq-dump --split-files --gzip --outdir raw_data/ "$accession"
done < scripts/SRR_Acc_List.txt
```

### Step 2 — Download reference genome

```bash
cd ~/rnaseq_project/reference

# Human genome FASTA (GRCh38, Ensembl 109)
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Gene annotation GTF
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

# Decompress
gunzip *.gz

# Build HISAT2 index (~60-90 minutes)
hisat2-build -p 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38_index
```

### Step 3 — Run the upstream pipeline

```bash
conda activate rnaseq_pipeline
bash ~/rnaseq_project/scripts/pipeline.sh
```

This runs FastQC → Trimmomatic → FastQC (post-trim) → HISAT2 → SAMtools → featureCounts for all 6 samples automatically.

### Step 4 — Run DESeq2 differential expression analysis

Open `scripts/deseq2_analysis.R` in RStudio or Posit Cloud and run the full script. This produces the count matrix, volcano plot, heatmap, PCA plot, and full results CSV.

---

## Key Results

### Alignment Rates

All samples achieved 70-73% alignment to GRCh38 — within normal range for leukocyte RNA-seq data.

### DESeq2 Summary

After filtering genes with < 10 counts in < 3 samples:
- **~16,000 genes tested**
- Significant at padj < 0.05 and |log2FC| > 1:
  - **Upregulated in COVID-19:** hundreds of genes
  - **Downregulated in COVID-19:** hundreds of genes

### Top Upregulated Genes in COVID-19

| Gene | log2FC | Fold Change | padj | Function |
|---|---|---|---|---|
| IFI27 | 11.15 | ~2,236x | 1.84×10⁻¹⁹ | Interferon-alpha inducible — antiviral defence |
| OTOF | 8.86 | ~462x | 1.61×10⁻⁹ | Immune signalling, elevated in severe COVID |
| IFI44L | 8.04 | ~264x | 1.52×10⁻⁸ | Interferon-induced antiviral protein |
| USP18 | 8.01 | ~256x | 4.57×10⁻³ | Regulates interferon signalling pathway |
| SIGLEC1 | 7.79 | ~220x | 1.64×10⁻⁷ | Macrophage activation marker |
| RSAD2 (Viperin) | 7.20 | ~147x | 4.78×10⁻⁹ | Blocks viral replication directly |
| OAS3 | 6.47 | ~89x | 8.67×10⁻¹⁴ | Destroys viral RNA |
| IFIT1 | 5.81 | ~56x | 2.56×10⁻¹² | Blocks viral mRNA translation |

### Top Downregulated Genes in COVID-19

| Gene | log2FC | Fold Change | padj | Function |
|---|---|---|---|---|
| HBA1 | -5.80 | ~56x lower | 1.96×10⁻⁴² | Haemoglobin alpha chain — oxygen transport |
| HBA2 | -6.34 | ~81x lower | 9.57×10⁻²⁴ | Haemoglobin alpha chain — oxygen transport |
| HPGD | -6.39 | ~84x lower | 8.10×10⁻¹⁴ | Degrades prostaglandins — disrupted in hyperinflammation |
| KLF14 | -7.57 | ~191x lower | 2.36×10⁻⁶ | Transcriptional regulator |

---

## Visualizations

### PCA Plot
The PCA plot shows PC1 (70% of variance) cleanly separating COVID-19 from non-COVID samples — confirming the COVID vs non-COVID difference is the dominant biological signal in the dataset.

![PCA Plot](results/counts/pca_plot.png)

### Volcano Plot
Volcano plot showing hundreds of significantly differentially expressed genes. Red = upregulated in COVID, Blue = downregulated, Grey = not significant. Dashed lines mark padj = 0.05 and |log2FC| = 1 thresholds.

![Volcano Plot](results/counts/volcano_plot.png)

### Heatmap
Heatmap of top 30 DE genes (15 up, 15 down) with hierarchical clustering. Clear separation of COVID and non-COVID sample groups.

![Heatmap](results/counts/heatmap.png)

---

## Biological Interpretation

### The Interferon-Stimulated Gene (ISG) Signature

All top upregulated genes belong to a single coordinated biological program — the **Type I Interferon antiviral response**:

```
SARS-CoV-2 infection
        │
        ▼
Viral dsRNA detected by pattern recognition receptors (TLR3, RIG-I)
        │
        ▼
Type I Interferons produced (IFN-α/β)
        │
        ▼
JAK-STAT signalling → ISGF3 complex binds ISREs
        │
        ▼
IFI27, IFI44L, RSAD2, OAS3, IFIT1, SIGLEC1, USP18
all massively upregulated
        │
        ▼
Antiviral state: block replication, destroy viral RNA,
kill infected cells, recruit immune cells
```

Each gene plays a specific antiviral role:
- **RSAD2 (Viperin)** — produces ddhCTP, a chain terminator for viral RNA polymerases; disrupts membrane composition to block viral budding
- **OAS3** — produces 2'-5'-oligoadenylates that activate RNase L, degrading viral RNA
- **IFIT1** — binds and sequesters viral RNA cap structures, blocking translation
- **IFI27** — promotes apoptosis of infected cells to prevent further viral spread
- **USP18** — negative regulator of interferon signalling; its induction prevents damaging over-activation of the immune response

### The Haematological Disruption

The dramatic downregulation of **HBA1** (padj = 1.96×10⁻⁴²) and **HBA2** (padj = 9.57×10⁻²⁴) provides molecular evidence for a well-documented clinical observation — COVID-19 impairs haemoglobin production and oxygen delivery through mechanisms beyond respiratory compromise alone. SARS-CoV-2 directly attacks erythroid precursor cells, disrupting haemopoiesis. The downregulation of **HPGD** (prostaglandin degradation enzyme) further supports the hyperinflammatory state characteristic of severe COVID-19.

### Validation Against Published Study

Results were compared against the full 126-sample published dataset. Directional consistency was **100%** across all examined genes — every gene upregulated in this analysis was upregulated in the full study, and vice versa. Fold change magnitudes were larger in this analysis due to the small sample size (3 per group vs 100 per group), which is expected and well-understood statistically — DESeq2's shrinkage estimation partially corrects for this.

---

## Limitations

1. **Small sample size (n=3 per group):** Inflates fold change estimates relative to the full dataset. DESeq2's shrinkage estimation mitigates but does not eliminate this effect.

2. **Control group:** Non-COVID samples are other hospitalized patients, not healthy volunteers. Any detected signal reflects COVID-19-specific biology over and above general critical illness.

3. **No gene name annotation in pipeline output:** Ensembl IDs (ENSG...) are used throughout. Conversion to HGNC symbols via biomaRt required for final interpretation.

4. **Sequencing depth variation:** SRR12544421 has ~23M reads vs ~40M for other COVID samples, creating secondary variation visible in PCA (PC2). DESeq2 normalization corrects for this in statistical testing.

---

## Tools and Versions

| Tool | Version | Purpose |
|---|---|---|
| FastQC | v0.11+ | Raw and post-trim quality control |
| Trimmomatic | v0.39 | Adapter and quality trimming |
| HISAT2 | v2.2.1 | Splice-aware alignment to GRCh38 |
| SAMtools | v1.17 | SAM/BAM conversion, sorting, indexing |
| featureCounts (Subread) | v2.0.1 | Read counting per gene |
| DESeq2 | v1.38+ | Differential expression analysis |
| ggplot2 | v3.4+ | Volcano and PCA plots |
| pheatmap | v1.0+ | Expression heatmap |
| R | v4.2+ | Statistical analysis environment |
| Ubuntu/WSL | 22.04 LTS | Computing environment |

---

## References

1. Overmyer, K.A. et al. (2021). Large-Scale Multi-omic Analysis of COVID-19 Severity. *Cell Systems*, 12(1):23-40. https://doi.org/10.1016/j.cels.2020.10.003

2. Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550. https://doi.org/10.1186/s13059-014-0550-8

3. Kim, D., Paggi, J.M., Park, C., Bennett, C. & Salzberg, S.L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37:907-915.

4. Liao, Y., Smyth, G.K. & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7):923-930.

5. Benjamini, Y. & Hochberg, Y. (1995). Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society B*, 57(1):289-300.

---

## Author

**Varshini**
B.E. Biotechnology — Sapthagiri College of Engineering, Bengaluru
MSc Computational Biology & Bioinformatics (Incoming, September 2026) — Ramaiah University of Applied Sciences, Bengaluru

*This project was completed independently as a pre-MSc computational biology portfolio project.*

---

*Dataset: GSE157103 | Reference genome: GRCh38 Ensembl 109 | Analysis: Ubuntu 22.04 WSL + R/Posit Cloud*
