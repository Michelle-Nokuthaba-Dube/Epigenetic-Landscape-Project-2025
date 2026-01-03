# Epigenetic Landscape of Stress-Related Genes in Mice Brain Tissue
# Abstract
Epigenetic regulation is central to stress-responsive transcriptional programs in the brain. This study characterizes promoter-associated histone modification patterns in mouse prefrontal cortex using ChIP-seq data for H3K4me3, aiming to identify regulatory regions and biological pathways associated with stress-related gene regulation.
 
# Experimental Context and Dataset
ChIP-seq data were obtained from Reshetnikov et al. (2021), which investigated the molecular effects of early-life stress and chronic social defeat stress (SDS) in mice. Data were retrieved from the NCBI Sequence Read Archive (SRA) under BioProject PRJNA610193.
The original experiment used male C57BL/6 mice subjected to:
- Maternal separation (MS): 3 hours daily from postnatal day 2–14
- Chronic social defeat stress (SDS): 15 consecutive days of 10-minute exposure to an aggressive CD-1 mouse in adulthood

Control animals were reared under identical housing conditions without stress exposure. Prefrontal cortex tissue was collected for ChIP-seq profiling of the H3K4me3 histone mark, a promoter-associated modification linked to active transcription.
For this bioinformatics workflow, two samples were analyzed:
- One mouse exposed to MS in early life and SDS in adulthood
- One unstressed control mouse (NS)
Due to the use of single samples per condition, all downstream analyses are exploratory and descriptive and do not support statistical inference.

# Bioinformatics Workflow
All analyses were performed in R (v4.5.1) and complemented with command-line bioinformatics tools. The pipeline consisted of read preprocessing, alignment, peak calling, peak annotation, functional enrichment analysis, and visualization. All scripts and outputs are available in this repository.

# Methodology
# 1. Read Retrieval and Preprocessing
Raw FASTQ files were downloaded from SRA using SRA Toolkit (v3.2.1). Quality control and adapter trimming were performed using FASTQ (v1.0.1) with default parameters to remove low-quality bases and sequencing adapters.

# 2a. Alignment
Trimmed reads were aligned to the mouse reference genome (GRCm38/mm10) using Bowtie2 (v2.5.4) with default settings optimized for ChIP-seq data. Alignment files were converted to sorted BAM format using SAMtools (v1.22.1). Mapping efficiency and alignment quality were assessed using samtools flagstat.

# 2b. Peak Calling
H3K4me3-enriched regions were identified using MACS2 via the Galaxy platform (v2.2.9.1 + Galaxy 0). Peak calling was performed in broad peak mode, appropriate for promoter-associated histone marks, with:
- q-value threshold: 0.05
- Default shift and extension parameters validated for H3K4me3

# 3. Peak Annotation
Peak annotation was performed using ChIPseeker (v1.46.1) in R. Peaks were annotated against the TxDb.Mmusculus.UCSC.mm10.knownGene reference and classified into promoters, exons, introns, and intergenic regions.
- Broad peaks were annotated using a ±3 kb window around transcription start sites (TSS) to capture extended regulatory regions
- Filtered promoter-associated peaks were restricted to a ±1 kb TSS window to focus on core promoter enrichment
Associated genes were extracted as Entrez Gene IDs for downstream analysis.

# 4. Functional Enrichment Analysis
Functional characterization was conducted using clusterProfiler (v4.18.2) and org.Mm.eg.db (v3.22.0). Gene Ontology Biological Process (GO:BP) enrichment analysis was performed using enrichGO with:
Multiple testing correction: Benjamini–Hochberg
- p-value cutoff: 0.05
- q-value cutoff: 0.2
Enrichment analyses were performed separately for genes associated with broad peaks and filtered promoter-associated peaks. Given the lack of biological replication, results are interpreted as descriptive trends. Stress-related pathways were identified through keyword-based filtering of GO terms containing the word “stress”.

# 5. Visualization
Data visualization was performed using ggplot2 (v4.0.1). Analyses included:
- Genomic feature distribution of H3K4me3 peaks
- Distance-to-TSS plots to assess promoter enrichment
- Bar plots comparing stress-related GO terms between broad and filtered peak sets

# Key Findings
- H3K4me3 enrichment was predominantly localized to promoter regions of stress-related genes
- Functional enrichment highlighted pathways associated with transcriptional regulation, neuronal function, and stress response
- Differences between broad and core promoter-associated peaks suggest layered epigenetic regulation

# Limitations
- Single-sample comparisons limit statistical inference
- Findings are exploratory and intended to highlight biological patterns rather than differential significance

# Repository Structure
├── data/          # Input and processed data
├── scripts/       # R scripts and command-line workflows
├── results/       # Peak files and enrichment outputs
├── figures/       # Publication-ready visualizations
└── README.md
