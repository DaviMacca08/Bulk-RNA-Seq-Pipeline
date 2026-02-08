Bulk RNA-seq – Upstream Analysis
BioProject: PRJNA1400483
Data source: NCBI SRA
Sequencing type: Illumina paired-end RNA-seq

============================================================
PROJECT OVERVIEW AND SCIENTIFIC SCOPE
============================================================

This study aimed to elucidate the spatiotemporal expression pattern of 
O-GlcNAc transferase (OGT) and its catalyzed O-GlcNAc modification in ovarian 
tissues, and to investigate their specific functions in ovarian granulosa cells.

To this end, high-throughput RNA sequencing (RNA-seq) was performed to analyze 
genome-wide transcriptional changes following OGT knockdown in the human ovarian 
granulosa cell line KGN.

Two specific siRNAs targeting OGT (siOGT1 and siOGT2) and a negative control 
siRNA (siNC) were designed. Each experimental group included two biological 
replicates, resulting in a total of six transcriptomic libraries.

Sequencing was carried out on the Illumina NovaSeq 6000 platform using paired-end 
150 bp reads.

By comparing siOGT-treated groups with the siNC control group:
- 300 differentially expressed genes (DEGs) were identified in the siOGT1 group
- 487 DEGs were identified in the siOGT2 group
- 176 DEGs overlapped between the two knockdown conditions

Pathway enrichment analyses indicated that OGT knockdown may affect biological 
processes related to immune defense response, cellular responses to external 
stimuli, and cytokine-mediated signaling pathways.

Overall, this dataset provides a comprehensive transcriptomic profile of human 
ovarian granulosa cells under OGT interference and represents a valuable 
resource for investigating the regulatory mechanisms of the OGT/O-GlcNAc axis 
in ovarian physiology and pathology.

============================================================
DATASET DESCRIPTION
============================================================

BioProject: PRJNA1400483
SRA accessions (paired-end):
- SRR36749378
- SRR36749379
- SRR36749380
- SRR36749381
- SRR36749382
- SRR36749383

Library layout: Paired-end
Sequencing platform: Illumina NovaSeq 6000
Read length: 2 x 150 bp

============================================================
UPSTREAM ANALYSIS WORKFLOW
============================================================

1. Data retrieval
   Raw sequencing data were downloaded from NCBI SRA using SRA Toolkit.

2. Initial quality control
   FastQC was applied to raw FASTQ files, followed by MultiQC aggregation.

3. Adapter and quality trimming
   Trimmomatic was used in paired-end mode to remove Illumina adapters and
   low-quality bases. Sliding window trimming and minimum length filtering
   were applied.

   NOTE:
   Poly(A) tails were intentionally NOT removed. This choice was made to:
   - Preserve biologically relevant transcript features
   - Remain consistent with the original experimental design
   - Avoid introducing undocumented preprocessing steps

4. Post-trimming quality control
   A second MultiQC report was generated to confirm successful trimming and
   quality improvement while retaining poly(A) features.

5. Reference genome preparation
   The reference genome and corresponding annotation were downloaded and
   indexed for alignment.

6. Alignment
   Trimmed reads were aligned to the reference genome using HISAT2.

7. BAM processing
   SAM files were converted to sorted BAM files using Samtools.

8. Gene-level quantification
   featureCounts was used to generate a gene-level count matrix based on the
   selected annotation.

============================================================
REPOSITORY STRUCTURE (UPSTREAM)
============================================================

Bulk-RNAseq/upstream/
├── metadata/
│   └── GSE316021_series_matrix.txt
├── report/
│   ├── Hisat2/
│   ├── Trimmomatic/
│   ├── featureCounts/
│   └── multiQC/
├── script/
│
└── README.txt
============================================================
FINAL OUTPUT OF UPSTREAM ANALYSIS
============================================================

The final product of the upstream analysis is a gene-level raw count matrix
(counts.txt), which serves as the input for downstream analyses.

============================================================

============================================================
DOWNSTREAM ANALYSIS WORKFLOW
============================================================

Downstream analyses were conducted in R using the DESeq2 framework and 
associated Bioconductor packages. The workflow includes:

1. **DESeq2 object construction**
   - The raw count matrix from upstream analysis was imported into R.
   - A DESeqDataSet object was created with the experimental design (~ condition) 
     and factor levels set according to the control group (siNC).

2. **Exploratory Data Analysis (EDA)**
   - Variance-stabilizing transformation (VST) was applied to normalize counts.
   - Sample QC included boxplots, hierarchical clustering, and PCA to evaluate 
     global expression patterns and sample relationships.
   - Top 500 most variable genes were used to highlight sample-specific variance.

3. **Differential Expression Analysis**
   - Pairwise comparisons were performed between each siOGT knockdown group 
     and the siNC control.
   - Adjusted p-value (<0.05) and |log2 fold change| (>0.5) thresholds were applied 
     to identify significantly differentially expressed genes (DEGs).
   - MA plots and volcano plots visualize gene-level changes.

4. **Gene Ontology (GO) and Pathway Enrichment**
   - DEGs common to both knockdown conditions were analyzed for enrichment in 
     biological processes (GO BP), KEGG pathways, and Reactome pathways.
   - Significant terms were filtered by adjusted p-value and visualized with barplots, 
     dotplots, and cnetplots for network representation.

============================================================
REPOSITORY STRUCTURE (DOWNSTREAM)
============================================================

Bulk-RNAseq/downstream/
├── R-scripts/
│   ├── 01_DESeq2_object.r
│   ├── 02_Exploratory_Data_Analysis.r
│   └── 03_Differentially_Expressed_Genes.r
├── counts.txt
└── images/
    ├── 01_QC/
    ├── 02_EDA/
    └── 03_DEG_GO_Enrichment/

============================================================
REPRODUCIBILITY AND RESULTS
============================================================

- All R scripts are fully reproducible using the provided count matrix.
- The downstream pipeline successfully replicates the original study’s results:
  - siOGT1 knockdown: 333 DEGs
  - siOGT2 knockdown: 556 DEGs
  - Overlapping DEGs between knockdowns: 191
- QC metrics, PCA, heatmaps, volcano plots, and GO/KEGG/Reactome enrichments 
  are included as PDF files in the `images/` folder.

============================================================
USAGE AND ACCESS
============================================================

- Scripts can be executed sequentially in R to reproduce the complete analysis.
- For private access to the raw count matrix or detailed logs, please contact 
  the repository owner.

============================================================
REFERENCES
============================================================

- PRJNA1400483 – NCBI BioProject
- Love MI, Huber W, Anders S. "Moderated estimation of fold change and dispersion 
  for RNA-seq data with DESeq2." Genome Biol. 2014;15(12):550.
- Bioconductor RNA-seq workflow: https://bioconductor.org/help/workflows/rnaseqGene/
- Yu G, Wang LG, Han Y, He QY. "clusterProfiler: an R package for comparing 
  biological themes among gene clusters." OMICS. 2012;16(5):284-7.
