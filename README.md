# Single-cell transcriptomic GRN analysis uncovers alteration of colonic mucosal immune cells in Nr4a1 KO mice

This repository contains the analysis code associated with the publication:

**"Single-cell transcriptomic gene regulatory network analysis uncovers the alteration of colonic mucosal immune cells in Nr4a1 knockout mice"**
*iScience (Cell Press) — manuscript ISCIENCE-D-25-13626R2*

Chapkin Lab · Texas A&M University

---

## Data availability

Raw sequencing data (h5 files) and processed count matrices are deposited at the Gene Expression Omnibus (GEO): **[GSE######]** *(accession to be updated upon acceptance).*
The final annotated Seurat RDS object (1.4 GB) is available from the corresponding author upon reasonable request.

---

## Repository structure

```
.
├── matlab_annotation/          # MATLAB scripts: raw h5 processing and initial cell annotation
├── seurat_R_processing/        # Core R analysis pipeline — run scripts in numbered order
├── extra_analysis/             # Supporting and validation analyses
├── exon_analysis_nr4a1.ipynb   # Python notebook: Nr4a1 exon detection across cell types
├── 2023_nr4a1_colon_meta_info.xlsx   # Sample-level metadata
├── Nr4a1_Exon_Detailed_Analysis.xlsx # Exon-level detection output
├── Nr4a1_Exon_Summary_Table.csv      # Exon summary table
├── metadata_subann4.xlsx             # Cell subannotation metadata
└── copy_sc_data.sh                   # Shell script for staging raw data files
```

---

## matlab_annotation/

MATLAB scripts for initial processing of raw 10x h5 files and cell type annotation. Run before the R pipeline.

| File | Description |
|------|-------------|
| `s1_processh5files.m` | Loads and processes raw 10x Genomics h5 files; outputs MAT files per sample |
| `s2_mergematfiles.m` | Merges per-sample MAT files into a single data matrix for import into R |
| `nr4a1_wholebody_ko_markers.txt` | Cell type marker genes used for cluster annotation |

> **Note on CCAT/SCENT potency scoring:** Differentiation potency scores were computed using the MATLAB SCENT toolbox ([Teschendorff et al.; github.com/aet21/SCENT](https://github.com/aet21/SCENT)) applied interactively via the MATLAB GUI. A scripted version of this step is not available. The exported potency scores are loaded directly in `seurat_R_processing/07_nr4a1_whole_body_KO_colonocytes_analysis_and_plots.R`. The SCENT toolbox is publicly available and the parameters used are described in the Methods section of the manuscript.

---

## seurat_R_processing/

Core numbered R analysis pipeline. Scripts are intended to be run in order (01 → 08). Each script saves output objects used by subsequent scripts.

| Script | Description |
|--------|-------------|
| `01_nr4a1_whole_body_KO_analysis_extraction.R` | Quality control, doublet removal, normalization, dimensionality reduction, and initial clustering of WT and Nr4a1 KO cells |
| `02_nr4a1_whole_body_KO_differential_expression.R` | Differential expression (DE) analysis between KO and WT for each annotated cell type |
| `03_nr4a1_whole_body_KO_differential_variability.R` | Differential variability (DV) analysis using scran; identifies genes with altered transcriptional noise in KO |
| `03_5_nr4a1_whole_body_KO_dv_de_enrichr.R` | GO enrichment analysis on DV+DE overlapping gene sets using Enrichr; generates Supplementary Table ST 3 |
| `04_nr4a1_whole_body_KO_cellchat.R` | Cell-cell communication analysis using CellChat 2.0; identifies altered ligand-receptor signaling in KO |
| `05_nr4a1_whole_body_KO_sctenifold_net_results.R` | Gene regulatory network (GRN) inference and differential regulation (DR) analysis using scTenifoldNet |
| `06_nr4a1_whole_body_KO_Tcells_analysis_and_plots.R` | Sub-clustering, analysis, and visualization of T cell populations |
| `07_nr4a1_whole_body_KO_colonocytes_analysis_and_plots.R` | Colonocyte sub-clustering, AUCell scoring, cell cycle scoring (Seurat CellCycleScoring), and loading of CCAT/SCENT potency scores exported from MATLAB |
| `08_nr4a1_whole_body_KO_immune_cells_proportions.R` | Quantification and visualization of immune cell type proportions across KO and WT conditions |
| `Genesets_cell_scores.xlsx` | Gene sets used for AUCell scoring and cell state annotation |

---

## extra_analysis/

Supporting analyses and validation scripts not part of the main numbered pipeline.

| File | Description |
|------|-------------|
| `nr4a1_whole_body_KO_analysis.R` | Original master analysis script; precursor to the numbered pipeline; retained for reference |
| `tcells_subannotation.R` | T cell subannotation script; assigns T cell subtypes used as input to script 06 |
| `plot_colonic_organoids.R` | Plots colonic organoid-forming efficiency data comparing KO and WT |
| `plot_qpcr_nr4a1_2.R` | Plots qPCR validation data for Nr4a1 expression |
| `nr4a1_qpcr.csv` | Raw qPCR data for Nr4a1 expression validation |

---

## Dependencies

**R packages:** Seurat (v4+), scran, CellChat (v2), scTenifoldNet, AUCell, Enrichr/enrichR, ggplot2, dplyr

**Python:** Jupyter, scanpy (for exon analysis notebook)

**MATLAB:** SCENT toolbox ([github.com/aet21/SCENT](https://github.com/aet21/SCENT)) — required only for potency scoring (see CCAT/SCENT note above)

---

## Citation

If you use this code, please cite:

> Romero et al. *Single-cell transcriptomic gene regulatory network analysis uncovers the alteration of colonic mucosal immune cells in Nr4a1 knockout mice.* iScience (2025). *(DOI to be updated upon acceptance)*

---

## Contact

Selim S. Romero — ssromerogon@tamu.edu
Chapkin Lab — [chapkinlab.tamu.edu](https://chapkinlab.tamu.edu)
