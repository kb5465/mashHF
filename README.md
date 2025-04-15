# mashHF: Multi-Trait Rare Variant Analysis for Heart Failure

This repository contains code for the thesis:

**Population Genetics for Heart Failure: A Multi-Trait Analysis of Rare, Protein-Coding Variants**  
Kiran J. Biddinger | Princeton University, 2025

The project uses exome-sequencing data from the UK Biobank to identify loss-of-function (LoF) variants associated with heart failure (HF), leveraging both single-trait burden tests and a novel multi-trait implementation of the mash framework.

---

## Repository Structure and Purpose

### 1. Phenotype Preparation and Summary Statistics
Location: `UKBB_data/`

- `01_curate_UKBB_phenotypes_with_LV_proteomics.R`  
  Constructs imaging and protein biomarker phenotypes.
- `02_curate_UKBB_phenotypes_phewas_traits.R`  
  Curates broader set of phenotypes for PheWAS scan.
- `03_generate_table1_summary_by_subset.R`  
  Generates descriptive statistics for table 1.

---

### 2. Single-Trait Association Testing (REGENIE)
Location: `single_trait_REGENIE/`

- `run_regenie_step1_QT.sh`, `run_regenie_step2_QT_chunked.sh`  
  Run REGENIE steps 1 and 2 for continuous traits.  
- `run_regenie_step1_BT.sh`, `run_regenie_step2_BT_chunked.sh`  
  Run REGENIE steps for binary traits (HF, DCM, HCM).  
- `filter_setinclusion_file_by_gene_and_tissue.R`  
  Filters REGENIE set files for selected tissues and genes.  
- `merge_canonical_LOF_results.R`, `extract_ncarriers_by_subset.R`  
  Post-process and summarize burden test results.

---

### 3. Multi-Trait Analysis Using `mashr`
Location: `multi_trait_mash/`

- `run_mashr_canonical_LOF_analysis.R`  
  Applies the mash empirical Bayes model to gene-trait effect estimates.
- `tables_and_figures/`  
  Contains visualization scripts:
  - `plot_genomic_inflation_lambda.R`
  - `plot_mash_nhits_per_trait.R`
  - `plot_significant_gene_trait_heatmap.R`
  - `visualize_mashr_zscores_and_qqplots.R`
  - `visualize_mash_covariance_structures.R`

---

### 4. Simulation Experiments for mash Validation
Location: `simulations/`

- `mashHF_simulations_final.R`  
  Main simulation used in the thesis to benchmark mash.
- `mashHF_simulations_noise.R`, `bias.R`, `domeffects.R`  
  Additional scenarios varying effect types and error structure.

---

### 5. GTEx Expression Enrichment
Location: `expression_GTEX/`

- `Thesis_GTEX_OT.ipynb`  
  Jupyter notebook for processing GTEx expression data and comparing expression across tissues.

---

### 6. eQTL Integration and Cauchy-Based Gene-Trait Aggregation
Location: `eQTL_integration/`

- `extract_eqtls_for_mashHF_genes.R`  
  Extracts GTEx-based eQTLs for nominated genes.
- `subset_gwas_for_mashHF_eqtls.R`  
  Matches eQTLs to GWAS results.
- `eqtl_gwas_cauchy_heatmap.R`  
  Aggregates p-values and visualizes tissue-specific signal.
- `intersect_rsids_across_gwas.R`  
  Identifies shared variants across summary datasets.

---

### 7. PheWAS Scan of Nominated Genes
Location: `pheWAS_scan/`

- `merge_singletrait_results_phewas.R`  
  Combines PheWAS output across traits.
- `canonical_LOF_gene_trait_heatmap_phewas.R`  
  Plots gene-phenotype association heatmap.
- `extract_phewas_summary_stats.sh`  
  Extracts relevant REGENIE results by gene.

---

## Thesis Reference

For background, methods, and interpretation, see:  
**Biddinger KJ. *Population Genetics for Heart Failure*. Princeton University, 2025.**

---