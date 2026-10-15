# Survival-GWAS

# Survival GWAS of Substance Use Initiation (ABCD) — REGENIE (T2E) + Meta-analysis

Mengman Wei\* and Qian Peng\*

This repository contains scripts to run **time-to-event (survival) GWAS** of substance use initiation in the ABCD Study using **REGENIE’s `--t2e` workflow**, followed by a **strand-aware meta-analysis** across ancestries (EUR/AFR/HISP) using fixed- and random-effects models.

> **Core scripts in this repo**
- `survivalGWAS.slurm` — per-ancestry survival GWAS (example: EUR), SLURM array over 4 traits  
- `meta_sud.R` — strand-aware, palindromic-safe meta-analysis across ancestries (EUR/AFR/HISP), chr-by-chr streaming

---

## 1) Overview of the workflow

### Step A — Per-ancestry survival GWAS (REGENIE `--t2e`)
For each ancestry group (e.g., EUR), and each trait:
1. Standardize chromosome coding in the genotype bedset (handles `chrX`, `chrY`, `chrMT`, etc.).
2. Build a **Step-1 SNP set** on **unrelated samples** via PLINK QC + pruning, excluding long-LD regions.
3. Run **REGENIE step 1** time-to-event model (`--t2e`) using `TIME` + `EVENT`.
   - Includes an auto-retry path if LOCO predictions are not produced.
4. Run **REGENIE step 2** GWAS scan with LOCO predictions and covariates.

### Step B — Strand-aware meta-analysis across ancestries
- Loads GWAS outputs for EUR/AFR/HISP **per chromosome** (1–22).
- Harmonizes alleles to a **reference ancestry** (default `EUR`) with:
  - allele match / swap
  - strand flip / strand+swap
  - palindromic filtering using allele frequency (safer defaults)
- Runs:
  - **Fixed-effect meta** (inverse-variance)
  - **Random-effect meta** (DerSimonian–Laird)
- Writes gzipped TSV outputs.

---

## 2) Requirements

### HPC / scheduler
- SLURM

### Software
- PLINK 1.9+
- R 4.3.0+ with `data.table`
- REGENIE (with `--t2e` support)

`survivalGWAS.slurm` currently uses modules:
```bash
module load plink/1.90b7.2
module load R/4.3.0

```

---

## 3) Citation

If you use this work, please cite:

Wei, M., & Peng, Q. (2026). *Multi-Ancestry Survival GWAS of Substance Use Initiation in the ABCD Study*. medRxiv. doi: https://doi.org/10.64898/2026.04.08.26350431.


---

## Data access and responsible use

This repository distributes analysis code and does not include participant-level research data. Users must obtain authorized access to the required ABCD Study data and comply with applicable data-use agreements and institutional requirements. The software license does not grant permission to redistribute external datasets or restricted derived outputs.

Generated files, including sample lists, phenotype and covariate tables, individual predictions, and some logs, may contain participant-level information. Do not publish these files without appropriate authorization.

For reproducibility, record the data release, accession or DOI where available, genome build, phenotype definitions, ancestry-group definitions, software versions, analysis settings, and exact code commit used.

## Validation status and updates

This repository contains research code that may include undetected errors or limitations. The current review covered source code, shell syntax, and limited synthetic checks; it did not include end-to-end execution or independent reproduction of the study results.

Source review identified issues in chromosome handling, unrelated-sample selection, prediction-file reuse, and variant-matching checks that require correction and validation. Their effects on the associated study results have not been established.

Confirmed issues and subsequent corrections will be documented in the commit history or release notes, including known effects on outputs. Users should assess suitability for their intended analyses.
