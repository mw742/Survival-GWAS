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



---

## 3) Citation

If you use or adapt this project's ideas and scripts, please cite the associated manuscript:

Multi-ancestry survival GWAS of substance use initiation in the ABCD Study and cross-trait comparison across alcohol, nicotine, cannabis, and any substance use

Mengman Wei\* and Qian Peng\* 
