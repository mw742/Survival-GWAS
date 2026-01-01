#!/usr/bin/env bash
#SBATCH --job-name=EUR_survGWAS
#SBATCH --output=EUR_survGWAS_%A_%a.log
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --partition=highmem

set -euo pipefail
IFS=$'\n\t'

diag(){ echo "[diag] $*" >&2; }
note(){ echo "[note]  $*" >&2; }
fatal(){ echo "[fatal] $*" >&2; exit 2; }

echo "[env] Host=$(hostname)  Date=$(date '+%F %T')"
echo "[env] SLURM_JOB_ID=${SLURM_JOB_ID:-NA}  ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-0}"

ROOT="/gpfs/group/ehlers/mwei/ABCD_2/ABCD"
WORK="${ROOT}/work_provided_assets"

EUR_DIR="${WORK}/ancestry/white_subsetID_20251102_201950"
EUR_PREFIX="${EUR_DIR}/ABCD_white_only"
EUR_KEEP_ALL="${EUR_DIR}/white.keep"                        # can be 1- or 2-col
EUR_UNREL_LOCAL="${EUR_DIR}/white_unrel.king.cutoff.out.id" # can be 1- or 2-col
CURATED_UNREL_1C="${ROOT}/genomics_sample03/ABCD_GRM/GENESIS/ABCD_202209.updated.nodups.curated_unrelateds.txt"

COVAR_TSV=""
for C in "${WORK}/covariates_final.fixed.tsv" "${WORK}/covariates_final.tsv"; do
  if [[ -s "$C" ]]; then COVAR_TSV="$C"; break; fi
done
[[ -n "${COVAR_TSV}" ]] || fatal "Could not find covariates_final.fixed.tsv or covariates_final.tsv under ${WORK}"
diag "Using covariates: ${COVAR_TSV}"

PHENO_ALC="/gpfs/group/ehlers/mwei/ABCD/init_phenos/alcohol_survival.tsv"
PHENO_NIC="/gpfs/group/ehlers/mwei/ABCD/init_phenos/nicotine_survival.tsv"
PHENO_CAN="/gpfs/group/ehlers/mwei/ABCD/init_phenos/cannabis_survival.tsv"
PHENO_ANY="/gpfs/group/ehlers/mwei/ABCD/init_phenos/any_survival.tsv"

INFO_INCLUDE_LIST=""
LONG_LD_BUILD="${LONG_LD_BUILD:-38}"

module load plink/1.90b7.2 2>/dev/null || true
module load R/4.3.0         2>/dev/null || true

# Prefer local REGENIE under WORK, fall back to PATH
REGENIE_BIN="${WORK}/mamba_env/bin/regenie"
if [[ ! -x "${REGENIE_BIN}" ]]; then
  REGENIE_BIN="$(command -v regenie || true)"
fi
[[ -x "${REGENIE_BIN}" ]] || fatal "REGENIE not found in PATH and ${WORK}/mamba_env/bin"
"${REGENIE_BIN}" --version || true

command -v plink   >/dev/null || fatal "plink not found in PATH"
command -v Rscript >/dev/null || fatal "Rscript not found in PATH"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export TMPDIR="${SLURM_TMPDIR:-/tmp}"
REGENIE_THREADS="${SLURM_CPUS_PER_TASK:-16}"

OUTTOP="${WORK}/surv_gwas_regenie_EUR_unrel"
mkdir -p "${OUTTOP}"
STD_PREFIX="${OUTTOP}/EUR_chrstd"
LOCKDIR="${OUTTOP}/.eur_chrstd.lockdir"

for ext in bed bim fam; do [[ -s "${EUR_PREFIX}.${ext}" ]] || fatal "Missing ${EUR_PREFIX}.${ext}"; done
[[ -s "${COVAR_TSV}" ]] || fatal "Missing covariates: ${COVAR_TSV}"

LR_FILE="${OUTTOP}/longLD_build${LONG_LD_BUILD}.ranges"
if [[ ! -s "${LR_FILE}" ]]; then
  if [[ "${LONG_LD_BUILD}" == "37" ]]; then
    cat > "${LR_FILE}" <<'EOF'
6 25000000 35000000 MHC
8 7000000 13000000 chr8inv
17 40000000 45000000 chr17q21
EOF
  else
    cat > "${LR_FILE}" <<'EOF'
6 28500000 33500000 MHC
8 7150000 12500000 chr8inv
17 43000000 48000000 chr17q21
EOF
  fi
fi

# ---------- Standardize chromosomes on the EUR bedset ----------
if [[ ! -s "${STD_PREFIX}.bed" ]]; then
  if mkdir "${LOCKDIR}" 2>/dev/null; then
    trap 'rmdir "${LOCKDIR}" 2>/dev/null || true' EXIT
    diag "Standardizing EUR chromosomes…"
    MAP_ALL="${OUTTOP}/eur_chrmap_all.txt"
    MAP_OK="${OUTTOP}/eur_chrmap_ok.txt"
    DROP_SNPS="${OUTTOP}/eur_drop_unknown_chr.txt"
    awk 'BEGIN{OFS="\t"}{
      chr=$1;snp=$2;gsub(/^chr/,"",chr);
      if(chr=="X") new=23; else if(chr=="Y") new=24; else if(chr=="XY") new=25;
      else if(chr=="MT"||chr=="M") new=26;
      else if(chr~/^[0-9]+$/&&chr>=1&&chr<=22) new=chr; else new="NA";
      print snp,new
    }' "${EUR_PREFIX}.bim" > "${MAP_ALL}"
    awk '$2!="NA"{print}' "${MAP_ALL}" > "${MAP_OK}"
    awk '$2=="NA"{print $1}' "${MAP_ALL}" > "${DROP_SNPS}"
    plink \
      --bfile "${EUR_PREFIX}" \
      --exclude "${DROP_SNPS}" \
      --update-chr "${MAP_OK}" 1 2 \
      --chr-set 26 \
      --make-bed \
      --out "${STD_PREFIX}"
  else
    note "Waiting for standardization done by another task…"
    while [[ ! -s "${STD_PREFIX}.bed" ]]; do sleep 3; done
  fi
else
  diag "Reusing standardized: ${STD_PREFIX}"
fi

# ---------- Select trait by array index ----------
case "${SLURM_ARRAY_TASK_ID:-0}" in
  0) TRAIT="alcohol";  PHENO="${PHENO_ALC}" ;;
  1) TRAIT="nicotine"; PHENO="${PHENO_NIC}" ;;
  2) TRAIT="cannabis"; PHENO="${PHENO_CAN}" ;;
  3) TRAIT="any";      PHENO="${PHENO_ANY}" ;;
  *) fatal "invalid array id" ;;
esac
export PHENO

OUTDIR="${OUTTOP}/${TRAIT}"
mkdir -p "${OUTDIR}"

PHENO_RG="${OUTDIR}/pheno.tsv"
COVAR_RG="${OUTDIR}/covar.tsv"
KEEP_IDS_ALL="${OUTDIR}/keep_all_EUR.txt"
KEEP_IDS_UNREL="${OUTDIR}/keep_unrel_EUR.txt"
QCBASE="${OUTDIR}/step1_qc"

# ---------- KEEP LISTS ----------
if [[ ! -s "${EUR_KEEP_ALL}" ]]; then
  note "EUR_KEEP_ALL missing/empty; deriving from ${EUR_PREFIX}.fam"
  awk 'BEGIN{OFS="\t"}{print $1,$2}' "${EUR_PREFIX}.fam" > "${EUR_KEEP_ALL}"
fi

if [[ ! -s "${EUR_UNREL_LOCAL}" ]]; then
  if [[ -s "${CURATED_UNREL_1C}" ]]; then
    note "EUR_UNREL_LOCAL missing; building from curated_unrelateds"
    awk 'BEGIN{OFS="\t"}{print $1,$1}' "${CURATED_UNREL_1C}" > "${EUR_UNREL_LOCAL}"
  else
    note "EUR_UNREL_LOCAL missing and no curated file; using EUR_KEEP_ALL as fallback"
    cp -f "${EUR_KEEP_ALL}" "${EUR_UNREL_LOCAL}"
  fi
fi

export COVAR_TSV PHENO PHENO_RG COVAR_RG KEEP_IDS_ALL KEEP_IDS_UNREL
export EUR_KEEP_ALL EUR_UNREL_LOCAL CURATED_UNREL_1C

# ---------- R PREP ----------
Rscript - <<'RS'
suppressWarnings(suppressMessages({ library(data.table) }))

get_keep_dt <- function(pth) {
  if (!nzchar(pth) || !file.exists(pth) || file.info(pth)$size == 0) return(NULL)
  dt <- tryCatch(fread(pth, header = FALSE), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  if (ncol(dt) == 1) {
    setnames(dt, "V1", "IID")
    dt[, FID := IID]
    setcolorder(dt, c("FID","IID"))
  } else {
    setnames(dt, c("FID","IID")[seq_len(min(2L, ncol(dt)))])
    dt <- dt[, .(FID,IID)]
  }
  unique(dt)
}

covar_tsv       <- Sys.getenv("COVAR_TSV")
keep_all_pth    <- Sys.getenv("EUR_KEEP_ALL")
unrel_local_pth <- Sys.getenv("EUR_UNREL_LOCAL")
curated_unrel   <- Sys.getenv("CURATED_UNREL_1C")
keep_all_out    <- Sys.getenv("KEEP_IDS_ALL")
keep_unrel_out  <- Sys.getenv("KEEP_IDS_UNREL")
covar_rg        <- Sys.getenv("COVAR_RG")
pheno_src       <- Sys.getenv("PHENO")
pheno_rg        <- Sys.getenv("PHENO_RG")

keep_all <- get_keep_dt(keep_all_pth)
unrel    <- get_keep_dt(unrel_local_pth)

if (nzchar(curated_unrel) && file.exists(curated_unrel) && file.info(curated_unrel)$size > 0) {
  cu <- fread(curated_unrel, header = FALSE)
  if (ncol(cu) == 1) {
    setnames(cu, "V1","IID"); cu[, FID := IID]
  } else {
    setnames(cu, c("FID","IID")[1:2])
  }
  cu <- cu[, .(FID,IID)]
  unrel <- unique(rbind(unrel, cu, fill = TRUE))
}

stopifnot(!is.null(keep_all) && nrow(keep_all) > 0)
stopifnot(!is.null(unrel)    && nrow(unrel)    > 0)

cov <- fread(covar_tsv)
setnames(cov, toupper(names(cov)))
if (!"IID" %in% names(cov) && "SUBJID" %in% names(cov)) setnames(cov, "SUBJID", "IID")
cov[, FID := IID]

if (!"SEX"  %in% names(cov) && "GENDER" %in% names(cov))  setnames(cov, "GENDER",  "SEX")
if (!"AGE"  %in% names(cov) && "AGE_YRS" %in% names(cov)) setnames(cov, "AGE_YRS", "AGE")
if (!"SITE" %in% names(cov) && "SITE_ID" %in% names(cov)) setnames(cov, "SITE_ID", "SITE")

req_cov <- c("SEX","AGE","SITE")
missing_cov <- setdiff(req_cov, names(cov))
if (length(missing_cov)) stop(sprintf("Missing covariates after normalization: %s", paste(missing_cov, collapse=",")))

cov[, AGE  := as.numeric(AGE)]
cov[, SEX  := as.numeric(SEX)]
cov[, SITE := as.character(SITE)]
pc_cols <- grep("^PC[0-9]+$", names(cov), value = TRUE)

cov <- cov[!is.na(AGE) & !is.na(SEX) & nzchar(SITE)]
cov <- cov[, c("FID","IID", req_cov, pc_cols), with = FALSE]

ph <- fread(pheno_src)
setnames(ph, toupper(names(ph)))
if (!"FID" %in% names(ph)) ph[, FID := IID]

reqp <- c("IID","TIME","EVENT")
if (!all(reqp %in% names(ph))) {
  stop(sprintf("Phenotype must have IID,TIME,EVENT (have: %s)", paste(names(ph), collapse=",")))
}

ph[, TIME  := as.numeric(TIME)]
ph[, EVENT := as.integer(EVENT)]
ph <- ph[is.finite(TIME) & !is.na(EVENT), .(FID,IID,TIME,EVENT)]

ids_all   <- merge(keep_all , cov[,.(FID,IID)], by = c("FID","IID"))
ids_unrel <- merge(unrel    , cov[,.(FID,IID)], by = c("FID","IID"))

fwrite(ids_all,   keep_all_out,   sep = "\t", col.names = FALSE)
fwrite(ids_unrel, keep_unrel_out, sep = "\t", col.names = FALSE)
fwrite(cov,       covar_rg,       sep = "\t", quote = FALSE, na = "NA")
fwrite(ph,        pheno_rg,       sep = "\t", quote = FALSE, na = "NA")

cat(sprintf("[diag] keep_all:  %d\n", nrow(ids_all)))
cat(sprintf("[diag] keep_unrel:%d\n", nrow(ids_unrel)))
cat(sprintf("[diag] cov rows:  %d; ph rows: %d\n", nrow(cov), nrow(ph)))
RS

diag "Preview keep_all (head):";   head -n 3 "${KEEP_IDS_ALL}"   || true
diag "Preview keep_unrel (head):"; head -n 3 "${KEEP_IDS_UNREL}" || true

# ---------- COVARIATE COLUMNS ----------
PC_LIST="$(
  awk 'NR==1{
    k=0;
    for(i=1;i<=NF;i++){
      if($i ~ /^PC[0-9][0-9]*$/){ v[++k]=$i }
    }
  }
  END{
    for(i=1;i<=k;i++){ printf "%s%s", v[i], (i<k?",":"\n") }
  }' "${COVAR_RG}" || true
)"
COVAR_COLS_FULL="SEX,AGE,SITE"
if [[ -n "${PC_LIST}" ]]; then COVAR_COLS_FULL="${COVAR_COLS_FULL},${PC_LIST}"; fi
HAS_BATCH=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="BATCH") f=1} END{print f?1:0}' "${COVAR_RG}")
CATCOV_FULL="SITE"
if [[ "${HAS_BATCH}" == "1" ]]; then CATCOV_FULL="SITE,BATCH"; fi

COVAR_COLS_LITE="SEX,AGE"
if [[ -n "${PC_LIST}" ]]; then COVAR_COLS_LITE="${COVAR_COLS_LITE},${PC_LIST}"; fi
CATCOV_LITE=""

diag "Covariates (full) : ${COVAR_COLS_FULL} ; cat=${CATCOV_FULL}"
diag "Covariates (lite) : ${COVAR_COLS_LITE} ; cat=${CATCOV_LITE}"

# ---------- STEP-1 SNP SET (PLINK QC ON UNRELATEDS) ----------
diag "[step1] Build ACGT, geno<=0.02, HWE 1e-6 midp, MAC>=10, exclude long-LD (UNRELATEDS)"
PLINK_BASE=(--bfile "${STD_PREFIX}" --keep "${KEEP_IDS_UNREL}" --chr-set 26 --snps-only just-acgt \
            --geno 0.02 --mac 10 --hwe 1e-6 midp --exclude range "${LR_FILE}")
if [[ -n "${INFO_INCLUDE_LIST}" && -s "${INFO_INCLUDE_LIST}" ]]; then
  plink "${PLINK_BASE[@]}" --extract "${INFO_INCLUDE_LIST}" --write-snplist --out "${QCBASE}"
else
  plink "${PLINK_BASE[@]}" --write-snplist --out "${QCBASE}"
fi

diag "[step1] LD prune r2=0.2 (UNRELATEDS)"
plink \
  --bfile "${STD_PREFIX}" \
  --keep "${KEEP_IDS_UNREL}" \
  --extract "${QCBASE}.snplist" \
  --chr-set 26 \
  --exclude range "${LR_FILE}" \
  --indep-pairwise 1000 100 0.2 \
  --out "${QCBASE}"

EXTRACT_SNPS="${QCBASE}.prune.in"
[[ -s "${EXTRACT_SNPS}" ]] || fatal "No SNPs left for Step-1 after filtering/pruning"

# ---------- REGENIE STEP-1 (runs on ALL EUR) ----------
STEP1_OUT="${OUTDIR}/fit_${TRAIT}"
diag "[regenie] STEP-1 (FULL covars, ALL EUR) -> ${STEP1_OUT}"
# Using SNP set built on unrelateds, but Step-1 fit on all EUR
"${REGENIE_BIN}" \
  --step 1 \
  --bed "${STD_PREFIX}" \
  --extract "${EXTRACT_SNPS}" \
  --covarFile "${COVAR_RG}" \
  --covarColList "${COVAR_COLS_FULL}" \
  --catCovarList "${CATCOV_FULL}" \
  --maxCatLevels 512 \
  --phenoFile "${PHENO_RG}" \
  --phenoColList TIME \
  --eventColList EVENT \
  --t2e \
  --keep "${KEEP_IDS_ALL}" \
  --bsize 1000 \
  --cv 5 \
  --threads "${REGENIE_THREADS}" \
  --out "${STEP1_OUT}" || true

STEP1_DEF_LIST="${STEP1_OUT}_pred.list"
NEED_RETRY=0
if [[ ! -s "${STEP1_DEF_LIST}" ]] || ! awk 'NF>=2{ok=1;exit} END{exit ok?0:1}' "${STEP1_DEF_LIST}"; then
  note "[warn] No usable LOCO predictions from default Step-1; building a larger balanced SNP set and retrying…"
  NEED_RETRY=1
fi

# ---------- AUTO-RETRY: PLINK QC STILL ON UNRELATEDS ----------
QCBASE_R="${OUTDIR}/step1_qc_retry"
if [[ "${NEED_RETRY}" == "1" ]]; then
  plink \
    --bfile "${STD_PREFIX}" \
    --keep "${KEEP_IDS_UNREL}" \
    --chr-set 26 \
    --snps-only just-acgt \
    --geno 0.02 \
    --mac 5 \
    --hwe 1e-6 midp \
    --exclude range "${LR_FILE}" \
    --write-snplist \
    --out "${QCBASE_R}"

  MIN_PER_CHR=800
  TGT_PER_CHR=2000
  R2=0.8
  WIN=1000
  STEP=100
  PRUNED_DIR="${OUTDIR}/rescue_perchr"
  mkdir -p "${PRUNED_DIR}"
  : > "${PRUNED_DIR}/extract_merged.snps"

  for CHR in $(seq 1 22); do
    diag "rescue chr${CHR}: pruning candidates (UNRELATEDS)…"
    awk -v c="${CHR}" 'NR==FNR{keep[$1]=1;next} $1==c && ($2 in keep){print $2}' \
        "${QCBASE_R}.snplist" "${STD_PREFIX}.bim" > "${PRUNED_DIR}/chr${CHR}.cand.snps"

    if [[ ! -s "${PRUNED_DIR}/chr${CHR}.cand.snps" ]]; then
      note "chr${CHR}: no candidates after first pass; relaxing to geno<=0.05, MAC>=3"
      plink \
        --bfile "${STD_PREFIX}" \
        --keep "${KEEP_IDS_UNREL}" \
        --chr ${CHR} \
        --snps-only just-acgt \
        --geno 0.05 \
        --mac 3 \
        --write-snplist \
        --out "${PRUNED_DIR}/chr${CHR}.loose" >/dev/null 2>&1 || true
      cp -f "${PRUNED_DIR}/chr${CHR}.loose.snplist" "${PRUNED_DIR}/chr${CHR}.cand.snps" 2>/dev/null || true
    fi
    [[ -s "${PRUNED_DIR}/chr${CHR}.cand.snps" ]] || { note "chr${CHR}: still empty, skipping"; continue; }

    plink \
      --bfile "${STD_PREFIX}" \
      --keep "${KEEP_IDS_UNREL}" \
      --chr ${CHR} \
      --extract "${PRUNED_DIR}/chr${CHR}.cand.snps" \
      --indep-pairwise ${WIN} ${STEP} ${R2} \
      --out "${PRUNED_DIR}/chr${CHR}" >/dev/null 2>&1 || true

    SRC="${PRUNED_DIR}/chr${CHR}.prune.in"
    [[ -s "${SRC}" ]] || { note "chr${CHR}: prune produced 0; skipping"; continue; }

    head -n "${TGT_PER_CHR}" "${SRC}" >> "${PRUNED_DIR}/extract_merged.snps"
  done

  awk '!seen[$0]++' "${PRUNED_DIR}/extract_merged.snps" > "${QCBASE_R}.balanced.preQC.snps"
  [[ -s "${QCBASE_R}.balanced.preQC.snps" ]] || fatal "Rescue: merged list empty"

  plink \
    --bfile "${STD_PREFIX}" \
    --keep "${KEEP_IDS_UNREL}" \
    --extract "${QCBASE_R}.balanced.preQC.snps" \
    --snps-only just-acgt \
    --geno 0.02 \
    --mac 5 \
    --hwe 1e-6 midp \
    --write-snplist \
    --out "${QCBASE_R}.balanced.qc" >/dev/null 2>&1 || true

  plink \
    --bfile "${STD_PREFIX}" \
    --keep "${KEEP_IDS_UNREL}" \
    --extract "${QCBASE_R}.balanced.qc.snplist" \
    --indep-pairwise ${WIN} ${STEP} ${R2} \
    --out "${QCBASE_R}.final" >/dev/null 2>&1 || true

  EXTRACT_R="${QCBASE_R}.final.prune.in"
  [[ -s "${EXTRACT_R}" ]] || fatal "Rescue: no SNPs after balanced final prune"

  awk 'NR==FNR{keep[$1]=1;next} ($2 in keep){c[$1]++} END{for(i=1;i<=22;i++) printf("chr%-2d : %d\n", i, (c[i]+0))}' \
      "${EXTRACT_R}" "${STD_PREFIX}.bim" >&2
  plink \
    --bfile "${STD_PREFIX}" \
    --keep "${KEEP_IDS_UNREL}" \
    --extract "${EXTRACT_R}" \
    --freq \
    --out "${OUTDIR}/step1_retry_freq" >/dev/null 2>&1 || true
  awk 'NR>1{n++; if($5==0 || $5==1) z++} END{printf("[diag] retry set: n=%d, zero/one MAF=%d\n", n, z)}' \
      "${OUTDIR}/step1_retry_freq.frq" || true

  STEP1_OUT_R="${OUTDIR}/fit_${TRAIT}_retryA"
  diag "[regenie] STEP-1 RETRY A (FULL covars, ALL EUR) -> ${STEP1_OUT_R}"
  "${REGENIE_BIN}" \
    --step 1 \
    --bed "${STD_PREFIX}" \
    --extract "${EXTRACT_R}" \
    --covarFile "${COVAR_RG}" \
    --covarColList "${COVAR_COLS_FULL}" \
    --catCovarList "${CATCOV_FULL}" \
    --maxCatLevels 512 \
    --phenoFile "${PHENO_RG}" \
    --phenoColList TIME \
    --eventColList EVENT \
    --t2e \
    --keep "${KEEP_IDS_ALL}" \
    --bsize 1000 \
    --cv 5 \
    --threads "${REGENIE_THREADS}" \
    --out "${STEP1_OUT_R}" || true

  STEP1_OUT_R2="${OUTDIR}/fit_${TRAIT}_retryB_lite"
  if [[ ! -s "${STEP1_OUT_R}_pred.list" ]] || ! awk 'NF>=2{ok=1;exit} END{exit ok?0:1}' "${STEP1_OUT_R}_pred.list"; then
    diag "[regenie] STEP-1 RETRY B (LITE covars, ALL EUR) -> ${STEP1_OUT_R2}"
    "${REGENIE_BIN}" \
      --step 1 \
      --bed "${STD_PREFIX}" \
      --extract "${EXTRACT_R}" \
      --covarFile "${COVAR_RG}" \
      --covarColList "${COVAR_COLS_LITE}" \
      --catCovarList "${CATCOV_LITE}" \
      --phenoFile "${PHENO_RG}" \
      --phenoColList TIME \
      --eventColList EVENT \
      --t2e \
      --keep "${KEEP_IDS_ALL}" \
      --bsize 1000 \
      --cv 5 \
      --threads "${REGENIE_THREADS}" \
      --out "${STEP1_OUT_R2}" || true
  fi
fi

BEST_LIST=""
for L in "${OUTDIR}/fit_${TRAIT}_retryB_lite_pred.list" \
         "${OUTDIR}/fit_${TRAIT}_retryA_pred.list" \
         "${OUTDIR}/fit_${TRAIT}_pred.list"; do
  if [[ -s "${L}" ]]; then
    if awk 'NF>=2{ok=1;exit} END{exit ok?0:1}' "${L}"; then
      BEST_LIST="${L}"
      break
    fi
  fi
done

[[ -n "${BEST_LIST}" ]] || fatal "No usable *_pred.list found in ${OUTDIR}. Check Step-1 logs."

ABS_PRED_LIST="${OUTDIR}/fit_${TRAIT}_pred.ABS.list"
awk -v outdir="${OUTDIR}" 'NF>=2{
  pheno=$1; file=$2;
  if (file !~ /^\//) file = outdir "/" file;
  print pheno, file;
}' "${BEST_LIST}" > "${ABS_PRED_LIST}"

N_GOOD=$(wc -l < "${ABS_PRED_LIST}" || echo 0)
diag "Usable LOCO prediction entries (phenotype,file): ${N_GOOD}"
(( N_GOOD >= 1 )) || fatal "Prediction list is empty after normalization."

STEP1_BEST="${ABS_PRED_LIST}"

# ---------- STEP-2 (ALL EUR, same covars) ----------
PC_LIST="$(
  awk 'NR==1{
    k=0;
    for(i=1;i<=NF;i++){
      if($i ~ /^PC[0-9][0-9]*$/){ v[++k]=$i }
    }
  }
  END{
    for(i=1;i<=k;i++){ printf "%s%s", v[i], (i<k?",":"\n") }
  }' "${COVAR_RG}" || true
)"
COVAR_COLS_STEP2="SEX,AGE,SITE"
if [[ -n "${PC_LIST}" ]]; then COVAR_COLS_STEP2="${COVAR_COLS_STEP2},${PC_LIST}"; fi
HAS_BATCH=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="BATCH") f=1} END{print f?1:0}' "${COVAR_RG}")
CATCOV_STEP2="SITE"
if [[ "${HAS_BATCH}" == "1" ]]; then CATCOV_STEP2="SITE,BATCH"; fi

STEP2_OUT="${OUTDIR}/gwas_${TRAIT}"
diag "[regenie] STEP-2 (ALL EUR) -> ${STEP2_OUT}"
"${REGENIE_BIN}" \
  --step 2 \
  --bed "${STD_PREFIX}" \
  --covarFile "${COVAR_RG}" \
  --covarColList "${COVAR_COLS_STEP2}" \
  --catCovarList "${CATCOV_STEP2}" \
  --maxCatLevels 512 \
  --phenoFile "${PHENO_RG}" \
  --phenoColList TIME \
  --eventColList EVENT \
  --t2e \
  --keep "${KEEP_IDS_ALL}" \
  --bsize 400 \
  --pred "${STEP1_BEST}" \
  --minMAC 20 \
  --threads "${REGENIE_THREADS}" \
  --out "${STEP2_OUT}"

echo "[done] Trait=${TRAIT} OUTDIR=${OUTDIR}"

