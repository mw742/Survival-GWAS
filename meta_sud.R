suppressPackageStartupMessages({
  library(data.table)
})

options(datatable.threads = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1"))))

# ---------------- CONFIG ----------------
root <- "/gpfs/group/ehlers/mwei/ABCD_2/ABCD/work_provided_assets"
ancestries <- c("EUR","AFR","HISP")
traits_all <- c("alcohol","any","cannabis","nicotine")

file_pattern <- file.path(
  root,
  "surv_gwas_regenie_%s_unrel",
  "%s",
  "gwas_%s_TIME.regenie"
)

out_dir <- file.path(root, "meta_surv_sud_strandaware")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Palindromic handling:
PAL_REF_DIST = 0.10  # drop palindromic SNPs with ref AF within 0.10 of 0.5 (ambiguous)
PAL_AF_TOL   = 0.15  # keep palindromic only if aligned AF matches ref AF within this tolerance
KEEP_PAL_WITHOUT_AF = FALSE  # safer FALSE; TRUE if accept risk to keep more SNPs

# ---------------- Args ----------------
args <- commandArgs(trailingOnly=TRUE)
trait <- if (length(args)>=1) args[1] else "ALL"
ref_anc <- if (length(args)>=2) args[2] else "EUR"

if (trait != "ALL" && !trait %in% traits_all) stop("Unknown trait: ", trait)
if (!ref_anc %in% ancestries) stop("ref_anc must be one of: ", paste(ancestries, collapse=", "))

# ---------------- Helpers ----------------
comp_base <- function(x) {
  # x is character vector of A/C/G/T
  # returns complement; others -> NA
  out <- x
  out[x=="A"] <- "T"
  out[x=="T"] <- "A"
  out[x=="C"] <- "G"
  out[x=="G"] <- "C"
  out[!(x %in% c("A","C","G","T"))] <- NA_character_
  out
}

is_pal_pair <- function(a0, a1) {
  # palindromic if set is {A,T} or {C,G}
  s <- paste0(pmin(a0,a1), pmax(a0,a1))
  s %in% c("AT","CG")
}

load_gwas_chr <- function(anc, trait, chr) {
  f <- sprintf(file_pattern, anc, trait, trait)
  if (!file.exists(f)) stop("[fatal] File not found: ", f)

  # Stream just this chromosome to reduce memory. Assumes CHROM is column 1 and header exists.
  cmd <- sprintf("awk 'NR==1 || $1==%d' %s", chr, shQuote(f))
  dt <- suppressWarnings(fread(cmd=cmd, showProgress=FALSE))

  # Standardize core names
  setnames(dt,
           old=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","BETA","SE","N","A1FREQ"),
           new=c("CHR","BP","SNP","A0","A1","BETA","SE","N","A1FREQ"),
           skip_absent=TRUE)

  # Some REGENIE outputs may call frequency differently
  if (!("A1FREQ" %in% names(dt))) {
    if ("AF" %in% names(dt)) setnames(dt,"AF","A1FREQ")
    if ("A1_FREQ" %in% names(dt)) setnames(dt,"A1_FREQ","A1FREQ")
  }

  needed <- c("CHR","BP","SNP","A0","A1","BETA","SE","N")
  miss <- setdiff(needed, names(dt))
  if (length(miss)>0) stop("[fatal] Missing columns in ", f, ": ", paste(miss, collapse=", "))

  keepcols <- c(needed, intersect("A1FREQ", names(dt)))
  dt <- dt[, ..keepcols]

  dt[, CHR := as.integer(CHR)]
  dt[, BP := as.integer(BP)]
  dt <- dt[CHR == chr]
  dt <- dt[CHR %in% 1:22]
  dt <- dt[is.finite(BETA) & is.finite(SE) & SE>0 & is.finite(N) & N>0]
  dt <- dt[!is.na(SNP) & SNP!=""]
  dt <- dt[!duplicated(SNP)]
  dt[, ANC := anc]
  dt[]
}

harmonize_to_ref_strandaware <- function(dt_all, ref_anc="EUR") {
  # Build reference table from ref ancestry
  ref <- dt_all[ANC==ref_anc]
  if (nrow(ref)==0) return(dt_all[0])

  # One row per SNP as reference
  if ("A1FREQ" %in% names(ref)) {
    ref_tbl <- ref[, .(SNP, CHR_ref=CHR, BP_ref=BP, A0_ref=A0, A1_ref=A1, AF_ref=A1FREQ)]
  } else {
    ref_tbl <- ref[, .(SNP, CHR_ref=CHR, BP_ref=BP, A0_ref=A0, A1_ref=A1, AF_ref=as.numeric(NA))]
  }
  ref_tbl <- ref_tbl[!duplicated(SNP)]

  dt <- merge(dt_all, ref_tbl, by="SNP", all=FALSE)

  # Determine allele relations
  dt[, A0c := comp_base(A0)]
  dt[, A1c := comp_base(A1)]

  same    <- dt[, A0==A0_ref & A1==A1_ref]
  swap    <- dt[, A0==A1_ref & A1==A0_ref]
  strand  <- dt[, A0c==A0_ref & A1c==A1_ref]
  strswap <- dt[, A0c==A1_ref & A1c==A0_ref]

  dt[, orient := fifelse(same, "same",
                 fifelse(swap, "swap",
                 fifelse(strand, "strand",
                 fifelse(strswap, "strswap", "bad"))))]

  # Drop unalignable
  dt <- dt[orient != "bad"]
  if (nrow(dt)==0) return(dt)

  # Apply flips for swap/strswap
  dt[orient %in% c("swap","strswap"), BETA := -BETA]

  # Align alleles to ref
  dt[, `:=`(CHR=CHR_ref, BP=BP_ref, A0=A0_ref, A1=A1_ref)]
  dt[, c("CHR_ref","BP_ref") := NULL]

  # Palindromic safety filter
  dt[, pal := is_pal_pair(A0_ref, A1_ref)]

  if (!("A1FREQ" %in% names(dt))) dt[, A1FREQ := as.numeric(NA)]
  # After alignment, if we flipped beta due to swap/strswap, the aligned A1 frequency should also flip:
  dt[, AF_aligned := A1FREQ]
  dt[orient %in% c("swap","strswap"), AF_aligned := 1 - AF_aligned]

  # Decide which pal SNPs are safe
  if (!KEEP_PAL_WITHOUT_AF) {
    # require both AF_ref and AF_aligned and ref not near 0.5
    dt[pal==TRUE & (!is.finite(AF_ref) | !is.finite(AF_aligned)), pal_keep := FALSE]
  } else {
    dt[pal==TRUE & (!is.finite(AF_ref) | !is.finite(AF_aligned)), pal_keep := TRUE]
  }
  dt[pal==TRUE & is.finite(AF_ref), pal_ref_ok := abs(AF_ref - 0.5) >= PAL_REF_DIST]
  dt[pal==TRUE & is.finite(AF_ref) & is.finite(AF_aligned), pal_af_ok := abs(AF_aligned - AF_ref) <= PAL_AF_TOL]

  # pal_keep logic
  dt[pal==FALSE, pal_keep := TRUE]
  dt[pal==TRUE & is.finite(AF_ref) & is.finite(AF_aligned),
     pal_keep := (pal_ref_ok & pal_af_ok)]
  dt[pal==TRUE & is.finite(AF_ref) & !is.finite(AF_aligned),
     pal_keep := FALSE]
  dt[pal==TRUE & !is.finite(AF_ref),
     pal_keep := ifelse(KEEP_PAL_WITHOUT_AF, TRUE, FALSE)]

  dt <- dt[pal_keep==TRUE]
  dt[, c("A0_ref","A1_ref","AF_ref","A0c","A1c","orient","pal","AF_aligned","pal_keep","pal_ref_ok","pal_af_ok") := NULL]
  dt[]
}

meta_fast <- function(dt) {
  if (nrow(dt)==0) return(dt[0])
  dt[, w := 1/(SE^2)]

  s <- dt[, .(
    CHR = CHR[1],
    BP  = BP[1],
    A0  = A0[1],
    A1  = A1[1],
    k   = .N,
    N_sum = sum(N),
    sumw = sum(w),
    sumw2 = sum(w^2),
    sumwb = sum(w * BETA),
    sumwB2 = sum(w * (BETA^2))
  ), by=SNP]

  s <- s[k >= 2]
  if (nrow(s)==0) return(s)

  # Fixed effect
  s[, beta_FE := sumwb / sumw]
  s[, se_FE := sqrt(1/sumw)]
  s[, z_FE := beta_FE / se_FE]
  s[, p_FE := 2*pnorm(-abs(z_FE))]

  # Heterogeneity
  s[, Q := sumwB2 - (sumwb^2 / sumw)]
  s[, df := pmax(1, k - 1)]
  s[, I2 := fifelse(Q > 0, pmax(0, (Q - df)/Q)*100, 0)]
  s[, C := sumw - (sumw2 / sumw)]
  s[, tau2_DL := fifelse(C > 0, pmax(0, (Q - df)/C), 0)]

  # Random effects (DL)
  dt2 <- merge(dt[, .(SNP, BETA, SE)], s[, .(SNP, tau2_DL)], by="SNP", all=FALSE)
  dt2[, w_re := 1/(SE^2 + tau2_DL)]
  re <- dt2[, .(sumw_re=sum(w_re), sumwb_re=sum(w_re*BETA)), by=SNP]
  s <- merge(s, re, by="SNP", all.x=TRUE)

  s[, beta_RE := sumwb_re / sumw_re]
  s[, se_RE := sqrt(1/sumw_re)]
  s[, z_RE := beta_RE / se_RE]
  s[, p_RE := 2*pnorm(-abs(z_RE))]

  s[, c("sumw2","sumwb","sumwB2","C","sumw_re","sumwb_re") := NULL]
  setcolorder(s, c("CHR","BP","SNP","A0","A1","k","N_sum",
                   "beta_FE","se_FE","z_FE","p_FE",
                   "beta_RE","se_RE","z_RE","p_RE",
                   "Q","I2","tau2_DL","df"))
  s[]
}

run_one_trait <- function(trait) {
  cat(sprintf("\n[diag] ===== Trait: %s (ref=%s) =====\n", trait, ref_anc))

  out_tsv <- file.path(out_dir, sprintf("meta_%s_EUR_AFR_HISP_surv_strandaware.tsv", trait))
  if (file.exists(out_tsv)) file.remove(out_tsv)
  if (file.exists(paste0(out_tsv, ".gz"))) file.remove(paste0(out_tsv, ".gz"))

  wrote_header <- FALSE
  total_snps <- 0L

  for (chr in 1:22) {
    cat(sprintf("[diag] Trait=%s CHR=%d loading...\n", trait, chr))
    lst <- lapply(ancestries, function(anc) load_gwas_chr(anc, trait, chr))
    dt <- rbindlist(lst, use.names=TRUE, fill=TRUE)
    cat(sprintf("[diag] Trait=%s CHR=%d rows loaded=%d\n", trait, chr, nrow(dt)))
    if (nrow(dt)==0) next

    before <- nrow(dt)
    dt <- harmonize_to_ref_strandaware(dt, ref_anc=ref_anc)
    cat(sprintf("[diag] Trait=%s CHR=%d rows after harmonization=%d (kept %.1f%%)\n",
                trait, chr, nrow(dt), 100*nrow(dt)/before))
    if (nrow(dt)==0) next

    res <- meta_fast(dt)
    cat(sprintf("[diag] Trait=%s CHR=%d SNPs meta-analyzed=%d\n", trait, chr, nrow(res)))
    if (nrow(res)==0) next

    fwrite(res, out_tsv, sep="\t", quote=FALSE, na="NA",
           col.names = !wrote_header, append = wrote_header)
    wrote_header <- TRUE
    total_snps <- total_snps + nrow(res)

    rm(dt,res,lst); gc()
  }

  if (!wrote_header) {
    warning("[warn] No results written for trait ", trait)
    return(invisible(NULL))
  }

  system(sprintf("gzip -f %s", shQuote(out_tsv)))
  cat(sprintf("[done] Trait %s complete. Total SNPs written=%d\n", trait, total_snps))
  invisible(TRUE)
}

cat("[diag] .libPaths():\n"); print(.libPaths())
cat(sprintf("[diag] SLURM_CPUS_PER_TASK=%s\n", Sys.getenv("SLURM_CPUS_PER_TASK","NA")))

if (trait=="ALL") {
  for (tr in traits_all) run_one_trait(tr)
} else {
  run_one_trait(trait)
}

cat("\n[done] Strand-aware meta-analysis finished.\n")

