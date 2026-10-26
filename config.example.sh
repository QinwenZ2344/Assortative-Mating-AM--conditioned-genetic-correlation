# Copy to config.sh and edit. config.sh is not committed.

# Clone of https://github.com/qlu-lab/GSUB (directory that contains GSUB.R)
GSUB_DIR=/path/to/GSUB

# All pipeline outputs go here
OUTDIR=/path/to/output

# ---------------------------------------------------------------------------
# Sumstats
# GSUB trait 1 is the factor to regress out.
# GSUB trait 2 is the outcome. The pipeline runs this twice (outcome A, outcome B).
# gamma2 is the outcome effect after removing the factor.
#
# Settings match a linear GWAS on a continuous trait:
#   se.logit=FALSE  OLS=TRUE  linprob=FALSE  sample.prev=NA  population.prev=NA
# For other designs see https://github.com/qlu-lab/GSUB/wiki/2.-Running-GSUB
# N: GWAS sample size. Use NA to take N from the sumstats file.
# ---------------------------------------------------------------------------

FACTOR_SUMSTATS=/path/to/factor.sumstats.gz
FACTOR_N=NA
FACTOR_SE_LOGIT=FALSE
FACTOR_OLS=TRUE
FACTOR_LINPROB=FALSE
FACTOR_SAMPLE_PREV=NA
FACTOR_POP_PREV=NA

TRAIT_A_SUMSTATS=/path/to/traitA.sumstats.gz
TRAIT_A_N=NA
TRAIT_A_SE_LOGIT=FALSE
TRAIT_A_OLS=TRUE
TRAIT_A_LINPROB=FALSE
TRAIT_A_SAMPLE_PREV=NA
TRAIT_A_POP_PREV=NA

TRAIT_B_SUMSTATS=/path/to/traitB.sumstats.gz
TRAIT_B_N=NA
TRAIT_B_SE_LOGIT=FALSE
TRAIT_B_OLS=TRUE
TRAIT_B_LINPROB=FALSE
TRAIT_B_SAMPLE_PREV=NA
TRAIT_B_POP_PREV=NA

# ---------------------------------------------------------------------------
# Reference files used by GSUB and LDSC
# LD_CHR / WLD_CHR are directories of per-chromosome LD scores
# (files named like eur_w_ld_chr/1.l2.ldscore.gz)
# ---------------------------------------------------------------------------
# Bundled with this repository (paths are relative to the repo root)
HM3=resources/w_hm3.snplist
REF=resources/reference.1000G.maf.0.005.txt.gz

# European LD scores from the LDSC wiki (eur_w_ld_chr)
LD_CHR=/path/to/eur_w_ld_chr
WLD_CHR=/path/to/eur_w_ld_chr

# ---------------------------------------------------------------------------
# Optional GSUB arguments. Defaults match a standard continuous-trait run.
# correction: standard | conserv | none
# ---------------------------------------------------------------------------
CORRECTION=standard
INFO_FILTER=0.9
MAF_FILTER=0.01
KEEP_INDEL=FALSE

# Clone of https://github.com/bulik/ldsc (directory that contains ldsc.py)
LDSC_DIR=/path/to/ldsc
MUNGE_CHUNKSIZE=50000
