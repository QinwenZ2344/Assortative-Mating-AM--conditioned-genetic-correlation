#!/usr/bin/env bash
# AM-conditioned genetic correlation.
# 1. GSUB: regress the factor out of trait A and out of trait B
# 2. Extract gamma2 sumstats
# 3. LDSC rg for the original GWAS and for the conditioned sumstats
# Results are written as a two-row table: unconditioned and conditioned.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${ROOT}/config.sh}"

if [[ ! -f "${CONFIG}" ]]; then
  echo "Config not found: ${CONFIG}" >&2
  echo "Copy config.example.sh to config.sh and edit it." >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG}"

: "${GSUB_DIR:?}"
: "${OUTDIR:?}"
: "${FACTOR_SUMSTATS:?}" "${FACTOR_N:?}"
: "${TRAIT_A_SUMSTATS:?}" "${TRAIT_A_N:?}"
: "${TRAIT_B_SUMSTATS:?}" "${TRAIT_B_N:?}"
: "${HM3:?}" "${LD_CHR:?}" "${WLD_CHR:?}" "${REF:?}"
: "${LDSC_DIR:?}"

# Paths in config.sh that do not start with / are relative to this repository.
resolve_path() {
  if [[ "$1" == /* ]]; then
    printf '%s\n' "$1"
  else
    printf '%s\n' "${ROOT}/$1"
  fi
}
HM3="$(resolve_path "${HM3}")"
REF="$(resolve_path "${REF}")"

CORRECTION="${CORRECTION:-standard}"
INFO_FILTER="${INFO_FILTER:-0.9}"
MAF_FILTER="${MAF_FILTER:-0.01}"
KEEP_INDEL="${KEEP_INDEL:-FALSE}"
MUNGE_CHUNKSIZE="${MUNGE_CHUNKSIZE:-50000}"
MUNGE_SUMSTATS="${LDSC_DIR}/munge_sumstats.py"
LDSC_PY="${LDSC_DIR}/ldsc.py"

for f in "${FACTOR_SUMSTATS}" "${TRAIT_A_SUMSTATS}" "${TRAIT_B_SUMSTATS}" \
         "${HM3}" "${REF}" "${GSUB_DIR}/GSUB.R" \
         "${MUNGE_SUMSTATS}" "${LDSC_PY}"; do
  if [[ ! -e "${f}" ]]; then
    echo "Missing input: ${f}" >&2
    exit 1
  fi
done
for d in "${LD_CHR}" "${WLD_CHR}"; do
  if [[ ! -d "${d}" ]]; then
    echo "Missing directory: ${d}" >&2
    exit 1
  fi
done
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Required command not on PATH: Rscript" >&2
  exit 1
fi

mkdir -p "${OUTDIR}/traitA" "${OUTDIR}/traitB" "${OUTDIR}/conditioned" "${OUTDIR}/ldsc"

run_gsub() {
  local outcome_sumstats="$1"
  local outcome_n="$2"
  local outcome_se_logit="$3"
  local outcome_ols="$4"
  local outcome_linprob="$5"
  local outcome_sample_prev="$6"
  local outcome_pop_prev="$7"
  local out_path="$8"

  (
    cd "${GSUB_DIR}"
    Rscript GSUB.R \
      --sumstat_files "${FACTOR_SUMSTATS},${outcome_sumstats}" \
      --output_path "${out_path}" \
      --correction "${CORRECTION}" \
      --N "${FACTOR_N},${outcome_n}" \
      --info.filter "${INFO_FILTER}" \
      --maf.filter "${MAF_FILTER}" \
      --sample.prev "${FACTOR_SAMPLE_PREV},${outcome_sample_prev}" \
      --population.prev "${FACTOR_POP_PREV},${outcome_pop_prev}" \
      --se.logit "${FACTOR_SE_LOGIT},${outcome_se_logit}" \
      --OLS "${FACTOR_OLS},${outcome_ols}" \
      --linprob "${FACTOR_LINPROB},${outcome_linprob}" \
      --keep.indel "${KEEP_INDEL}" \
      --hm3 "${HM3}" \
      --ld "${LD_CHR}" \
      --wld "${WLD_CHR}" \
      --ref "${REF}"
  )
}

echo "Step 1: GSUB for trait A and trait B"
run_gsub "${TRAIT_A_SUMSTATS}" "${TRAIT_A_N}" \
  "${TRAIT_A_SE_LOGIT}" "${TRAIT_A_OLS}" "${TRAIT_A_LINPROB}" \
  "${TRAIT_A_SAMPLE_PREV}" "${TRAIT_A_POP_PREV}" \
  "${OUTDIR}/traitA" &
pid_a=$!
run_gsub "${TRAIT_B_SUMSTATS}" "${TRAIT_B_N}" \
  "${TRAIT_B_SE_LOGIT}" "${TRAIT_B_OLS}" "${TRAIT_B_LINPROB}" \
  "${TRAIT_B_SAMPLE_PREV}" "${TRAIT_B_POP_PREV}" \
  "${OUTDIR}/traitB" &
pid_b=$!
wait "${pid_a}"
wait "${pid_b}"

echo "Step 2: extract gamma2 sumstats"
Rscript "${ROOT}/scripts/02_extract_conditioned.R" \
  "${OUTDIR}/traitA/gsub_analytical_all.txt.gz" \
  "${OUTDIR}/conditioned/traitA.gamma2.txt.gz"
Rscript "${ROOT}/scripts/02_extract_conditioned.R" \
  "${OUTDIR}/traitB/gsub_analytical_all.txt.gz" \
  "${OUTDIR}/conditioned/traitB.gamma2.txt.gz"

munge_one() {
  # --ignore BETA: sumstats contain both BETA and Z; keep Z as the signed column.
  python "${MUNGE_SUMSTATS}" \
    --merge-alleles "${HM3}" \
    --sumstats "$1" \
    --ignore BETA \
    --chunksize "${MUNGE_CHUNKSIZE}" \
    --out "$2"
}

echo "Step 3: munge original and conditioned sumstats"
munge_one "${TRAIT_A_SUMSTATS}" "${OUTDIR}/ldsc/uncond_traitA.munged" &
munge_one "${TRAIT_B_SUMSTATS}" "${OUTDIR}/ldsc/uncond_traitB.munged" &
munge_one "${OUTDIR}/conditioned/traitA.gamma2.txt.gz" "${OUTDIR}/ldsc/cond_traitA.munged" &
munge_one "${OUTDIR}/conditioned/traitB.gamma2.txt.gz" "${OUTDIR}/ldsc/cond_traitB.munged" &
wait

echo "Step 3: LDSC genetic correlations"
python "${LDSC_PY}" \
  --rg "${OUTDIR}/ldsc/uncond_traitA.munged.sumstats.gz,${OUTDIR}/ldsc/uncond_traitB.munged.sumstats.gz" \
  --ref-ld-chr "${LD_CHR}/" \
  --w-ld-chr "${WLD_CHR}/" \
  --out "${OUTDIR}/ldsc/unconditioned_rg"

python "${LDSC_PY}" \
  --rg "${OUTDIR}/ldsc/cond_traitA.munged.sumstats.gz,${OUTDIR}/ldsc/cond_traitB.munged.sumstats.gz" \
  --ref-ld-chr "${LD_CHR}/" \
  --w-ld-chr "${WLD_CHR}/" \
  --out "${OUTDIR}/ldsc/conditioned_rg"

echo "Collect unconditioned and conditioned results"
Rscript "${ROOT}/scripts/04_collect_results.R" \
  "${OUTDIR}/ldsc/unconditioned_rg.log" \
  "${OUTDIR}/ldsc/conditioned_rg.log" \
  "${OUTDIR}/rg_results.tsv"

echo "Done. Results: ${OUTDIR}/rg_results.tsv"
