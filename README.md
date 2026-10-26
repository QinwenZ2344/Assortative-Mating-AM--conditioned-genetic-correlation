# Assortative Mating (AM)-conditioned genetic correlation

Pipeline for the genetic correlation estimation between two GWAS after removing partner choice index component using our previous software [GSUB](https://github.com/qlu-lab/GSUB).


## Dependencies

- [GSUB](https://github.com/qlu-lab/GSUB) and the R packages it loads (`tidyverse`, `data.table`, `optparse`, `GenomicSEM`)
- [LDSC](https://github.com/bulik/ldsc) (`munge_sumstats.py`, `ldsc.py`), called from `LDSC_DIR`
- R with `data.table`
- Python environment with the packages in the LDSC `environment.yml` (numpy, pandas, scipy)
- European LD scores (not included here; download them from LDSC)
- HapMap3 SNP list and 1000 Genomes MAF reference bundled in `resources/`

Summary statistics must follow the [GenomicSEM sumstats format](https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-Resources-and-Key-Information) (`SNP`, `A1`, `A2`, `BETA`, `P`, and usually `Z` and `N`).

## Setup

Clone GSUB and LDSC, then install the LDSC Python environment. The pipeline calls `ldsc.py` and `munge_sumstats.py` from `LDSC_DIR`. It does not use an environment module.

Download the European LD scores from the [LDSC wiki](https://github.com/bulik/ldsc) (`eur_w_ld_chr`). This repository already includes the two reference files GSUB and munge need:

- `resources/w_hm3.snplist` — HapMap3 SNPs used to align alleles
- `resources/reference.1000G.maf.0.005.txt.gz` — 1000 Genomes SNPs with MAF > 0.005 (`SNP`, `MAF`, `A1`, `A2`)


```bash
git clone https://github.com/qlu-lab/GSUB.git
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create -f environment.yml
conda activate ldsc
cd ..
cp config.example.sh config.sh
```

In `config.sh`, set `GSUB_DIR` and `LDSC_DIR` to those clones. Activate the `ldsc` conda environment before `run_pipeline.sh` so `python` can import the LDSC packages.

The other values a user normally changes are the three sumstats paths, the GSUB settings for each file, and `OUTDIR`.

For a linear GWAS of a continuous trait, the settings used by this pipeline are:

| argument | value |
| --- | --- |
| `se.logit` | `FALSE` |
| `OLS` | `TRUE` |
| `linprob` | `FALSE` |
| `sample.prev` | `NA` |
| `population.prev` | `NA` |
| `N` | GWAS sample size, or `NA` to read the `N` column |

Other designs (logistic, linear probability) are documented in the [GSUB running guide](https://github.com/qlu-lab/GSUB/wiki/2.-Running-GSUB). Optional GSUB arguments (`correction`, `info.filter`, `maf.filter`, `keep.indel`) are also in `config.sh`.

## Run

```bash
bash run_pipeline.sh config.sh
```

`config.sh` is the default setting assuming two GWAS are quantitative traits. User should download reference files and specify their path.

## Output

```
$OUTDIR/
  traitA/                         GSUB output for trait A
  traitB/                         GSUB output for trait B
  conditioned/traitA.gamma2.txt.gz
  conditioned/traitB.gamma2.txt.gz
  ldsc/unconditioned_rg.log
  ldsc/conditioned_rg.log
  rg_results.tsv
```

`rg_results.tsv` has one row per analysis:

| column | meaning |
| --- | --- |
| `analysis` | `unconditioned` (original GWAS) or `conditioned` (gamma2) |
| `rg` | genetic correlation |
| `rg_se` | standard error of `rg` |
| `p` | p-value for `rg` |

`N` on the conditioned sumstats is the GSUB effective sample size for gamma2 (`n.gamma2`), not the original GWAS `N`. GSUB sets that value to missing when MAF is outside 0.1–0.4, and munge drops those SNPs.

## What this repository does not include

GSUB, LDSC, the LD score files, and the GWAS summary statistics stay outside this repository. Point `config.sh` at local copies. The HapMap3 SNP list and the 1000 Genomes reference are included under `resources/`.
