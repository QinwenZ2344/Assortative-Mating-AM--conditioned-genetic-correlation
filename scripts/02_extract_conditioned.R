#!/usr/bin/env Rscript
# Pull GSUB gamma2 (outcome after regressing out the factor) into an LDSC-ready sumstats file.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 02_extract_conditioned.R <gsub_analytical_all.txt.gz> <out.txt.gz>")
}

infile <- args[[1]]
outfile <- args[[2]]

suppressPackageStartupMessages({
  library(data.table)
})

need <- c(
  "SNP", "CHR", "BP", "MAF", "A1", "A2",
  "gamma2", "se.gamma2", "z.gamma2", "p.gamma2", "n.gamma2"
)
dt <- fread(infile)
missing <- setdiff(need, names(dt))
if (length(missing)) {
  stop("Missing columns in ", infile, ": ", paste(missing, collapse = ", "))
}

out <- dt[, .(
  SNP, CHR, BP, MAF, A1, A2,
  BETA = gamma2,
  SE = se.gamma2,
  Z = z.gamma2,
  P = p.gamma2,
  N = n.gamma2
)]
fwrite(out, outfile, sep = "\t")
message("Wrote ", outfile, " (", nrow(out), " SNPs; ",
        sum(!is.na(out$N)), " with N)")
