#!/usr/bin/env Rscript
# Read the unconditioned and conditioned LDSC --rg logs into one data.frame.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 04_collect_results.R <unconditioned.log> <conditioned.log> <out.tsv>")
}

parse_est <- function(line) {
  m <- regexec(":\\s*([-0-9.eE+]+)\\s*\\(([-0-9.eE+]+)\\)", line)
  hit <- regmatches(line, m)[[1]]
  if (length(hit) != 3) {
    return(c(est = NA_real_, se = NA_real_))
  }
  c(est = as.numeric(hit[2]), se = as.numeric(hit[3]))
}

parse_rg_log <- function(path, analysis) {
  if (!file.exists(path)) {
    stop("LDSC log not found: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  if (any(grepl("^ERROR", lines))) {
    stop("LDSC reported an error in ", path)
  }

  section <- function(start_pat, end_pat) {
    i <- grep(start_pat, lines)
    if (!length(i)) {
      stop("Could not find '", start_pat, "' in ", path)
    }
    i <- i[1]
    j <- grep(end_pat, lines)
    j <- j[j > i]
    end <- if (length(j)) j[1] - 1L else length(lines)
    lines[i:end]
  }

  rg <- section("^Genetic Correlation$", "^Summary of Genetic Correlation Results")

  line_matching <- function(block, pat) {
    hit <- grep(pat, block, value = TRUE)
    if (!length(hit)) stop("Missing '", pat, "' in ", path)
    hit[1]
  }

  rg_est <- parse_est(line_matching(rg, "^Genetic Correlation:"))
  p <- as.numeric(sub("^P:\\s*", "", line_matching(rg, "^P:")))

  data.frame(
    analysis = analysis,
    rg = unname(rg_est["est"]),
    rg_se = unname(rg_est["se"]),
    p = p,
    stringsAsFactors = FALSE
  )
}

res <- rbind(
  parse_rg_log(args[[1]], "unconditioned"),
  parse_rg_log(args[[2]], "conditioned")
)
write.table(res, args[[3]], sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", args[[3]])
print(res, row.names = FALSE)
