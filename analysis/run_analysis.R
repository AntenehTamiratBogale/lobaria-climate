args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script) != 1L) stop("Run with Rscript analysis/run_analysis.R")
root <- dirname(dirname(normalizePath(script)))
dir.create(file.path(root, "results"), showWarnings = FALSE)
setwd(file.path(root, "data"))
set.seed(123)
pdf(file.path(root, "results", "analysis_plots.pdf"), width = 12, height = 9)
tryCatch(
  source(file.path(root, "analysis", "lichen_data_analysis.R"), echo = FALSE),
  finally = {
    dev.off()
    capture.output(sessionInfo(), file = file.path(root, "results", "sessionInfo.txt"))
  }
)
