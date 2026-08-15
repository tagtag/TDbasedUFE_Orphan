## Record the software environment used for a reproducible analysis run.
source("R/config.R")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
capture.output(sessionInfo(), file = file.path(RESULTS_DIR, "sessionInfo.txt"))
cat("sessionInfo written to ", file.path(RESULTS_DIR, "sessionInfo.txt"), "\n", sep = "")
