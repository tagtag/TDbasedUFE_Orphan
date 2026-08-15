## ---------------------------------------------------------------------------
## 05_figure2.R
## Reproduce manuscript Figure 2 from the standardized-expression analysis.
## Default: Table 2 with zero-abundance transcripts INCLUDED, which gives the
## 21/3, 26/27, 19/51, and 69/0 counts shown in the manuscript figure.
##
## Optional command line:
##   Rscript R/05_figure2.R table4 excluded figures/figure_rank.pdf
## ---------------------------------------------------------------------------

source("R/config.R")

args <- commandArgs(trailingOnly = TRUE)
TABLE <- if (length(args) >= 1L) args[1] else "table2"
ZERO  <- if (length(args) >= 2L) args[2] else "included"
OUT   <- if (length(args) >= 3L) args[3] else file.path(FIGURE_DIR, "figure2.pdf")

if (!TABLE %in% c("table2", "table4")) stop("TABLE must be table2 or table4")
if (!ZERO %in% c("included", "excluded")) stop("ZERO must be included or excluded")

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

COL_UP   <- "#C1454F"
COL_DOWN <- "#3C6FB0"
COL_NONE <- "#B8B8B8"

load_res <- function(cohort) {
  f <- file.path(RESULTS_DIR,
                 sprintf("perpatient_%s_%s_zero_%s.csv", TABLE, cohort, ZERO))
  if (!file.exists(f)) stop("run R/02_patient_level_analysis.R first: ", f)
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}

xlabel <- if (TABLE == "table2") {
  "Mean standardized expression difference (normal - disease)"
} else {
  "Mean within-sample rank difference (normal - disease)"
}

pdf(OUT, width = 9, height = 11, useDingbats = FALSE)
layout(matrix(seq_len(2L * length(COHORTS)), ncol = 2, byrow = TRUE),
       widths = c(1.15, 1))

for (cohort in names(COHORTS)) {
  res <- load_res(cohort)
  d <- res$mean_d
  col <- ifelse(res$direction == "up", COL_UP,
         ifelse(res$direction == "down", COL_DOWN, COL_NONE))

  par(mar = c(4.4, 1.2, 2.6, 1.0), cex = 0.85)
  plot(d, jitter(rep(0, length(d)), amount = 0.8),
       xlim = range(pretty(d), finite = TRUE), ylim = c(-2, 2),
       pch = 19, cex = 0.9, col = col,
       xlab = xlabel, ylab = "", yaxt = "n", bty = "n",
       main = cohort)
  abline(v = 0, col = "red", lwd = 1.2)

  cnt <- c(
    up   = sum(res$direction == "up"),
    none = sum(res$direction == "none"),
    down = sum(res$direction == "down")
  )
  n <- sum(cnt)

  par(mar = c(4.4, 1.0, 2.6, 1.0), cex = 0.85)
  barplot(matrix(cnt, ncol = 1), horiz = TRUE, beside = FALSE,
          col = c(COL_UP, COL_NONE, COL_DOWN), border = NA,
          xlim = c(0, n), xlab = "Number of patients", axes = TRUE)
  pos <- cumsum(cnt) - cnt / 2
  for (k in seq_along(cnt)) {
    if (cnt[k] > 0) {
      txt_col <- if (names(cnt)[k] == "none") "black" else "white"
      text(pos[k], 0.7, sprintf("%d\n(%.1f%%)", cnt[k], 100 * cnt[k] / n),
           cex = 0.75, col = txt_col, font = 2)
    }
  }
}
dev.off()
cat("figure written to ", OUT, "\n", sep = "")
