## ---------------------------------------------------------------------------
## 07_validate_manuscript_counts.R
## Optional final check: compare generated direction counts with those reported
## in the current manuscript. This script fails loudly if a change in metadata,
## reference, or code alters the reported Tables 2-5.
## ---------------------------------------------------------------------------

source("R/config.R")

read_required <- function(path) {
  if (!file.exists(path)) stop("missing result file: ", path, "\nRun the analysis first.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

## Expected (up, down, none) for included/excluded zero rules.
expected <- list(
  table2 = list(
    GSE244679 = c(included_up=21, included_down=3, included_none=0,
                  excluded_up=21, excluded_down=3, excluded_none=0),
    GSE127165 = c(included_up=26, included_down=27, included_none=4,
                  excluded_up=23, excluded_down=23, excluded_none=11),
    GSE144269 = c(included_up=19, included_down=51, included_none=0,
                  excluded_up=19, excluded_down=50, excluded_none=1),
    GSE40419 = c(included_up=69, included_down=0, included_none=0,
                 excluded_up=69, excluded_down=0, excluded_none=0)
  ),
  table3 = list(
    GSE244679 = c(included_up=19, included_down=5, included_none=0,
                  excluded_up=19, excluded_down=5, excluded_none=0),
    GSE127165 = c(included_up=25, included_down=31, included_none=1,
                  excluded_up=25, excluded_down=28, excluded_none=4),
    GSE144269 = c(included_up=19, included_down=51, included_none=0,
                  excluded_up=19, excluded_down=51, excluded_none=0),
    GSE40419 = c(included_up=69, included_down=0, included_none=0,
                 excluded_up=69, excluded_down=0, excluded_none=0)
  ),
  table4 = list(
    GSE244679 = c(included_up=16, included_down=6, included_none=2,
                  excluded_up=16, excluded_down=6, excluded_none=2),
    GSE127165 = c(included_up=5, included_down=14, included_none=38,
                  excluded_up=3, excluded_down=13, excluded_none=41),
    GSE144269 = c(included_up=15, included_down=30, included_none=25,
                  excluded_up=12, excluded_down=23, excluded_none=35),
    GSE40419 = c(included_up=1, included_down=68, included_none=0,
                 excluded_up=0, excluded_down=69, excluded_none=0)
  )
)

long <- read_required(file.path(RESULTS_DIR, "tables_2_3_4_counts_long.csv"))
problems <- character()
for (tb in names(expected)) {
  for (cohort in names(expected[[tb]])) {
    got <- long[long$table == tb & long$cohort == cohort, ]
    if (nrow(got) != 2L) {
      problems <- c(problems, paste(tb, cohort, "missing included/excluded result"))
      next
    }
    inc <- got[got$zero == "included", ]
    exc <- got[got$zero == "excluded", ]
    gv <- c(included_up=inc$up, included_down=inc$down, included_none=inc$none,
            excluded_up=exc$up, excluded_down=exc$down, excluded_none=exc$none)
    ev <- expected[[tb]][[cohort]]
    if (!identical(as.numeric(gv), as.numeric(ev))) {
      problems <- c(problems,
                    paste(tb, cohort, "expected", paste(ev, collapse="/"),
                          "got", paste(gv, collapse="/")))
    }
  }
}

if (length(problems)) {
  cat(paste("-", problems), sep = "\n")
  stop("\nGenerated counts differ from the manuscript.")
}
cat("Tables 2-4 match the manuscript counts.\n")

## Table 5 is summarized as strings; verify the manuscript summary exactly.
t5 <- read_required(file.path(RESULTS_DIR, "table5.csv"))
expected_t5 <- data.frame(
  GEO_ID = rep(c("GSE244679","GSE127165","GSE144269","GSE40419"), each=2),
  transcript_set = rep(c("Orphan","Control"), 4),
  disease_gt_normal = c("16","2 (1-3)","5","3 (1-4)","15","4.5 (3-6)","1","68 (67-68)"),
  disease_lt_normal = c("6","1 (0-3)","14","1 (0-4)","30","18 (14-23)","68","0 (0-0)"),
  no_direction = c("2","21 (19-23)","38","53 (50-55)","25","47 (42-51)","0","1 (1-2)"),
  stringsAsFactors = FALSE
)
for (nm in colnames(expected_t5)) t5[[nm]] <- as.character(t5[[nm]])
if (!identical(t5[, colnames(expected_t5)], expected_t5)) {
  warning("Table 5 summary differs from the current manuscript; inspect random seeds and metadata.")
} else {
  cat("Table 5 matches the manuscript summary.\n")
}
