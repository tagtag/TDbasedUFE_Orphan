source("R/00_functions.R")

cat("1. rank transformation / ties ... ")
r <- transform_expression(c(0, 0, 2, 5), "rank")
stopifnot(identical(as.numeric(r), c(1.5, 1.5, 3, 4)))
cat("ok\n")

cat("2. standardization ... ")
x <- transform_expression(1:10, "scale")
stopifnot(abs(mean(x)) < 1e-12, abs(var(x) - 1) < 1e-12)
cat("ok\n")

cat("3. sign convention ... ")
## Three target features move upward in disease relative to the background.
n <- c(1, 2, 3, 10, 20, 30)
d <- c(4, 5, 6, 10, 20, 30)
target <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
ans <- paired_test(n, d, target, transform = "rank", zero_rule = "none")
stopifnot(ans$mean_d < 0)
cat("ok\n")

cat("4. OR zero rule ... ")
n <- c(0, 1, 2, 3, 4, 5)
d <- c(1, 0, 2, 3, 4, 5)
ans <- paired_test(n, d, rep(TRUE, 6), transform = "rank", zero_rule = "or")
stopifnot(ans$n_kept == 6)
cat("ok\n")

cat("5. matched control size ... ")
set.seed(1)
n <- c(rep(0, 20), seq_len(80))
d <- c(rep(0, 20), seq_len(80) + rep(c(-1, 1), 40))
is_orphan <- rep(FALSE, 100)
is_orphan[c(1:5, 21:30)] <- TRUE
ctrl <- draw_control_set(n, d, is_orphan, n_strata = 5)
stopifnot(sum(ctrl) == sum(is_orphan), !any(ctrl & is_orphan))
cat("ok\n")

cat("All smoke tests passed.\n")
