#!/usr/bin/env Rscript
## 抽出ばらつきが回帰の残差のどれだけを占めるかを、既存の CSV から出す。
## diag_slope.R を回し直す必要はない（slope_<COHORT>.csv に全部入っている）。
##
##   Rscript check_draw_noise.R              # カレントディレクトリを見る
##   Rscript check_draw_noise.R <ディレクトリ>  # CSV のある場所を指定
##   NSET=200 Rscript check_draw_noise.R     # NSET は実行時の値に合わせる
##
## CSV は diag_slope.R を走らせたディレクトリ（= run_verify.sh を起動した
## リポジトリのルート。logs/ の隣）に書き出されている。

args <- commandArgs(trailingOnly = TRUE)
DIR  <- if (length(args) >= 1) args[1] else "."
NSET <- as.integer(Sys.getenv("NSET", "200"))
cos  <- c("GSE244679", "GSE127165", "GSE144269")

cat(sprintf("探す場所: %s   （NSET = %d）\n\n", normalizePath(DIR, mustWork = FALSE), NSET))

cat(sprintf("%-12s %10s %10s %10s %10s %13s\n",
            "コホート", "残差sd", "抽出SE", "残差比", "分散比", "切片の回帰SE"))
found <- 0L
for (CO in cos) {
  f <- file.path(DIR, sprintf("slope_%s.csv", CO))
  if (!file.exists(f)) { cat(sprintf("%-12s （%s がありません）\n", CO, basename(f))); next }
  found <- found + 1L
  d <- read.csv(f)
  need <- c("d_zr_all", "d_zr_orphan", "d_zr_control", "d_zr_control_sd")
  if (!all(need %in% names(d))) {
    cat(sprintf("%-12s （列が足りません: %s）\n", CO,
                paste(setdiff(need, names(d)), collapse = ", ")))
    next
  }
  d$excess <- d$d_zr_orphan - d$d_zr_control
  fe <- lm(excess ~ d_zr_all, data = d)
  s_resid <- summary(fe)$sigma
  s_draw  <- mean(d$d_zr_control_sd) / sqrt(NSET)
  se_int  <- summary(fe)$coefficients[1, 2]
  cat(sprintf("%-12s %10.5f %10.6f %9.2f%% %9.3f%% %13.6f\n",
              CO, s_resid, s_draw, 100*s_draw/s_resid,
              100*(s_draw/s_resid)^2, se_int))
}

if (found == 0L) {
  cat("\n")
  cat("CSV が 1 つも見つかりません。diag_slope.R はカレントディレクトリに\n")
  cat("slope_<COHORT>.csv を書き出すので、run_verify.sh を起動した\n")
  cat("ディレクトリ（logs/ の隣）で実行するか、その場所を引数で渡してください。\n")
  cat("  例: Rscript check_draw_noise.R ~/orphan-repo\n")
  cat("  探す: find ~ -name 'slope_GSE*.csv' 2>/dev/null\n")
  quit(save = "no", status = 1)
}

cat("\n")
cat("読み方\n")
cat("  抽出ノイズは excess の測定誤差なので、回帰の残差にそのまま入り、\n")
cat("  切片の標準誤差に既に伝播している。**別途足し込むのは二重計上。**\n")
cat("  見るべきは「分散比」= 残差分散のうち抽出ノイズが占める割合。\n")
cat("  数 % 以下なら NSET は足りている。5 % を超えるなら NSET を増やす。\n")
cat("  「残差比」は参考（sd どうしの比）。判断には分散比を使う。\n")
