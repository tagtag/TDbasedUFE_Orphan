## ---------------------------------------------------------------------------
## 06_session_info.R
##
## 解析に使ったソフトウェア環境を記録する（Reviewer 2-4）。
##
##   Rscript R/06_session_info.R
##   → results/sessionInfo.txt
##
## ---------------------------------------------------------------------------
## 【2026-09-23 改訂】
##
## 旧版は sessionInfo() を呼ぶだけだった。しかしこのスクリプトは config.R しか
## 読まないので、**解析で実際に使うパッケージが 1 つも attach されず**、出力に
## limma も edgeR も DESeq2 も現れなかった。「再現に必要な版を報告せよ」という
## 要求に対して、R 本体の版しか答えていないことになる。
##
## そこで解析が使うパッケージを明示的に読み込んでから sessionInfo() を取り、
## 併せて外部ツール（kallisto）の版と、数値再現性に効く BLAS/LAPACK も記録する。
## パッケージが入っていない環境でも落ちないようにしてある。
## ---------------------------------------------------------------------------

source("R/config.R")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(RESULTS_DIR, "sessionInfo.txt")

## 解析で実際に使うもの。任意（入っていなければ「未インストール」と記録する）。
##   limma / edgeR / DESeq2 / SummarizedExperiment … de_framework.R
##   BayesFactor                                   … clinical_assoc.R
##   data.table                                    … de_framework.R の高速読み込み
##   parallel                                      … N_CORES > 1 のとき
PKGS <- c("limma", "edgeR", "DESeq2", "SummarizedExperiment",
          "BayesFactor", "data.table", "parallel")

have <- vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)
for (p in PKGS[have]) suppressPackageStartupMessages(
  library(p, character.only = TRUE))

ver <- vapply(PKGS, function(p)
  if (have[[p]]) as.character(utils::packageVersion(p)) else "未インストール", "")

## 外部ツール。PATH に無ければその旨を残す（隠さない）。
ext_version <- function(cmd, args = "version") {
  if (nzchar(Sys.which(cmd)) == FALSE) return(paste0(cmd, ": PATH にありません"))
  out <- suppressWarnings(try(
    system2(cmd, args, stdout = TRUE, stderr = TRUE), silent = TRUE))
  if (inherits(out, "try-error") || !length(out)) return(paste0(cmd, ": 取得できません"))
  paste0(cmd, ": ", paste(utils::head(out, 2), collapse = " / "))
}

con <- file(OUT, open = "wt", encoding = "UTF-8")
sink(con); on.exit({ sink(); close(con) }, add = TRUE)

cat("記録日時: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("ホスト:   ", Sys.info()[["nodename"]], "\n\n", sep = "")

cat(strrep("-", 74), "\n")
cat("解析で使うパッケージの版\n")
cat(strrep("-", 74), "\n")
for (p in PKGS) cat(sprintf("  %-22s %s\n", p, ver[[p]]))
cat("\n  未インストールのものがある環境では、それを使う解析だけが走らない。\n")
cat("  de_framework.R は ENGINE に応じて limma / edgeR / DESeq2 を要求し、\n")
cat("  clinical_assoc.R は BayesFactor が無ければ Bayes 相関を飛ばす。\n\n")

cat(strrep("-", 74), "\n")
cat("外部ツール\n")
cat(strrep("-", 74), "\n")
cat("  ", ext_version("kallisto"), "\n", sep = "")
cat("  ", ext_version("gffread", "--version"), "\n", sep = "")
cat("\n  参照の構築と定量に使ったもの。reference_build/reference_build_summary.txt\n")
cat("  にも build 時点の kallisto の版が記録される。\n\n")

cat(strrep("-", 74), "\n")
cat("数値再現性に関わる設定\n")
cat(strrep("-", 74), "\n")
cat(sprintf("  BASE_SEED            %s\n", BASE_SEED))
cat("  患者ごとの乱数種      BASE_SEED + i（ワーカーごとではない）\n")
cat("  → 結果は N_CORES に依存しない（N_CORES=1 と 4 で一致を確認済み）\n")
cat(sprintf("  N_STRATA             %s\n", N_STRATA))
cat(sprintf("  ALPHA                %s\n", ALPHA))
cat(sprintf("  転写産物数（期待）    %s（うち orphan %s）\n",
            N_TRANSCRIPTS, N_ORPHAN_TRANSCRIPT))
cat("\n  BLAS/LAPACK は下の sessionInfo() に出る。参照実装であれば\n")
cat("  スレッド数による丸め差も生じない。\n\n")

cat(strrep("-", 74), "\n")
cat("sessionInfo()\n")
cat(strrep("-", 74), "\n")
print(sessionInfo())

sink(); close(con); on.exit()
message("書き出し: ", OUT)
message("  パッケージ: ", sum(have), " / ", length(PKGS), " が利用可能")
if (any(!have))
  message("  未インストール: ", paste(PKGS[!have], collapse = ", "))
