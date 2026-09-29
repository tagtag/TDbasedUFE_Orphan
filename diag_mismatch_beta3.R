#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## diag_mismatch_beta3.R
##
## diag_mismatch.R は orph_excess（ゼロ率の差の差）についてペアを崩した
## 陰性対照をかけ、3コホートすべてで sd 比 ≈ 1 を得た。しかし beta3 と
## orph_excess の相関はコホートで違う（-0.95 / -0.72 / -0.38）ので、
## beta3 そのものに同じ検定をかける必要がある。
##
##   実ペア    : normal(i) と disease(i) から beta3 を計算
##   崩したペア: normal(i) と disease(j), j != i から beta3 を計算
##
## マッチング変数（external2）は normal 側とコホート全体の normal から
## 作るので、disease を入れ替えても変わらない。したがって崩したペアでも
## 対照集合の構成は実ペアと同一であり、変わるのは disease 側だけである。
## 乱数種も患者ごとに固定し、実ペアと崩したペアで同じ抽出列を使う。
##
## 使い方
##   N_CORES=12 N_PERM=200 TRANSFORMS=rank NMIS=5 \
##   KALLISTO_ROOT=Revised/GSE244679 METADATA_DIR=Revised/metadata \
##   COHORT=GSE244679 Rscript diag_mismatch_beta3.R
##
## 判定
##   sd 比 ≈ 1 かつ 崩したペアでも up/down が同様に割れる
##      → beta3 の患者間変動は患者固有ではない。患者レベルの主張は不可。
##   sd 比 > 1 かつ 実ペアの方が極端
##      → その超過が患者固有。患者レベルの主張が部分的に残る。
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

CO      <- Sys.getenv("COHORT", names(COHORTS)[1])
N_PERM  <- as.integer(Sys.getenv("N_PERM", "200"))
TF      <- strsplit(Sys.getenv("TRANSFORMS", "rank"), ",")[[1]][1]
NMIS    <- as.integer(Sys.getenv("NMIS", "5"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "1"))
N_STRATA2 <- as.integer(Sys.getenv("N_STRATA2", "4"))
ZERO_RULE <- "or"

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

## 本体から統計部分をそのまま読み込む（実装を二重に持たない）
SRC <- if (file.exists("R/07_interaction_permutation.R"))
         "R/07_interaction_permutation.R" else "07_interaction_permutation.R"
src <- readLines(SRC)
eval(parse(text = paste(src[grep("^cell_coefs <- function", src):
                            (grep("^run_cohort_interaction <- function", src) - 1L)],
                        collapse = "\n")))

sheet <- read_sheet(CO)
dat   <- load_cohort_tpm(sheet, KALLISTO_ROOT, "abundance.tsv",
                         expected_features = N_TRANSCRIPTS)
is_o  <- orphan_index(dat$id, ORPHAN_ID_FILE, expected_n = N_ORPHAN_TRANSCRIPT)

pts <- unique(sheet$patient_id)
lb  <- function(pid, cond)
  sheet$library_id[sheet$patient_id == pid & sheet$condition == cond][1]
ln <- vapply(pts, lb, "", cond = "normal")
ld <- vapply(pts, lb, "", cond = "disease")
ok <- !is.na(ln) & !is.na(ld)
pts <- pts[ok]; ln <- ln[ok]; ld <- ld[ok]
np <- length(pts)
stopifnot(np >= 6L)

normal_mat   <- dat$tpm[, ln, drop = FALSE]
sum_normal   <- rowSums(normal_mat)
count_normal <- rowSums(normal_mat > 0)

## external2 の層別変数。normal 側だけで決まるので disease を替えても不変。
mvars <- function(i) {
  xn <- dat$tpm[, ln[i]]
  list((sum_normal - xn) / max(np - 1L, 1L),
       (count_normal - (xn > 0)) / max(np - 1L, 1L))
}

## 1 ペア分の beta3 と p
one <- function(i, j) {
  out <- run_patient(dat$tpm[, ln[i]], dat$tpm[, ld[j]], mvars(i), is_o,
                     TF, N_PERM, c(N_STRATA, N_STRATA2),
                     BASE_SEED + i, "external2")   # 種は i で固定
  if (is.null(out)) return(NULL)
  data.frame(patient_id = pts[i], partner = pts[j], kind = if (i == j) "real" else "mismatch",
             beta3 = out$beta3, z_mad = out$z_mad,
             p_up = out$p_up, p_down = out$p_down,
             zero_gap_asym = out$zero_gap_asym,
             stringsAsFactors = FALSE)
}

set.seed(BASE_SEED)
jobs <- do.call(rbind, lapply(seq_len(np), function(i) {
  j <- c(i, sample(setdiff(seq_len(np), i), min(NMIS, np - 1L)))
  data.frame(i = i, j = j)
}))

msg(CO, " / 患者 ", np, " 人 / 崩したペア 各 ", NMIS,
    " 通り / 計 ", nrow(jobs), " 回 / N_PERM = ", N_PERM,
    " / transform = ", TF, " / cores = ", N_CORES)

f <- function(k) one(jobs$i[k], jobs$j[k])
lst <- if (N_CORES > 1L) {
  parallel::mclapply(seq_len(nrow(jobs)), f, mc.cores = N_CORES)
} else {
  lapply(seq_len(nrow(jobs)), f)
}
bad <- which(!vapply(lst, function(x) is.null(x) || is.data.frame(x), logical(1)))
if (length(bad)) {
  msg("失敗: ", paste(utils::head(as.character(lst[[bad[1]]]), 1), collapse = " "))
  stop(length(bad), " 件でエラー。N_CORES=1 で原因を確認してください。", call. = FALSE)
}
res <- do.call(rbind, Filter(Negate(is.null), lst))

R <- res[res$kind == "real", ]
M <- res[res$kind == "mismatch", ]
stopifnot(nrow(R) >= 6L, nrow(M) >= 6L)

q <- function(v) sprintf("中位 %+.4g  四分位 %+.4g / %+.4g  範囲 %+.4g 〜 %+.4g",
                         median(v), quantile(v, .25), quantile(v, .75), min(v), max(v))

hr(); msg("beta3 の分布（", TF, "）")
msg("  実ペア      ", q(R$beta3))
msg("  崩したペア  ", q(M$beta3))
msg("")
msg(sprintf("  sd   実 %.4g   崩し %.4g   比 %.2f", sd(R$beta3), sd(M$beta3),
            sd(R$beta3) / sd(M$beta3)))
msg(sprintf("  IQR  実 %.4g   崩し %.4g   比 %.2f", IQR(R$beta3), IQR(M$beta3),
            IQR(R$beta3) / IQR(M$beta3)))

hr(); msg("向きの内訳（生の片側 p < 0.05 で判定。BH は患者数が違うので使わない）")
dir_of <- function(d) ifelse(d$p_up < 0.05, "up", ifelse(d$p_down < 0.05, "down", "none"))
tr <- table(factor(dir_of(R), levels = c("up", "none", "down")))
tm <- table(factor(dir_of(M), levels = c("up", "none", "down")))
msg(sprintf("  実ペア      up %d / none %d / down %d   （n = %d）",
            tr["up"], tr["none"], tr["down"], nrow(R)))
msg(sprintf("  崩したペア  up %d / none %d / down %d   （n = %d）",
            tm["up"], tm["none"], tm["down"], nrow(M)))
msg(sprintf("  有意割合  実 %.3f   崩し %.3f", 1 - tr["none"]/nrow(R), 1 - tm["none"]/nrow(M)))
msg(sprintf("  up の割合（有意な中で） 実 %.3f   崩し %.3f",
            tr["up"]/max(tr["up"]+tr["down"], 1), tm["up"]/max(tm["up"]+tm["down"], 1)))

hr(); msg("患者ごとの経験的 p（自分の崩したペア分布に対して）")
pv <- vapply(seq_len(np), function(i) {
  r <- R$beta3[R$patient_id == pts[i]]
  m <- M$beta3[M$patient_id == pts[i]]
  if (!length(r) || length(m) < 2L) return(NA_real_)
  min((1 + sum(m >= r)) / (length(m) + 1), (1 + sum(m <= r)) / (length(m) + 1))
}, 0)
pv <- pv[is.finite(pv)]
msg(sprintf("  最小可能 p = %.3f（崩したペア %d 通りなので）", 1/(NMIS+1), NMIS))
msg(sprintf("  最小値を取った患者 %d / %d  （帰無なら約 %.1f 人）",
            sum(pv <= 1/(NMIS+1) + 1e-9), length(pv), length(pv) * 2/(NMIS+1)))
msg(sprintf("  p の中位 %.3f", median(pv)))

hr(); msg("参考: beta3 と zero_gap_asym の相関（実ペア）")
msg(sprintf("  Spearman %.3f", cor(R$beta3, R$zero_gap_asym, method = "spearman")))

out <- paste0("mismatch_beta3_", CO, "_", TF, ".csv")
write.csv(res, out, row.names = FALSE)
hr(); msg("書き出し: ", out)
msg("")
msg("読み方")
msg("  sd 比 ≈ 1 で、崩したペアでも有意割合と up/down の比が実ペアと同様なら、")
msg("  beta3 の患者間変動は患者固有ではない。観測されているのはコホート全体の")
msg("  条件差と、ライブラリ個体のばらつきである。患者レベルの主張はできない。")
msg("")
msg("  sd 比が 1 より明確に大きく、実ペアの有意割合が崩したペアを上回るなら、")
msg("  その超過が患者固有の成分である。")
msg("")
msg("  崩したペア数が ", NMIS, " なので患者ごとの p は ", sprintf("%.3f", 1/(NMIS+1)),
    " より小さくならない。")
msg("  全体の分布（sd 比・有意割合）で判断すること。")
