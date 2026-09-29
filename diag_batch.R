#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## diag_batch.R
##
## normal と disease が別のバッチに分かれていないかを調べる。
## 分かれていれば、疾患差とバッチ差は原理的に区別できない。
##
## GSE40419 では normal が ERX1359 系、disease が ERX140 系で相互に0本と
## 判明している。それが技術指標（深さ・マップ率・ゼロ率）の系統差として
## 出ているかを確認し、他コホートと比較する。
##
## 使い方（リポジトリのルートで、コホートごとに）
##   KALLISTO_ROOT=Revised/GSE40419 METADATA_DIR=Revised/metadata \
##   COHORT=GSE40419 Rscript diag_batch.R
##
## 判定
##   accession の範囲が normal と disease で重なっていない  → バッチ交絡
##   マップ率やゼロ率が条件間で系統的に違う                 → バッチ差の実体
##   交絡しているコホートでは、疾患差の主張はできない。
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")
CO <- Sys.getenv("COHORT", names(COHORTS)[1])

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

sheet <- read_sheet(CO)
dat   <- load_cohort_tpm(sheet, KALLISTO_ROOT, "abundance.tsv",
                         expected_features = N_TRANSCRIPTS)
is_o  <- orphan_index(dat$id, ORPHAN_ID_FILE, expected_n = N_ORPHAN_TRANSCRIPT)

gj <- function(dir, key) {
  p <- file.path(KALLISTO_ROOT, dir, "run_info.json")
  if (!file.exists(p)) return(NA_real_)
  s <- paste(readLines(p, warn = FALSE), collapse = " ")
  m <- regmatches(s, regexpr(paste0('"', key, '":\\s*[0-9.eE+-]+'), s))
  if (!length(m)) NA_real_ else as.numeric(sub('.*:\\s*', '', m))
}

sheet$depth  <- vapply(sheet$directory, gj, 0, key = "n_processed")
sheet$pmap   <- vapply(sheet$directory, gj, 0, key = "p_pseudoaligned")
sheet$zr     <- vapply(sheet$library_id, function(l) mean(dat$tpm[, l] == 0), 0)
sheet$zr_o   <- vapply(sheet$library_id, function(l) mean(dat$tpm[is_o, l] == 0), 0)
## accession の番号部分
sheet$acc_n  <- suppressWarnings(as.numeric(gsub("[^0-9]", "", sheet$library_id)))

N <- sheet[sheet$condition == "normal", ]
D <- sheet[sheet$condition == "disease", ]

hr(); msg(CO, " / normal ", nrow(N), " 本 / disease ", nrow(D), " 本")

## ---------------------------------------------------------------------------
## (1) accession の範囲が重なるか
## ---------------------------------------------------------------------------
hr(); msg("(1) run accession の番号範囲")
msg(sprintf("  normal   %.0f 〜 %.0f", min(N$acc_n), max(N$acc_n)))
msg(sprintf("  disease  %.0f 〜 %.0f", min(D$acc_n), max(D$acc_n)))
ov <- max(0, min(max(N$acc_n), max(D$acc_n)) - max(min(N$acc_n), min(D$acc_n)))
n_in_d <- sum(N$acc_n >= min(D$acc_n) & N$acc_n <= max(D$acc_n))
d_in_n <- sum(D$acc_n >= min(N$acc_n) & D$acc_n <= max(N$acc_n))
msg(sprintf("  disease の範囲に入る normal: %d / %d", n_in_d, nrow(N)))
msg(sprintf("  normal の範囲に入る disease: %d / %d", d_in_n, nrow(D)))
if (n_in_d == 0L && d_in_n == 0L) {
  msg("  → **完全に分離。バッチ交絡の疑いが強い。**")
} else if (n_in_d > 0.8 * nrow(N)) {
  msg("  → 入り混じっている。番号の上では交絡していない。")
} else {
  msg("  → 部分的に分離。内訳を確認すること。")
}

## ---------------------------------------------------------------------------
## (2) 技術指標の条件間比較
## ---------------------------------------------------------------------------
hr(); msg("(2) 技術指標（normal / disease）")
cmp <- function(lab, vn, vd, fmt = "%.4g") {
  kn <- is.finite(vn); kd <- is.finite(vd)
  if (sum(kn) < 3L || sum(kd) < 3L) { msg(sprintf("  %-16s 取得不可", lab)); return(invisible()) }
  w <- suppressWarnings(wilcox.test(vn[kn], vd[kd]))
  msg(sprintf(paste0("  %-16s N 中位 ", fmt, "  D 中位 ", fmt,
                     "   Wilcoxon p = %.3g%s"),
              lab, median(vn[kn]), median(vd[kd]), w$p.value,
              if (w$p.value < 0.01) "  **" else ""))
}
cmp("深さ",           N$depth, D$depth)
cmp("マップ率 %",     N$pmap,  D$pmap, "%.1f")
cmp("ゼロ率 全体",    N$zr,    D$zr,   "%.4f")
cmp("ゼロ率 orphan",  N$zr_o,  D$zr_o, "%.4f")

## ---------------------------------------------------------------------------
## (3) 極端に浅いライブラリ
## ---------------------------------------------------------------------------
hr(); msg("(3) 深さの外れ値")
thr <- median(sheet$depth, na.rm = TRUE) / 5
bad <- sheet[is.finite(sheet$depth) & sheet$depth < thr, ]
if (nrow(bad)) {
  msg(sprintf("  中位の 1/5 未満のライブラリ %d 本（閾値 %.3g）:", nrow(bad), thr))
  b <- bad[order(bad$depth), ]
  for (i in seq_len(min(nrow(b), 8L)))
    msg(sprintf("    %-14s %-8s patient %-8s depth %.3g  マップ率 %.1f%%  ゼロ率 %.3f",
                b$library_id[i], b$condition[i], b$patient_id[i],
                b$depth[i], b$pmap[i], b$zr[i]))
  if (nrow(b) > 8L) msg("    ... 他 ", nrow(b) - 8L, " 本")
  msg("  該当患者は解析から除くか、深さを揃えたサブサンプリングが必要。")
} else {
  msg("  なし（全ライブラリが中位の 1/5 以上）")
}

## ---------------------------------------------------------------------------
## (4) 条件を混ぜた対照実験 — normal どうしを疑似ペアにする
## ---------------------------------------------------------------------------
## 同一条件（normal のみ）の中で患者を任意に2つずつ組にして「疾患効果なし」の
## 対照とする。ここで有意な患者が多数出れば、検定はバッチや個体差に反応して
## いるということになり、疾患差の主張は成り立たない。
hr(); msg("(4) normal どうしの疑似ペア（疾患効果ゼロの陰性対照）")
if (nrow(N) >= 6L) {
  set.seed(BASE_SEED)
  o <- sample(nrow(N))
  k <- floor(length(o) / 2)
  a <- N$library_id[o[seq_len(k)]]; b <- N$library_id[o[k + seq_len(k)]]
  zr_gap <- vapply(seq_len(k), function(i) {
    vA <- dat$tpm[, a[i]]; vB <- dat$tpm[, b[i]]
    (mean(vB[is_o] == 0) - mean(vB[!is_o] == 0)) -
      (mean(vA[is_o] == 0) - mean(vA[!is_o] == 0))
  }, 0)
  msg(sprintf("  %d 組。orph_excess 相当の中位 %+.4f   範囲 %+.4f 〜 %+.4f",
              k, median(zr_gap), min(zr_gap), max(zr_gap)))
  msg(sprintf("  |値| > 0.05 の組: %d / %d", sum(abs(zr_gap) > 0.05), k))
  msg("  実データの orph_excess と同程度の広がりがあれば、その量は疾患に")
  msg("  固有ではなく個体間・バッチ間のばらつきを測っている。")
} else {
  msg("  normal が少なすぎて実施できません。")
}

hr()
msg("この診断の結論の使い方")
msg("  (1) が完全分離なら、そのコホートで疾患差を主張してはならない。")
msg("  (2) でマップ率やゼロ率に系統差があれば、それがバッチ差の実体。")
msg("  (4) は疾患効果ゼロの陰性対照。実データと同程度なら現象は疾患由来でない。")
