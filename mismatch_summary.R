#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## mismatch_summary.R
##
## `diag_mismatch_slope.R` が各コホートに書き出した mismatch_slope_<CO>.csv を
## 読んで、**判定に必要な 1 枚の表**にまとめる。
##
## 見るべきものは 2 つ。
##
##   (1) コホートごと。崩したペアの切片が 0 付近か、本物と同程度か。
##   (2) **プールした切片。** 論文の主張は 3 コホートをプールした
##       +0.0544 ± 0.0059（z = 9.2）である。崩したペアで同じプールを作り、
##       その分布に本物を当てる。**これが主張そのものに対する検定になる。**
##
##   プールは反復ごとに行う。反復 r について 3 コホートの切片を逆分散重みで
##   平均し、それを NMIS 個集めて帰無分布とする。各コホートの置換は独立なので、
##   どの組み合わせも正当な帰無標本である。
##
## 使い方（mismatch_slope_*.csv のあるディレクトリで）
##   Rscript mismatch_summary.R
##   COHORTS=GSE244679,GSE127165,GSE144269 Rscript mismatch_summary.R
##
##   POOL_COHORTS  既定 GSE244679,GSE127165,GSE144269
##                 （**GSE40419 はプールに入れない。** 除外したコホートなので）
##   DIR           既定 .
## ---------------------------------------------------------------------------

DIR  <- Sys.getenv("DIR", ".")
ALL  <- trimws(strsplit(Sys.getenv("COHORTS",
          "GSE244679,GSE127165,GSE144269,GSE40419"), ",")[[1]])
POOL <- trimws(strsplit(Sys.getenv("POOL_COHORTS",
          "GSE244679,GSE127165,GSE144269"), ",")[[1]])
ALL  <- ALL[nzchar(ALL)]; POOL <- POOL[nzchar(POOL)]

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 78), "\n")

d <- list()
for (CO in ALL) {
  f <- file.path(DIR, sprintf("mismatch_slope_%s.csv", CO))
  if (!file.exists(f)) { msg("（", basename(f), " がありません。飛ばします）"); next }
  x <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("kind", "intercept", "intercept_se") %in% names(x))) {
    msg("（", basename(f), " の列が足りません。飛ばします）"); next
  }
  d[[CO]] <- x
}
if (!length(d)) stop("読める CSV がありません。", call. = FALSE)

## ---------------------------------------------------------------- (1) 個別
hr(); msg("(1) コホートごと ― 切片")
msg("")
## 見出しは半角で書く（全角は端末の表示幅が合わず列がずれる）
msg(sprintf("%-11s %4s %11s %11s %9s %9s %9s",
            "cohort", "NMIS", "true", "mis_median", "mis_sd", "ratio", "p_emp"))
msg("  （mis_median = 崩したペアの中位、ratio = 崩し中位 / 本物、p_emp = 両側）")
rows <- list()
for (CO in names(d)) {
  x  <- d[[CO]]
  tr <- x$intercept[x$kind == "true"][1]
  mm <- x$intercept[x$kind == "mismatched"]
  if (!length(mm)) next
  p  <- (1 + sum(abs(mm) >= abs(tr))) / (length(mm) + 1)
  msg(sprintf("%-11s %4d %+11.5f %+11.5f %9.5f %9.3f %9.4g",
              CO, length(mm), tr, median(mm), sd(mm), median(mm) / tr, p))
  rows[[CO]] <- data.frame(cohort = CO, nmis = length(mm), true = tr,
                           mis_median = median(mm), mis_sd = sd(mm),
                           ratio = median(mm) / tr, p_emp = p,
                           stringsAsFactors = FALSE)
}
msg("")
msg("  崩し/本物 が 0 付近  → 患者内の対応が効いている。**主張が通る。**")
msg("  崩し/本物 が 1 付近  → 条件づけの産物。**主張は崩れる。**")
msg("  両側 p の下限は 1/(NMIS+1)。NMIS=20 なら 0.048。")

## ---------------------------------------------------------------- (2) プール
pool_avail <- intersect(POOL, names(d))
hr(); msg("(2) プールした切片 ― **これが主張そのものに対する検定**")
msg("")
if (length(pool_avail) < 2L) {
  msg("  プールできるコホートが ", length(pool_avail), " 個しかありません。")
} else {
  msg("  プール対象: ", paste(pool_avail, collapse = ", "),
      "（GSE40419 は除外したコホートなので入れない）")

  wmean <- function(est, se) {
    w <- 1 / se^2
    m <- sum(w * est) / sum(w)
    c(est = m, se = sqrt(1 / sum(w)),
      Q = sum(w * (est - m)^2), df = length(est) - 1L)
  }
  tr_e  <- vapply(pool_avail, function(CO)
    d[[CO]]$intercept[d[[CO]]$kind == "true"][1], 0)
  tr_s  <- vapply(pool_avail, function(CO)
    d[[CO]]$intercept_se[d[[CO]]$kind == "true"][1], 0)
  TP <- wmean(tr_e, tr_s)
  msg("")
  msg(sprintf("  本物      %+.5f ± %.5f   z = %.2f   Q = %.2f (df %d, p = %.3f)",
              TP["est"], TP["se"], TP["est"] / TP["se"], TP["Q"], TP["df"],
              stats::pchisq(TP["Q"], TP["df"], lower.tail = FALSE)))

  nr <- min(vapply(pool_avail, function(CO) sum(d[[CO]]$kind == "mismatched"), 0L))
  if (nr < 2L) {
    msg("  崩しの反復が足りません。")
  } else {
    pooled <- t(vapply(seq_len(nr), function(r) {
      e <- vapply(pool_avail, function(CO) {
        y <- d[[CO]][d[[CO]]$kind == "mismatched", ]; y$intercept[r] }, 0)
      s <- vapply(pool_avail, function(CO) {
        y <- d[[CO]][d[[CO]]$kind == "mismatched", ]; y$intercept_se[r] }, 0)
      wmean(e, s)
    }, numeric(4)))
    pe <- pooled[, "est"]; pz <- pooled[, "est"] / pooled[, "se"]
    msg(sprintf("  崩したペア 中位 %+.5f   sd %.5f   範囲 %+.5f 〜 %+.5f （%d 反復）",
                median(pe), sd(pe), min(pe), max(pe), nr))
    msg(sprintf("             z の中位 %.2f   範囲 %.2f 〜 %.2f",
                median(pz), min(pz), max(pz)))
    p_pool <- (1 + sum(abs(pe) >= abs(TP["est"]))) / (nr + 1)
    msg("")
    msg(sprintf("  **両側の経験的 p = %.4g**   （下限 %.4g）", p_pool, 1 / (nr + 1)))
    msg(sprintf("  比（崩し中位 / 本物） = %.3f", median(pe) / TP["est"]))
    write.csv(data.frame(rep = seq_len(nr), pooled),
              "mismatch_pooled.csv", row.names = FALSE)
    msg("  書き出し: mismatch_pooled.csv")
  }
}

## ------------------------------------------------- (3) GSE40419 との対比
hr(); msg("(3) 除外したコホートとの対比")
msg("")
if ("GSE40419" %in% names(d)) {
  x  <- d[["GSE40419"]]
  tr <- x$intercept[x$kind == "true"][1]
  mm <- x$intercept[x$kind == "mismatched"]
  msg(sprintf("  GSE40419   本物 %+.5f   崩し中位 %+.5f   比 %.3f",
              tr, median(mm), median(mm) / tr))
  msg("")
  msg("  **GSE40419 は条件 = バッチなので、崩しても切片が残るはずである。**")
  msg("  崩しても残る（比が 1 に近い）なら、それは疾患ではなくコホート全体の")
  msg("  条件差を測っている証拠であり、除外の正しさの確認になる。")
  msg("")
  msg("  残す 3 コホートで比が 0 付近、GSE40419 で比が 1 付近、という対比が")
  msg("  **「残した 3 コホートに交絡がない」ことの最も直接的な証拠**である。")
  msg("  この 2 行を並べて Methods の妥当性確認に書く。")
} else {
  msg("  GSE40419 の結果がありません。**対比が作れないので必ず回すこと。**")
}

if (length(rows)) {
  write.csv(do.call(rbind, rows), "mismatch_summary.csv", row.names = FALSE)
  hr(); msg("書き出し: mismatch_summary.csv")
}

hr(); msg("判定")
msg("")
msg("  (2) のプールが決め手である。論文の主張は 3 コホートをプールした")
msg("  +0.0544 ± 0.0059（z = 9.2）なので、崩したペアで同じプールを作って")
msg("  そこに本物を当てるのが、主張そのものに対する検定になる。")
msg("")
msg("  経験的 p が下限に張り付き、比が 0 付近 → **通った。**")
msg("    Methods の妥当性確認に書き、Limitations にも 1 文置く。")
msg("  比が 0.5 を超える → **主張の立て方を見直す。** Results を書き始めない。")
msg("    条件づけとマッチングが効果を作っている可能性を潰すまで先に進まない。")
