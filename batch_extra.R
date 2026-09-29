#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## batch_extra.R
##
## `diag_batch.R` が出していない 3 点を補う。**diag_batch.R は改変しない**
## （既に記録した数値の再現性を壊さないため）。
##
##   (A) 浅いライブラリの**全件**と条件別の内訳
##       diag_batch.R は先頭 8 本しか印字しない。GSE40419 では 12 本あり、
##       方針書 R3-1 の英文「seven libraries with fewer than 4e6 reads are
##       all on the normal side」が**誤り**（ERR164584 が disease で 3.86e6）
##       なので、正しい内訳を出す。
##
##   (B) 陰性対照を normal どうし **と** disease どうしの両方で、
##       ランダムな組み分けを NREP 回繰り返して安定させる。
##       diag_batch.R は normal どうし 1 回（34 組）だけ。
##
##   (C) **本物のペアの orph_excess を並べて出す。** diag_batch.R は陰性対照
##       しか出さないので、「実データと同程度の広がりか」を判定できなかった。
##       ここで同じ画面に並べる。
##
## ---------------------------------------------------------------------------
## 【量の定義 — 分母に注意】
##
##   f(lib)      = ゼロ率(orphan, lib) − ゼロ率(非 orphan, lib)
##   orph_excess = f(disease) − f(normal)
##
##   **分母は全転写産物（198,507 本）である。** diag_batch.R の (4) と同じで、
##   `diag_slope.R` の `d_zr_*`（`(normal>0)|(disease>0)` で絞った後）とは
##   分母が違う。**数値を混ぜてはならない。** 方針書 §0 (D) の注意と同じ。
##   ここで見たいのは「本物のペアと任意のペアで広がりが違うか」なので、
##   分母が揃っていれば足りる。
##
## 【使い方】リポジトリのルートで
##   KALLISTO_ROOT=Revised/GSE40419 METADATA_DIR=Revised/metadata \
##   COHORT=GSE40419 Rscript batch_extra.R
##
##   NREP    既定 200   ランダム組み分けの反復回数
##   SHALLOW 既定 4e6   「浅い」の絶対閾値（中位の 1/5 も併せて出す）
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

CO      <- Sys.getenv("COHORT", names(COHORTS)[1])
NREP    <- as.integer(Sys.getenv("NREP", "200"))
SHALLOW <- as.numeric(Sys.getenv("SHALLOW", "4e6"))

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

sheet <- read_sheet(CO)
dat   <- load_cohort_tpm(sheet, KALLISTO_ROOT, "abundance.tsv",
                         expected_features = N_TRANSCRIPTS)
is_o  <- orphan_index(dat$id, ORPHAN_ID_FILE, expected_n = N_ORPHAN_TRANSCRIPT)

## run_info.json の読み出しは diag_batch.R と同じ実装
gj <- function(dir, key) {
  p <- file.path(KALLISTO_ROOT, dir, "run_info.json")
  if (!file.exists(p)) return(NA_real_)
  s <- paste(readLines(p, warn = FALSE), collapse = " ")
  m <- regmatches(s, regexpr(paste0('"', key, '":\\s*[0-9.eE+-]+'), s))
  if (!length(m)) NA_real_ else as.numeric(sub('.*:\\s*', '', m))
}
sheet$depth <- vapply(sheet$directory, gj, 0, key = "n_processed")
sheet$pmap  <- vapply(sheet$directory, gj, 0, key = "p_pseudoaligned")
sheet$zr    <- vapply(sheet$library_id, function(l) mean(dat$tpm[, l] == 0), 0)
sheet$zr_o  <- vapply(sheet$library_id, function(l) mean(dat$tpm[is_o, l] == 0), 0)

hr(); msg(CO, " / ライブラリ ", nrow(sheet), " 本")
msg("  分母は全転写産物 ", nrow(dat$tpm), " 本。diag_slope.R の d_zr_* とは分母が違う。")

## ===========================================================================
## (A) 浅いライブラリの全件と条件別内訳
## ===========================================================================
hr(); msg("(A) 浅いライブラリの全件")
med <- median(sheet$depth, na.rm = TRUE)
thr_rel <- med / 5
msg(sprintf("  深さの中位 %.4g   中位の1/5 = %.4g   絶対閾値 = %.4g",
            med, thr_rel, SHALLOW))

report_thr <- function(thr, lab) {
  b <- sheet[is.finite(sheet$depth) & sheet$depth < thr, ]
  msg("")
  msg(sprintf("  【%s（< %.4g）】 該当 %d 本", lab, thr, nrow(b)))
  if (!nrow(b)) { msg("    なし"); return(invisible()) }
  tb <- table(b$condition)
  msg("    条件別: ", paste(sprintf("%s %d 本", names(tb), as.integer(tb)),
                             collapse = " / "))
  b <- b[order(b$depth), ]
  msg(sprintf("    %-14s %-9s %-10s %10s %8s %8s",
              "library_id", "condition", "patient", "depth", "map%", "zr"))
  for (i in seq_len(nrow(b)))
    msg(sprintf("    %-14s %-9s %-10s %10.3g %7.1f%% %8.3f",
                b$library_id[i], b$condition[i], b$patient_id[i],
                b$depth[i], b$pmap[i], b$zr[i]))
  ## 「最も浅い k 本がすべて同一条件」の k を出す。回答文に書けるのはこの形。
  runlen <- 1L
  while (runlen < nrow(b) && b$condition[runlen + 1L] == b$condition[1L])
    runlen <- runlen + 1L
  if (runlen == nrow(b)) {
    msg(sprintf("    → 該当する %d 本すべてが %s 側", runlen, b$condition[1L]))
  } else {
    msg(sprintf("    → 最も浅い %d 本がすべて %s 側（%d 本目 %s が %s）",
                runlen, b$condition[1L], runlen + 1L,
                b$library_id[runlen + 1L], b$condition[runlen + 1L]))
  }
}
report_thr(SHALLOW, "絶対閾値")
report_thr(thr_rel, "中位の1/5")

## ===========================================================================
## (B)(C) 陰性対照 ― 本物のペアと任意のペアを並べる
## ===========================================================================
## f(lib) はライブラリ 1 本で決まる量なので、任意の 2 本の差が取れる。
f <- vapply(sheet$library_id,
            function(l) mean(dat$tpm[is_o, l] == 0) - mean(dat$tpm[!is_o, l] == 0),
            0)
names(f) <- sheet$library_id

lb <- function(pid, cond)
  sheet$library_id[sheet$patient_id == pid & sheet$condition == cond][1]
pts <- unique(sheet$patient_id)
ln  <- vapply(pts, lb, "", cond = "normal")
ld  <- vapply(pts, lb, "", cond = "disease")
ok  <- !is.na(ln) & !is.na(ld)
pts <- pts[ok]; ln <- ln[ok]; ld <- ld[ok]

true_ex <- f[ld] - f[ln]          # 本物のペア（= 疾患効果 ＋ もしあればバッチ）

## 同一条件内のランダム組み分けを NREP 回。疾患効果ゼロの帰無。
## 反復ごとの中位も保持する。**本物の中位を帰無に当てるときは、帰無側も
## 中位でなければならない。** 個別値の分布に中位を当てると、帰無の裾が
## 広すぎて何でも「埋もれる」と出る（初版の誤り）。
set.seed(BASE_SEED)
pseudo <- function(libs) {
  if (length(libs) < 6L) return(list(vals = numeric(0), meds = numeric(0), k = 0L))
  k <- floor(length(libs) / 2)
  v <- vector("list", NREP); m <- numeric(NREP)
  for (r in seq_len(NREP)) {
    o <- sample(length(libs))
    d <- f[libs[o[k + seq_len(k)]]] - f[libs[o[seq_len(k)]]]
    v[[r]] <- unname(d); m[r] <- median(d)
  }
  list(vals = unlist(v, use.names = FALSE), meds = m, k = k)
}
PN <- pseudo(sheet$library_id[sheet$condition == "normal"])
PD <- pseudo(sheet$library_id[sheet$condition == "disease"])
pn <- PN$vals; pd <- PD$vals

summ <- function(v, lab) {
  if (!length(v)) { msg(sprintf("  %-26s （計算できません）", lab)); return(invisible()) }
  msg(sprintf("  %-26s n=%6d  中位 %+.4f  sd %.4f  IQR %.4f  |x|>0.05 %5.1f%%",
              lab, length(v), median(v), sd(v), IQR(v), 100 * mean(abs(v) > 0.05)))
}
hr(); msg("(B)(C) orph_excess ― 本物のペア vs 同一条件のランダムな組")
msg("  orph_excess = [ゼロ率(orphan) − ゼロ率(非orphan)] の disease − normal")
msg("  **正 = 疾患側で orphan が非 orphan より余分にゼロになる**")
msg("")
summ(true_ex, sprintf("本物のペア (%d 組)", length(true_ex)))
summ(pn, sprintf("normal どうし (%d 組 × %d)", PN$k, NREP))
summ(pd, sprintf("disease どうし (%d 組 × %d)", PD$k, NREP))

## --- 中位どうしの比較（次元を揃える）
meds <- c(PN$meds, PD$meds)
if (length(meds)) {
  obs   <- median(true_ex)
  p_emp <- (1 + sum(abs(meds - median(meds)) >= abs(obs - median(meds)))) /
           (length(meds) + 1)
  msg("")
  msg(sprintf("  本物の中位 %+.5f", obs))
  msg(sprintf("  帰無の中位の分布（反復 %d 回 × 2 条件）  中位 %+.5f  sd %.5f",
              NREP, median(meds), sd(meds)))
  msg(sprintf("  両側の経験的 p = %.4g   （下限 %.4g）",
              p_emp, 1 / (length(meds) + 1)))
  msg(sprintf("  ※ 帰無側は 1 反復あたり %d／%d 組、本物は %d 組。組数が少ない分",
              PN$k, PD$k, length(true_ex)))
  msg("     帰無の中位は本物より広く散るので、この p は保守的である。")
}
if (length(pn) && length(pd)) {
  null_all <- c(pn, pd)
  msg("")
  msg(sprintf("  sd 比（本物 / 帰無、いずれも個別値） = %.2f",
              sd(true_ex) / sd(null_all)))
}
msg("")
msg("  読み方")
msg("    本物の中位が帰無の中位の分布に埋もれる → その量は疾患に固有ではない。")
msg("    sd 比が 1 付近 → 患者ごとのばらつきは個体差で説明でき、")
msg("                     患者固有の疾患応答ではない（§0 (B)(C) と同じ結論）。")
msg("    **コホート = バッチのコホートでは、本物の中位が大きく出ても")
msg("    それは疾患効果ではなくバッチ効果である**（GSE40419）。")
msg("    したがって GSE40419 でこの値が大きいことは、除外の正しさの確認であって")
msg("    疾患効果の証拠ではない。3 コホートとの対比に使う。")

## 患者ごとの値を書き出す（後で 3 コホートを並べるため）
out <- data.frame(cohort = CO, patient_id = pts,
                  lib_normal = ln, lib_disease = ld,
                  f_normal = unname(f[ln]), f_disease = unname(f[ld]),
                  orph_excess = unname(true_ex), stringsAsFactors = FALSE)
write.csv(out, sprintf("batch_extra_%s.csv", CO), row.names = FALSE)
nullout <- data.frame(cohort = CO,
                      kind = c(rep("normal_normal", length(pn)),
                               rep("disease_disease", length(pd))),
                      value = c(pn, pd), stringsAsFactors = FALSE)
write.csv(nullout, sprintf("batch_extra_null_%s.csv", CO), row.names = FALSE)
hr(); msg("書き出し: batch_extra_", CO, ".csv / batch_extra_null_", CO, ".csv")
