#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## diag_mismatch.R
##
## 患者内の対応が本当に効いているかを、ペアを崩した陰性対照で確かめる。
##
## diag_batch.R の (4)（normal どうしの疑似ペア）は陰性対照として不適切だった。
## 別の患者同士を比べるので個体間変動を含み、対応のある設計が除いている分を
## 混ぜてしまう。比較対象として厳しすぎる。
##
## ここでは条件差は保ったまま患者内の対応だけを壊す。
##   実ペア    : normal(i) と disease(i)
##   崩したペア: normal(i) と disease(j), j != i
## コホート全体の条件差（およびバッチ差）は崩したペアにも同じだけ入る。
## したがって実ペアが崩したペアより極端であれば、その分が患者固有である。
##
## 使い方
##   KALLISTO_ROOT=Revised/GSE244679 METADATA_DIR=Revised/metadata \
##   COHORT=GSE244679 Rscript diag_mismatch.R [崩したペアの本数 既定 20]
##
## ---------------------------------------------------------------------------
## 【2026-09-23 訂正】帰無期待値の印字が誤っていた。**統計は変えていない。**
##
## 患者ごとの p は 2 つの片側経験的 p の**小さい方**であり、2 倍していない。
## したがって閾値 0.05 で切ったときの帰無での大きさは 0.05 ではない。
##
##   実ペアの値と崩したペア NMIS 個、計 NMIS+1 個は、患者固有成分が無ければ
##   交換可能である。実ペアの順位は一様なので
##       P(p < 0.05) = 2 * #{ r : (1+r)/(NMIS+1) < 0.05 } / (NMIS+1)
##   NMIS = 20 なら (1+r)/21 < 0.05 を満たすのは r = 0 だけで、
##       P = 2/21 = 0.0952        期待人数 = 0.0952 * np
##   すなわち「実ペアが 21 個の中で片側の端に来る」場合だけが該当する。
##
## 旧版は `0.05 * np` を期待値として印字していた（24 人なら 1.2 人）。
## **正しくは約 2.3 人。**2 倍の取り違えである。
## 以前の実行結果（2/24、4/57、9/70）はそのまま有効で、解釈だけが変わる。
## 旧: 2 が 1.2 を上回る → わずかに患者固有成分がありそう
## 新: 2 が 2.3 を下回る → 患者固有成分は見えない（本文の結論と整合）
##
## 同じ取り違えが test_2d.R にもあった。diag_mismatch_beta3.R は正しかった。
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")
CO   <- Sys.getenv("COHORT", names(COHORTS)[1])
a    <- commandArgs(TRUE)
NMIS <- if (length(a) >= 1) as.integer(a[1]) else 20L

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

## 閾値 alpha で切ったときの、この経験的 p の帰無での大きさ。
## p は 2 つの片側の min なので 0.05 にはならない。両側 2 倍ぶんを数える。
null_size <- function(nmis, alpha = 0.05) {
  r <- 0:nmis                       # 実ペアより極端な崩しの本数
  min(2 * sum((1 + r) / (nmis + 1) < alpha) / (nmis + 1), 1)
}

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

## orph_excess = (orphan のゼロ率変化) - (非 orphan のゼロ率変化)
## マッチングも層別化も使わない。
oe <- function(libN, libD) {
  vN <- dat$tpm[, libN]; vD <- dat$tpm[, libD]
  (mean(vD[is_o] == 0) - mean(vN[is_o] == 0)) -
    (mean(vD[!is_o] == 0) - mean(vN[!is_o] == 0))
}

real <- vapply(seq_len(np), function(i) oe(ln[i], ld[i]), 0)

set.seed(BASE_SEED)
mis <- unlist(lapply(seq_len(np), function(i) {
  j <- sample(setdiff(seq_len(np), i), min(NMIS, np - 1L))
  vapply(j, function(jj) oe(ln[i], ld[jj]), 0)
}))

q <- function(v) sprintf("中位 %+.4f  四分位 %+.4f / %+.4f  範囲 %+.4f 〜 %+.4f",
                         median(v), quantile(v, .25), quantile(v, .75),
                         min(v), max(v))

hr(); msg(CO, " / 患者 ", np, " 人 / 崩したペア ", length(mis), " 通り")
msg("")
msg("orph_excess の分布")
msg("  実ペア      ", q(real))
msg("  崩したペア  ", q(mis))
msg("")
msg(sprintf("  ばらつき  実ペア sd %.4f  IQR %.4f", sd(real), IQR(real)))
msg(sprintf("            崩し   sd %.4f  IQR %.4f", sd(mis),  IQR(mis)))
msg(sprintf("  sd 比 (実/崩し) = %.2f", sd(real) / sd(mis)))

## 患者ごとの経験的 p。実ペアが崩したペアの分布のどこに来るか。
## **2 つの片側の min であり、2 倍していない。**下の印字を参照。
pv <- vapply(seq_len(np), function(i) {
  j <- setdiff(seq_len(np), i)
  jj <- j[seq_len(min(NMIS, length(j)))]
  m  <- vapply(jj, function(x) oe(ln[i], ld[x]), 0)
  min((1 + sum(m >= real[i])) / (length(m) + 1),
      (1 + sum(m <= real[i])) / (length(m) + 1))
}, 0)

a05 <- null_size(NMIS, 0.05)
a10 <- null_size(NMIS, 0.10)
hr(); msg("患者ごとの経験的 p（自分の normal を他患者の disease と組にした分布に対して）")
msg("  **2 つの片側 p の小さい方であり、2 倍していない。**")
msg(sprintf("  取りうる最小値 %.4f = 1/(NMIS+1)", 1 / (NMIS + 1)))
msg(sprintf("  閾値 0.05 で切ったときの帰無での大きさ = %.4f（0.05 ではない）", a05))
msg(sprintf("  p < 0.05 の患者 %d / %d  （**帰無なら約 %.1f 人**）",
            sum(pv < 0.05), np, a05 * np))
msg(sprintf("  p < 0.10 の患者 %d / %d  （帰無なら約 %.1f 人）",
            sum(pv < 0.10), np, a10 * np))
msg(sprintf("  p の中位 %.3f", median(pv)))
msg("  観測が期待を**下回る**なら、患者固有成分は見えていない。")

## 条件差の全体成分（崩したペアの中位）を引いた残差
hr(); msg("コホート全体の条件差を引いた残差")
base <- median(mis)
msg(sprintf("  崩したペアの中位 = %+.4f  ← コホート全体の条件差（バッチ含む）", base))
msg(sprintf("  実ペア - その中位: %s", q(real - base)))
msg(sprintf("  符号が + の患者 %d / %d", sum(real - base > 0), np))

hr(); msg("読み方")
msg("  sd 比が 1 付近で p<0.05 の患者が上の期待人数と同程度なら、患者内の")
msg("  対応は効いておらず、観測されているのはコホート全体の条件差")
msg("  （またはバッチ差）である。")
msg("  sd 比が 1 より大きく p<0.05 が期待人数を超えるなら、その超過分が患者固有。")
msg("")
msg("  **期待人数は 0.05 * np ではない。**上の a05 を使うこと。旧版はここを")
msg("  0.05 * np と印字しており、期待値を半分に見積もっていた（2026-09-23 訂正）。")
msg("")
msg("  注意: 崩したペアにもコホート全体の条件差は入る。したがってこの検定は")
msg("  「疾患効果があるか」ではなく「患者ごとに違うか」を見ている。論文の")
msg("  主張が患者レベルの異質性なので、問うべきはこちらである。")

write.csv(data.frame(patient_id = pts, orph_excess = real, p_mismatch = pv),
          paste0("mismatch_", CO, ".csv"), row.names = FALSE)
msg("")
msg("書き出し: ", paste0("mismatch_", CO, ".csv"))
