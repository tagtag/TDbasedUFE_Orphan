#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## rho_compare.R — R2-1／R3-2 の決め手。camera との食い違いを裁定する。
##
## 2 つのことを同時に測る。どちらも患者内の対数比行列
##   r[t, i] = log2( (TPM_disease[t,i] + PRIOR) / (TPM_normal[t,i] + PRIOR) )
## から出る。
##
## ── (1) 集合内の平均対相関 rho_bar ─────────────────────────────────
##
## camera は転写産物間相関 rho を推定し、分散を VIF = 1 + (m-1)*rho 倍する。
## 実測で rho = 0.185 / 0.020 / 0.072、m = 1538 / 611 / 649 なので
## VIF = 285 / 13 / 48、実効的な独立単位は 5 / 46 / 14 本にまで落ちる。
## これが camera p が 0.11〜0.21 にとどまる理由である。
##
## 一方、本研究の置換検定は「同一ライブラリから引いたマッチ対照集合」を帰無に
## するので、対照集合が持つ相関はすでに帰無分布に入っている。
##
## **ここで測る rho_bar は仮定の検査ではない。置換検定の妥当性の条件でもない。**
##
## 置換検定の帰無仮説は「orphan 集合はマッチ非 orphan 集合と交換可能である」で
## ある。orphan の転写産物間相関が対照より強ければ、**その時点で帰無仮説は既に
## 破れている**。したがって棄却は正しい棄却であり、第一種の過誤ではない。
## 相関が一致していることを確かめる必要はない。
##
## そのうえ rho_bar が大きいことは「雑音が信号に見える」状況ではない。2190 本が
## 揃って動くとは orphan 固有の共通因子が実在するということで、それ自体が所見で
## ある。
##
## 残るのは 2 点だけで、どちらも分散の話ではない。
##
##  (1) 解釈の一意性。棄却は「交換可能でない」を意味し、平均のずれでも相関構造の
##      違いでも起こる。だから平均差と rho_bar の両方を記述として出す。これは
##      p 値の妥当性の問題ではなく結論文の書き方の問題である。
##  (2) 効果量の外挿。共通因子が強ければコホートを変えたときの振れ幅は対照集合の
##      再抽出から出る SE より大きい。効果の存在と向きには影響しないが、大きさの
##      信頼区間には影響する。これに答えるのは分散補正ではなく**反復**であり、
##      3 コホートで向きが一致し各コホート内で 18/24、43/57、57/70〔**ゼロ率**の
##      超過。このスクリプトが出す log2 比の Delta<0 の数（18/24、35/57、54/70）
##      とは別の量なので取り違えないこと〕という形で既に答えが出ている。
##
## したがって rho_bar は **R2-1 が求めた記述**として測る。転写産物単位の置換 p は
## 妥当なので引用してよい。camera が保守的になる理由の説明にもなる。
##
## **本当の残余リスクは分散ではなく交絡である。** orphan の共通因子が技術的な
## もので（短い、反復配列に多い等）、その技術変数が腫瘍と正常で系統的に違えば、
## 3 コホートで向きが揃っても転写の下方制御とは言えない。手当てはマッチング
## （発現水準・検出頻度）、深さとマップ率での調整、ペアを崩した陰性対照であり、
## 相関の大きさではなくこの経路を潰すことが要点である。
##
##   rho_bar = ( sum_{t,u} cor(r[t,], r[u,]) - m ) / ( m (m-1) )
##
## m x m 行列は作らない。各行を標準化した Z について
##   sum_{t,u} cor = ||colSums(Z)||^2 / (np - 1)
## なので O(m * np) で終わる。cor() と完全一致することを検算済み。
##
## ── (2) 患者を単位とした log2 比の検定 ─────────────────────────────
##
## camera への反論は「主張の単位は転写産物ではなく患者である」だが、camera が
## 検定したのは log2FC で、rev.18 の 118/151 はゼロ率である。**別の量で反論して
## いることになる**（しかも 118/151 はコホートをプールした値で、異質性のため
## rev.27・rev.35 で撤回済み）。そこで **camera と同一の量（log2 比）で、単位
## だけ患者に変えた検定**を行う。患者ごとに
##
##   Delta_i = mean(r[orphan, i]) - mean_k( mean(r[control_k, i]) )
##
## を計算し、**1 標本 t（両側）を主**、Wilcoxon と Delta<0 の数を頑健性として
## 併記する。転写産物間相関は Delta_i の精度を落とすが、患者 i と患者 j の
## 独立性は壊さないので、この検定は camera が適用する分散膨張を受けない。
##
## **検定は両側である（rev.36 で片側から変更）。** 以前の版は
## alternative = "less" / "greater" で片側を印字しており、論文の表（両側）と
## 食い違っていた。事前に向きを決めていないので両側が正しい。
##
## de_logfc_<COHORT>_voom.csv が手元にあれば、**DE 解析が実際に使った転写産物に
## 限定した版**も併せて計算する。これで camera と土俵が完全に揃う。
##
## 使い方
##   NSET=20 N_CORES=12 KALLISTO_ROOT=Revised/GSE244679 \
##   METADATA_DIR=Revised/metadata COHORT=GSE244679 Rscript rho_compare.R
##
## 環境変数
##   NSET     既定 20   マッチ対照集合の組数
##   PRIOR    既定 1    対数比の下駄（TPM 単位）。0 は不可（log が発散する）
##   DE_FILE  既定 de_logfc_<COHORT>_voom.csv   無ければ全集合のみ
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

CO        <- Sys.getenv("COHORT", names(COHORTS)[1])
NSET      <- as.integer(Sys.getenv("NSET", "20"))
N_CORES   <- as.integer(Sys.getenv("N_CORES", "1"))
N_STRATA2 <- as.integer(Sys.getenv("N_STRATA2", "4"))
PRIOR     <- as.numeric(Sys.getenv("PRIOR", "1"))
DE_FILE   <- Sys.getenv("DE_FILE", sprintf("de_logfc_%s_voom.csv", CO))
stopifnot(is.finite(PRIOR), PRIOR > 0, NSET >= 2L)
## rho_bar の経験的 p の下限は 1/(NSET+1)。NSET が小さいと判定の閾値 0.2 に
## 届かず、常に「orphan の方が強く相関」と出てしまう。計算は軽いので大きく取る。
if (NSET < 20L)
  warning(sprintf(paste0("NSET = %d では経験的 p の下限が %.3f で、判定が常に",
                         "「orphan の方が強く相関」に倒れます。NSET >= 50 を推奨。"),
                  NSET, 1/(NSET+1)), call. = FALSE)

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 76), "\n")

## 本体から層別化と抽出を借りる（実装を二重に持たない）
SRC <- if (file.exists("R/07_interaction_permutation.R"))
         "R/07_interaction_permutation.R" else "07_interaction_permutation.R"
src <- readLines(SRC)
eval(parse(text = paste(src[grep("^make_strata <- function", src):
                            (grep("^run_patient <- function", src) - 1L)],
                        collapse = "\n")))
stopifnot(exists("make_strata_multi"), exists("draw_matched_sets"))

sheet <- read_sheet(CO)
tab <- table(sheet$patient_id, sheet$condition)
okp <- rownames(tab)[tab[, "normal"] >= 1 & tab[, "disease"] >= 1]
sheet <- sheet[sheet$patient_id %in% okp, ]
pts <- unique(sheet$patient_id)
lb  <- function(pid, cond)
  sheet$library_id[sheet$patient_id == pid & sheet$condition == cond][1]
ln <- vapply(pts, lb, "", cond = "normal")
ld <- vapply(pts, lb, "", cond = "disease")
np <- length(pts)
stopifnot(np >= 6L)

dat  <- load_cohort_tpm(sheet, KALLISTO_ROOT, "abundance.tsv",
                        expected_features = N_TRANSCRIPTS)
is_o <- orphan_index(dat$id, ORPHAN_ID_FILE, expected_n = N_ORPHAN_TRANSCRIPT)

## 患者内の対数比。TPM は既に総和 1e6 に揃っているのでそのまま使える。
R <- log2((dat$tpm[, ld, drop = FALSE] + PRIOR) /
          (dat$tpm[, ln, drop = FALSE] + PRIOR))
colnames(R) <- pts

## 両条件でゼロの転写産物は対数比が恒等的に 0 で分散を持たず相関が定義できない。
## 分散ゼロの行も落とす。本体と同じ絞り込みを先にかける。
nm   <- dat$tpm[, ln, drop = FALSE]
keep <- rowSums(dat$tpm[, c(ln, ld), drop = FALSE] > 0) > 0
sdv  <- apply(R, 1, sd)
keep <- keep & is.finite(sdv) & sdv > 0

tid    <- dat$id[keep]
R      <- R[keep, , drop = FALSE]
orph   <- is_o[keep]

## 層別変数。**leave-one-out ではない（rev.36 で明記）。**
##
## 本体と diag_slope.R は患者ごとに対照集合を引き直すので、当該患者の normal を
## 抜いて層を作る（抜かないと対照がその患者に寄り、検定が甘くなる）。
## ここは違う。**コホート水準の 1 組の対照集合を全患者に共通に使う**ので、
## 特定の患者に寄りようがなく leave-one-out の必要がない。全 normal の平均と
## 検出頻度で層を作る。de_framework.R も同じ理由で leave-one-out を使わない
## （ただし向こうは CPM、こちらは TPM）。
##
## **Methods に書き分けること。** 書かないと「マッチングが解析ごとに違う」と
## 読まれる。方針書 §0 (E) の表と文案を参照。
mean_n <- rowMeans(nm)[keep]
freq_n <- rowMeans(nm > 0)[keep]

## ---------------------------------------------------------------------------
## 1 つの転写産物部分集合について、rho_bar と患者ごとの Delta を計算する。
## sel は R の行に対する論理ベクトル。
## ---------------------------------------------------------------------------
rho_bar <- function(M) {
  m <- nrow(M)
  if (m < 2L) return(NA_real_)
  Z <- (M - rowMeans(M)) / apply(M, 1, sd)
  s <- sum(colSums(Z)^2) / (ncol(M) - 1)   # = sum_{t,u} cor(t,u)
  (s - m) / (m * (m - 1))
}

analyze <- function(sel, label) {
  Rs <- R[sel, , drop = FALSE]; os <- orph[sel]
  if (sum(os) < 20L || sum(!os) < 20L) {
    msg("  ", label, ": orphan ", sum(os), " 本 / 非 orphan ", sum(!os),
        " 本しかないので飛ばします。")
    return(NULL)
  }
  st <- make_strata_multi(list(mean_n[sel], freq_n[sel]), c(N_STRATA, N_STRATA2))

  set.seed(BASE_SEED)
  ctrl <- vector("list", NSET); short <- 0L
  for (k in seq_len(NSET)) {
    s1 <- draw_matched_sets(st, os, n_sets = 1L)
    short <- short + attr(s1, "shortfall")
    ctrl[[k]] <- s1[[1]]
  }

  f <- function(k) {
    Mc <- Rs[ctrl[[k]], , drop = FALSE]
    list(rho = rho_bar(Mc), mean_by_patient = colMeans(Mc))
  }
  lst <- if (N_CORES > 1L) {
    parallel::mclapply(seq_len(NSET), f, mc.cores = min(N_CORES, NSET))
  } else {
    lapply(seq_len(NSET), f)
  }
  bad <- which(!vapply(lst, is.list, logical(1)))
  if (length(bad))
    stop(label, ": ", length(bad), " 組でエラー。N_CORES=1 で確認してください。",
         call. = FALSE)

  ## 患者ごとの Delta。orphan の患者内平均 log2 比 - 対照の同じ量（NSET 組の平均）。
  o_by_pt <- colMeans(Rs[os, , drop = FALSE])
  c_by_pt <- rowMeans(do.call(cbind, lapply(lst, `[[`, "mean_by_patient")))
  list(label = label, m = sum(os), n_kept = nrow(Rs), shortfall = short,
       rho_o = rho_bar(Rs[os, , drop = FALSE]),
       rho_c = vapply(lst, `[[`, 0, "rho"),
       delta = o_by_pt - c_by_pt,
       o_by_pt = o_by_pt, c_by_pt = c_by_pt)
}

report <- function(a) {
  if (is.null(a)) return(invisible())
  hr(); msg("【", a$label, "】 解析対象 ", a$n_kept, " 転写産物 / orphan ", a$m, " 本",
            " / shortfall ", a$shortfall)
  if (a$shortfall > 0)
    msg("  ** shortfall があるので対照が揃っていません。N_STRATA を下げること。**")

  ## ---- (1) rho_bar
  rc <- a$rho_c; ro <- a$rho_o; m <- a$m
  vf <- function(r) 1 + (m - 1) * max(r, 0)
  pe <- (1 + sum(rc >= ro)) / (length(rc) + 1)
  msg("")
  msg("(1) 集合内の平均対相関 rho_bar（患者をまたいだ対数比の共変動）")
  msg(sprintf("    orphan        %+.4f", ro))
  msg(sprintf("    マッチ対照    中位 %+.4f   範囲 %+.4f 〜 %+.4f   （%d 組）",
              median(rc), min(rc), max(rc), length(rc)))
  msg(sprintf("    経験的片側 p（対照より orphan が大きい） = %.3f  ← 最小 %.3f",
              pe, 1/(length(rc)+1)))
  msg(sprintf("    VIF  orphan %8.1f（実効 %6.1f 本）   対照中位 %8.1f（実効 %6.1f 本）",
              vf(ro), m/vf(ro), vf(median(rc)), m/vf(median(rc))))
  if (pe > 0.2) {
    msg("    → 対照集合の分布内。orphan に固有の共通因子は検出されない。")
  } else {
    msg("    → orphan の方が強く相関している。**これは所見であって不都合ではない。**")
    msg("      orphan 固有の共通因子が実在するということであり、その時点で帰無仮説")
    msg("      （マッチ対照と交換可能）は破れている。置換検定の棄却は正しい棄却。")
  }
  msg("    この量は R2-1 が求めた記述である。置換検定の妥当性の条件ではないので、")
  msg("    どちらの結果でも転写産物単位の p は引用できる。結論文では平均差と")
  msg("    この値の両方を出し、棄却が平均のずれによることを明示する。")
  msg("    なお論文に載せるのは camera が推定した相関（de_framework.R の")
  msg("    Correlation 列。0.185 / 0.020 / 0.072）であって、この自前の値ではない。")

  ## ---- (2) 患者を単位とした log2 比の検定
  d <- a$delta
  npp <- length(d); k <- sum(d < 0)
  ## **すべて両側（rev.36 で片側から変更）。** 事前に向きを決めていないので
  ## 両側が正しく、論文の表の値も両側である。以前の版は alternative = "less"
  ## / "greater" を渡しており、再実行すると論文と違う値が出ていた。
  bt <- binom.test(k, npp, 0.5)
  wt <- suppressWarnings(wilcox.test(d))
  tt <- t.test(d)
  se <- sd(d) / sqrt(npp)
  ## 発現量に直した 95% CI。Delta が負なら conf.int[2] 側が「減少が小さい」端。
  pc <- function(x) 100 * abs(1 - 2^x)
  msg("")
  msg("(2) 患者を単位とした検定（camera と同一の量 = log2 比。単位だけ患者に変更）")
  msg(sprintf("    orphan の患者内平均 log2 比   中位 %+.4f", median(a$o_by_pt)))
  msg(sprintf("    対照の同じ量                  中位 %+.4f", median(a$c_by_pt)))
  msg(sprintf("    Delta = orphan - 対照         平均 %+.5f ± %.5f (SE)   中位 %+.5f",
              mean(d), se, median(d)))
  msg(sprintf("    95%% CI（両側 t）              %+.5f 〜 %+.5f",
              tt$conf.int[1], tt$conf.int[2]))
  msg(sprintf("    発現量に直すと                %.1f%% %s（95%% CI %.1f 〜 %.1f%%）",
              pc(mean(d)), if (mean(d) < 0) "減" else "増",
              pc(tt$conf.int[2]), pc(tt$conf.int[1])))
  msg(sprintf("    **1 標本 t p（両側） = %.3g**   ← 主たる検定", tt$p.value))
  msg(sprintf("    Wilcoxon 符号付き順位 p（両側） = %.3g", wt$p.value))
  msg(sprintf("    Delta < 0 の患者  %d / %d   符号検定 p（両側） = %.3g",
              k, npp, bt$p.value))
  msg("    Delta_i は 611〜1538 本の log2 比の平均なのでほぼ正規である。")
  msg("    したがって 1 標本 t を主とし、Wilcoxon と Delta<0 の数は頑健性として")
  msg("    全コホートに一律に併記する（rev.27 で検定の事後選択を撤回した）。")
  msg("    転写産物間相関は Delta_i の精度を落とすが、患者 i と患者 j の独立性は")
  msg("    壊さない。したがってこの検定は camera の分散膨張を受けない。")
}

hr(); msg(CO, " / 患者 ", np, " 人 / 対数比の下駄 ", PRIOR, " TPM / 対照 ", NSET, " 組")

A <- analyze(rep(TRUE, nrow(R)), "全集合（コホート内のどこかで非ゼロな転写産物すべて）")
report(A)

B <- NULL
if (file.exists(DE_FILE)) {
  de <- read.csv(DE_FILE, stringsAsFactors = FALSE)
  stopifnot("transcript" %in% names(de))
  sel <- tid %in% de$transcript
  msg("")
  msg("camera 同一集合（DE = differential expression、差次的発現解析の",
      " filterByExpr 通過分）でも計算します: ", DE_FILE, "、", sum(sel), " 本一致")
  B <- analyze(sel, "camera 同一集合（DE 解析の filterByExpr 通過分。camera が検定した集合そのもの）")
  report(B)
} else {
  msg("")
  msg("※ ", DE_FILE, " が無いので camera 同一集合の版は計算していません。")
  msg("  de_framework.R を走らせた後に再実行すると、camera と同じ転写産物集合で")
  msg("  患者単位の検定が出ます（camera への反論を同一の土俵で行うため）。")
}

out <- do.call(rbind, lapply(Filter(Negate(is.null), list(A, B)), function(a)
  data.frame(cohort = CO, subset = a$label, patient_id = names(a$delta),
             orphan_log2 = a$o_by_pt, control_log2 = a$c_by_pt,
             delta = a$delta, row.names = NULL, stringsAsFactors = FALSE)))
write.csv(out, paste0("rho_delta_", CO, ".csv"), row.names = FALSE)
out2 <- do.call(rbind, lapply(Filter(Negate(is.null), list(A, B)), function(a)
  data.frame(cohort = CO, subset = a$label, m = a$m, shortfall = a$shortfall,
             rho_orphan = a$rho_o, rho_control = a$rho_c,
             set = seq_along(a$rho_c), row.names = NULL, stringsAsFactors = FALSE)))
write.csv(out2, paste0("rho_", CO, ".csv"), row.names = FALSE)
hr(); msg("書き出し: ", paste0("rho_", CO, ".csv"), " と ",
          paste0("rho_delta_", CO, ".csv"))
msg("")
msg("**3 コホートをプールしないこと（rev.27・rev.35 で撤回）。**")
msg("効果量の異質性の Q 検定が p = 0.001 で、コホート間で 6 倍違う")
msg("（-0.117 / -0.019 / -0.046）。固定効果の重み付き平均も、Delta<0 の患者数を")
msg("合計した符号検定も使わない。**3 コホートを個別に報告する。**")
msg("各コホートで既に有意なので、プールは何も足さない。")
msg("以前このファイルは「3コホート揃ったら合計して符号検定をかけること」と")
msg("書いていたが、それが撤回した方針である。")
msg("")
msg("camera 同一集合の行が camera と同じ転写産物・同じ量での回答になる。")
msg("それが camera と同一の量・患者単位での主張の根拠になる。論文の主表には")
msg("この行を使い、全集合の行は Supplementary に回す。")
msg("")
msg("(1) の結果で主張が変わることはない。相関が対照より強ければ、その時点で")
msg("帰無仮説は破れており棄却は正しい。(1) は記述として報告する。")
msg("")
msg("注意すべきは分散ではなく交絡である。orphan の共通因子が技術的で、その")
msg("技術変数が条件間で違うなら、向きが揃っても転写の下方制御とは言えない。")
msg("マッチング・深さ調整・ペアを崩した陰性対照がその経路への手当てである。")
msg("主張は測定水準（マッチ対照に対する検出の低下）にとどめること。")
