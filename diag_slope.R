#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## diag_slope.R
##
## 生き残っている主張（コホート水準）を 1 コホート分、数値で出す。
##
##   患者ごとに、疾患側 - 正常側 のゼロ率変化を 3 つ計算する。
##     d_zr_all     全転写産物             … 検出力の全体的な移動
##     d_zr_orphan  orphan 転写産物
##     d_zr_control 発現量マッチ対照転写産物（external2 の層で orphan と同数）
##
##   そして 3 本の回帰を当てる。コホート内の患者が 1 点である（n = 24 / 57 / 70）。
##     fit_o : d_zr_orphan  ~ d_zr_all      → 「orphan の傾き」
##     fit_c : d_zr_control ~ d_zr_all      → 「対照の傾き」
##     fit_e : excess       ~ d_zr_all      （excess = d_zr_orphan - d_zr_control）
##
##   論文の表の各列は次の通り。
##     orphan の傾き … fit_o の傾き
##     対照の傾き   … fit_c の傾き
##     増幅率       … 上の 2 つの比（d_zr_all の分母の取り方に不変）
##     傾きの差     … **fit_e の傾き**
##     切片         … **fit_e の切片**   ← 傾きの差と同じ 1 本の回帰から出る
##
##   報告する ± は回帰の標準誤差である。**対照集合を NSET 組引くばらつきは
##   これに既に含まれている。** excess = d_zr_orphan - d_zr_control なので、
##   抽出ノイズは y の測定誤差として回帰の残差に入り、切片の SE に伝播する
##   （患者ごとに独立に引いているため iid 誤差になる）。二乗和で別途足し込むのは
##   二重計上であり、やってはならない。下の確認で見るのは「残差のうち抽出
##   ノイズが占める割合」で、これが小さければ NSET は十分という意味である。
##
## 置換検定は不要なので速い。マッチングは本体と同じ external2。
##
## 使い方（リポジトリのルートで）
##   NSET=200 N_CORES=12 KALLISTO_ROOT=Revised/GSE244679 METADATA_DIR=Revised/metadata \
##   COHORT=GSE244679 Rscript diag_slope.R
##   ※ 論文の数値は NSET=200 で出している。既定の 20 では小数第 3 位が動く。
##
## 読み方（rev.36 で更新。**実測を踏まえた現行の方針はこちら**）
##   対照の傾きが 1 付近なら、対照は全体の検出移動をそのまま写しているだけ。
##   orphan の傾きがそれを超えれば、orphan は全体の移動を増幅している。
##
##   **切片が 0 と区別できるかどうかが要点である。**
##   実測では 3 コホートで切片 +0.0250 / +0.0553 / +0.0659、逆分散重み付き平均
##   +0.0544 ± 0.0059（z = 9.2）、異質性の Q 検定 p = 0.163（異質性なし）。
##   **つまり切片は 0 ではない。** 全体の検出移動がゼロでも orphan は対照より
##   ゼロ率が 5.4 ポイント多く増える。これが論文の中核の主張である。
##
##   増幅率（傾きの比）はコホート間で 3.53 / 1.52 / 1.27 と揃わず（Q 検定
##   p = 5.0e-09）、しかも検出されている orphan 数と完全に順位相関する
##   （rho = 1）ので、ライブラリの検出感度との交絡が切れない。
##   **したがって切片を主、傾きの差を副として報告する。増幅率をプールしたり、
##   「3.5 倍」を見出しにしたりしてはならない。**
##
##   （この方針は rev.15 で確定した。それ以前のこのファイルのコメントは
##   「切片は 0 になるはずで、主張は増幅である」と書いていたが、3 コホートの
##   実測でそうならなかったため撤回した。）
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

CO        <- Sys.getenv("COHORT", names(COHORTS)[1])
N_STRATA2 <- as.integer(Sys.getenv("N_STRATA2", "4"))
NSET      <- as.integer(Sys.getenv("NSET", "20"))   # 対照集合の反復数（平均を取る）
N_CORES   <- as.integer(Sys.getenv("N_CORES", "1"))
## EXCLUDE_LIBS: 除外するライブラリ ID（カンマ区切り）。そのライブラリを持つ
## 患者をペアごと落とす。浅いライブラリの感度分析（R2-4）用。
## 例: EXCLUDE_LIBS=SRR8631680
EXCL <- trimws(strsplit(Sys.getenv("EXCLUDE_LIBS", ""), ",")[[1]])
EXCL <- EXCL[nzchar(EXCL)]
## 乱数種は one() の中で BASE_SEED + i と患者ごとに固定してあるので、
## 並列数を変えても結果は bit 単位で一致する。
## 本体は N_PERM = 1000 回引いた平均を zero_control_* として出す。ここは NSET 回。
## 対照集合の平均は回数に対してすぐ収束するので 20 でほぼ一致するが、小数第4位
## あたりは動く。手計算との比較では傾きの小数第2位までを見ること。
## 大きく違う場合は回数ではなく集合の構成が違っているので、n_kept と n_orphan を
## 本体の出力と突き合わせる。

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

## 本体から層別化と抽出だけ借りる（実装を二重に持たない）
SRC <- if (file.exists("R/07_interaction_permutation.R"))
         "R/07_interaction_permutation.R" else "07_interaction_permutation.R"
src <- readLines(SRC)
eval(parse(text = paste(src[grep("^make_strata <- function", src):
                            (grep("^run_patient <- function", src) - 1L)],
                        collapse = "\n")))
stopifnot(exists("make_strata_multi"), exists("draw_matched_sets"))

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

if (length(EXCL)) {
  unknown <- setdiff(EXCL, sheet$library_id)
  if (length(unknown))
    stop("EXCLUDE_LIBS にこのコホートに無いライブラリがあります: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  drop <- ln %in% EXCL | ld %in% EXCL
  message("  EXCLUDE_LIBS により患者 ", sum(drop), " 人を除外: ",
          paste(pts[drop], collapse = ", "))
  pts <- pts[!drop]; ln <- ln[!drop]; ld <- ld[!drop]
}

np  <- length(pts)
stopifnot(np >= 6L)

normal_mat   <- dat$tpm[, ln, drop = FALSE]
sum_normal   <- rowSums(normal_mat)
count_normal <- rowSums(normal_mat > 0)

one <- function(i) {
  vN <- dat$tpm[, ln[i]]; vD <- dat$tpm[, ld[i]]

  ## 解析対象の集合を本体 run_patient と厳密に揃える。本体は
  ##   keep <- (normal > 0) | (disease > 0)
  ## で絞った後にゼロ率も層別化も行う。ここで絞らないと、両条件でゼロの
  ## 転写産物が分母に入ってゼロ率の差が薄まり、本体の出力（zero_orphan_normal
  ## 等）と突き合わせられなくなる。恒等式 excess = -zero_gap_asym の検算も
  ## 揃えていることが前提である。
  ##
  ## なお全転写産物を分母にすると d_zr_all だけが小さくなるので傾きが一律に
  ## 約 1.2 倍になる（増幅率は不変）。論文の数値はすべてこの絞り込み後の
  ## 分母である。分母を混ぜて報告しないこと。
  keep <- (vN > 0) | (vD > 0)
  vN <- vN[keep]; vD <- vD[keep]
  orph <- is_o[keep]

  ## external2 の層別変数。leave-one-out なので自分の normal を抜く。
  ## 絞った後の添字で作る（本体も keep 後に make_strata_multi を呼ぶ）。
  ##
  ## ここで leave-one-out が要るのは、**患者ごとに対照集合を引き直す**からで
  ## ある。当該患者の normal で層を作ると対照がその患者に寄り、検定が甘くなる。
  ## コホート水準で 1 組の対照集合を全患者に共通に使う解析（de_framework.R、
  ## rho_compare.R）では leave-one-out を行っていない。Methods で書き分けること。
  ext  <- ((sum_normal - dat$tpm[, ln[i]]) / max(np - 1L, 1L))[keep]
  freq <- ((count_normal - (dat$tpm[, ln[i]] > 0)) / max(np - 1L, 1L))[keep]
  st   <- make_strata_multi(list(ext, freq), c(N_STRATA, N_STRATA2))

  ## 対照集合を NSET 回引く。**独立に 1 組ずつ**引くこと。
  ##
  ## draw_matched_sets(..., n_sets = NSET) は互いに素な NSET 組を一度に作るので、
  ## 層ごとに orphan 数の NSET 倍の候補を同時に要求する。本体は 1 組 + 互いに素な
  ## 2 組（3 倍）しか要求しないので足りるが、NSET = 20 では層が枯れて shortfall が
  ## 出る（GSE127165 で 209 件）。枯れた層は対照に入らないので平均が偏る。
  ##
  ## ここで欲しいのは対照のゼロ率変化の期待値だけで、組の間の独立性は不要
  ## （互いに素が必要なのは置換検定の帰無分布を作るときだけ）。1 組ずつ引けば
  ## 不偏で、要求は本体と同じ 1 倍になる。
  set.seed(BASE_SEED + i)
  short <- 0L
  cs <- numeric(NSET)

  zr <- function(idx) c(mean(vN[idx] == 0), mean(vD[idx] == 0))
  for (k in seq_len(NSET)) {
    one_set <- draw_matched_sets(st, orph, n_sets = 1L)
    short <- short + attr(one_set, "shortfall")
    cs[k] <- diff(zr(one_set[[1]]))
  }
  a <- zr(seq_along(vN))
  o <- zr(which(orph))

  data.frame(patient_id = pts[i], n_kept = length(vN), n_orphan = sum(orph),
             ## excess = -zero_gap_asym という恒等式が成り立つので、本体の
             ## per-patient 出力の zero_gap_asym と符号を反転して一致するはず。
             ## 一致しなければ集合の構成が本体とずれている。
             zr_all_normal = a[1], zr_all_disease = a[2], d_zr_all = diff(a),
             zr_orph_normal = o[1], zr_orph_disease = o[2], d_zr_orphan = diff(o),
             d_zr_control = mean(cs), d_zr_control_sd = sd(cs),
             shortfall = short, stringsAsFactors = FALSE)
}

## トップレベルでは else を行頭に置けない（構文エラーになる）。波括弧で包む。
lst <- if (N_CORES > 1L) {
  parallel::mclapply(seq_len(np), one, mc.cores = N_CORES)
} else {
  lapply(seq_len(np), one)
}
## mclapply は子プロセスのエラーを例外ではなく try-error オブジェクトとして
## 返すので、rbind すると意味不明な失敗になる。ここで捕まえる。
bad <- which(!vapply(lst, is.data.frame, logical(1)))
if (length(bad)) {
  msg("失敗した患者 ", length(bad), " 人。最初の内容:")
  msg("  ", paste(utils::head(as.character(lst[[bad[1]]]), 1), collapse = " "))
  stop("N_CORES=1 で走らせると原因が見えます。", call. = FALSE)
}
d <- do.call(rbind, lst)
d$excess <- d$d_zr_orphan - d$d_zr_control

hr(); msg(CO, " / 患者 ", np, " 人 / 対照集合 ", NSET, " 反復の平均",
          if (length(EXCL)) paste0(" / 除外 ", paste(EXCL, collapse = ",")) else "")
msg(sprintf("  解析対象 %d 転写産物（中位、(normal>0)|(disease>0) で絞った後）/ orphan %d 本",
            round(median(d$n_kept)), round(median(d$n_orphan))))
msg("  ゼロ率はすべてこの絞り込み後の集合で計算している。全転写産物を分母に")
msg("  すると d_zr_all だけが小さくなり、傾きが一律に大きく出る（比は不変）。")
if (max(d$shortfall) > 0) {
  msg("  注意: shortfall 最大 ", max(d$shortfall),
      " 件（", NSET, " 回の合計）。いずれかの層で非 orphan の候補が足りていません。")
  msg("  N_STRATA を下げて再実行し、増幅率が動かないことを確認すること。")
} else {
  msg("  shortfall 0。すべての層で対照を orphan と同数引けている。")
}

q <- function(v) sprintf("中位 %+.4f  範囲 %+.4f 〜 %+.4f", median(v), min(v), max(v))
hr(); msg("ゼロ率の変化（疾患 - 正常）")
msg("  全転写産物  ", q(d$d_zr_all))
msg("  orphan      ", q(d$d_zr_orphan))
msg("  マッチ対照  ", q(d$d_zr_control))
msg("  差 (orphan - 対照)  ", q(d$excess))

## ------------------------------------------------------------------ 回帰
fit_o <- lm(d_zr_orphan  ~ d_zr_all, data = d)
fit_c <- lm(d_zr_control ~ d_zr_all, data = d)
sl <- function(f) {
  s <- summary(f)$coefficients
  c(int = s[1, 1], int_se = s[1, 2], int_p = s[1, 4],
    slope = s[2, 1], slope_se = s[2, 2], slope_p = s[2, 4],
    r2 = summary(f)$r.squared)
}
so <- sl(fit_o); sc <- sl(fit_c)

hr(); msg("回帰 （応答 ~ 全体のゼロ率変化 d_zr_all。患者 ", np, " 人が 1 点）")
msg(sprintf("  orphan   傾き %.3f ± %.3f   切片 %+.4f ± %.4f (p = %.3g)   R2 %.3f",
            so["slope"], so["slope_se"], so["int"], so["int_se"], so["int_p"], so["r2"]))
msg(sprintf("  対照     傾き %.3f ± %.3f   切片 %+.4f ± %.4f (p = %.3g)   R2 %.3f",
            sc["slope"], sc["slope_se"], sc["int"], sc["int_se"], sc["int_p"], sc["r2"]))
msg("  ※ 論文の表に載せる切片は下の fit_e のものであって、この 2 つではない。")

## ---------------------------------------------------------------- 傾きの差
## これが主たる統計量である。患者ごとの超過 (orphan - 対照) を全体の移動に
## 回帰する。両者は同一患者・同一ライブラリから出るので対応がとれており、
## 差をとった時点で患者固有の水準は消えている。
##
## **論文の表の「傾きの差」と「切片」はどちらもこの 1 本の回帰から出る。**
## 表の切片の列が 1 つしかないのはそのためである。
##
## 生の超過の平均（下の対応のある t 検定）ではない。患者ごとに d_zr_all が
## 大きく違うため、条件づけないとその分散が差を埋める。
fit_e <- lm(excess ~ d_zr_all, data = d)
ce <- summary(fit_e)$coefficients
hr(); msg("傾きの差と切片（**論文の表に載るのはこの回帰**）")
msg(sprintf("  傾き (orphan - 対照) = %.3f ± %.3f   p = %.3g",
            ce[2, 1], ce[2, 2], ce[2, 4]))
msg(sprintf("  切片                 = %+.4f ± %.4f   p = %.3g",
            ce[1, 1], ce[1, 2], ce[1, 4]))
msg(sprintf("  R2 = %.3f", summary(fit_e)$r.squared))
msg("  切片は「全体の検出移動が 0 のときの超過」である。3 コホートの実測では")
msg("  これが 0 と区別でき（プールで z = 9.2、異質性 Q p = 0.163）、しかも")
msg("  コホート間で一定だった。**これが論文の中核の主張である。**")
msg("")
msg(sprintf("  増幅率 = orphan の傾き / 対照の傾き = %.3f / %.3f = %.2f",
            so["slope"], sc["slope"], so["slope"] / sc["slope"]))
msg("  増幅率は d_zr_all の分母の取り方に不変なので、コホート間の比較は")
msg("  この値で行う。ただし増幅率は検出されている orphan 数と完全に順位相関")
msg("  する（3 コホートで rho = 1）ため、ライブラリの検出感度との交絡が")
msg("  切れない。**増幅率をプールしたり見出しにしたりしないこと。**")

## -------------------------------------- 切片の解釈に必要な確認（rev.36 で追加）
## (1) 切片は d_zr_all = 0 における超過である。患者の d_zr_all が 0 を挟んで
##     いなければ、データのない点への外挿を主張の中核に据えていることになる。
## (2) 対照集合の抽出ばらつきが、回帰の残差のどれだけを占めるか。
##     **この成分は切片の ± に既に含まれている**（上のヘッダ参照）ので、
##     足し込む必要はない。NSET が足りているかの確認として見る。
rng <- range(d$d_zr_all)
inside <- rng[1] < 0 && rng[2] > 0
hr(); msg("切片の解釈に必要な確認")
msg(sprintf("  d_zr_all の範囲  %+.5f 〜 %+.5f", rng[1], rng[2]))
msg(sprintf("  切片を取る点 (d_zr_all = 0) はデータ範囲の%s",
            if (inside) "**内側**。内挿である。" else "**外側**。外挿である。"))
if (!inside) {
  msg(sprintf("    0 に最も近い患者の d_zr_all = %+.5f",
              d$d_zr_all[which.min(abs(d$d_zr_all))]))
  msg("    ** Limitations に「切片は観測範囲外への外挿である」と明記すること。**")
  msg("    査読者に先に指摘されると主張の中核が揺らぐので、自分で書く。")
}
se_draw <- mean(d$d_zr_control_sd) / sqrt(NSET)
s_resid <- summary(fit_e)$sigma
msg(sprintf("  対照集合の抽出による d_zr_control の SE = %.6f（%d 組の平均）",
            se_draw, NSET))
msg(sprintf("  回帰の残差 sd %.5f に対して %.2f%%（残差分散の %.3f%%）",
            s_resid, 100 * se_draw / s_resid, 100 * (se_draw / s_resid)^2))
msg("  **この成分は切片の ± に既に含まれている。** excess の測定誤差として")
msg("  残差に入り、切片の標準誤差に伝播しているので、二乗和で足し込むのは")
msg("  二重計上である。参考: 切片の回帰 SE = ", sprintf("%.6f", ce[1, 2]))
if ((se_draw / s_resid)^2 > 0.05) {
  msg("    ** 残差分散の 5% 超を抽出ノイズが占める。NSET を増やすこと。**")
} else {
  msg("    残差のほとんどは患者間の真のばらつきである。NSET は足りている。")
}

## ------------------------------------------------------------- 対応のある検定
## ------------------------------------------- 生の差（参考。条件づけない形）
hr(); msg("生の超過の平均（参考。d_zr_all で条件づけない形）")
tt <- t.test(d$d_zr_orphan, d$d_zr_control, paired = TRUE)
wt <- suppressWarnings(wilcox.test(d$d_zr_orphan, d$d_zr_control, paired = TRUE))
msg(sprintf("  平均差 %+.5f   95%% CI %+.5f 〜 %+.5f   対応のある t 検定 p = %.3g",
            unname(tt$estimate), tt$conf.int[1], tt$conf.int[2], tt$p.value))
msg(sprintf("  Wilcoxon 符号付き順位 p = %.3g", wt$p.value))
msg(sprintf("  差が + の患者 %d / %d", sum(d$excess > 0), np))
msg("  分布は裾が重いので t 検定ではなく Wilcoxon と符号で判断すること。")
msg("")
msg("  **コホートをまたいでプールしないこと（rev.27・rev.35 で撤回）。**")
msg("  実測は 3 コホートで 18/24、43/57、57/70 が正だが、効果量に異質性が")
msg("  あるため（患者単位 log2 比の Q 検定 p = 0.001）、合計して符号検定を")
msg("  かけた値（118/151、p = 1.0e-12）は**撤回済み**である。使わないこと。")
msg("  論文では 3 コホートを個別に報告する。各コホートで既に有意なので、")
msg("  プールは何も足さない。")

## 比（全体の移動で割った増幅率）。全体の移動が 0 に近い患者は除く。
use <- abs(d$d_zr_all) > 1e-4
if (sum(use) >= 6L) {
  ro <- d$d_zr_orphan[use]  / d$d_zr_all[use]
  rc <- d$d_zr_control[use] / d$d_zr_all[use]
  tr <- t.test(ro, rc, paired = TRUE)
  msg("")
  msg(sprintf("  増幅率（応答 / 全体）  orphan %.3f ± %.3f   対照 %.3f ± %.3f   （%d 人）",
              mean(ro), sd(ro)/sqrt(length(ro)), mean(rc), sd(rc)/sqrt(length(rc)), sum(use)))
  msg(sprintf("  対応のある差 %.3f ± %.3f   p = %.3g",
              unname(tr$estimate), tr$stderr, tr$p.value))
} else {
  msg("")
  msg("  全体のゼロ率変化が小さい患者が多く、増幅率は計算していません。")
}

## 除外つきの結果が既定の結果を上書きしないようにする。
out <- paste0("slope_", CO, if (length(EXCL)) "_excl" else "", ".csv")
write.csv(d, out, row.names = FALSE)
hr(); msg("書き出し: ", out)
msg("")
msg("読み方（rev.36 で更新）")
msg("  対照の傾きが 1 付近 → 対照は全体の検出移動をそのまま写しているだけ。")
msg("  orphan の傾きがそれを超える → orphan は全体の移動を増幅している。")
msg("  **切片が 0 と区別できる → 全体が動かないときでも orphan だけが下がる。**")
msg("")
msg("  実測では切片が 3 コホートで一定かつ 0 でなかった（プール +0.0544、")
msg("  z = 9.2、異質性 Q p = 0.163）。したがって論文の中核の主張は")
msg("  **「全体の検出移動がゼロでも、orphan は発現水準と検出頻度を揃えた")
msg("  対照より系統的に検出を失う」**である。増幅は副次的な所見として書く。")
msg("")
msg("  増幅率（傾きの比）はコホート間で揃わず（Q p = 5.0e-09）、検出されて")
msg("  いる orphan 数と完全に順位相関する。疾患差として解釈してはならない。")
msg("  傾きの差は 1 患者の影響を強く受ける（GSE144269 で浅いライブラリ 1 人を")
msg("  抜くと 63% 動く）。切片は −1%〜+13% で相対的に安定している。")
msg("")
msg("  主張は測定水準にとどめること。「orphan の転写が疾患で下がる」と書くと")
msg("  配列依存的な定量バイアス（長さ・GC）の可能性に対して脆くなる。")
msg("  マッチングは発現水準と検出頻度の 2 次元のみで、配列上の性質は揃えて")
msg("  いない。Abstract と Discussion でこの線を越えないこと。")
