#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## diag_mismatch_slope.R
##
## **中核の主張（切片）に対するペア崩し陰性対照。**
##
## 方針書 §0 (B)(C) のペア崩し陰性対照は β3 と生の orph_excess に対して行った。
## **切片に対してはまだ行っていない。** ところが rev.38 以降、論文の主張は
## 切片 +0.0544 ± 0.0059 に一本化されている。**主張そのものに陰性対照が
## ついていない状態である。**
##
## しかも 2026-09-20 の監査で、残す 3 コホートの「生の」orph_excess は
## 同一条件のランダムな組と区別できないことが分かった
## （経験的 p = 0.125 / 0.085 / 0.244、sd 比 1.00 / 0.89 / 0.98）。
## つまり効果は **d_zr_all で条件づけ、マッチ対照を引いて初めて現れる**。
## 条件づけが効果を作っていないことを示す必要がある。
##
## 【やること】
##   患者 i の normal と患者 perm(i) の disease を組にして（perm は不動点なしの
##   置換）、diag_slope.R と**完全に同じ手続き**で切片を計算する。これを NMIS 回
##   繰り返し、本物の切片と比べる。
##
##   崩したペアでも +0.05 前後が出る → **主張は崩れる。** 条件づけの産物である。
##   崩したペアでは 0 付近      → これまでで最も強い証拠になる。
##
## 【設計上の注意】
##   - normal 側は患者 i 自身なので、層別変数の leave-one-out は **i** を抜く
##     （diag_slope.R と同じ）。disease を入れ替えても対照集合の作り方は変わらない。
##     したがって変わるのは disease 側だけであり、§0 (B) の設計と同じ論理である。
##   - コホート全体の条件差（あれば）は崩したペアにも同じだけ入る。したがって
##     この検定が問うのは「疾患効果があるか」ではなく
##     **「患者内の対応が効いているか」**である。
##     **コホート = バッチのコホート（GSE40419）では崩しても大きく出るはずで、
##     それは正しい挙動である。** 3 コホートとの対比に使う。
##   - **本物の切片も同じ NSET で計算し直して比べる。** 論文の値は NSET=200、
##     ここの既定は NSET=20。NSET が違う値どうしを比べてはならない。
##
## 【使い方】リポジトリのルートで
##   NSET=20 NMIS=20 N_CORES=12 KALLISTO_ROOT=Revised/GSE244679 \
##   METADATA_DIR=Revised/metadata COHORT=GSE244679 Rscript diag_mismatch_slope.R
##
##   所要は diag_slope.R の約 (NMIS+1) 倍。NSET=20 なら 1 コホート数分程度。
##   NSET=200 にすると 10 倍かかるので、まず既定で傾向を見ること。
##
## 【出力】
##   mismatch_slope_<COHORT>.csv   反復ごとの切片・傾きの差・増幅率
##   標準出力に本物との比較
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

CO        <- Sys.getenv("COHORT", names(COHORTS)[1])
N_STRATA2 <- as.integer(Sys.getenv("N_STRATA2", "4"))
NSET      <- as.integer(Sys.getenv("NSET", "20"))
NMIS      <- as.integer(Sys.getenv("NMIS", "20"))
N_CORES   <- as.integer(Sys.getenv("N_CORES", "1"))

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

## 本体から層別化と抽出を借りる（diag_slope.R と同じ実装を使う）
SRC <- if (file.exists("R/07_interaction_permutation.R"))
         "R/07_interaction_permutation.R" else "07_interaction_permutation.R"
src <- readLines(SRC)
eval(parse(text = paste(src[grep("^make_strata <- function", src):
                            (grep("^run_patient <- function", src) - 1L)],
                        collapse = "\n")))
stopifnot(exists("make_strata_multi"), exists("draw_matched_sets"))
## R/config.R から来る定数。無いと make_strata_multi の中で分かりにくく落ちる。
if (!exists("N_STRATA"))
  stop("N_STRATA がありません。R/config.R が読めているか確認してください。",
       call. = FALSE)

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
np  <- length(pts)
stopifnot(np >= 6L)

normal_mat   <- dat$tpm[, ln, drop = FALSE]
sum_normal   <- rowSums(normal_mat)
count_normal <- rowSums(normal_mat > 0)

## ---------------------------------------------------------------------------
## 1 患者分。j が disease 側の患者番号（j == i なら本物のペア）。
## diag_slope.R の one() と同じ手続きで、disease 側だけ差し替える。
## ---------------------------------------------------------------------------
one <- function(i, j, seed_off) {
  vN <- dat$tpm[, ln[i]]; vD <- dat$tpm[, ld[j]]
  keep <- (vN > 0) | (vD > 0)
  vN <- vN[keep]; vD <- vD[keep]
  orph <- is_o[keep]

  ## 層別変数は normal 側だけから作る。normal は患者 i 自身なので i を抜く。
  ext  <- ((sum_normal - dat$tpm[, ln[i]]) / max(np - 1L, 1L))[keep]
  freq <- ((count_normal - (dat$tpm[, ln[i]] > 0)) / max(np - 1L, 1L))[keep]
  st   <- make_strata_multi(list(ext, freq), c(N_STRATA, N_STRATA2))

  set.seed(BASE_SEED + seed_off * 100003L + i)
  zr <- function(idx) c(mean(vN[idx] == 0), mean(vD[idx] == 0))
  cs <- numeric(NSET)
  for (k in seq_len(NSET)) cs[k] <- diff(zr(draw_matched_sets(st, orph, n_sets = 1L)[[1]]))

  a <- zr(seq_along(vN)); o <- zr(which(orph))
  c(d_zr_all = diff(a), d_zr_orphan = diff(o), d_zr_control = mean(cs))
}

## 1 反復分。perm は disease 側の並べ替え。
fit_one_rep <- function(perm, seed_off) {
  lst <- if (N_CORES > 1L) {
    parallel::mclapply(seq_len(np), function(i) one(i, perm[i], seed_off),
                       mc.cores = N_CORES)
  } else lapply(seq_len(np), function(i) one(i, perm[i], seed_off))
  bad <- which(!vapply(lst, is.numeric, logical(1)))
  if (length(bad)) return(NULL)
  d <- as.data.frame(do.call(rbind, lst))
  d$excess <- d$d_zr_orphan - d$d_zr_control
  fe <- lm(excess ~ d_zr_all, data = d)
  fo <- lm(d_zr_orphan  ~ d_zr_all, data = d)
  fc <- lm(d_zr_control ~ d_zr_all, data = d)
  se <- summary(fe)$coefficients
  ## unname を忘れると c(slope_orphan = coef(fo)[2]) が
  ## "slope_orphan.d_zr_all" という列名になる（R の名前付き結合の仕様）。
  c(intercept = se[1, 1], intercept_se = se[1, 2],
    slope_diff = se[2, 1], slope_diff_se = se[2, 2],
    slope_orphan = unname(coef(fo)[2]), slope_control = unname(coef(fc)[2]),
    amp = unname(coef(fo)[2] / coef(fc)[2]),
    median_excess = median(d$excess))
}

## 不動点のない置換（derangement）。単純な棄却法で十分。
derange <- function(n) {
  repeat { p <- sample(n); if (!any(p == seq_len(n))) return(p) }
}

hr(); msg(CO, " / 患者 ", np, " 人 / NSET = ", NSET, " / 崩し ", NMIS, " 反復")
msg("  **NSET を論文の 200 から下げている場合、本物の値も同じ NSET で計算し直す。**")
msg("  NSET が違う値どうしを比べてはならない。")

## --- 本物（恒等置換）
msg(""); msg("本物のペアを計算中...")
truth <- fit_one_rep(seq_len(np), 0L)
if (is.null(truth)) stop("本物のペアの計算に失敗しました。N_CORES=1 で確認してください。",
                         call. = FALSE)
msg(sprintf("  切片 %+.5f ± %.5f   傾きの差 %+.4f ± %.4f   増幅率 %.3f",
            truth["intercept"], truth["intercept_se"],
            truth["slope_diff"], truth["slope_diff_se"], truth["amp"]))

## --- 崩したペア
msg(""); msg("崩したペアを計算中（", NMIS, " 反復）...")
set.seed(BASE_SEED)
perms <- lapply(seq_len(NMIS), function(r) derange(np))
res <- vector("list", NMIS)
for (r in seq_len(NMIS)) {
  res[[r]] <- fit_one_rep(perms[[r]], r)
  cat(sprintf("\r  %d / %d", r, NMIS)); flush.console()
}
cat("\n")
res <- res[!vapply(res, is.null, logical(1))]
if (!length(res)) stop("崩したペアがすべて失敗しました。", call. = FALSE)
M <- as.data.frame(do.call(rbind, res))

## ---------------------------------------------------------------------------
hr(); msg("切片 ― 本物 vs 崩したペア")
ti <- unname(truth["intercept"]); mi <- M$intercept
msg(sprintf("  本物          %+.5f", ti))
msg(sprintf("  崩したペア    中位 %+.5f   sd %.5f   範囲 %+.5f 〜 %+.5f   （%d 反復）",
            median(mi), sd(mi), min(mi), max(mi), length(mi)))
p_emp <- (1 + sum(abs(mi) >= abs(ti))) / (length(mi) + 1)
msg(sprintf("  両側の経験的 p = %.4g   （下限 %.4g）", p_emp, 1 / (length(mi) + 1)))
msg(sprintf("  比（崩し中位 / 本物） = %.3f", median(mi) / ti))

hr(); msg("傾きの差 ― 本物 vs 崩したペア")
td <- unname(truth["slope_diff"]); md <- M$slope_diff
msg(sprintf("  本物          %+.4f", td))
msg(sprintf("  崩したペア    中位 %+.4f   sd %.4f   範囲 %+.4f 〜 %+.4f",
            median(md), sd(md), min(md), max(md)))
msg(sprintf("  両側の経験的 p = %.4g",
            (1 + sum(abs(md) >= abs(td))) / (length(md) + 1)))

hr(); msg("増幅率 ― 本物 vs 崩したペア")
msg(sprintf("  本物 %.3f   崩したペア 中位 %.3f（%.3f 〜 %.3f）",
            truth["amp"], median(M$amp), min(M$amp), max(M$amp)))

out <- rbind(data.frame(rep = 0L, kind = "true",  t(truth)),
             data.frame(rep = seq_len(nrow(M)), kind = "mismatched", M))
write.csv(out, sprintf("mismatch_slope_%s.csv", CO), row.names = FALSE)
hr(); msg("書き出し: ", sprintf("mismatch_slope_%s.csv", CO))

msg("")
msg("読み方")
msg("")
msg("  **崩したペアの切片が本物と同程度なら、主張は崩れる。**")
msg("  その場合、切片は患者内の対応ではなくコホート全体の条件差を測っている。")
msg("  d_zr_all で条件づけてマッチ対照を引く手続きが、効果を作っていることになる。")
msg("")
msg("  崩したペアの切片が 0 付近なら、**中核の主張に対する陰性対照が通った**")
msg("  ことになる。方針書 §0 (B)(C) の陰性対照は β3 と生の orph_excess に")
msg("  対するもので、切片に対しては存在しなかった。これがその穴を埋める。")
msg("")
msg("  **GSE40419 は条件 = バッチなので、崩しても大きく出るはずである。**")
msg("  それは想定どおりの挙動であって、このコホートを救う材料ではない。")
msg("  3 コホートで崩すと消え、GSE40419 では消えない、という対比が")
msg("  「残した 3 コホートには交絡が無い」ことの最も直接的な証拠になる。")
msg("")
msg("  結果は Results ではなく **Methods の妥当性確認**として書き、")
msg("  Limitations にも 1 文置く。査読者に先に問われる前に自分で出す。")
