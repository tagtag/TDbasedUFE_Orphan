## ---------------------------------------------------------------------------
## 07_interaction_permutation.R
##
## 患者ごとの disease x transcript class 交互作用と、
## 発現量マッチしたランダム集合のアンサンブルから作る経験的 p 値。
##
##   Reviewer 1-1 : 交互作用項 beta3 を患者ごとに推定する
##   Reviewer 1-2 : var(beta3) を患者間分散として報告する
##   Reviewer 3-2 : 名目 p 値ではなくランダム集合の帰無分布で判定する
##
## 既存のリポジトリ構成に乗せる。R/config.R と R/00_functions.R から
##   read_sheet(), load_cohort_tpm(), orphan_index(), transform_expression()
##   KALLISTO_ROOT, RESULTS_DIR, ORPHAN_ID_FILE, ALPHA, N_STRATA, BASE_SEED,
##   N_TRANSCRIPTS, N_ORPHAN_TRANSCRIPT, COHORTS
## を利用する。
##
## 実行:  Rscript R/07_interaction_permutation.R
##
## 環境変数
##   N_PERM      既定 1000        患者あたりの抽出回数
##   TRANSFORMS  既定 "scale,rank"
##   MATCH_ON    既定 "external2" 下の「マッチング変数」を参照。
##                              カンマ区切りで複数指定すると1回の実行で
##                              すべてを計算し比較表を出す。
##                              例: MATCH_ON="external,pair_mean"
##   COHORT      既定 空          走らせるコホート。カンマ区切り。空なら全部。
##                              例: COHORT="GSE244679"
##   N_CORES     既定 1           parallel::mclapply の並列数（Windows では 1）
##
## config.R が環境変数で受ける KALLISTO_ROOT / METADATA_DIR / RESULTS_DIR も
## そのまま効く。手元の配置に合わせる場合は make_revised_layout.R を参照。
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

N_PERM     <- as.integer(Sys.getenv("N_PERM", "1000"))
TRANSFORMS <- strsplit(Sys.getenv("TRANSFORMS", "scale,rank"), ",")[[1]]
MATCH_ONS  <- trimws(strsplit(Sys.getenv("MATCH_ON", "external2"), ",")[[1]])
MATCH_ONS  <- MATCH_ONS[nzchar(MATCH_ONS)]
N_CORES    <- as.integer(Sys.getenv("N_CORES", "1"))
## external2 の第 2 次元（他患者での検出頻度）の層数。4 で shortfall が
## 出ず、ゼロ率の不一致がほぼ解消することを合成データで確認している。
N_STRATA2  <- as.integer(Sys.getenv("N_STRATA2", "4"))
ZERO_RULE  <- "or"          # 既存の Table 4/5 と同じ（どちらかが非ゼロなら残す）
OUT_PREFIX <- file.path(RESULTS_DIR, "interaction")

## 一部のコホートだけ走らせる場合。手元に kallisto 出力が揃っていない
## コホートを黙って飛ばすことはせず、明示的に指定させる。
COHORT_SEL <- trimws(strsplit(Sys.getenv("COHORT", ""), ",")[[1]])
COHORT_SEL <- COHORT_SEL[nzchar(COHORT_SEL)]
if (length(COHORT_SEL)) {
  unknown <- setdiff(COHORT_SEL, names(COHORTS))
  if (length(unknown))
    stop("COHORT に config.R の COHORTS にない名前があります: ",
         paste(unknown, collapse = ", "), "\n  使える名前: ",
         paste(names(COHORTS), collapse = ", "), call. = FALSE)
  COHORTS <- COHORTS[COHORT_SEL]
}

stopifnot(all(MATCH_ONS %in% c("external2", "external", "pair_mean", "normal")))
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

## 経験的 p の下限は 1/(N_PERM + 1)。これが ALPHA に近いと、BH 後の q が
## どの患者も ALPHA を割れず、direction が全員 "none" になる。効果が無いのか
## 置換回数が足りないのか区別がつかなくなるので、ここで止める。
.p_floor <- 1 / (N_PERM + 1)
if (.p_floor > ALPHA / 5) {
  warning(sprintf(
    paste0("N_PERM = %d では p の下限が %.4g、ALPHA = %.3g に近すぎる。\n",
           "  direction が全員 none になっても効果が無い証拠にはならない。\n",
           "  N_PERM >= %d を推奨。この実行では medZ / frac05 で判断すること。"),
    N_PERM, .p_floor, ALPHA, ceiling(5 / ALPHA)), call. = FALSE)
}

## ---------------------------------------------------------------------------
## 符号の約束
##
## beta3 は cond ダミーを disease = 1 として定義する。すなわち
##
##     beta3 = [orphan の disease - normal] - [control の disease - normal]
##
## 正の beta3 は「対照に比べ orphan が疾患側で上昇」を意味する。
## 既存の run_cohort() が返す mean_d は normal - disease なので符号が逆である。
## 出力には beta3 と beta3_normal_minus_disease の両方を出す。
##
## ---------------------------------------------------------------------------
## マッチング変数（MATCH_ON）
##
## 合成データでの検証結果（患者12人、各シナリオ、p<0.05 となる患者の割合）:
##
##   シナリオ                     pair_mean   normal   external
##   疾患効果なし（偽陽性）            0.00     0.58       0.00
##   組成シフトのみ（偽陽性）          0.08     0.33       0.08
##   低発現全体が上昇（偽陽性）        1.00     0.08       0.08
##   orphan 固有の効果（検出力）       1.00     1.00       1.00
##
##  pair_mean : その患者の対平均 (normal + disease)/2。既存の draw_control_set()
##              と同じ。疾患効果を受けた転写産物は対平均が押し上げられて別の層へ
##              移るため、低発現域全体が疾患側で動く状況では対照が揃わず、
##              偽陽性が 100% に達する。
##  normal    : その患者の normal のみ。対照だけがその患者の雑音で条件付けられる
##              ため disease 側で平均回帰が起き、疾患効果ゼロでも偽陽性が出る。
##  external  : その患者を除いた同一コホートの normal の平均（leave-one-out）。
##              その患者の疾患状態にも測定雑音にも依存しないので、効果を吸収せず
##              選択による平均回帰も起こさない。ただし平均だけでは足りない
##              （下記）。
##  external2 : external に「他患者での検出頻度」を加えた 2 次元の層化。推奨。
##
## なぜ external だけでは足りないか（実データで観測）:
##   平均が同じでも「少数の患者で高く大半でゼロ」の転写産物と「どの患者でも
##   中程度」の転写産物は同じ層に入る。orphan は前者に偏るため、その患者の
##   normal で orphan のゼロ率が対照より高く出る。GSE244679 の 1 患者では
##   normal で orphan 0.531 / 対照 0.216、disease では 0.075 / 0.069 と
##   ほぼ一致していた。この非対称性がそのまま差の差 beta3 を作る。
##
##   検出頻度も揃えた合成データでの検証（患者24人、疾患効果ゼロ、
##   disease 側の検出が飽和する条件）:
##     マッチング              normal のゼロ率差
##     external  (1 次元)      -0.108
##     external2 (2 次元)      -0.008
##   実データの -0.0087 -> +0.00028 と同じ向き・同じ桁の改善である。
##   （2026-09-23 に test_2d.R を再実行して確認。以前ここに書いていた
##   「偽陽性 4/24 対 1/24」は再現せず、4 行すべて 1/24 だった。さらに
##   この検定の p は較正されていないので人数を根拠に使わないこと。）
##   shortfall はどちらも 0。層数は 20 -> 21 / 74。
##
## 既存の Table 5 と揃えた比較が要る場合は MATCH_ON に両方を並べて
##     MATCH_ON="external,pair_mean" Rscript R/07_interaction_permutation.R
## と 1 回で走らせる。出力ファイル名に match_on が入るので上書きは起きず、
## results/interaction_match_on_comparison.csv に並べた比較表が出る。
## ---------------------------------------------------------------------------

## 4つのセル平均から交互作用係数を閉形式で求める。
## lm(y ~ class * cond) の係数と厳密に一致する（verify_against_lm で確認）。
cell_coefs <- function(yN, yD, idx_o, idx_c) {
  c_N <- mean(yN[idx_c]); c_D <- mean(yD[idx_c])
  o_N <- mean(yN[idx_o]); o_D <- mean(yD[idx_o])
  c(b0 = c_N,                       # 対照・normal
    b1 = o_N - c_N,                 # class 主効果（マッチが効いていれば約 0）
    b2 = c_D - c_N,                 # 対照の疾患応答
    b3 = (o_D - o_N) - (c_D - c_N)) # 交互作用
}

## 閉形式が lm と一致することを一度だけ確認し、結果をファイルに残す。
verify_against_lm <- function(yN, yD, idx_o, idx_c, path) {
  long <- data.frame(
    y     = c(yN[idx_o], yD[idx_o], yN[idx_c], yD[idx_c]),
    class = factor(rep(c("orphan", "control"),
                       c(2 * length(idx_o), 2 * length(idx_c))),
                   levels = c("control", "orphan")),
    cond  = factor(c(rep(c("normal", "disease"), each = length(idx_o)),
                     rep(c("normal", "disease"), each = length(idx_c))),
                   levels = c("normal", "disease"))
  )
  fit  <- lm(y ~ class * cond, data = long)
  cf   <- coef(fit)
  mine <- cell_coefs(yN, yD, idx_o, idx_c)
  writeLines(c(
    "lm(y ~ class * cond) と閉形式の照合",
    "",
    capture.output(summary(fit)),
    "",
    "閉形式:",
    sprintf("  %-4s %+.10f", names(mine), mine),
    "",
    sprintf("交互作用の差: %.3e", abs(unname(cf[4]) - unname(mine["b3"])))
  ), path)
  stopifnot(abs(unname(cf[4]) - unname(mine["b3"])) < 1e-8)
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## 発現量マッチした層別抽出
##
## 既存の draw_control_set() と同じ層別（厳密ゼロを第1層、正の値を n_strata
## 分位で分割）を再現したうえで、1回の sample() で複数の互いに素な集合を
## 引けるようにしたもの。帰無分布には互いに素な2集合が要るため必要。
## ---------------------------------------------------------------------------

make_strata <- function(m, n_strata) {
  s <- integer(length(m))
  pos <- m > 0
  s[!pos] <- 1L
  if (any(pos)) {
    br <- unique(quantile(m[pos], probs = seq(0, 1, length.out = n_strata + 1),
                          na.rm = TRUE, type = 7))
    ## 変数が縮退している（正の値がすべて同一など）と unique() が 1 点に
    ## なる。cut() はそれを「区間の個数」と解釈して落ちるので、まとめて
    ## 1 層にする。第 2 次元の検出頻度で実際に起こりうる。
    if (length(br) < 2L) { s[pos] <- 2L; return(s) }
    br[1] <- -Inf; br[length(br)] <- Inf
    s[pos] <- as.integer(cut(m[pos], breaks = br, include.lowest = TRUE)) + 1L
  }
  s
}

## 各層で orphan と同数を非 orphan から復元なしで抽出し、
## 互いに素な n_sets 個の集合を返す。
draw_matched_sets <- function(strata, is_orphan, n_sets = 1L) {
  out <- replicate(n_sets, integer(0), simplify = FALSE)
  short <- 0L
  for (s in sort(unique(strata))) {
    need <- sum(strata == s & is_orphan)
    if (need == 0L) next
    pool <- which(strata == s & !is_orphan)
    take <- need * n_sets
    if (length(pool) < take) { short <- short + (take - length(pool)); take <- length(pool) }
    if (take == 0L) next
    picked <- pool[sample.int(length(pool), take)]
    cuts <- split(picked, rep(seq_len(n_sets), length.out = take))
    for (k in seq_along(cuts)) out[[k]] <- c(out[[k]], cuts[[k]])
  }
  attr(out, "shortfall") <- short
  out
}

## ---------------------------------------------------------------------------
## 患者 1 人分
##
##  観測: orphan 対 マッチ集合 A_r        (r = 1..n_perm)
##  帰無: マッチ集合 B_r 対 マッチ集合 C_r (B_r と C_r は互いに素)
##
## 抽出は transform に依存しないので 1 度だけ行い、両方の表現で使い回す。
## ---------------------------------------------------------------------------

## 複数のマッチング変数の層を掛け合わせ、1..K に詰め直す。
make_strata_multi <- function(vars, n_strata_vec) {
  ids <- lapply(seq_along(vars), function(j) make_strata(vars[[j]], n_strata_vec[j]))
  out <- ids[[1]]
  for (j in seq_along(ids)[-1]) out <- (out - 1L) * max(ids[[j]]) + ids[[j]]
  match(out, sort(unique(out)))
}

run_patient <- function(tpm_n, tpm_d, match_value, is_orphan, transforms,
                        n_perm, n_strata, seed, match_on) {
  keep <- (tpm_n > 0) | (tpm_d > 0)               # zero_rule = "or"
  if (!any(keep)) return(NULL)

  xn <- tpm_n[keep]; xd <- tpm_d[keep]
  mvs <- if (is.list(match_value)) lapply(match_value, `[`, keep)
         else list(match_value[keep])
  orph <- is_orphan[keep]
  if (sum(orph) < 2L || sum(!orph) < 2L) return(NULL)

  strata <- make_strata_multi(mvs, rep_len(n_strata, length(mvs)))

  set.seed(seed)
  idx_obs  <- vector("list", n_perm)
  idx_null <- vector("list", n_perm)
  short <- 0L
  for (r in seq_len(n_perm)) {
    a  <- draw_matched_sets(strata, orph, 1L)
    bc <- draw_matched_sets(strata, orph, 2L)
    short <- short + attr(a, "shortfall") + attr(bc, "shortfall")
    idx_obs[[r]] <- a[[1]]; idx_null[[r]] <- bc
  }

  idx_o <- which(orph)

  ## マッチングの達成度。normal 側と disease 側でともに約 1 であること。
  ## 両者が 1 から反対向きに外れる場合、マッチングが疾患効果を吸収している。
  ratio_n <- mean(vapply(idx_obs, function(a) mean(xn[a]), 0)) / mean(xn[idx_o])
  ratio_d <- mean(vapply(idx_obs, function(a) mean(xd[a]), 0)) / mean(xd[idx_o])

  ## ゼロ率の一致。平均 TPM の比は裾の数本に支配されるので単独では
  ## 判断材料にならない。マッチングの失敗はまずここに出る。
  ## とくに zero_gap_normal と zero_gap_disease が違う値になっていたら、
  ## その非対称性がそのまま beta3 を作るので beta3 は信用できない。
  z0_on <- mean(xn[idx_o] == 0); z0_od <- mean(xd[idx_o] == 0)
  z0_cn <- mean(vapply(idx_obs, function(a) mean(xn[a] == 0), 0))
  z0_cd <- mean(vapply(idx_obs, function(a) mean(xd[a] == 0), 0))

  res <- list()
  for (tf in transforms) {
    yN <- transform_expression(xn, tf)
    yD <- transform_expression(xd, tf)
    if (all(is.na(yN)) || all(is.na(yD))) next

    obs <- nul <- b1 <- b2 <- numeric(n_perm)
    for (r in seq_len(n_perm)) {
      cf <- cell_coefs(yN, yD, idx_o, idx_obs[[r]])
      obs[r] <- cf["b3"]; b1[r] <- cf["b1"]; b2[r] <- cf["b2"]
      nul[r] <- cell_coefs(yN, yD, idx_null[[r]][[1]], idx_null[[r]][[2]])["b3"]
    }

    t_obs <- mean(obs)                 # 抽出ノイズを平均した観測統計量
    ## 経験的 p 値。+1 補正により下限は 1/(n_perm + 1)。
    p_up   <- (1 + sum(nul >= t_obs)) / (n_perm + 1)
    p_down <- (1 + sum(nul <= t_obs)) / (n_perm + 1)
    p_two  <- (1 + sum(abs(nul - mean(nul)) >= abs(t_obs - mean(nul)))) / (n_perm + 1)

    res[[tf]] <- data.frame(
      transform = tf, zero_rule = ZERO_RULE, match_on = match_on,
      n_kept = length(xn), n_orphan = length(idx_o), n_control = length(idx_obs[[1]]),
      shortfall = short,
      match_ratio_normal = ratio_n, match_ratio_disease = ratio_d,
      n_strata_used = length(unique(strata)),
      zero_orphan_normal = z0_on, zero_control_normal = z0_cn,
      zero_orphan_disease = z0_od, zero_control_disease = z0_cd,
      zero_gap_normal  = z0_cn - z0_on,
      zero_gap_disease = z0_cd - z0_od,
      zero_gap_asym    = (z0_cd - z0_od) - (z0_cn - z0_on),
      b1_mean = mean(b1), b2_mean = mean(b2),
      beta3 = t_obs, beta3_sd = sd(obs),
      beta3_normal_minus_disease = -t_obs,
      null_mean = mean(nul), null_sd = sd(nul),
      ## 帰無分布から何 SD 離れているか。p 値と違って 1/(n_perm+1) で床打ちしない。
      ## 感度分析で match_on 間を比べるときはこちらを見る。
      z_vs_null = (t_obs - mean(nul)) / sd(nul),
      ## 帰無分布は外れ値で裾が重く（z 化した実データで sd/mad の中位 28、最大 214 を観測）、
      ## sd 基準の z は過小に出る。中位数と mad 基準のこちらを使う。
      z_mad = (t_obs - median(nul)) / max(mad(nul), 1e-300),
      null_sd_over_mad = sd(nul) / max(mad(nul), 1e-300),
      p_up = p_up, p_down = p_down, p_two = p_two,
      stringsAsFactors = FALSE
    )
  }
  if (!length(res)) return(NULL)
  do.call(rbind, res)
}

## ---------------------------------------------------------------------------
## コホート 1 つ分
## ---------------------------------------------------------------------------

run_cohort_interaction <- function(cohort, match_on, verify_path = NULL) {
  message("== ", cohort, " / match_on = ", match_on, " ==")
  sheet <- read_sheet(cohort)
  dat   <- load_cohort_tpm(sheet, KALLISTO_ROOT, "abundance.tsv",
                           expected_features = N_TRANSCRIPTS)
  is_orphan <- orphan_index(dat$id, ORPHAN_ID_FILE,
                            expected_n = N_ORPHAN_TRANSCRIPT)

  patients <- unique(sheet$patient_id)
  lib <- function(pid, cond) {
    v <- sheet$library_id[sheet$patient_id == pid & sheet$condition == cond]
    if (length(v) == 1L) v else NA_character_
  }
  ln_all <- vapply(patients, lib, "", cond = "normal")
  ld_all <- vapply(patients, lib, "", cond = "disease")
  ok <- !is.na(ln_all) & !is.na(ld_all)
  if (any(!ok)) warning(cohort, ": ペアが 1 対 1 でない患者を除外 -> ",
                        paste(patients[!ok], collapse = ", "))
  patients <- patients[ok]; ln_all <- ln_all[ok]; ld_all <- ld_all[ok]

  normal_mat <- dat$tpm[, ln_all, drop = FALSE]
  sum_normal <- rowSums(normal_mat)
  ## 他患者での検出頻度（leave-one-out）。external2 の第 2 次元。
  ## 発現量の平均だけで層化すると、平均が同じでも「少数の患者で高く
  ## 大半でゼロ」の転写産物と「どの患者でも中程度」の転写産物が同じ層に
  ## 入り、その患者でのゼロ率が揃わない。orphan は前者に偏るため、
  ## normal 側で orphan のゼロ率が対照より高く出る。検出頻度も揃えると
  ## この不一致がほぼ解消する（合成データで確認済み）。
  count_normal <- rowSums(normal_mat > 0)
  np <- length(patients)

  match_vars <- function(xn, xd) {
    ext <- (sum_normal - xn) / max(np - 1L, 1L)
    switch(match_on,
      external2 = list(ext, (count_normal - (xn > 0)) / max(np - 1L, 1L)),
      external  = list(ext),
      pair_mean = list((xn + xd) / 2),
      normal    = list(xn))
  }

  one <- function(i) {
    xn <- dat$tpm[, ln_all[i]]; xd <- dat$tpm[, ld_all[i]]
    mv <- match_vars(xn, xd)
    out <- run_patient(xn, xd, mv, is_orphan, TRANSFORMS,
                       N_PERM, c(N_STRATA, N_STRATA2), BASE_SEED + i, match_on)
    if (is.null(out)) return(NULL)
    cbind(GEO_ID = cohort, patient_id = patients[i], out, stringsAsFactors = FALSE)
  }

  if (!is.null(verify_path)) {
    xn <- dat$tpm[, ln_all[1]]; xd <- dat$tpm[, ld_all[1]]
    mv <- match_vars(xn, xd)
    keep <- (xn > 0) | (xd > 0)
    xn <- xn[keep]; xd <- xd[keep]; orph <- is_orphan[keep]
    set.seed(BASE_SEED)
    st0 <- make_strata_multi(lapply(mv, `[`, keep),
                             rep_len(c(N_STRATA, N_STRATA2), length(mv)))
    a <- draw_matched_sets(st0, orph, 1L)[[1]]
    verify_against_lm(transform_expression(xn, TRANSFORMS[1]),
                      transform_expression(xd, TRANSFORMS[1]),
                      which(orph), a, verify_path)
    message("  lm との一致を確認: ", verify_path)
  }

  lst <- if (N_CORES > 1L) parallel::mclapply(seq_along(patients), one, mc.cores = N_CORES)
         else lapply(seq_along(patients), one)

  ## mclapply は子プロセスのエラーを停止させず try-error として返す。
  ## NULL だけ弾くと try-error が rbind に流れ込み、原因の分からない
  ## エラーになる。あるいは患者が黙って減る。必ず明示的に検査する。
  bad <- which(!vapply(lst, function(x) is.null(x) || is.data.frame(x), logical(1)))
  if (length(bad)) {
    for (i in utils::head(bad, 5))
      message("  失敗: 患者 ", patients[i], " -> ",
              paste(utils::head(as.character(lst[[i]]), 1), collapse = " "))
    stop(cohort, ": ", length(bad), " 人でエラーが発生しました。",
         "N_CORES=1 で走らせると原因が見えます。", call. = FALSE)
  }
  skipped <- which(vapply(lst, is.null, logical(1)))
  if (length(skipped))
    warning(cohort, ": orphan か対照が 2 本未満で除外した患者 ", length(skipped),
            " 人: ", paste(utils::head(patients[skipped], 5), collapse = ", "),
            call. = FALSE)

  res <- do.call(rbind, Filter(Negate(is.null), lst))
  if (is.null(res) || !nrow(res)) stop(cohort, ": 結果が空です。", call. = FALSE)

  ## BH 補正は既存 run_cohort() と同じくコホート内・表現ごとに行う
  res$q_up <- NA_real_; res$q_down <- NA_real_
  for (tf in unique(res$transform)) {
    k <- res$transform == tf
    res$q_up[k]   <- p.adjust(res$p_up[k],   method = "BH")
    res$q_down[k] <- p.adjust(res$p_down[k], method = "BH")
  }
  res$direction <- ifelse(res$q_up   < ALPHA, "up",
                   ifelse(res$q_down < ALPHA, "down", "none"))
  res
}

## ---------------------------------------------------------------------------
## 実行
## ---------------------------------------------------------------------------

message("N_PERM = ", N_PERM, " / transforms = ", paste(TRANSFORMS, collapse = ", "),
        " / match_on = ", paste(MATCH_ONS, collapse = ", "),
        " / cores = ", N_CORES)

## match_on ごとに独立に全コホートを回す。出力ファイル名には必ず match_on を
## 入れる。入れないと 2 つ目の match_on が 1 つ目を上書きし、比較という
## 本来の目的が壊れる。
res_by_mo <- list(); summ_by_mo <- list()

for (mo in MATCH_ONS) {

  all_res <- list(); first <- TRUE
  for (cohort in names(COHORTS)) {
    vp <- if (first) sprintf("%s_lm_check_%s.txt", OUT_PREFIX, mo) else NULL
    r <- run_cohort_interaction(cohort, match_on = mo, verify_path = vp)
    first <- FALSE
    all_res[[cohort]] <- r
    for (tf in unique(r$transform))
      write.csv(r[r$transform == tf, ],
                sprintf("%s_perpatient_%s_%s_%s.csv", OUT_PREFIX, cohort, tf, mo),
                row.names = FALSE)
  }
  res <- do.call(rbind, all_res)
  rownames(res) <- NULL
  write.csv(res, sprintf("%s_perpatient_all_%s.csv", OUT_PREFIX, mo),
            row.names = FALSE)

  ## -------------------------------------------------------------------------
  ## コホート要約
  ## -------------------------------------------------------------------------
  summ <- do.call(rbind, lapply(
    split(res, list(res$GEO_ID, res$transform), drop = TRUE),
    function(d) data.frame(
      GEO_ID = d$GEO_ID[1], transform = d$transform[1],
      match_on = d$match_on[1],
      n_patients        = nrow(d),
      disease_gt_normal = sum(d$direction == "up"),
      disease_lt_normal = sum(d$direction == "down"),
      no_direction      = sum(d$direction == "none"),
      frac_positive     = mean(d$beta3 > 0),
      var_beta3         = var(d$beta3),   # 患者間分散 = Reviewer 1-2 (a) への回答
      beta3_mean        = mean(d$beta3),  # 診断用。主張には使わない
      ## 閾値で床打ちしない効果量。感度分析はこれで比べる。
      med_abs_zmad      = median(abs(d$z_mad)),      # 外れ値に強い効果量
      med_abs_z         = median(abs(d$z_vs_null)),  # 参考（sd 基準。過小に出る）
      med_sd_over_mad   = median(d$null_sd_over_mad),
      frac_p05          = mean(pmin(d$p_up, d$p_down) < 0.05),
      ## マッチングの合否。gap が 0 から離れ、normal と disease で
      ## 値が違う（asym != 0）場合、beta3 はマッチングの失敗を測っている。
      zero_gap_normal   = mean(d$zero_gap_normal),
      zero_gap_disease  = mean(d$zero_gap_disease),
      zero_gap_asym     = mean(d$zero_gap_asym),
      n_strata_used     = median(d$n_strata_used),
      b1_mean           = mean(d$b1_mean),        # 約 0 ならマッチは妥当
      match_ratio_normal  = mean(d$match_ratio_normal),
      match_ratio_disease = mean(d$match_ratio_disease),
      shortfall_max     = max(d$shortfall),
      stringsAsFactors = FALSE)))
  rownames(summ) <- NULL
  write.csv(summ, sprintf("%s_summary_%s.csv", OUT_PREFIX, mo), row.names = FALSE)

  res_by_mo[[mo]]  <- res
  summ_by_mo[[mo]] <- summ
}

## 全 match_on を縦に積んだもの。これが報告用の一次ファイル。
res_all  <- do.call(rbind, res_by_mo);  rownames(res_all)  <- NULL
summ_all <- do.call(rbind, summ_by_mo); rownames(summ_all) <- NULL
write.csv(res_all,  paste0(OUT_PREFIX, "_perpatient_all.csv"), row.names = FALSE)
write.csv(summ_all, paste0(OUT_PREFIX, "_summary.csv"),        row.names = FALSE)
print(summ_all)

## ---------------------------------------------------------------------------
## match_on 間の比較（感度分析）
## MATCH_ON に 2 つ以上指定したときだけ意味がある。
## ---------------------------------------------------------------------------
if (length(MATCH_ONS) > 1L) {
  key <- paste(summ_all$GEO_ID, summ_all$transform, sep = " / ")
  cmp <- data.frame(cohort_transform = unique(key), stringsAsFactors = FALSE)
  for (mo in MATCH_ONS) {
    s <- summ_all[summ_all$match_on == mo, ]
    k <- paste(s$GEO_ID, s$transform, sep = " / ")
    idx <- match(cmp$cohort_transform, k)
    cmp[[paste0("up_",   mo)]] <- s$disease_gt_normal[idx]
    cmp[[paste0("down_", mo)]] <- s$disease_lt_normal[idx]
    ## q < ALPHA の内訳だけでは match_on 間の差が閾値で潰れて見えなくなる。
    ## 床打ちしない medZ と、緩い閾値での frac05 を必ず併記する。
    cmp[[paste0("medZmad_", mo)]] <- signif(s$med_abs_zmad[idx], 3)
    cmp[[paste0("zgapN_",  mo)]] <- signif(s$zero_gap_normal[idx], 3)
    cmp[[paste0("zgapAsym_",mo)]] <- signif(s$zero_gap_asym[idx], 3)
    cmp[[paste0("frac05_",mo)]] <- signif(s$frac_p05[idx], 3)
    cmp[[paste0("b1_",   mo)]] <- signif(s$b1_mean[idx], 3)
    cmp[[paste0("mrN_",  mo)]] <- signif(s$match_ratio_normal[idx], 4)
  }
  write.csv(cmp, paste0(OUT_PREFIX, "_match_on_comparison.csv"), row.names = FALSE)
  cat("\n--- match_on 間の比較 ---\n"); print(cmp)
  cat("\n読み方（上から順に見る）:\n",
      "  zgapN_*    … normal 側のゼロ率の差（対照 - orphan）。0 に近いこと。\n",
      "                マッチングの失敗はまずここに出る。平均 TPM の比\n",
      "                (mrN_*) は裾の数本に支配されるので単独では使えない。\n",
      "  zgapAsym_* … normal と disease でゼロ率の差が違う量。ここが 0 から\n",
      "                離れていると、その非対称性がそのまま beta3 を作る。\n",
      "                beta3 を疾患応答として読めるのは、ここが約 0 のときだけ。\n",
      "  medZmad_*  … |beta3 - median(null)| / mad(null) の患者中央値。帰無分布\n",
      "                は裾が重く（z 化で sd/mad の中位 28）sd 基準の z は過小に\n",
      "                出るので、効果量はこちらで比べる。\n",
      "  frac05_*   … 生の両側 p < 0.05 の患者割合。q < ALPHA より緩い読み。\n",
      "  b1_*       … normal 側の class 主効果。scale では約 0 が望ましい。\n",
      "  up/down    … match_on 間で食い違う場合、結論はマッチング変数の選択に\n",
      "                依存している。zgapN と zgapAsym が 0 に近い列を採る。\n",
      "                既定の external2 が最も揃う（冒頭の検証表を参照）。\n",
      sep = "")
}

## ---------------------------------------------------------------------------
## 結果を読むときの注意（Methods と Response letter に反映すること）
##
##  1. beta3 をコホート内で平均してはならない。符号が患者間で割れていれば平均は
##     構造的にゼロになり、現象を検定したことにならない。報告するのは
##     direction の内訳、frac_positive、var_beta3 である。
##     summ の beta3_mean は診断目的でのみ出している。
##
##  2. p 値の下限は 1/(N_PERM + 1)。N_PERM = 1000 なら 0.000999。
##     BH 後の q がこれを下回ることはない。
##
##  3. 診断値の読み方。
##       match_ratio_normal  … 約 1 であること。ここが 1 から外れていれば
##                             発現量マッチ自体が効いていない。
##       match_ratio_disease … 1 から外れてよい。orphan に本物の疾患応答が
##                             あれば対照は disease 側で揃わなくなる。これは
##                             効果の証拠であって不具合ではない。
##       b1_mean             … normal 側の class 主効果。約 0 であること。
##     match_on=pair_mean では normal 側の比も崩れる（疾患効果を吸収するため）。
##     その場合 beta3 は信用できない。
##
##  4. shortfall_max が 0 でない患者は、いずれかの層で非 orphan の候補が
##     不足している。その患者では対照集合が orphan より小さい。
##
##  5. 既存 Table 5 と揃えた比較が要る場合は
##       MATCH_ON="external,pair_mean" Rscript R/07_interaction_permutation.R
##     と 1 回で両方走らせ、_match_on_comparison.csv を感度分析として併記する。
##     ただし冒頭の検証表の通り、pair_mean は低発現域全体が疾患側で動く状況で
##     偽陽性を出す。主たる結果は external、pair_mean は既存結果との対応を
##     示すためだけに載せる。
## ---------------------------------------------------------------------------

message("完了: ", paste0(OUT_PREFIX, "_summary.csv"),
        " (match_on = ", paste(MATCH_ONS, collapse = ", "), ")")
