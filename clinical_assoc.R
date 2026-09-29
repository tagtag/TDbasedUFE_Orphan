#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## clinical_assoc.R
##
## 患者ごとのオーファン特異的効果と臨床変数の関連を評価する。
## 旧原稿の Discussion D2（PASI との Bayes 相関）と D3（多変量回帰）の差し替え。
##
## 既存の metadata_output/GEO_patient_metadata_combined.tsv をそのまま使い、
## disease_run_accession で結合する（従来の手順と同じ結合キー）。
##
## ---------------------------------------------------------------------------
## 【なぜ差し替えが必要か】
##   旧原稿は相関・回帰の左辺に rank(P^disease_j) を置いていた。これは §3.4 の
##   「z 標準化した発現量に片側対応 t 検定をかけて得た患者ごとの p 値」であり、
##   改訂で撤回する量そのものである。撤回した量と臨床変数の相関は報告できない。
##
##   従来の手順にはもう一つ問題がある。片側検定を 2 本持っていたため、
##   コホートごとに左辺が違っていた（GSE127165 は TTEST_greater、
##   GSE40419 は TTEST_less）。どちらを使うかがデータを見た後の選択になる。
##   両側の単一量に変えることでこの自由度が消える。**これは改善なので
##   response letter に書く価値がある。**
##
## 【新しい患者水準量】
##   excess_j = d_zr_orphan_j - d_zr_control_j
##     その患者で、オーファンが発現量・検出頻度を揃えた対照よりどれだけ
##     余分にゼロ率を増やしたか。**正 = 疾患側でオーファンが余分に検出を失う。**
##
##   ただし excess は全体の検出移動 d_zr_all に依存する（それが fit_e の傾き）。
##   患者間の比較に使うべきは、その依存を除いた残差である。
##
##   e_adj_j = excess_j - (a_hat + b_hat * d_zr_all_j)      ... fit_e の残差
##
##   構成上 mean(e_adj) = 0 になるが、切片は全患者に共通の定数なので
##   患者間の比較には影響しない。
##
## 【推定の方針】
##   (1) Bayes 相関  … 旧原稿と同じ correlationBF(rscale = 1/3) + posterior(50000)。
##       **手法・事前分布・rank 変換をすべて据え置き、量だけ e_adj に差し替える。**
##       correlationBF は共変量を取れないので、ここは残差方式でなければならない。
##       rank 版を主（旧 Table 7 の直系の後継）、生の値を感度分析とする。
##   (2) 多変量線形モデル … excess ~ <臨床変数> + d_zr_all
##       従来は lm(rank(LHS) ~ ...) だったが、excess は比率の差という
##       意味のある物理量なので生のまま使う。係数が解釈できる形になる。
##       d_zr_all を共変量に入れるので残差を作る必要がない（二段階推定より
##       標準誤差が正しい）。
##   (3) Spearman … 分布の仮定を置かない確認。
##
## 【旧原稿からの変更点で、本文に書く必要があるもの】
##   - GSE40419 は除外（バッチ完全交絡）。旧原稿の第 2 式は消える。
##   - 片側検定はすべて両側に改めた（左辺の選択の自由度が消える）。
##   - 臨床係数に BH 補正をかける（旧原稿は無補正で「有意なし」と述べていた）。
##   - 影響点の診断（Cook 距離）を付ける。GSE144269 で 1 人が傾きを 63% 動かした
##     前例があるので、臨床モデルでも同じ確認をする。
##
## ---------------------------------------------------------------------------
## 【使い方】
##   slope_<COHORT>.csv のあるディレクトリ（run_verify.sh を起動した場所）で:
##
##     METADATA_DIR=Revised/metadata Rscript clinical_assoc.R
##
##     METADATA_DIR=Revised/metadata \
##     META_FILE=../metadata_output/GEO_patient_metadata_combined.tsv \
##       Rscript clinical_assoc.R
##
##     COHORTS=GSE244679,GSE127165 Rscript clinical_assoc.R
##     VARS_GSE127165=stage,alcohol,age Rscript clinical_assoc.R
##
##   環境変数
##     METADATA_DIR  read_sheet() が読むサンプルシートの場所。
##                   **run_verify.sh が全スクリプトに渡しているのと同じもの。**
##                   これを忘れると患者と run accession の対応が作れない。
##     META_FILE     臨床メタデータ。既定は
##                   metadata_output/GEO_patient_metadata_combined.tsv
##     META_KEY      結合キーの列名。既定は disease_run_accession を自動検出
##     SLOPE_DIR     slope_<COHORT>.csv の場所。既定はカレント
##     VARS_<COHORT> 使う臨床変数（カンマ区切り）。未指定なら自動検出
##
##   臨床変数を指定しなければ、コホートごとに**使える列を自動検出**して
##   一覧を出す。GSE144269 に何が入っているかもこれで分かる
##   （旧原稿は "missing biological labels" としているが未確認）。
##
## 【患者と疾患 run accession の対応】
##   サンプルシートの library_id が run accession そのものである
##   （run_verify.sh の EXCLUDE_LIBS が SRR8631680 / SRR10969201 を
##   library_id として渡していることから確認できる）。したがって
##
##     read_sheet(CO) の condition == "disease" の行の library_id
##       = メタデータの disease_run_accession
##
##   が対応になる。diag_slope.R / rho_compare.R / de_framework.R が
##   ln, ld を作っているのと同じ経路で、**再実行は一切不要**。
##
##   ただし read_sheet() は METADATA_DIR を見るので、
##   **run_verify.sh と同じように METADATA_DIR を渡すこと。**
##
##     METADATA_DIR=Revised/metadata Rscript clinical_assoc.R
##
##   探す順序:
##     1. slope_<COHORT>.csv に lib_disease / disease_run_accession 列がある
##     2. read_sheet(CO) から作る（上記。既定の経路）
##     3. map_<COHORT>.csv（patient_id, disease_run_accession）
##
## 【出力】
##   clinical_assoc_coefficients.csv   全係数（BH 補正済み）
##   clinical_assoc_posterior.csv      Bayes 相関の事後分布の要約
##   clinical_assoc_coverage.csv       コホート別の臨床変数の埋まり具合
## ---------------------------------------------------------------------------

SLOPE_DIR <- Sys.getenv("SLOPE_DIR", ".")
META_FILE <- Sys.getenv("META_FILE",
                        "metadata_output/GEO_patient_metadata_combined.tsv")
SEED      <- as.integer(Sys.getenv("SEED", "20260920"))
ITER      <- as.integer(Sys.getenv("BF_ITER", "50000"))
RSCALE    <- as.numeric(Sys.getenv("BF_RSCALE", "0.3333333333333333"))
MIN_COV   <- as.numeric(Sys.getenv("MIN_COVERAGE", "0.8"))

COHORTS <- trimws(strsplit(
  Sys.getenv("COHORTS", "GSE244679,GSE127165,GSE144269"), ",")[[1]])
COHORTS <- COHORTS[nzchar(COHORTS)]

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")

## 臨床変数として使わない列。識別子・管理情報・配列統計。
## 自動検出のときだけ効く。VARS_<COHORT> で名指しすれば通る。
DROP_RE <- paste0(
  "(accession|run$|_run|^run|sample|patient|subject|gsm|gse|srr|srx|srs|",
  "srp|biosample|library|^id$|_id$|title|file|path|fastq|url|link|",
  "condition|tissue$|group$|batch|lane|read|bases|spots|bytes|md5|",
  "date|submit|platform|instrument|layout|strategy|source)")

## --------------------------------------------------------------- 読み込み
if (!dir.exists(SLOPE_DIR))
  stop("SLOPE_DIR がありません: ", SLOPE_DIR, call. = FALSE)

have_bf <- requireNamespace("BayesFactor", quietly = TRUE)
if (!have_bf)
  msg("注意: BayesFactor が無いので Bayes 相関は飛ばします。",
      "install.packages(\"BayesFactor\")\n")

if (!file.exists(META_FILE))
  stop("臨床メタデータがありません: ", META_FILE, "\n",
       "  META_FILE=... で場所を指定してください。", call. = FALSE)

sepc <- if (grepl("\\.tsv$|\\.txt$", META_FILE)) "\t" else ","
meta <- read.csv(META_FILE, sep = sepc, stringsAsFactors = FALSE,
                 check.names = FALSE)

## 結合キーの列名を決める
KEY <- Sys.getenv("META_KEY", "")
if (!nzchar(KEY)) {
  cand <- c("disease_run_accession", "run_accession", "disease_run", "run", "srr")
  KEY  <- cand[cand %in% names(meta)][1]
}
if (is.na(KEY) || !nzchar(KEY) || !KEY %in% names(meta))
  stop("結合キーの列が見つかりません。META_KEY=<列名> で指定してください。\n",
       "  ", basename(META_FILE), " の列: ",
       paste(utils::head(names(meta), 40), collapse = ", "), call. = FALSE)

msg("臨床メタデータ: ", META_FILE)
msg("  ", nrow(meta), " 行 / 結合キー: ", KEY)
if (!nzchar(Sys.getenv("METADATA_DIR")))
  msg("  METADATA_DIR は未設定。sample_sheet_<COHORT>.csv は自分で探します。\n",
      "  （run_verify.sh と同じ値を渡せば探索は省けます）")

## -------------------------------------------- patient_id -> 疾患 run accession
## diag_slope.R / rho_compare.R が ln, ld を作っているのと同じ経路。再実行は不要。
## サンプルシートの library_id が run accession そのものである。
##
## 失敗したら **理由をそのまま返す**。握り潰すと「対応が作れません」としか
## 出ず、METADATA_DIR の渡し忘れなのか R/ が無いのかが分からなくなる。
emsg <- function(e) {
  cc <- attr(e, "condition")
  if (is.null(cc)) paste(as.character(e), collapse = " ") else conditionMessage(cc)
}

sheet_map <- function(CO) {
  if (!exists("read_sheet")) {
    for (f in c("R/config.R", "R/00_functions.R")) {
      if (!file.exists(f))
        return(list(map = NULL, err = paste0(f, " がありません（カレントが",
                                             "リポジトリのルートではない）")))
      e <- try(source(f), silent = TRUE)
      if (inherits(e, "try-error"))
        return(list(map = NULL, err = paste0(f, " の読み込みに失敗: ", emsg(e))))
    }
  }
  if (!exists("read_sheet"))
    return(list(map = NULL, err = "R/00_functions.R に read_sheet() がありません"))

  sh   <- try(read_sheet(CO), silent = TRUE)
  from <- "read_sheet()"

  ## read_sheet() は <METADATA_DIR>/sample_sheet_<COHORT>.csv を見る。
  ## METADATA_DIR の渡し忘れで落ちることが多いので、自分で探して直接読む。
  ## サンプルシートは patient_id / library_id / condition を持つ単なる CSV なので
  ## read_sheet() を通さなくても同じものが得られる。
  if (inherits(sh, "try-error")) {
    e1   <- emsg(sh)
    base <- sprintf("sample_sheet_%s.csv", CO)
    dirs <- c(Sys.getenv("METADATA_DIR"),
              "metadata", "Revised/metadata", "metadata_output",
              "../metadata", "../Revised/metadata", "Revised", ".", "..")
    dirs <- dirs[nzchar(dirs)]
    cand <- file.path(dirs, base)
    hit  <- cand[file.exists(cand)][1]
    if (is.na(hit)) {
      ## 決め打ちで見つからなければ、深さを限って探す（大きなリポジトリでも軽い）
      found <- list.files(".", pattern = paste0("^", base, "$"),
                          recursive = TRUE, full.names = TRUE)
      found <- found[lengths(strsplit(found, "/")) <= 5L]
      hit <- found[1]
    }
    if (is.na(hit) || !length(hit))
      return(list(map = NULL,
                  err = paste0("read_sheet(\"", CO, "\") が失敗: ", e1,
                               "\n        ", base, " も見つかりませんでした")))
    sh <- try(read.csv(hit, stringsAsFactors = FALSE), silent = TRUE)
    if (inherits(sh, "try-error"))
      return(list(map = NULL,
                  err = paste0(hit, " の読み込みに失敗: ", emsg(sh))))
    from <- paste0("自動で見つけた ", hit)
  }

  if (!all(c("patient_id", "library_id", "condition") %in% names(sh)))
    return(list(map = NULL,
                err = paste0("サンプルシートに必要な列がありません。列: ",
                             paste(names(sh), collapse = ", "))))
  dd <- sh[sh$condition == "disease", c("patient_id", "library_id")]
  if (!nrow(dd))
    return(list(map = NULL,
                err = paste0("condition == \"disease\" の行がありません。",
                             "値: ", paste(unique(sh$condition), collapse = ", "))))
  names(dd)[2] <- "lib_disease"
  list(map = dd[!duplicated(dd$patient_id), , drop = FALSE], err = NA_character_,
       from = from)
}

coef_rows <- list(); post_rows <- list(); cov_rows <- list()
auto_used <- character(0)

## 旧原稿が実際に使っていた変数。これを「事前指定」の集合とみなす。
PRESPEC <- c("PASI", "pasi", "stage", "alcohol", "age", "smoking")

## ------------------------------------------------------- コホートごとの処理
for (CO in COHORTS) {

  f <- file.path(SLOPE_DIR, sprintf("slope_%s.csv", CO))
  if (!file.exists(f)) {
    hr(); msg(CO, " : ", basename(f), " がありません。飛ばします。"); next
  }
  d <- read.csv(f, stringsAsFactors = FALSE)
  need <- c("patient_id", "d_zr_all", "d_zr_orphan", "d_zr_control")
  if (!all(need %in% names(d))) {
    hr(); msg(CO, " : 列が足りません（",
              paste(setdiff(need, names(d)), collapse = ", "), "）。飛ばします。"); next
  }
  d$excess <- d$d_zr_orphan - d$d_zr_control
  fit_e    <- lm(excess ~ d_zr_all, data = d)
  d$e_adj  <- residuals(fit_e)

  hr(); msg(CO, " / slope に患者 ", nrow(d), " 人")

  ## --- 疾患 run accession の列を用意する
  src <- NA_character_; sheet_err <- NA_character_
  if ("lib_disease" %in% names(d)) {
    src <- "slope_*.csv の lib_disease 列"
  } else if ("disease_run_accession" %in% names(d)) {
    d$lib_disease <- d$disease_run_accession
    src <- "slope_*.csv の disease_run_accession 列"
  } else {
    sm <- sheet_map(CO)
    if (!is.null(sm$map)) {
      d <- merge(d, sm$map, by = "patient_id", all.x = TRUE)
      src <- paste0(sm$from, "（library_id = run accession）")
    } else {
      sheet_err <- sm$err
      mf <- file.path(SLOPE_DIR, sprintf("map_%s.csv", CO))
      if (file.exists(mf)) {
        mm <- read.csv(mf, stringsAsFactors = FALSE)
        nm <- intersect(c("lib_disease", "disease_run_accession"), names(mm))[1]
        if (!is.na(nm)) {
          mm$lib_disease <- mm[[nm]]
          d <- merge(d, mm[, c("patient_id", "lib_disease")], by = "patient_id",
                     all.x = TRUE)
          src <- basename(mf)
        }
      }
    }
  }
  if (!"lib_disease" %in% names(d)) {
    msg("  patient_id と疾患 run accession の対応が作れません。")
    if (exists("sheet_err") && !is.na(sheet_err)) msg("  理由: ", sheet_err)
    msg("")
    msg("  対応はサンプルシートから作れます（library_id が run accession）。")
    msg("  run_verify.sh と同じように METADATA_DIR を渡し、リポジトリの")
    msg("  ルート（R/config.R のある場所）で実行してください:")
    msg("")
    msg("    METADATA_DIR=Revised/metadata Rscript clinical_assoc.R")
    msg("")
    msg("  それでも駄目なら map_", CO, ".csv を置いてください。2 列だけです:")
    msg("    patient_id,disease_run_accession")
    msg("  次の 1 行で作れます:")
    msg("    METADATA_DIR=Revised/metadata Rscript -e '",
        "source(\"R/config.R\");source(\"R/00_functions.R\");",
        "s<-read_sheet(\"", CO, "\");s<-s[s$condition==\"disease\",];",
        "write.csv(data.frame(patient_id=s$patient_id,",
        "disease_run_accession=s$library_id),\"map_", CO, ".csv\",row.names=FALSE)'")
    next
  }
  msg("  対応の出所: ", src)

  ## --- 結合
  d$.key <- sub("_kallisto$", "", as.character(d$lib_disease))
  mk     <- sub("_kallisto$", "", as.character(meta[[KEY]]))
  hit    <- match(d$.key, mk)
  msg(sprintf("  メタデータに一致 %d / %d 人", sum(!is.na(hit)), nrow(d)))
  if (sum(!is.na(hit)) < 8L) {
    msg("  一致が少なすぎます。結合キーを確認してください。")
    msg("    slope 側の例: ", paste(utils::head(d$.key, 3), collapse = ", "))
    msg("    meta  側の例: ", paste(utils::head(mk, 3), collapse = ", "))
    next
  }
  dm <- d[!is.na(hit), , drop = FALSE]
  mm <- meta[hit[!is.na(hit)], , drop = FALSE]

  ## --- 使える臨床変数を決める
  vars <- trimws(strsplit(Sys.getenv(paste0("VARS_", CO), ""), ",")[[1]])
  vars <- vars[nzchar(vars)]
  auto <- !length(vars)

  cand <- setdiff(names(mm), KEY)
  rep_rows <- list()
  usable <- character(0)
  for (v in cand) {
    x  <- mm[[v]]
    nn <- sum(!is.na(x) & trimws(as.character(x)) != "")
    cv <- nn / nrow(mm)
    nu <- length(unique(x[!is.na(x) & trimws(as.character(x)) != ""]))
    idlike <- grepl(DROP_RE, v, ignore.case = TRUE)
    why <- if (idlike) "識別子・管理情報とみなして除外"
           else if (cv < MIN_COV) sprintf("埋まり %.0f%% < %.0f%%", 100*cv, 100*MIN_COV)
           else if (nu < 2L) "値が 1 種類しかない"
           else if (nu > nrow(mm) * 0.9 && !is.numeric(x)) "ほぼ全行で異なる（識別子か）"
           else ""
    if (!nzchar(why)) usable <- c(usable, v)
    rep_rows[[length(rep_rows)+1L]] <- data.frame(
      cohort = CO, column = v, n_nonmissing = nn, coverage = cv,
      n_distinct = nu, usable = !nzchar(why),
      note = if (nzchar(why)) why else "使える", stringsAsFactors = FALSE)
  }
  cov_rows <- c(cov_rows, rep_rows)

  msg("")
  msg("  臨床変数の埋まり具合")
  for (r in rep_rows) {
    if (!r$usable && r$note == "識別子・管理情報とみなして除外") next
    msg(sprintf("    %-22s 非欠測 %3d/%3d (%3.0f%%)  値 %3d 種  %s",
                r$column, r$n_nonmissing, nrow(mm), 100*r$coverage,
                r$n_distinct, if (r$usable) "**使える**" else r$note))
  }

  if (auto) {
    vars <- usable
    auto_used <- c(auto_used, CO)
    msg("  → 自動検出: ", if (length(vars)) paste(vars, collapse = ", ") else "なし")
    extra <- setdiff(vars, PRESPEC)
    if (length(extra)) {
      msg("    うち事前指定の外: ", paste(extra, collapse = ", "))
      msg("    これらは探索的。確認的解析は VARS_", CO, " で旧原稿の変数に")
      msg("    絞ること（下の注意を参照）。")
    }
  } else {
    miss <- setdiff(vars, names(mm))
    if (length(miss)) msg("  指定された列が無い: ", paste(miss, collapse = ", "))
    vars <- intersect(vars, names(mm))
  }
  if (!length(vars)) {
    msg("  使える臨床変数がありません。飛ばします。"); next
  }

  for (v in vars) dm[[v]] <- mm[[v]]

  ## 欠測を落とす（モデル間で n を揃えるため一括）
  blank <- function(x) is.na(x) | trimws(as.character(x)) == ""
  cc <- !Reduce(`|`, lapply(vars, function(v) blank(dm[[v]])))
  if (any(!cc)) { msg("  欠測により ", sum(!cc), " 人を除外"); dm <- dm[cc, , drop = FALSE] }
  if (nrow(dm) < 8L) { msg("  n = ", nrow(dm), " は少なすぎます。飛ばします。"); next }

  ## 型を決めて報告
  msg("")
  for (v in vars) {
    x <- dm[[v]]
    if (!is.numeric(x)) {
      xn <- suppressWarnings(as.numeric(as.character(x)))
      dm[[v]] <- if (!any(is.na(xn))) xn else factor(as.character(x))
    }
    msg(sprintf("    %-18s %s", v,
        if (is.factor(dm[[v]]))
          paste0("因子 [", paste(levels(dm[[v]]), collapse = " / "), "]")
        else sprintf("数値 中位 %.4g 範囲 %.4g 〜 %.4g",
                     median(dm[[v]]), min(dm[[v]]), max(dm[[v]]))))
  }
  msg(sprintf("    %-18s 中位 %+.5f 範囲 %+.5f 〜 %+.5f", "excess",
              median(dm$excess), min(dm$excess), max(dm$excess)))
  msg("    符号: excess > 0 = 疾患側でオーファンが対照より余分に検出を失う")

  ## -------------------------------------------- (1) 多変量線形モデル
  fml <- stats::as.formula(paste("excess ~",
                                 paste(c(vars, "d_zr_all"), collapse = " + ")))
  fit <- lm(fml, data = dm)
  s   <- summary(fit)$coefficients
  ci  <- suppressWarnings(stats::confint(fit))

  msg(""); msg("  多変量線形モデル  ", deparse(fml))
  msg(sprintf("    n = %d   残差自由度 %d   R2 = %.3f   調整 R2 = %.3f",
              nrow(dm), fit$df.residual, summary(fit)$r.squared,
              summary(fit)$adj.r.squared))
  for (nm in rownames(s)) {
    if (nm == "(Intercept)") next
    msg(sprintf("    %-22s %+.5f ± %.5f   95%%CI %+.5f 〜 %+.5f   p = %.3g%s",
                nm, s[nm,1], s[nm,2], ci[nm,1], ci[nm,2], s[nm,4],
                if (nm == "d_zr_all") "  [共変量]" else ""))
    if (nm != "d_zr_all")
      coef_rows[[length(coef_rows)+1L]] <- data.frame(
        cohort = CO, model = "lm(excess ~ vars + d_zr_all)", term = nm,
        n = nrow(dm), estimate = s[nm,1], se = s[nm,2],
        ci_lo = ci[nm,1], ci_hi = ci[nm,2], p_raw = s[nm,4],
        stringsAsFactors = FALSE)
  }

  ## 影響点
  cd <- stats::cooks.distance(fit); thr <- 4/nrow(dm); bad <- which(cd > thr)
  msg(""); msg(sprintf("    Cook 距離 最大 %.3f（患者 %s、閾値 4/n = %.3f）",
                       max(cd), dm$patient_id[which.max(cd)], thr))
  if (length(bad)) {
    msg(sprintf("    閾値超え %d 人: %s", length(bad),
                paste(dm$patient_id[bad], collapse = ", ")))
    s2 <- summary(lm(fml, data = dm[-bad, , drop = FALSE]))$coefficients
    msg("    これらを除いたとき:")
    for (nm in rownames(s2)) {
      if (nm %in% c("(Intercept)", "d_zr_all")) next
      sh <- if (!nm %in% rownames(s)) ""
            else if (abs(s[nm,1]) > s[nm,2])
              sprintf("  （%+.0f%%）", 100*(s2[nm,1]-s[nm,1])/abs(s[nm,1]))
            else "  （元の係数が 0 と区別できないので変化率は出しません）"
      msg(sprintf("      %-20s %+.5f ± %.5f   p = %.3g%s",
                  nm, s2[nm,1], s2[nm,2], s2[nm,4], sh))
    }
    msg("    係数が大きく動く場合は本文に影響点の存在を明記すること。")
  } else msg("    閾値を超える影響点なし。")

  ## -------------------------------------------- (2) Spearman
  msg(""); msg("  Spearman 順位相関（両側）")
  for (v in vars) {
    if (is.factor(dm[[v]])) { msg(sprintf("    %-18s 因子なので取りません", v)); next }
    for (q in c("excess", "e_adj")) {
      ct <- suppressWarnings(stats::cor.test(dm[[q]], dm[[v]], method = "spearman"))
      msg(sprintf("    %-18s vs %-7s rho = %+.3f   p = %.3g",
                  v, q, unname(ct$estimate), ct$p.value))
      coef_rows[[length(coef_rows)+1L]] <- data.frame(
        cohort = CO, model = paste0("spearman(", q, ")"), term = v,
        n = nrow(dm), estimate = unname(ct$estimate), se = NA_real_,
        ci_lo = NA_real_, ci_hi = NA_real_, p_raw = ct$p.value,
        stringsAsFactors = FALSE)
    }
  }

  ## -------------------------------------------- (3) Bayes 相関（旧 Table 7）
  ## 旧: correlationBF(rank(TTEST1_less), rank(y$PASI[index]), rscale = 1/3)
  ## 新: 量だけ e_adj に差し替える。rank 変換・rscale・反復数は据え置き。
  if (have_bf) for (v in vars) {
    if (is.factor(dm[[v]])) next
    for (mode in c("rank", "raw")) {
      xx <- if (mode == "rank") rank(dm[[v]])  else dm[[v]]
      yy <- if (mode == "rank") rank(dm$e_adj) else dm$e_adj
      set.seed(SEED)
      bf <- try(BayesFactor::correlationBF(y = yy, x = xx, rscale = RSCALE),
                silent = TRUE)
      if (inherits(bf, "try-error")) {
        msg("    ", v, " (", mode, ") : correlationBF が失敗しました"); next
      }
      po <- NULL
      invisible(utils::capture.output(suppressWarnings(suppressMessages({
        po <- try(BayesFactor::posterior(bf, iterations = ITER, progress = FALSE),
                  silent = TRUE)
        if (inherits(po, "try-error"))
          po <- try(BayesFactor::posterior(bf, iterations = ITER), silent = TRUE)
      }))))
      bf10 <- as.numeric(BayesFactor::extractBF(bf)$bf)
      if (is.null(po) || inherits(po, "try-error")) {
        msg(sprintf("    Bayes 相関 e_adj vs %s (%s): posterior 失敗。BF10 = %.3f",
                    v, mode, bf10)); next
      }
      cn  <- colnames(po)
      rho <- as.numeric(po[, if ("rho" %in% cn) "rho" else cn[1]])
      qs  <- stats::quantile(rho, c(0.025, 0.5, 0.975))
      msg("")
      msg(sprintf("  Bayes 相関  e_adj vs %s  [%s]  (rscale = %.4f, %d 反復)%s",
                  v, mode, RSCALE, ITER,
                  if (mode == "rank") "   ← 旧 Table 7 と同じ変換" else "   ← 感度分析"))
      msg(sprintf("    事後中央値 rho = %+.4f   95%% 確信区間 %+.4f 〜 %+.4f",
                  qs[2], qs[1], qs[3]))
      msg(sprintf("    BF10 = %.3f   P(rho > 0) = %.3f", bf10, mean(rho > 0)))
      post_rows[[length(post_rows)+1L]] <- data.frame(
        cohort = CO, variable = v, quantity = "e_adj", transform = mode,
        n = nrow(dm), rscale = RSCALE, q2.5 = qs[1], median = qs[2],
        q97.5 = qs[3], bf10 = bf10, p_rho_gt0 = mean(rho > 0),
        stringsAsFactors = FALSE)
    }
  }
}

## ------------------------------------------------------------------ まとめ
if (length(cov_rows))
  write.csv(do.call(rbind, cov_rows), "clinical_assoc_coverage.csv",
            row.names = FALSE)

if (!length(coef_rows)) {
  hr(); msg("係数が 1 つも得られませんでした。")
  if (length(cov_rows)) msg("clinical_assoc_coverage.csv に列の埋まり具合を出しました。")
  quit(save = "no", status = 1L)
}

cf  <- do.call(rbind, coef_rows)
## BH は線形モデルの係数のみを族とする。
## Spearman と Bayes 相関は同じ仮説の感度分析なので二重に数えない。
isl <- cf$model == "lm(excess ~ vars + d_zr_all)"
cf$p_bh <- NA_real_
cf$p_bh[isl] <- stats::p.adjust(cf$p_raw[isl], method = "BH")

hr(); msg("臨床係数のまとめ（BH は線形モデルの係数 ", sum(isl), " 個を族とする）"); msg("")
msg(sprintf("%-11s %-30s %-16s %4s %10s %10s %10s",
            "cohort", "model", "term", "n", "estimate", "p", "p(BH)"))
for (i in seq_len(nrow(cf)))
  msg(sprintf("%-11s %-30s %-16s %4d %+10.5f %10.3g %10s",
              cf$cohort[i], cf$model[i], cf$term[i], cf$n[i], cf$estimate[i],
              cf$p_raw[i], if (is.na(cf$p_bh[i])) "-" else sprintf("%.3g", cf$p_bh[i])))

write.csv(cf, "clinical_assoc_coefficients.csv", row.names = FALSE)
msg(""); msg("書き出し: clinical_assoc_coefficients.csv")
if (length(post_rows)) {
  write.csv(do.call(rbind, post_rows), "clinical_assoc_posterior.csv",
            row.names = FALSE)
  msg("書き出し: clinical_assoc_posterior.csv")
}
if (length(cov_rows)) msg("書き出し: clinical_assoc_coverage.csv")

hr(); msg("読み方"); msg("")
msg("  左辺は excess = d_zr_orphan - d_zr_control（ゼロ率。表現に依存しない）。")
msg("  **正 = 疾患側でオーファンが発現量・検出頻度を揃えた対照より余分に")
msg("  検出を失う。** 旧原稿の rank(P^disease_j) とは別の量なので、旧 Table 7")
msg("  の 0.309 とは比較できない。「0.309 が X になった」と書かないこと。")
msg("")
msg("  d_zr_all は共変量として必ず入れてある。入れないと、全体の検出移動")
msg("  （ライブラリの深さ・腫瘍純度などと相関しうる）が臨床変数との見かけの")
msg("  関連を作る。これを外した結果を主結果にしてはならない。")
msg("")
msg("  Bayes 相関は rank 版を主とする。旧 Table 7 と同じ rank 変換・同じ")
msg("  rscale = 1/3・同じ 50000 反復なので、**変わったのは量だけ**と言える。")
msg("  生の値の版は感度分析として併記する。")
msg("")
msg("  BH 後に有意でなくても結果として十分である。旧原稿も L544 で")
msg("  「No significant associations were detected」と報告しており、")
msg("  撤回後の慎重な論調にはむしろ合う。null の場合の書き方は")
msg("  「臨床変数では説明できない患者間異質性が残る」であって、")
msg("  「臨床的に無意味」ではない。検出力の限界を必ず併記すること。")
msg("")
msg("  Cook 距離が閾値を超えた患者がいた場合、除外前後の両方を本文に書く。")
msg("  GSE144269 で 1 人が傾きを 63% 動かした前例がある。")

if (length(auto_used)) {
  msg("")
  hr()
  msg("** 自動検出を使ったコホート: ", paste(auto_used, collapse = ", "), " **")
  msg("")
  msg("  自動検出は「メタデータに何が入っているか」を見るための機能である。")
  msg("  見つかった列をすべて回帰に入れて BH をかけると、検定の族が")
  msg("  データを見た後に決まったことになり、査読で必ず突かれる。")
  msg("")
  msg("  **論文に載せる確認的解析では VARS_<COHORT> で変数を明示すること。**")
  msg("  旧原稿が使っていた事前指定の集合は次のとおり:")
  msg("    GSE244679 … PASI")
  msg("    GSE127165 … stage, alcohol（旧原稿では age はコメントアウト済み）")
  msg("    GSE40419  … stage, age, smoking（**このコホートは除外するので消える**）")
  msg("")
  msg("    VARS_GSE244679=PASI VARS_GSE127165=stage,alcohol Rscript clinical_assoc.R")
  msg("")
  msg("  それ以外の列（sex など）を使う場合は、本文で探索的と明記し、")
  msg("  確認的解析の BH の族には入れないこと。")
}
