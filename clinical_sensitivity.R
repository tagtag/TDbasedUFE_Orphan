#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## clinical_sensitivity.R
##
## clinical_assoc.R が出した臨床変数の関連を、論文に載せる前に潰す。
## 現時点の対象は GSE127165 の stage（生 p = 0.024、BH = 0.072、BF10 = 2.7）。
##
## **この所見を本文に書くかどうかは、下の (F) 技術交絡の結果で決まる。**
## stage が検出深度やライブラリの性質と相関しているなら、
## 「進行期ほどオーファンの検出損失が小さい」は生物ではなく技術の話になる。
##
## 走らせる確認:
##   (A) 基準モデルの再掲
##   (B) stage を因子として扱う（線形トレンドの仮定を外す）
##   (C) 利用可能な共変量を 1 つずつ追加（smoking, age など）
##   (D) d_zr_all を外したらどうなるか（共変量の効き方を見る）
##   (E) 1 人ずつ抜く（leave-one-out）
##   (F) **技術交絡** — stage と検出深度・転写産物数・ゼロ率の相関
##   (G) 並べ替え検定（分布の仮定を置かない p）
##
## 使い方（clinical_assoc.R と同じ場所で）
##   METADATA_DIR=Revised/metadata COHORT=GSE127165 FOCUS=stage \
##     Rscript clinical_sensitivity.R
##
##   COHORT   既定 GSE127165
##   FOCUS    既定 stage        注目する臨床変数
##   BASE     既定 stage,alcohol  基準モデルに入れる変数（カンマ区切り）
##   NPERM    既定 10000       並べ替え回数
## ---------------------------------------------------------------------------

CO        <- Sys.getenv("COHORT", "GSE127165")
FOCUS     <- Sys.getenv("FOCUS", "stage")
BASE      <- trimws(strsplit(Sys.getenv("BASE", "stage,alcohol"), ",")[[1]])
BASE      <- BASE[nzchar(BASE)]
NPERM     <- as.integer(Sys.getenv("NPERM", "10000"))
SEED      <- as.integer(Sys.getenv("SEED", "20260920"))
SLOPE_DIR <- Sys.getenv("SLOPE_DIR", ".")
META_FILE <- Sys.getenv("META_FILE",
                        "metadata_output/GEO_patient_metadata_combined.tsv")

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 74), "\n")
emsg <- function(e) { cc <- attr(e, "condition")
  if (is.null(cc)) paste(as.character(e), collapse = " ") else conditionMessage(cc) }

## ------------------------------------------------------------------ 読み込み
f <- file.path(SLOPE_DIR, sprintf("slope_%s.csv", CO))
if (!file.exists(f)) stop(f, " がありません。", call. = FALSE)
d <- read.csv(f, stringsAsFactors = FALSE)
d$excess <- d$d_zr_orphan - d$d_zr_control
d$e_adj  <- residuals(lm(excess ~ d_zr_all, data = d))

sepc <- if (grepl("\\.tsv$|\\.txt$", META_FILE)) "\t" else ","
meta <- read.csv(META_FILE, sep = sepc, stringsAsFactors = FALSE, check.names = FALSE)
KEY  <- Sys.getenv("META_KEY", "disease_run_accession")
stopifnot(KEY %in% names(meta))

## patient_id -> 疾患 run accession（clinical_assoc.R と同じ経路）
sh <- NULL
if (!exists("read_sheet")) {
  for (p in c("R/config.R", "R/00_functions.R"))
    if (file.exists(p)) try(source(p), silent = TRUE)
}
if (exists("read_sheet")) { s <- try(read_sheet(CO), silent = TRUE)
  if (!inherits(s, "try-error")) sh <- s }
if (is.null(sh)) {
  base <- sprintf("sample_sheet_%s.csv", CO)
  dirs <- c(Sys.getenv("METADATA_DIR"), "metadata", "Revised/metadata",
            "metadata_output", "../metadata", "../Revised/metadata", ".", "..")
  dirs <- dirs[nzchar(dirs)]
  hit  <- file.path(dirs, base); hit <- hit[file.exists(hit)][1]
  if (is.na(hit)) stop("サンプルシートが見つかりません: ", base, call. = FALSE)
  sh <- read.csv(hit, stringsAsFactors = FALSE)
}
dd <- sh[sh$condition == "disease", c("patient_id", "library_id")]
dd <- dd[!duplicated(dd$patient_id), ]
names(dd)[2] <- "lib_disease"
d  <- merge(d, dd, by = "patient_id", all.x = TRUE)

hit <- match(sub("_kallisto$", "", d$lib_disease),
             sub("_kallisto$", "", meta[[KEY]]))
d   <- d[!is.na(hit), , drop = FALSE]
mm  <- meta[hit[!is.na(hit)], , drop = FALSE]

## 使える臨床列（欠測ゼロでないもの）を拾っておく
## 識別子・管理情報は共変量にしない。clinical_assoc.R と同じ除外パターン。
## （これを写し忘れていたため patient_id / normal_GSM / disease_GSM が
##   共変量に紛れ込み、完全共線で NaN を出していた。）
DROP_RE <- paste0(
  "(accession|run$|_run|^run|sample|patient|subject|gsm|gse|srr|srx|srs|",
  "srp|biosample|library|^id$|_id$|title|file|path|fastq|url|link|",
  "condition|tissue$|group$|batch|lane|read|bases|spots|bytes|md5|",
  "date|submit|platform|instrument|layout|strategy|source|",
  "paired_complete|disease_name|orphan_delta)")

blank  <- function(x) is.na(x) | trimws(as.character(x)) == ""
avail  <- setdiff(names(mm), KEY)
dropped <- avail[grepl(DROP_RE, avail, ignore.case = TRUE)]
avail  <- avail[vapply(avail, function(v) {
  x <- mm[[v]]
  mean(!blank(x)) >= 0.8 && length(unique(x[!blank(x)])) >= 2L &&
    !grepl(DROP_RE, v, ignore.case = TRUE)
}, logical(1))]
for (v in avail) d[[v]] <- mm[[v]]

num_or_fac <- function(x) {
  if (is.numeric(x)) return(x)
  xn <- suppressWarnings(as.numeric(as.character(x)))
  if (!any(is.na(xn))) xn else factor(as.character(x))
}
for (v in avail) d[[v]] <- num_or_fac(d[[v]])

BASE  <- intersect(BASE, avail)
if (!FOCUS %in% avail)
  stop(FOCUS, " がメタデータにありません。使えるのは: ",
       paste(avail, collapse = ", "), call. = FALSE)
keep <- !Reduce(`|`, lapply(unique(c(BASE, FOCUS)), function(v) blank(d[[v]])))
d    <- d[keep, , drop = FALSE]

hr(); msg(CO, " / 注目する変数: ", FOCUS, " / n = ", nrow(d))
msg("  メタデータで使える列: ", paste(avail, collapse = ", "))
if (length(dropped))
  msg("  識別子として除外した列: ", paste(dropped, collapse = ", "))
msg("  基準モデルの変数: ", paste(BASE, collapse = ", "), " + d_zr_all")

## ------------------------------------------------------------------- 道具
co <- function(fit, term) {
  s <- summary(fit)$coefficients
  r <- grep(paste0("^", term), rownames(s))
  if (!length(r)) return(c(NA, NA, NA))
  s[r[1], c(1, 2, 4)]
}
fml <- function(vs) stats::as.formula(paste("excess ~",
                                            paste(c(vs, "d_zr_all"), collapse = " + ")))

## ------------------------------------------------------- (A) 基準モデル
fit0 <- lm(fml(BASE), data = d)
c0   <- co(fit0, FOCUS)
hr(); msg("(A) 基準モデル  ", deparse(fml(BASE)))
msg(sprintf("    %s  %+.5f ± %.5f   p = %.4g   （n = %d, R2 = %.3f）",
            FOCUS, c0[1], c0[2], c0[3], nrow(d), summary(fit0)$r.squared))

## -------------------------------------------- (B) 因子として扱う
if (is.numeric(d[[FOCUS]]) && length(unique(d[[FOCUS]])) <= 8L) {
  d$.fac <- factor(d[[FOCUS]])
  ff <- lm(fml(c(setdiff(BASE, FOCUS), ".fac")), data = d)
  an <- anova(fit0, ff)
  hr(); msg("(B) ", FOCUS, " を因子として扱う（線形トレンドの仮定を外す）")
  tb <- tapply(d$e_adj, d$.fac, function(x) c(n = length(x), mean = mean(x)))
  for (lv in names(tb))
    msg(sprintf("    %s = %-4s  n = %2d   e_adj 平均 %+.5f",
                FOCUS, lv, tb[[lv]]["n"], tb[[lv]]["mean"]))
  msg(sprintf("    線形 vs 因子の F 検定  p = %.4g", an$`Pr(>F)`[2]))
  msg("    p が大きければ線形トレンドで十分。小さければ単調でない。")
  msg("    水準ごとの平均が単調に減っているかを目で確認すること。")
  d$.fac <- NULL
}

## -------------------------------------------- (C) 共変量を 1 つずつ足す
## 臨床の列だけでなく、**技術変数（検出の広さ）も足してみる**。
## (F) で相関が出たときに効いてくるのはこちら。
add <- c(setdiff(avail, unique(c(BASE, FOCUS))),
         intersect(c("n_kept", "n_orphan"), names(d)))
hr(); msg("(C) 共変量を 1 つずつ追加（", FOCUS, " の係数が生き残るか）")
msg("    n_kept / n_orphan は技術変数。(F) で相関が出たらこの行が効く。")
if (!length(add)) msg("    追加できる列がありません。")
for (v in add) {
  if (all(blank(d[[v]]))) next
  ok <- try(lm(fml(c(BASE, v)), data = d), silent = TRUE)
  if (inherits(ok, "try-error")) { msg("    + ", v, " : 失敗 (", emsg(ok), ")"); next }
  ## 完全共線で落とされた共変量は「足しても何も変わらない」ように見えるので
  ## 必ず名指しする。黙って通すと基準比 +0% を頑健性と読み違える。
  ali <- names(which(is.na(coef(ok))))
  if (length(ali)) {
    msg(sprintf("    + %-14s ** %s が共線のため lm に落とされました（この行は無効）**",
                v, paste(ali, collapse = ", ")))
    msg(sprintf("      %s と既存の変数の対応表を確認すること: table(%s, %s)",
                v, v, paste(BASE, collapse = ", ")))
    next
  }
  cc <- co(ok, FOCUS)
  msg(sprintf("    + %-14s %s %+.5f ± %.5f   p = %.4g   （基準比 %+.0f%%）",
              v, FOCUS, cc[1], cc[2], cc[3], 100 * (cc[1] - c0[1]) / abs(c0[1])))
}

## -------------------------------------------- (D) d_zr_all を外す
fitD <- lm(stats::as.formula(paste("excess ~", paste(BASE, collapse = " + "))), data = d)
cD   <- co(fitD, FOCUS)
hr(); msg("(D) d_zr_all を共変量から外した場合")
msg(sprintf("    %s  %+.5f ± %.5f   p = %.4g   （基準比 %+.0f%%）",
            FOCUS, cD[1], cD[2], cD[3], 100 * (cD[1] - c0[1]) / abs(c0[1])))
msg("    大きく動くなら、全体の検出移動が関連の一部を担っている。")
msg("    **主結果は d_zr_all を入れた (A) の方である。**")

## -------------------------------------------- (E) leave-one-out
hr(); msg("(E) 1 人ずつ抜く（leave-one-out、", nrow(d), " 通り）")
loo <- t(vapply(seq_len(nrow(d)), function(i) co(lm(fml(BASE), data = d[-i, ]), FOCUS),
                numeric(3)))
w   <- which.max(abs(loo[, 1] - c0[1]))
msg(sprintf("    係数の範囲  %+.5f 〜 %+.5f   （基準 %+.5f）",
            min(loo[,1]), max(loo[,1]), c0[1]))
msg(sprintf("    p    の範囲  %.4g 〜 %.4g", min(loo[,3]), max(loo[,3])))
msg(sprintf("    最も動かす患者  %s   係数 %+.5f（%+.0f%%）  p = %.4g",
            d$patient_id[w], loo[w,1], 100*(loo[w,1]-c0[1])/abs(c0[1]), loo[w,3]))
msg(sprintf("    p > 0.05 になる抜き方  %d / %d", sum(loo[,3] > 0.05), nrow(d)))
if (max(loo[,3]) > 0.05) {
  msg("    ** 1 人抜くと有意でなくなる。本文に必ず書くこと。**")
} else {
  msg("    どの 1 人を抜いても p < 0.05 は保たれる。")
}

## -------------------------------------------- (F) 技術交絡  ★ここが要点
hr(); msg("(F) 技術交絡の確認  ★ この所見を書くかどうかはここで決まる")
msg("    ", FOCUS, " が検出深度やライブラリの性質と相関していれば、")
msg("    「進行期ほど検出損失が小さい」は生物ではなく技術の話になる。")
msg("")
## **3 群に分けて出す。** 以前はまとめて出していたため、応答変数の構成要素
## （d_zr_orphan など）が「要注意」と表示され、交絡と読み違える原因になった。
## stage と d_zr_orphan が相関するのは所見そのものであって交絡ではない。
GRP <- list(
  "技術変数（★ 交絡の判定はここだけで行う）" = c("n_kept", "n_orphan"),
  "ゼロ率の水準（正常側との相関は交絡の疑い）" =
    c("zr_all_normal", "zr_all_disease", "zr_orph_normal", "zr_orph_disease"),
  "応答変数の構成要素（相関は所見そのもの。交絡ではない）" =
    c("d_zr_all", "d_zr_orphan", "d_zr_control"))
okcol <- function(v) v %in% names(d) && is.numeric(d[[v]]) &&
                     is.finite(sd(d[[v]])) && sd(d[[v]]) > 0
GRP  <- lapply(GRP, function(vs) vs[vapply(vs, okcol, logical(1))])
tech <- unlist(GRP, use.names = FALSE)
fx <- d[[FOCUS]]
for (g in names(GRP)) {
  if (!length(GRP[[g]])) next
  msg("    [", g, "]")
  for (v in GRP[[g]]) {
    if (is.factor(fx)) {
      kt <- suppressWarnings(kruskal.test(d[[v]], fx))
      est <- NA_real_; pv2 <- kt$p.value; lab <- "Kruskal-Wallis"
    } else {
      ct <- suppressWarnings(cor.test(fx, d[[v]], method = "spearman"))
      est <- unname(ct$estimate); pv2 <- ct$p.value; lab <- "Spearman"
    }
    ## 警告を出すのは交絡になりうる 2 群だけ。応答の構成要素には出さない。
    flag <- if (!is.na(pv2) && pv2 < 0.05 && !grepl("^応答", g)) "  ** 要注意 **" else ""
    msg(sprintf("      %-16s %s rho = %s p = %.4g%s", v, lab,
                if (is.na(est)) "   -   " else sprintf("%+.3f ", est), pv2, flag))
  }
}
msg("")
msg("    判定は **技術変数の群だけ** で行う。")
msg("      n_kept   … その患者で解析対象になった転写産物の本数（検出の広さ）")
msg("      n_orphan … うちオーファンの本数")
msg("    ここが有意なら深度の交絡を疑い、(C) の + n_kept / + n_orphan の行で")
msg("    係数が生き残るかを見る。生き残る（むしろ強くなる）なら交絡ではない。")
msg("")
msg("    ゼロ率の水準の群では、**正常側**（zr_*_normal）との相関に注意する。")
msg("    正常組織は疾患の進行を知らないはずなので、そこが相関していれば")
msg("    患者側かライブラリ側の何かが一緒に動いている。Limitations に書く。")
msg("")
msg("    応答変数の構成要素の群は参考。excess = d_zr_orphan - d_zr_control なので")
msg("    相関するのは当然であり、**交絡ではなく所見の内訳である**。")
msg("    ここの相関を根拠に所見を取り下げてはならない。")

## -------------------------------------------- (G) 並べ替え検定
hr(); msg("(G) 並べ替え検定（", FOCUS, " のラベルだけ入れ替え、", NPERM, " 回）")
set.seed(SEED)
obs <- c0[1]
pv  <- replicate(NPERM, {
  dd2 <- d; dd2[[FOCUS]] <- sample(dd2[[FOCUS]])
  co(lm(fml(BASE), data = dd2), FOCUS)[1]
})
pperm <- (1 + sum(abs(pv) >= abs(obs))) / (NPERM + 1)
msg(sprintf("    観測 %+.5f   帰無分布 sd %.5f", obs, sd(pv, na.rm = TRUE)))
msg(sprintf("    両側の並べ替え p = %.4g   （下限 %.5f）", pperm, 1/(NPERM+1)))
msg("    正規性や等分散の仮定を置かない p。(A) の p と大きく違うなら")
msg("    こちらを報告すること。")

hr(); msg("まとめて判断するときの順序")
msg("  1. (F) で n_kept / n_orphan と強く相関していないか。していたら技術交絡。")
msg("  2. (E) で 1 人抜くと消えないか。消えるなら所見として弱い。")
msg("  3. (B) で単調か。単調でなければ「進行期ほど」とは書けない。")
msg("  4. (C) で他の共変量を入れても残るか。")
msg("  5. そのうえで、BH 後に 0.05 を超えることを本文に明記する。")
msg("     生 p・BH 後 p・BF10 を全部書き、確認的ではなく所見として提示する。")
