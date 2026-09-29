#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## de_framework.R — R3-3（DESeq2／edgeR と比較せよ）と R3-2（転写産物間相関）への回答
##
## 主張がコホート水準になったので、「患者ごとに N=1 ペアなので DE 枠組みは走らない」
## という反論は成立しない。患者をブロック因子にした標準的な対応のある設計である。
## そこで実際に走らせる。
##
##   design = ~ patient + condition        （edgeR/limma-voom、必要なら DESeq2）
##
## ただし DE 枠組みの出力は転写産物ごとの検定結果で、本研究の主張は
## **集合対集合の対比**である。したがって DE は代替ではなく入力として使う。
## log2FC を取ったうえで、orphan 集合とマッチ対照集合を competitive に比較する。
## 使うのは limma の camera で、**転写産物間相関を明示的に補正する**ので
## R3-2（相関で独立性が破れる）にも同時に答える。
##
## 比較の作り方が重要である。camera は「集合 対 それ以外の全遺伝子」を比べるので、
## そのままでは発現水準の違いを拾ってしまう。そこで解析対象を
## **orphan ∪ マッチ対照** に限定してから camera をかける。これで対比は
## 本研究の他の解析と同じ「orphan 対 発現水準・検出頻度マッチ対照」になる。
##
## そして DE 枠組みの限界を必ず報告する。filterByExpr は低カウント転写産物を
## 落とすので、解析対象は「検出されている orphan」の部分集合である。
## 本研究の信号はゼロ率に乗っているので、**DE が見ているのは狭い問い**である。
## 何本の orphan が filtering を生き残ったかを出力の先頭に出す。
##
## 使い方
##   # まず配管だけ確認（Bioconductor 不要）
##   DRY_RUN=1 KALLISTO_ROOT=Revised/GSE244679 METADATA_DIR=Revised/metadata \
##   COHORT=GSE244679 Rscript de_framework.R
##
##   # 本番
##   NSET=20 ENGINE=voom KALLISTO_ROOT=Revised/GSE244679 \
##   METADATA_DIR=Revised/metadata COHORT=GSE244679 Rscript de_framework.R
##
## 環境変数
##   ENGINE    既定 "voom"   "voom" | "edger" | "deseq2"
##                           voom を推奨。患者数が 57/70 でブロックが多いので
##                           limma-voom が最も安定して速い。DESeq2 は査読者が
##                           名指ししたので 1 コホートで確認用に走らせる。
##                           DESeq2 はブロック係数が多いと非常に遅いので、
##                           患者 24 人の GSE244679 で確認すること。
##                           70 人のコホートでは数時間かかりうる。
##   NSET      既定 20       マッチ対照集合の反復数。camera を各回かけて分布を見る。
##   MIN_COUNT 既定 10       filterByExpr の min.count。既定は edgeR の既定値。
##   DRY_RUN   既定 0        1 なら Bioconductor を使わず配管だけ確認する。
##   N_CORES   既定 1        abundance.tsv の読み込みの並列数。読み込みが
##                           所要時間のほぼ全部なので、まずこれを上げる。
##   REBUILD   既定 0        1 でカウント行列のキャッシュを無視して読み直す。
##
## 所要時間はほぼ abundance.tsv の読み込み（198,507 行 × ライブラリ数）で決まる。
## カウント行列は counts_<COHORT>.rds にキャッシュするので、DRY_RUN の確認と
## 本番で読み直さない。data.table があれば fread を使い数倍速くなる。
## ---------------------------------------------------------------------------

source("R/config.R")
source("R/00_functions.R")

CO        <- Sys.getenv("COHORT", names(COHORTS)[1])
ENGINE    <- Sys.getenv("ENGINE", "voom")
NSET      <- as.integer(Sys.getenv("NSET", "20"))
MIN_COUNT <- as.numeric(Sys.getenv("MIN_COUNT", "10"))
N_STRATA2 <- as.integer(Sys.getenv("N_STRATA2", "4"))
DRY_RUN   <- nzchar(Sys.getenv("DRY_RUN", "")) && Sys.getenv("DRY_RUN") != "0"
N_CORES   <- as.integer(Sys.getenv("N_CORES", "1"))
REBUILD   <- nzchar(Sys.getenv("REBUILD", "")) && Sys.getenv("REBUILD") != "0"
stopifnot(ENGINE %in% c("voom", "edger", "deseq2"))

msg <- function(...) cat(..., "\n", sep = "")
hr  <- function() cat(strrep("-", 76), "\n")

## 段階ごとの所要時間。どこが遅いのか分からないまま待つことをなくす。
.T0 <- proc.time()[3]; .TL <- .T0
tick <- function(label) {
  now <- proc.time()[3]
  message(sprintf("  [%6.1f 秒] %s（この段階 %.1f 秒）", now - .T0, label, now - .TL))
  .TL <<- now
}

## N_CORES が効く範囲を明示する。
## 読み込み          … mclapply で並列（ただしキャッシュがあると飛ぶ）
## filterByExpr/TMM  … edgeR 内部。単一スレッド
## voom / lmFit      … limma 内部。**単一スレッド。N_CORES では変わらない**
##                     行列演算は BLAS 依存なので、速くしたいなら R に
##                     マルチスレッド BLAS（OpenBLAS 等）を入れる話になる
## camera × NSET     … mclapply で並列（下記）

## ---------------------------------------------------------------------------
## 依存パッケージ。足りないものは名前を挙げて止める（黙って別の道に逃げない）。
## ---------------------------------------------------------------------------
need <- switch(ENGINE, voom = c("limma", "edgeR"), edger = c("edgeR", "limma"),
               deseq2 = c("DESeq2", "limma", "edgeR", "SummarizedExperiment"))
if (!DRY_RUN) {
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("必要なパッケージがありません: ", paste(miss, collapse = ", "), "\n",
         '  BiocManager::install(c("', paste(miss, collapse = '", "'), '"))\n',
         "  配管だけ確認するなら DRY_RUN=1 を付けてください。", call. = FALSE)
}

## ---------------------------------------------------------------------------
## カウント行列。kallisto の est_counts をそのまま使う。
##
## tximport があれば使う（査読者向けには「標準パイプライン」と言える）。
## なければ est_counts を直接読む。DE 枠組みは整数カウントを期待するので
## 丸めるが、丸めによる違いは effective length 補正を入れないことより小さい。
## ---------------------------------------------------------------------------
## data.table があれば fread。なければ read.delim で必要 2 列だけ読む。
.HAVE_DT <- requireNamespace("data.table", quietly = TRUE)
read_one <- function(path) {
  if (.HAVE_DT) {
    d <- data.table::fread(path, select = c("target_id", "est_counts"),
                           showProgress = FALSE, data.table = FALSE)
  } else {
    ## colClasses で不要列を捨てると read.delim も大幅に速くなる。
    hd <- strsplit(readLines(path, n = 1L), "\t")[[1]]
    cc <- ifelse(hd == "target_id", "character",
                 ifelse(hd == "est_counts", "numeric", "NULL"))
    d <- utils::read.delim(path, colClasses = cc, stringsAsFactors = FALSE)
  }
  stopifnot(all(c("target_id", "est_counts") %in% names(d)))
  d
}

build_counts <- function(sheet) {
  libs  <- sheet$library_id
  paths <- vapply(seq_along(libs), function(k)
    resolve_abundance_path(KALLISTO_ROOT, sheet$directory[k], "abundance.tsv"), "")

  message("  ", length(libs), " 本を読み込みます（cores = ", N_CORES,
          if (.HAVE_DT) ", data.table::fread" else ", read.delim（2 列のみ）", "）")
  lst <- if (N_CORES > 1L) {
    parallel::mclapply(paths, read_one, mc.cores = N_CORES)
  } else {
    lapply(paths, read_one)
  }
  bad <- which(!vapply(lst, is.data.frame, logical(1)))
  if (length(bad))
    stop(length(bad), " 本の読み込みに失敗しました。最初: ", paths[bad[1]],
         "\n  ", paste(utils::head(as.character(lst[[bad[1]]]), 1), collapse = " "),
         call. = FALSE)

  ids <- lst[[1]]$target_id
  m <- matrix(0, length(ids), length(libs), dimnames = list(ids, libs))
  for (k in seq_along(libs)) {
    d <- lst[[k]]
    if (!identical(d$target_id, ids)) {
      ## 行順が違うだけなら並べ替える。集合が違うなら止める。
      if (!setequal(d$target_id, ids))
        stop("転写産物集合がライブラリ間で一致しません: ", libs[k], call. = FALSE)
      d <- d[match(ids, d$target_id), ]
    }
    m[, k] <- d$est_counts
  }
  round(m)
}

sheet <- read_sheet(CO)
## 対応のとれた患者だけ使う。片方しかない患者はブロック因子が推定できない。
tab <- table(sheet$patient_id, sheet$condition)
okp <- rownames(tab)[tab[, "normal"] >= 1 & tab[, "disease"] >= 1]
sheet <- sheet[sheet$patient_id %in% okp, ]
sheet <- sheet[order(sheet$patient_id, match(sheet$condition, c("normal", "disease"))), ]
np <- length(unique(sheet$patient_id))
stopifnot(np >= 6L)

## カウント行列をキャッシュする。DRY_RUN の確認と本番で読み直さないため。
## キャッシュの鍵はコホート名とライブラリ集合。サンプルシートが変わったら
## 作り直す（黙って古い行列を使うと原因不明の不一致になる）。
CACHE <- sprintf("counts_%s.rds", CO)
cnt <- NULL
if (!REBUILD && file.exists(CACHE)) {
  cc <- try(readRDS(CACHE), silent = TRUE)
  if (!inherits(cc, "try-error") && is.list(cc) && is.matrix(cc$m) &&
      identical(cc$libs, sheet$library_id)) {
    cnt <- cc$m
    message("  キャッシュを再利用: ", CACHE,
            "（読み直すには REBUILD=1）")
  } else {
    message("  キャッシュ ", CACHE, " はライブラリ集合が違うので作り直します。")
  }
}
if (is.null(cnt)) {
  t0 <- proc.time()[3]
  cnt <- build_counts(sheet)
  saveRDS(list(m = cnt, libs = sheet$library_id), CACHE)
  message(sprintf("  読み込み %.0f 秒。キャッシュ: %s", proc.time()[3] - t0, CACHE))
}
is_o <- orphan_index(rownames(cnt), ORPHAN_ID_FILE, expected_n = N_ORPHAN_TRANSCRIPT)

patient   <- factor(sheet$patient_id)
condition <- factor(sheet$condition, levels = c("normal", "disease"))
design    <- model.matrix(~ patient + condition)
COEF      <- "conditiondisease"
stopifnot(COEF %in% colnames(design))

hr(); msg(CO, " / 患者 ", np, " 人 / ライブラリ ", ncol(cnt), " 本 / engine = ", ENGINE)
msg("  design = ~ patient + condition   （患者をブロック因子にした対応のある設計）")
msg("  転写産物 ", nrow(cnt), " 本 / うち orphan ", sum(is_o), " 本")

## ---------------------------------------------------------------------------
## filterByExpr。**これが DE 枠組みの限界の本体である。**
## 低カウント転写産物が落ちるので、orphan の何割が生き残るかを必ず報告する。
## ---------------------------------------------------------------------------
if (!DRY_RUN) {
  y    <- edgeR::DGEList(cnt)
  keep <- edgeR::filterByExpr(y, design, min.count = MIN_COUNT)
} else {
  ## filterByExpr の簡易版（CPM 換算で min.count 相当を最小群サイズ本で満たすか）。
  ## 配管確認用。数値は本番と一致しないが桁は合う。
  cpm0 <- t(t(cnt) / pmax(colSums(cnt), 1)) * 1e6
  thr  <- MIN_COUNT / (median(colSums(cnt)) / 1e6)
  keep <- rowSums(cpm0 >= thr) >= min(table(condition))
}

hr(); msg("filterByExpr の通過状況（DE 枠組みの限界の本体）")
msg(sprintf("  全転写産物  %6d → %6d 本 通過（%.1f%%）",
            nrow(cnt), sum(keep), 100 * mean(keep)))
msg(sprintf("  orphan      %6d → %6d 本 通過（%.1f%%）",
            sum(is_o), sum(keep & is_o), 100 * sum(keep & is_o) / sum(is_o)))
msg(sprintf("  非 orphan   %6d → %6d 本 通過（%.1f%%）",
            sum(!is_o), sum(keep & !is_o), 100 * sum(keep & !is_o) / sum(!is_o)))
msg("  orphan の通過率が非 orphan より低ければ、それ自体が DE 枠組みでは")
msg("  この現象を十分に見られないことの定量である。Limitations に書く。")
if (sum(keep & is_o) < 20L)
  msg("  ** orphan が 20 本未満しか残らない。集合レベルの検定は意味を持たない。**")

## ---------------------------------------------------------------------------
## マッチ対照集合。コホート水準なので全 normal の平均と検出頻度で層化する
## （患者ごとの leave-one-out は患者ごとの解析でのみ必要）。
## 層化は filter を通過した転写産物の中で行う。DE が見る集合の中で対比するため。
## ---------------------------------------------------------------------------
SRC <- if (file.exists("R/07_interaction_permutation.R"))
         "R/07_interaction_permutation.R" else "07_interaction_permutation.R"
src <- readLines(SRC)
eval(parse(text = paste(src[grep("^make_strata <- function", src):
                            (grep("^run_patient <- function", src) - 1L)],
                        collapse = "\n")))
stopifnot(exists("make_strata_multi"), exists("draw_matched_sets"))

nlibs <- sheet$library_id[sheet$condition == "normal"]
cpm_n <- t(t(cnt[, nlibs, drop = FALSE]) / pmax(colSums(cnt[, nlibs, drop = FALSE]), 1)) * 1e6
mean_n <- rowMeans(cpm_n)
freq_n <- rowMeans(cpm_n > 0)

ki   <- which(keep)
st   <- make_strata_multi(list(mean_n[ki], freq_n[ki]), c(N_STRATA, N_STRATA2))
oi   <- is_o[ki]
stopifnot(sum(oi) >= 2L, sum(!oi) >= 2L)

set.seed(BASE_SEED)
ctrl_sets <- vector("list", NSET); short <- 0L
for (k in seq_len(NSET)) {
  s1 <- draw_matched_sets(st, oi, n_sets = 1L)
  short <- short + attr(s1, "shortfall")
  ctrl_sets[[k]] <- s1[[1]]
}
hr(); msg("マッチ対照集合（filter 通過集合の中で、CPM 平均 20 分位 × 検出頻度 4 分位）")
msg(sprintf("  層 %d / orphan %d 本 / 対照 %d 本 × %d 組 / shortfall %d",
            length(unique(st)), sum(oi), length(ctrl_sets[[1]]), NSET, short))
if (short > 0)
  msg("  ** shortfall があるので対照が揃っていない。N_STRATA を下げて再実行。**")

if (DRY_RUN) {
  hr(); msg("DRY_RUN=1 なので DE 本体は走らせません。ここまでの配管は健全です。")
  msg("本番は DRY_RUN を外して実行してください。")
  quit(save = "no", status = 0)
}

## ---------------------------------------------------------------------------
## DE 本体
## ---------------------------------------------------------------------------
tick("ここまで（読み込み・filter・対照集合）")
yk <- edgeR::calcNormFactors(edgeR::DGEList(cnt[ki, , drop = FALSE]), method = "TMM")
tick("TMM 正規化")

if (ENGINE == "voom") {
  v   <- limma::voom(yk, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  tt  <- limma::topTable(fit, coef = COEF, number = Inf, sort.by = "none")
  lfc <- tt$logFC
  obj <- v
} else if (ENGINE == "edger") {
  yk  <- edgeR::estimateDisp(yk, design)
  fit <- edgeR::glmQLFit(yk, design)
  qlf <- edgeR::glmQLFTest(fit, coef = COEF)
  lfc <- qlf$table$logFC
  obj <- edgeR::cpm(yk, log = TRUE)
} else {
  dds <- DESeq2::DESeqDataSetFromMatrix(cnt[ki, , drop = FALSE],
           data.frame(patient = patient, condition = condition), ~ patient + condition)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  ## DESeq2 は係数名を付け替える。model.matrix の "conditiondisease" ではなく
  ## "condition_disease_vs_normal" になるので、resultsNames から引く。
  ## 決め打ちにすると「couldn't find results」で落ちる。
  rn <- DESeq2::resultsNames(dds)
  cn <- grep("^condition", rn, value = TRUE)
  if (length(cn) != 1L)
    stop("DESeq2 の condition 係数が一意に決まりません。resultsNames: ",
         paste(rn, collapse = ", "), call. = FALSE)
  message("  DESeq2 の係数名: ", cn)
  rr  <- DESeq2::results(dds, name = cn, independentFiltering = FALSE)
  lfc <- rr$log2FoldChange
  obj <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
}
stopifnot(length(lfc) == length(ki))
tick(paste0(ENGINE, " の当てはめ"))

## ---------------------------------------------------------------------------
## 集合対集合の対比。camera を orphan ∪ 対照 に限定してかける。
## camera は転写産物間相関を補正するので R3-2 にも同時に答える。
## ---------------------------------------------------------------------------
hr(); msg("log2FC の分布（", COEF, "。負が疾患側で低い）")
qq <- function(v) sprintf("中位 %+.4f  平均 %+.4f  四分位 %+.4f / %+.4f",
                          median(v, na.rm = TRUE), mean(v, na.rm = TRUE),
                          quantile(v, .25, na.rm = TRUE), quantile(v, .75, na.rm = TRUE))
msg("  orphan      ", qq(lfc[oi]))
msg("  対照（1組）  ", qq(lfc[ctrl_sets[[1]]]))

## voom の出力は EList で行列ではないので [i, , drop = FALSE] が通らない。
## EList は [i, ] で行を取れる。行列とどちらでも動く形にする。
sub_rows <- function(x, i) if (is.matrix(x)) x[i, , drop = FALSE] else x[i, ]

one_set_test <- function(k) {
  sub <- c(which(oi), ctrl_sets[[k]])
  sub <- sub[!duplicated(sub)]
  idx <- list(orphan = which(sub %in% which(oi)))
  ## inter.gene.cor = NA で camera に相関を推定させる（既定の 0.01 は決め打ち）。
  ## 推定させると Correlation 列が返り、それが R2-1／R3-2 で問われた量になる。
  ca <- try(limma::camera(sub_rows(obj, sub), idx, design,
                          contrast = COEF, inter.gene.cor = NA), silent = TRUE)
  if (inherits(ca, "try-error") && k == 1L)
    message("  camera が失敗: ", conditionMessage(attr(ca, "condition")))
  wt <- suppressWarnings(wilcox.test(lfc[oi], lfc[ctrl_sets[[k]]]))
  data.frame(set = k,
             mean_orphan  = mean(lfc[oi], na.rm = TRUE),
             mean_control = mean(lfc[ctrl_sets[[k]]], na.rm = TRUE),
             camera_dir = if (inherits(ca, "try-error")) NA_character_ else ca$Direction[1],
             camera_p   = if (inherits(ca, "try-error")) NA_real_ else ca$PValue[1],
             camera_cor = if (inherits(ca, "try-error") || is.null(ca$Correlation))
                            NA_real_ else ca$Correlation[1],
             wilcox_p   = wt$p.value, stringsAsFactors = FALSE)
}

## camera の 20 組は互いに独立なので、ここは並列にできる。
## fork なので obj のコピーは起きない（copy-on-write）。
lst2 <- if (N_CORES > 1L) {
  parallel::mclapply(seq_len(NSET), one_set_test, mc.cores = min(N_CORES, NSET))
} else {
  lapply(seq_len(NSET), one_set_test)
}
bad2 <- which(!vapply(lst2, is.data.frame, logical(1)))
if (length(bad2))
  stop(length(bad2), " 組でエラー。最初: ",
       paste(utils::head(as.character(lst2[[bad2[1]]]), 1), collapse = " "),
       "\n  N_CORES=1 で走らせると原因が見えます。", call. = FALSE)
res <- do.call(rbind, lst2)
tick(paste0("camera × ", NSET, " 組"))
res$diff <- res$mean_orphan - res$mean_control

hr(); msg("competitive 集合検定（orphan 対 マッチ対照、", NSET, " 組）")
msg(sprintf("  log2FC 平均差 (orphan - 対照)  中位 %+.4f   範囲 %+.4f 〜 %+.4f",
            median(res$diff), min(res$diff), max(res$diff)))
msg(sprintf("  差が負（orphan が低い）の組  %d / %d", sum(res$diff < 0), NSET))
if (any(is.finite(res$camera_p))) {
  msg(sprintf("  camera p  中位 %.3g   範囲 %.3g 〜 %.3g",
              median(res$camera_p, na.rm = TRUE), min(res$camera_p, na.rm = TRUE),
              max(res$camera_p, na.rm = TRUE)))
  msg(sprintf("  camera Direction  Down %d / Up %d",
              sum(res$camera_dir == "Down", na.rm = TRUE),
              sum(res$camera_dir == "Up",   na.rm = TRUE)))
  msg(sprintf("  camera が推定した転写産物間相関  中位 %.4f",
              median(res$camera_cor, na.rm = TRUE)))
  msg("  ← この値が R2-1／R3-2 で問われた相関そのものである。camera はこれを")
  msg("    分散に織り込むので、相関を仮定せず補正した p が上の値である。")
} else {
  msg("  camera が全組で失敗しました。上の Wilcoxon で代替してください。")
}
msg(sprintf("  Wilcoxon p  中位 %.3g", median(res$wilcox_p)))

## MIN_COUNT を出力名に入れる。入れないと感度分析が既定の結果を上書きする
## （MATCH_ON で一度踏んだのと同じ失敗）。既定の 10 のときは付けない。
.mc  <- if (MIN_COUNT != 10) sprintf("_mc%g", MIN_COUNT) else ""
out  <- sprintf("de_%s_%s%s.csv", CO, ENGINE, .mc)
out2 <- sprintf("de_logfc_%s_%s%s.csv", CO, ENGINE, .mc)
write.csv(res, out, row.names = FALSE)
write.csv(data.frame(transcript = rownames(cnt)[ki], is_orphan = oi, logFC = lfc),
          out2, row.names = FALSE)
hr(); msg("書き出し: ", out, " と ", out2)

msg("")
msg("読み方")
msg("  camera が Down で有意なら、DE 枠組みでも同じ向きの結論が出る。R3-3 に")
msg("  そのまま答えられ、R3-2 の相関補正も同じ検定で済む。")
msg("  有意でないなら、それは主張の反証ではなく **DE 枠組みの感度の問題** で")
msg("  ある可能性が高い。上の filterByExpr 通過率を根拠として示すこと。")
msg("  どちらの結果でも、通過率と camera の推定相関は本文に載せる。")
