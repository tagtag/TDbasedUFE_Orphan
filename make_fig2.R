#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## make_fig2.R  —  Figure 2 を作る
##
## 患者ごとの超過 E_i = d_zr_orphan - d_zr_control を、全体の検出移動
## d_zr_all に対して散布し、回帰直線と**切片**を描く。4 パネル。
##   A/B/C  残す 3 コホート
##   D      GSE40419（除外。交絡が作る値の較正）
##
## この図で本文から落とせるもの:
##   - 「切片は内挿であって外挿ではない」… x = 0 の破線が点群の内側にあることが見える
##   - 「GSE40419 では 13 倍」          … D パネルに A-C の y 範囲を矩形で重ねる
##   - 「患者 18 のレバレッジ」          … C パネルで別記号 + ラベル
##
## ---------------------------------------------------------------------------
## 【使い方】slope_*.csv のあるディレクトリで
##   Rscript make_fig2.R
##   OUT=Fig2.pdf W=7.2 H=6.4 Rscript make_fig2.R
##   COHORTS=GSE244679,GSE127165,GSE144269,GSE40419 Rscript make_fig2.R
##
## 【環境変数】
##   DIR       既定 .        slope_<COHORT>.csv の場所
##   OUT       既定 Fig2.pdf 出力（.pdf 推奨。.png も可）
##   W, H      既定 7.2, 6.4 インチ。MDPI の \linewidth に収まる比率
##   SHARE_Y   既定 1        A-C の y 軸を共通にする（0 で各パネル自由）
##   PNG       既定 1        確認用に同名 .png も出す
##
## 【依存】base R のみ。外部パッケージ不要。
##
## 【配色】dataviz の検証済みパレットから 2 色だけ使う。
##   青 #2a78d6（患者）／橙 #eb6834（切片）
##   validate_palette.js で全項目 PASS（CVD ΔE 24.7、通常視 33.6）。
##   **色に加えて形と位置でも区別させてあるので白黒印刷でも読める。**
##   除外コホートは色ではなく**地のトーンとラベル**で示す（色は量を表し、
##   採否を表さない）。
## ---------------------------------------------------------------------------

DIR     <- Sys.getenv("DIR", ".")
OUT     <- Sys.getenv("OUT", "Fig2.pdf")
## 既定を MDPI の 1 段幅（約 16 cm = 6.3 in）に近づけてある。\linewidth に
## 入れたときの縮小率が小さいほど、パネル内の文字が小さくなりすぎない。
W       <- as.numeric(Sys.getenv("W", "6.7"))
H       <- as.numeric(Sys.getenv("H", "6.0"))
SHARE_Y <- Sys.getenv("SHARE_Y", "1") == "1"
PNG     <- Sys.getenv("PNG", "1") == "1"
COS     <- trimws(strsplit(Sys.getenv("COHORTS",
             "GSE244679,GSE127165,GSE144269,GSE40419"), ",")[[1]])
COS     <- COS[nzchar(COS)]

LAB <- c(GSE244679 = "Psoriasis", GSE127165 = "LSCC",
         GSE144269 = "HCC",       GSE40419  = "LAC")
EXCLUDED <- "GSE40419"

## ---- 配色（検証済み。2 色だけ） -------------------------------------------
SERIES  <- "#2a78d6"   # 患者
ACCENT  <- "#eb6834"   # 切片
INK     <- "#0b0b0b"
MUTED   <- "#52514e"
GRID    <- "#e6e5e2"
BAND    <- "#d8d7d3"
TINT    <- "#f2f1ee"   # 除外パネルの地
SURFACE <- "#ffffff"

msg <- function(...) cat(..., "\n", sep = "")

## ---- 読み込みと当てはめ ----------------------------------------------------
fits <- list()
for (CO in COS) {
  f <- file.path(DIR, sprintf("slope_%s.csv", CO))
  if (!file.exists(f)) { msg("（", basename(f), " がありません。飛ばします）"); next }
  d <- read.csv(f, stringsAsFactors = FALSE)
  need <- c("d_zr_all", "d_zr_orphan", "d_zr_control")
  if (!all(need %in% names(d))) { msg("（", basename(f), " の列が足りません）"); next }
  d$excess <- d$d_zr_orphan - d$d_zr_control
  fit <- lm(excess ~ d_zr_all, data = d)
  s   <- summary(fit)$coefficients
  d$cook <- cooks.distance(fit)
  fits[[CO]] <- list(d = d, fit = fit,
                     a = s[1, 1], a_se = s[1, 2],
                     b = s[2, 1], b_se = s[2, 2])
}
if (!length(fits)) stop("読める CSV がありません。", call. = FALSE)

keep <- setdiff(names(fits), EXCLUDED)

## A-C の共通 y 範囲（D に矩形として重ねる）
ylim_keep <- if (length(keep)) {
  r <- range(unlist(lapply(fits[keep], function(z) z$d$excess)), finite = TRUE)
  r + c(-1, 1) * 0.06 * diff(r)
} else NULL

## ---- 1 パネル --------------------------------------------------------------
panel <- function(CO, letter) {
  z <- fits[[CO]]; d <- z$d
  excl <- identical(CO, EXCLUDED)

  xr <- range(d$d_zr_all, finite = TRUE)
  xr <- xr + c(-1, 1) * 0.08 * diff(xr)
  yr <- if (!excl && SHARE_Y && !is.null(ylim_keep)) ylim_keep else {
    b <- c(d$excess, z$a - 2*z$a_se, z$a + 2*z$a_se)
    ## **除外パネルは A-C の帯を同じ軸に載せる。**これをしないと矩形が
    ## 表示域の外に落ちて「13 倍」が図から読めなくなる。
    if (excl && !is.null(ylim_keep)) b <- c(b, ylim_keep)
    r <- range(b, finite = TRUE)
    r + c(-1, 1) * 0.08 * diff(r)
  }

  plot(NA, xlim = xr, ylim = yr, xlab = "", ylab = "", axes = FALSE,
       xaxs = "i", yaxs = "i")
  ## 地（除外パネルだけトーンを敷く。色ではなく明度で示す）
  rect(xr[1], yr[1], xr[2], yr[2], col = if (excl) TINT else SURFACE, border = NA)

  ## 目盛と recessive なグリッド
  xt <- pretty(xr, 5); yt <- pretty(yr, 5)
  abline(v = xt, col = GRID, lwd = 0.6)
  abline(h = yt, col = GRID, lwd = 0.6)

  ## 信頼帯と回帰直線
  nx <- seq(xr[1], xr[2], length.out = 200)
  pr <- predict(z$fit, newdata = data.frame(d_zr_all = nx), interval = "confidence")
  polygon(c(nx, rev(nx)), c(pr[, "lwr"], rev(pr[, "upr"])),
          col = adjustcolor(BAND, alpha.f = 0.55), border = NA)
  lines(nx, pr[, "fit"], col = INK, lwd = 1.8)

  ## 参照線。x = 0 が切片を取る位置
  abline(h = 0, col = MUTED, lty = 3, lwd = 0.9)
  abline(v = 0, col = MUTED, lty = 2, lwd = 1.1)

  ## 患者。白のリングで重なりを読めるようにする
  hi <- which.max(d$cook)
  big <- !excl && d$cook[hi] > 1        # 影響点を別記号にする閾値
  pchv <- rep(21, nrow(d)); if (big) pchv[hi] <- 24
  points(d$d_zr_all, d$excess, pch = pchv,
         bg = adjustcolor(SERIES, alpha.f = 0.75), col = SURFACE,
         cex = if (nrow(d) > 60) 0.85 else 1.0, lwd = 0.7)

  ## 切片。橙の菱形＋誤差棒。色・形・位置の 3 つで区別される
  ## （誤差棒が画面上で潰れる場合 arrows は角度不定で落ちるので、線に落とす）
  if (2 * z$a_se > 0.004 * diff(yr)) {
    arrows(0, z$a - z$a_se, 0, z$a + z$a_se, code = 3, angle = 90,
           length = 0.03, col = ACCENT, lwd = 1.6)
  } else {
    segments(0, z$a - z$a_se, 0, z$a + z$a_se, col = ACCENT, lwd = 1.6)
  }
  points(0, z$a, pch = 23, bg = ACCENT, col = SURFACE, cex = 1.5, lwd = 0.9)

  ## D パネルに A-C の y 範囲を重ねる。
  ## **rect() は使わない。**左右の辺がパネルの枠線と重なって見えないので、
  ## 「破線の矩形」と説明しても読者にはそう見えない。実際に見えるのは
  ## 水平な破線 2 本なので、最初からそう描いてキャプションもそう書く。
  if (excl && !is.null(ylim_keep)) {
    abline(h = ylim_keep, col = MUTED, lty = 2, lwd = 1.0)
    text(xr[1] + 0.03 * diff(xr), ylim_keep[2], adj = c(0, 1.45),
         labels = "range of panels A-C", col = MUTED, cex = 0.62)
  }

  ## 影響点のラベル
  if (big) text(d$d_zr_all[hi], d$excess[hi], adj = c(1.25, 0.4),
                labels = sprintf("Cook's D = %.2f", d$cook[hi]),
                col = MUTED, cex = 0.62)

  box(col = MUTED, lwd = 0.8)
  axis(1, at = xt, col = MUTED, col.axis = MUTED, cex.axis = 0.72,
       tck = -0.018, mgp = c(2, 0.35, 0))
  axis(2, at = yt, col = MUTED, col.axis = MUTED, cex.axis = 0.72,
       tck = -0.018, las = 1, mgp = c(2, 0.45, 0))

  ## 見出しと数値（直接ラベル。点ごとの数値は打たない）
  mtext(sprintf("%s  %s", letter, CO), side = 3, line = 0.55, adj = 0,
        cex = 0.78, font = 2, col = INK)
  nm <- if (CO %in% names(LAB)) LAB[[CO]] else CO
  mtext(sprintf("%s, n = %d%s", nm, nrow(d),
                if (excl) "  (excluded)" else ""),
        side = 3, line = -0.1, adj = 0, cex = 0.66, col = MUTED)
  ## 直接ラベル。**legend() は使わない。**ラベル列と数値列を別々に打って
  ## 桁を揃え、切片の行の頭に橙の菱形を置いて「色 = 切片」を紐づける。
  ## ± は plotmath の %+-% で出す。**文字列に非 ASCII を入れると device に
  ## よっては ".." に化ける。**
  ## 置き場所は**点群と回帰線を避けて自動で選ぶ。**固定の topleft だと
  ## 除外パネルのように点が左上に寄るコホートで文字が data に重なる。
  ## ブロックの幅は**実測する。**決め打ちにすると右寄せの候補で枠から出る。
  CEXK <- 0.66
  v1 <- as.expression(bquote(.(sprintf("%+.4f", z$a)) %+-% .(sprintf("%.4f", z$a_se))))
  v2 <- as.expression(bquote(.(sprintf("%+.3f", z$b)) %+-% .(sprintf("%.3f", z$b_se))))
  wpad <- 0.045 * diff(xr)                                   # 菱形の分
  wlab <- max(strwidth(c("intercept", "slope"), cex = CEXK)) * 1.18
  wval <- max(strwidth(v1, cex = CEXK), strwidth(v2, cex = CEXK))
  BWu  <- wpad + wlab + wval + 0.03 * diff(xr)               # user 単位
  BW   <- BWu / diff(xr)                                     # 画面比
  BH   <- 0.135
  fx <- (d$d_zr_all - xr[1]) / diff(xr)
  fy <- (d$excess   - yr[1]) / diff(yr)
  ## 候補は上・下・中の左右 6 か所。**順番が優先順位**で、同点なら先頭が勝つ。
  xL <- 0.035; xR <- max(xL, 0.975 - BW)
  cand <- list(c(xL, 0.975), c(xR, 0.975), c(xL, 0.175),
               c(xR, 0.175), c(xL, 0.560), c(xR, 0.560))
  ## 除外パネルで引く水平線（A-C の帯の上下端）。文字を重ねない
  hline <- if (excl && !is.null(ylim_keep))
    (ylim_keep - yr[1]) / diff(yr) else numeric(0)
  wrngf <- strwidth("range of panels A-C", cex = 0.62) / diff(xr)
  cost <- vapply(cand, function(cc) {
    x0 <- cc[1]; x1 <- x0 + BW; y1 <- cc[2]; y0 <- y1 - BH
    n  <- sum(fx >= x0 - 0.02 & fx <= x1 + 0.02 & fy >= y0 - 0.02 & fy <= y1 + 0.02)
    ## 回帰線がブロックを横切るか
    xe <- xr[1] + c(x0, x1) * diff(xr)
    ye <- (predict(z$fit, data.frame(d_zr_all = xe)) - yr[1]) / diff(yr)
    if (max(ye) >= y0 - 0.02 && min(ye) <= y1 + 0.02) n <- n + 6
    ## 帯の枠線そのもの
    for (h in hline) if (h >= y0 - 0.02 && h <= y1 + 0.02) n <- n + 6
    ## 帯の上端の直下に置く "range of panels A-C" の行。**幅は実測して判定する。**
    ## ここを「左半分なら」と決め打ちにすると、ブロックが広いコホートで
    ## 逃げ場がなくなって全候補が同点になり、結局重なる。
    if (length(hline)) {
      h <- max(hline)
      if (x0 <= 0.03 + wrngf && x1 >= 0.03 && y1 >= h - 0.10 && y0 <= h)
        n <- n + 6
    }
    n
  }, 0)
  cc <- cand[[which.min(cost)]]
  X0 <- xr[1] + cc[1] * diff(xr)
  gx <- X0 + 0.016 * diff(xr)                # 記号
  lx <- X0 + wpad                            # ラベル
  vx <- lx + wlab                            # 数値
  ly <- yr[1] + (cc[2] - c(0.035, 0.105)) * diff(yr)
  ## 薄い地を敷く。候補選択で避けきれなかった罫線（x = 0 の破線など）の上に
  ## 数字が乗っても読めるようにするため。枠線は付けない
  rect(X0, ly[2] - 0.042 * diff(yr), X0 + BWu, ly[1] + 0.042 * diff(yr),
       col = adjustcolor(if (excl) TINT else SURFACE, alpha.f = 0.82),
       border = NA)
  points(gx, ly[1], pch = 23, bg = ACCENT, col = SURFACE, cex = 1.0, lwd = 0.7)
  text(lx, ly[1], "intercept", adj = c(0, 0.5), cex = CEXK, col = INK)
  text(lx, ly[2], "slope",     adj = c(0, 0.5), cex = CEXK, col = MUTED)
  text(vx, ly[1], labels = v1, adj = c(0, 0.5), cex = CEXK, col = INK)
  text(vx, ly[2], labels = v2, adj = c(0, 0.5), cex = CEXK, col = MUTED)
}

## ---- 描画 ------------------------------------------------------------------
draw <- function() {
  op <- par(mfrow = c(2, 2), mar = c(2.6, 3.0, 2.2, 0.8),
            oma = c(2.2, 2.0, 0.4, 0.4), family = "sans", bg = SURFACE)
  on.exit(par(op))
  ord <- c(keep, intersect(EXCLUDED, names(fits)))
  for (i in seq_along(ord)) panel(ord[i], LETTERS[i])
  ## 軸名は **1 回の mtext** で出す。2 回に分けると同じ line に重ね打ちされる。
  mtext(expression("Shift in overall zero rate, " * Delta^{"all"} *
                   " (disease " - " normal)"),
        side = 1, outer = TRUE, line = 0.7, cex = 0.78, adj = 0.5, col = INK)
  mtext(expression("Orphan-specific excess, " * italic(E)),
        side = 2, outer = TRUE, line = 0.4, cex = 0.78, adj = 0.5, col = INK)
}

## cairo が無い環境（--with-x なしの R など）では既定の device に落とす
open_png <- function(path, res) {
  ok <- try(png(path, width = W, height = H, units = "in", res = res,
                type = "cairo"), silent = TRUE)
  if (inherits(ok, "try-error"))
    png(path, width = W, height = H, units = "in", res = res)
  invisible(TRUE)
}

## **PDF は cairo_pdf を優先する。**素の pdf() は base-14 を埋め込まないので、
## Symbol（Delta）が環境によって欠落し、投稿先のフォント埋め込み要件も満たさない。
## 出力後に `pdffonts Fig2.pdf` で emb が yes になっていることを確認すること。
open_pdf <- function(path) {
  ok <- try(cairo_pdf(path, width = W, height = H), silent = TRUE)
  if (inherits(ok, "try-error")) {
    msg("（cairo_pdf が使えないので pdf() に落とします。",
        "**フォントが埋め込まれないので投稿前に確認すること**）")
    pdf(path, width = W, height = H, useDingbats = FALSE)
  }
  invisible(TRUE)
}

ext <- tolower(tools::file_ext(OUT))
if (ext == "png") {
  open_png(OUT, 600); draw(); dev.off()
} else {
  open_pdf(OUT)
  draw(); dev.off()
  if (PNG) {
    p <- sub("\\.[Pp][Dd][Ff]$", ".png", OUT)
    ok <- try({ open_png(p, 200); draw(); dev.off() }, silent = TRUE)
    if (!inherits(ok, "try-error")) msg("確認用: ", p)
  }
}
msg("書き出し: ", OUT)

## ---- 図の中身の要約（キャプションの数値を照合するため） --------------------
cat("\n")
cat(sprintf("%-11s %4s %11s %11s %9s %9s\n",
            "cohort", "n", "intercept", "se", "slope", "x range"))
for (CO in c(keep, intersect(EXCLUDED, names(fits)))) {
  z <- fits[[CO]]
  cat(sprintf("%-11s %4d %+11.5f %11.5f %+9.3f  %+.3f..%+.3f\n",
              CO, nrow(z$d), z$a, z$a_se, z$b,
              min(z$d$d_zr_all), max(z$d$d_zr_all)))
}
cat("\n")
cat("本文・キャプションに写す前に、上の値が Table 2 と一致することを確認すること。\n")
cat("一致しない場合は slope_*.csv が NSET=200 のものか確かめる（論文の値は 200）。\n")
