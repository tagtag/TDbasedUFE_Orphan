#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## make_fig1.R — Figure 1（解析の流れ）を作る。旧 overview.pdf の置き換え。
##
## 旧 Figure 1 は「Analysis 1--5」の構成で、その 5 本のうち 3 本が撤回された
## ので描き直しになった。新しい図が示すのは次の 4 つだけである。
##
##   A  コホートの選別（4 → 3。GSE40419 をなぜ外したか）
##   B  患者 1 人あたりの量（K_i → 3 つのゼロ率 → E_i）
##   C  患者をまたぐ推定（E を Δ^all に回帰 → **切片**が主張）
##   D  二次解析（β3、陰性対照、解析単位、感度解析、臨床変数）がどの表になるか
##
## **この図はデータを持たない。**すべて構造の説明なので、数値は主張の中核
## （切片のプール値と異質性）だけを入れてある。それ以外の数字は表に任せる。
##
## 【使い方】
##   Rscript make_fig1.R
##   OUT=Fig1.pdf W=6.7 H=7.6 Rscript make_fig1.R
##
## 【依存】base R のみ。
##
## 【配色】Figure 2 と同じ。青 #2a78d6（主たる流れ）／橙 #eb6834（切片）。
##   validate_palette.js で全項目 PASS。**色・形・位置の 3 つで区別**するので
##   白黒印刷でも読める。除外は色ではなく地のトーンで示す。
##
## 【PDF】cairo_pdf を使う。素の pdf() はフォントを埋め込まない（方針書 §9）。
##
## 【D パネルのテーブル番号は \ref の対象外である（重要）】
##   図の中に焼き込んだ "Table 3" 等は LaTeX の相互参照を通らない。本文の表を
##   増減・移動したら**必ずここを手で合わせること**。現行の対応は
##     Table 1 Table:samples      Table 5 Table:units
##     Table 2 Table:slopes       Table 6 Table:benchmark
##     Table 3 Table:interaction  Table 7 Table:clinical
##     Table 4 Table:excess-null  Table S1 Table:synthetic（巻末。カウンタを戻す）
##   Figure 1 = Fig:ovverview、Figure 2 = Fig:slopes。
##
## 【改訂】
##   2026-09-26  D パネルを 4 箱 → 5 箱。Table 6（Table:benchmark、表現・正規化・
##               モデル・マッチングの感度解析）が抜けていた。Results は 4 小節
##               あり 3.4 がこの表の節なので、二次解析の一覧から落ちているのは
##               不整合だった。幅と字の大きさを詰め、Table 6 は 4 行なので箱の
##               高さを 8.0 → 9.0 にした。
##   2026-09-26  B パネルの orphan の箱を "2,190" → "(2,190 in the reference)"。
##               2,190 は参照配列全体の orphan 数であって、K_i に絞った後の数
##               ではない（絞った後は 2,166 / 2,008 / 1,981、Table S3）。箱が
##               K_i の下流にあるので、無条件に 2,190 と書くと別の量に読める。
## ---------------------------------------------------------------------------

OUT <- Sys.getenv("OUT", "Fig1.pdf")
W   <- as.numeric(Sys.getenv("W", "6.7"))
H   <- as.numeric(Sys.getenv("H", "7.6"))
PNG <- Sys.getenv("PNG", "1") == "1"

SERIES  <- "#2a78d6"; ACCENT <- "#eb6834"
INK     <- "#0b0b0b"; MUTED  <- "#52514e"
LINE    <- "#b9b8b4"; FILL   <- "#eef4fc"   # 青の淡い地
AFILL   <- "#fdeee7"                        # 橙の淡い地
TINT    <- "#f2f1ee"; SURFACE <- "#ffffff"

## --- 部品 -------------------------------------------------------------------
## 箱。lab は複数行可（"\n" で改行）。
bx <- function(x, y, w, h, lab, fill = SURFACE, border = LINE,
               col = INK, cex = 0.62, font = 1, lwd = 0.8) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = border, lwd = lwd)
  text(x, y, lab, cex = cex, col = col, font = font)
}
## 縦の矢印
va <- function(x, y0, y1, col = SERIES, lwd = 1.3)
  arrows(x, y0, x, y1, length = 0.055, angle = 22, col = col, lwd = lwd)
## 横の矢印
ha <- function(x0, x1, y, col = SERIES, lwd = 1.3)
  arrows(x0, y, x1, y, length = 0.055, angle = 22, col = col, lwd = lwd)
## パネルの見出し
pl <- function(x, y, letter, title) {
  text(x, y, letter, adj = c(0, 0.5), cex = 0.80, font = 2, col = INK)
  text(x + 3.0, y, title, adj = c(0, 0.5), cex = 0.72, col = MUTED)
}

draw <- function() {
  op <- par(mar = c(0, 0, 0, 0), family = "sans", bg = SURFACE)
  on.exit(par(op))
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")

  ## ===================== A  コホートの選別 =================================
  pl(2, 98, "A", "Cohort selection")
  co <- c("GSE244679\npsoriasis, 24", "GSE127165\nLSCC, 57",
          "GSE144269\nHCC, 70", "GSE40419\nLAC, 69")
  xs <- c(16, 38, 60, 84)
  for (i in 1:4)
    bx(xs[i], 91.5, 20, 6.6, co[i], fill = if (i == 4) TINT else FILL, cex = 0.56)
  for (i in 1:4) va(xs[i], 88.0, 85.6)
  bx(50, 83.0, 88, 5.2,
     "Library-level diagnostics: accession order, sequencing depth, pseudoalignment rate",
     fill = SURFACE, cex = 0.58)
  va(35, 80.2, 77.8); va(84, 80.2, 77.8, col = MUTED)
  bx(35, 75.0, 56, 5.4, "Retained: 3 cohorts, 151 matched pairs",
     fill = FILL, cex = 0.60, font = 2)
  bx(84, 75.0, 28, 5.4, "Excluded (calibration)", fill = TINT, cex = 0.58)
  text(84, 71.2, "depth and mapping differ,\nsubmissions disjoint",
       cex = 0.50, col = MUTED)

  ## ===================== B  患者 1 人あたりの量 =============================
  pl(2, 66.5, "B", "Per patient: the quantity")
  bx(24, 61.0, 26, 5.0, "Normal library", fill = SURFACE, cex = 0.60)
  bx(56, 61.0, 26, 5.0, "Disease library", fill = SURFACE, cex = 0.60)
  va(24, 58.3, 55.9); va(56, 58.3, 55.9)
  bx(40, 53.2, 58, 5.0,
     expression("Retained set " * italic(K)[italic(i)] *
                ": detected in either library"), fill = FILL, cex = 0.60)
  for (x in c(15, 40, 65)) va(x, 50.6, 48.2)
  segments(15, 50.6, 65, 50.6, col = SERIES, lwd = 1.3)
  bx(15, 45.4, 26, 5.4, expression("All of " * italic(K)[italic(i)]),
     fill = SURFACE, cex = 0.58)
  ## 2,190 は参照配列全体の数。K_i に絞った後の数ではないので但し書きを付ける。
  bx(40, 45.4, 22, 5.4, "Orphan\n(2,190 in the reference)",
     fill = SURFACE, cex = 0.50)
  bx(65, 45.4, 26, 5.4, "Matched control: abundance\nand detection frequency",
     fill = SURFACE, cex = 0.53)
  text(15, 41.2, expression(Delta^{"all"}), cex = 0.66, col = MUTED)
  text(40, 41.2, expression(Delta^{"orph"}), cex = 0.66, col = MUTED)
  text(65, 41.2, expression(Delta^{"ctrl"}), cex = 0.66, col = MUTED)
  text(80, 41.2, "zero rate,\ndisease minus normal",
       adj = c(0, 0.5), cex = 0.50, col = MUTED)
  ## Delta^orph と Delta^ctrl だけが E に入る。Delta^all はパネル C の横軸。
  va(40, 39.4, 36.7); va(65, 39.4, 36.7)
  bx(52.5, 34.2, 32, 5.0,
     expression(italic(E)[italic(i)] == Delta^{"orph"} - Delta^{"ctrl"}),
     fill = AFILL, border = ACCENT, cex = 0.62)
  text(15, 37.0, expression(Delta^{"all"} * " -> panel C"),
       cex = 0.50, col = MUTED)

  ## ===================== C  患者をまたぐ推定 ================================
  pl(2, 28.0, "C", "Across patients: the estimate")
  bx(30, 22.8, 54, 5.2,
     expression("Regress " * italic(E) * " on " * Delta^{"all"} *
                ", one point per patient"), fill = FILL, cex = 0.60)
  ha(57.5, 63.5, 22.8, col = ACCENT)
  bx(81, 22.8, 34, 5.2, "Intercept: the excess\nwhen detection does not move",
     fill = AFILL, border = ACCENT, cex = 0.56, font = 2)
  points(65.5, 18.8, pch = 23, bg = ACCENT, col = SURFACE, cex = 1.2, lwd = 0.8)
  text(67.8, 18.8, adj = c(0, 0.5), cex = 0.60, col = INK,
       labels = expression(bold("+0.0544") %+-% bold("0.0059") *
                           "   pooled,  " * italic(Q) * " test " *
                           italic(P) == 0.163))
  text(30, 18.8, "Table 2, Figure 2", cex = 0.54, col = MUTED)

  ## ===================== D  二次解析 ========================================
  ## **箱を増減したら xd / wd の合計を確認すること。**
  ## 合計 15+18+19+20+14 = 86、間隔 4 × 1.5 = 6、合わせて 92。x = 4〜96 に収まる。
  pl(2, 12.6, "D", "Secondary analyses, and where they are reported")
  d <- c("Interaction\ncoefficient\nTable 3",
         "Negative controls\non the pairing\nTables 3-4",
         "Transcript versus\npatient as the unit\nTable 5",
         "Representation,\nnormalisation,\nmodel, matching\nTable 6",
         "Clinical\nvariables\nTable 7")
  xd <- c(11.5, 29.5, 49.5, 70.5, 89.0)
  wd <- c(15, 18, 19, 20, 14)
  for (i in 1:5) bx(xd[i], 6.4, wd[i], 9.0, d[i], fill = SURFACE, cex = 0.50)

  ## パネルを分ける薄い罫線
  for (y in c(69.5, 30.8, 15.4)) segments(2, y, 98, y, col = "#e6e5e2", lwd = 0.7)
}

open_pdf <- function(p) {
  ok <- try(cairo_pdf(p, width = W, height = H), silent = TRUE)
  if (inherits(ok, "try-error")) {
    message("（cairo_pdf が使えないので pdf() に落とします。フォント埋め込みを確認すること）")
    pdf(p, width = W, height = H, useDingbats = FALSE)
  }
}
open_png <- function(p, r) {
  ok <- try(png(p, width = W, height = H, units = "in", res = r, type = "cairo"),
            silent = TRUE)
  if (inherits(ok, "try-error")) png(p, width = W, height = H, units = "in", res = r)
}

if (tolower(tools::file_ext(OUT)) == "png") { open_png(OUT, 600); draw(); dev.off() } else {
  open_pdf(OUT); draw(); dev.off()
  if (PNG) { p <- sub("\\.[Pp][Dd][Ff]$", ".png", OUT)
             open_png(p, 200); draw(); dev.off(); cat("確認用:", p, "\n") }
}
cat("書き出し:", OUT, "\n")
cat("投稿前に pdffonts", OUT, "で emb が yes であることを確認すること。\n")
