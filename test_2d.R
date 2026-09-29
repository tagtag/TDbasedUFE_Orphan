## ---------------------------------------------------------------------------
## test_2d.R
##
## external2（発現量 × 検出頻度の 2 次元マッチング）が 1 次元より良いことを、
## 疾患効果ゼロの合成データで確かめる。
##
##   1D      : 他患者の normal 平均を 20 分位
##   2D_glob : それに他患者での検出頻度 4 分位を掛ける（= external2）
##   2D_orph : 層の切り方を orphan の分布で決める版（参考。既定の走査には入れない）
##
## 疾患側は DS 倍してから床 0.02 で切るので、疾患効果はゼロのまま**検出だけ**が
## 押し上がる。orphan は検出確率 p_expr を (0.02, 0.40)、非 orphan を
## (0.02, 1.00) から引き、存在するときの発現量を base/p_expr にしてある。
## **平均発現量は p_expr に依らず同じで、違うのは検出頻度だけ**になる。
## これが「発現量だけ揃えても解けない」配置である。
##
## 使い方（リポジトリのルートで）
##   Rscript test_2d.R
##   DS=2.0 Rscript test_2d.R
##
## ---------------------------------------------------------------------------
## 【2026-09-23 訂正 2 件。**統計は 1 行も変えていない。**】
##
## (1) 07_interaction_permutation.R の読み込みに R/ のフォールバックが無く、
##     リポジトリのルートで走らせると落ちた。他の診断スクリプトに合わせた。
##
## (2) 「(名目 1.2)」の印字をやめた。**この p は較正されていない。**
##
##     t <- mean(obs) は 150 回の抽出を平均した量で、比較先の nul は 1 回ごとの
##     生の値の分布である。obs は orphan 集合を固定して対照 1 組を引くのに対し、
##     nul は対照 2 組を引くので、nul の方が構造的に広い。そこへ平均を当てるので
##     **この検定は保守側に偏る。**帰無での大きさは 0.05 でも 0.10 でもなく、
##     解析的には出ない。
##
##     旧版の「名目 1.2」（= 0.05 * 24）も、それを 2 倍した 2.4 も、どちらも
##     根拠がない。**期待人数は印字しない。**方式間の比較（1D 対 2D）は同じ
##     規則で計算しているので有効であり、そちらを見ること。
##
##     判定に使うのは **normal 側のゼロ率差**である。これは生の TPM から
##     計算するので表現にも検定の較正にも依存しない。1D で −0.108、
##     2D_glob で −0.008（DS=1.6、seed 11）。実データの −0.0087 → +0.00028 と
##     同じ向き・同じ桁の改善であり、external2 を採る根拠はこちらに置く。
## ---------------------------------------------------------------------------

transform_expression <- function(v, transform=c("rank","scale")) {
  transform <- match.arg(transform)
  if (transform=="rank") return(rank(v, ties.method="average"))
  s <- sd(v); if(!is.finite(s)||s==0) return(rep(NA_real_,length(v))); (v-mean(v))/s
}
## R/ 配下を先に試す。他の診断スクリプトと同じ流儀にした（2026-09-23）。
SRC <- if (file.exists("R/07_interaction_permutation.R"))
         "R/07_interaction_permutation.R" else "07_interaction_permutation.R"
if (!file.exists(SRC))
  stop("07_interaction_permutation.R が見つかりません。リポジトリのルートか ",
       "そのファイルのあるディレクトリで実行してください。", call. = FALSE)
src <- readLines(SRC)
eval(parse(text=paste(src[grep("^cell_coefs <- function",src):
                          (grep("^run_cohort_interaction <- function",src)-1L)],collapse="\n")))

## 層を「orphan の分布」で切る。全体の分位で切ると orphan が住む領域が
## 1 つの箱に潰れ、マッチングが効かない。
strata_on_orphan <- function(v, is_orphan, n) {
  br <- unique(quantile(v[is_orphan], probs = seq(0, 1, length.out = n + 1),
                        na.rm = TRUE, type = 7))
  if (length(br) < 2L) return(rep(1L, length(v)))
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(v, breaks = br, include.lowest = TRUE))
}
combine <- function(...) { l <- list(...); out <- l[[1]]
  for (j in seq_along(l)[-1]) out <- (out - 1L) * max(l[[j]]) + l[[j]]
  match(out, sort(unique(out))) }

set.seed(11)
NT <- 30000; NO <- 1500; NP <- 24
base <- rlnorm(NT, -3, 2.2)
is_o <- rep(FALSE, NT); is_o[sample(NT, NO)] <- TRUE
## 非 orphan も散発的なものを含む（orphan と重なる範囲を持つ）
p_expr <- ifelse(is_o, runif(NT, 0.02, 0.40), runif(NT, 0.02, 1.00))
amp <- base / p_expr

sim_lib <- function(ds) {
  on <- runif(NT) < p_expr
  v <- ifelse(on, amp * rlnorm(NT, 0, .5), 0) * ds
  v[v < 0.02] <- 0
  1e6 * v / max(sum(v), 1e-9)
}
DS <- as.numeric(Sys.getenv("DS","1.6"))
N <- replicate(NP, sim_lib(1.0)); D <- replicate(NP, sim_lib(DS))   # 疾患効果なし

one <- function(i, mode, tf = "scale", nperm = 150) {
  xn <- N[,i]; xd <- D[,i]; oth <- N[,-i,drop=FALSE]
  mu <- rowMeans(oth); df <- rowMeans(oth > 0)
  keep <- (xn>0)|(xd>0); xn<-xn[keep]; xd<-xd[keep]; orph<-is_o[keep]
  st <- switch(mode,
    "1D"      = make_strata(mu[keep], 20L),
    "2D_glob" = combine(make_strata(mu[keep],20L), make_strata(df[keep],4L)),
    "2D_orph" = combine(strata_on_orphan(mu[keep],orph,20L),
                        strata_on_orphan(df[keep],orph,5L)))
  set.seed(100+i)
  yN <- transform_expression(xn,tf); yD <- transform_expression(xd,tf)
  io <- which(orph); obs<-nul<-numeric(nperm); zo<-zc<-dzo<-dzc<-numeric(nperm); shrt<-0L
  for (r in seq_len(nperm)) {
    a <- draw_matched_sets(st,orph,1L); bc <- draw_matched_sets(st,orph,2L)
    shrt <- shrt + attr(a,"shortfall") + attr(bc,"shortfall")
    obs[r] <- cell_coefs(yN,yD,io,a[[1]])["b3"]
    nul[r] <- cell_coefs(yN,yD,bc[[1]],bc[[2]])["b3"]
    zo[r] <- mean(xn[io]==0); zc[r] <- mean(xn[a[[1]]]==0)
    dzo[r] <- mean(xd[io]==0); dzc[r] <- mean(xd[a[[1]]]==0)
  }
  t <- mean(obs)
  c(p=min((1+sum(nul>=t))/(nperm+1),(1+sum(nul<=t))/(nperm+1)),
    zo=mean(zo), zc=mean(zc), dzo=mean(dzo), dzc=mean(dzc),
    sdmad=sd(nul)/max(mad(nul),1e-300),
    K=length(unique(st)), sh=shrt)
}

cat(sprintf("DS = %.2f   患者 %d 人   転写産物 %d 本（orphan %d）  seed 11\n",
            DS, NP, NT, NO))
cat("** p < 0.05 の人数は方式間の比較にだけ使う。帰無での期待人数は出せない。\n")
cat("   t = mean(obs) を 1 回ごとの nul の分布に当てているので検定は保守側に\n")
cat("   偏る。判定は右端の『normal ゼロ率差』で行う（0 に近いほど良い）。\n\n")
for (tf in c("scale","rank")) for (m in c("1D","2D_glob")) {
  r <- t(sapply(1:NP, one, mode=m, tf=tf))
  cat(sprintf("%-5s %-8s p<0.05 %2d/%d  sd/mad %5.2f  **normal ゼロ率差 %+.3f**  非対称 %+.3f  層 %3d sh %d\n",
    tf, m, sum(r[,"p"]<0.05), NP, mean(r[,"sdmad"]),
    mean(r[,"zc"])-mean(r[,"zo"]),
    (mean(r[,"dzc"])-mean(r[,"dzo"])) - (mean(r[,"zc"])-mean(r[,"zo"])),
    round(mean(r[,"K"])), sum(r[,"sh"])))
}
cat("\n読み方\n")
cat("  normal ゼロ率差 … 疾患効果が存在しえない側での、対照集合と orphan 集合の\n")
cat("    ゼロ率の差。0 から離れているほど対照が揃っていない。**生の TPM から\n")
cat("    計算するので rank と scale で同一の値になる。**\n")
cat("  sh … shortfall。0 でなければ層のどこかで非 orphan が足りていない。\n")
cat("  p<0.05 の人数 … 1D と 2D の比較にだけ使う。上記のとおり較正されていない。\n")
