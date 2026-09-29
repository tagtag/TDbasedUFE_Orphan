## マッチング変数の候補を比較する。
##  pair_mean : 対平均 (n+d)/2            … 論文の方式
##  normal    : その患者の normal のみ
##  external  : 他の患者の normal の平均   … その患者の疾患状態にも雑音にも依存しない
transform_expression <- function(v,transform=c("rank","scale")){
  transform <- match.arg(transform)
  if(transform=="rank") return(rank(v,ties.method="average"))
  s<-sd(v); if(!is.finite(s)||s==0) return(rep(NA_real_,length(v))); (v-mean(v))/s }
src <- readLines("07_interaction_permutation.R")
eval(parse(text=paste(src[grep("^cell_coefs <- function",src):
   (grep("^run_cohort_interaction <- function",src)-1L)],collapse="\n")))

NT<-20000; NO<-500; NP<-12   # 転写産物, orphan, 患者数
tp <- function(x) 1e6*x/sum(x)

make_cohort <- function(mode, seed){
  set.seed(seed)
  base <- rlnorm(NT,0,2.5); ord <- order(base)
  o <- rep(FALSE,NT); o[ord[sample(8000,NO)]] <- TRUE
  N <- D <- matrix(0,NT,NP)
  for(p in 1:NP){
    pat <- base*rlnorm(NT,0,.35)
    n <- pat*rlnorm(NT,0,.3); d <- pat*rlnorm(NT,0,.3)
    if(mode=="composition") d[ord[(NT-199):NT]] <- d[ord[(NT-199):NT]]*3
    if(mode=="low_wide"){ low <- base<quantile(base,.4); d[low] <- d[low]*exp(.5) }
    if(mode=="orphan_only") d[o] <- d[o]*exp(.5)
    N[,p] <- tp(n); D[,p] <- tp(d) }
  list(N=N, D=D, o=o) }

run <- function(co, p, match_on, nperm=150, seed=1){
  xn <- co$N[,p]; xd <- co$D[,p]
  keep <- (xn>0)|(xd>0); xn<-xn[keep]; xd<-xd[keep]; orph<-co$o[keep]
  mv <- switch(match_on,
    pair_mean = (xn+xd)/2,
    normal    = xn,
    external  = rowMeans(co$N[keep,-p,drop=FALSE]))
  st <- make_strata(mv,20); set.seed(seed)
  yN <- transform_expression(xn,"scale"); yD <- transform_expression(xd,"scale")
  io <- which(orph); obs<-nul<-numeric(nperm)
  for(r in 1:nperm){ a<-draw_matched_sets(st,orph,1L)[[1]]; bc<-draw_matched_sets(st,orph,2L)
    obs[r]<-cell_coefs(yN,yD,io,a)["b3"]; nul[r]<-cell_coefs(yN,yD,bc[[1]],bc[[2]])["b3"] }
  t<-mean(obs); (1+sum(abs(nul-mean(nul))>=abs(t-mean(nul))))/(nperm+1) }

cat(sprintf("%-13s %12s %12s %12s\n","シナリオ","pair_mean","normal","external"))
cat(strrep("-",53),"\n")
for(m in c("null","composition","low_wide","orphan_only")){
  co <- make_cohort(m, 900)
  r <- sapply(c("pair_mean","normal","external"), function(mo)
        mean(sapply(1:NP, function(p) run(co,p,mo,150,900+p)) < .05))
  cat(sprintf("%-13s %12.2f %12.2f %12.2f\n", m, r[1],r[2],r[3])) }
cat("\n上3行は偽陽性率(低いほど良い)、最終行は検出力(高いほど良い)\n")
