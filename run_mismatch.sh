#!/usr/bin/env bash
## ---------------------------------------------------------------------------
## run_mismatch.sh
##
## 中核の主張（切片）に対するペア崩し陰性対照を 4 コホートで回す。
##
##   1. diag_mismatch_slope.R × 4 コホート（**GSE40419 も含める**）
##   2. mismatch_summary.R で 1 枚の表にまとめ、プールした切片で判定する
##
## **GSE40419 を必ず入れること。** 条件 = バッチのコホートでは崩しても切片が
## 残るはずで、残す 3 コホートで消えることの対照になる。3 コホートだけ回すと
## 「消えた」が正しい挙動なのか検出力不足なのか区別できない。
##
## 【使い方】リポジトリのルート（R/ の隣）で
##   ./run_mismatch.sh
##   NSET=20 NMIS=20 N_CORES=12 ./run_mismatch.sh
##   NMIS=5 ./run_mismatch.sh        # まず傾向だけ見る（p の下限は 1/6 = 0.17）
##   DRY=1 ./run_mismatch.sh
##
## 【所要時間】
##   diag_slope.R（NSET=200、3 コホート）が約 14 分。ここは NSET=20 なので
##   1 コホート 1 反復あたりその 1/10 強。(NMIS+1) 反復 × 4 コホートで、
##   既定（NSET=20, NMIS=20）なら **30〜60 分**程度を見込む。
##   N_CORES を上げると患者方向に並列化される。
##
##   **NSET を論文の 200 に上げるのは、既定で傾向を確認してからにすること。**
##   10 倍かかる。判定に必要な精度は 20 で足りる（本物も同じ 20 で計算し直す）。
##
## 【環境変数】
##   REVISED_ROOT  既定 Revised
##   METADATA_DIR  既定 Revised/metadata
##   COHORTS       既定 4 コホート
##   NSET          既定 20    対照集合の組数（本物も同じ値で再計算される）
##   NMIS          既定 20    崩しの反復数。経験的 p の下限は 1/(NMIS+1)
##   N_CORES       既定 4
##   DRY           1 でコマンド表示のみ
## ---------------------------------------------------------------------------
set -u

REVISED_ROOT="${REVISED_ROOT:-Revised}"
METADATA_DIR="${METADATA_DIR:-${REVISED_ROOT}/metadata}"
COHORTS="${COHORTS:-GSE244679 GSE127165 GSE144269 GSE40419}"
NSET="${NSET:-20}"
NMIS="${NMIS:-20}"
N_CORES="${N_CORES:-4}"
DRY="${DRY:-0}"

LOGDIR="logs_mismatch"
SUMMARY="mismatch_report.txt"
mkdir -p "$LOGDIR"

hr() { printf '%.0s-' {1..74}; echo; }
FAILED=()

for CO in $COHORTS; do
  if [ ! -d "${REVISED_ROOT}/${CO}" ]; then
    echo "**${REVISED_ROOT}/${CO} がありません。REVISED_ROOT を確認してください。**"
    exit 1
  fi
done

if ! echo "$COHORTS" | grep -q GSE40419; then
  echo
  echo "**注意: COHORTS に GSE40419 が入っていません。**"
  echo "  条件 = バッチのコホートで崩しても切片が残ることが、"
  echo "  残す 3 コホートで消えることの対照になります。入れてください。"
  echo
fi

echo
hr
echo "切片に対するペア崩し陰性対照   NSET=${NSET}  NMIS=${NMIS}  N_CORES=${N_CORES}"
echo "経験的 p の下限 = 1/(NMIS+1) = $(awk -v n="$NMIS" 'BEGIN{printf "%.4f", 1/(n+1)}')"
hr

START=$(date +%s)
for CO in $COHORTS; do
  LOG="${LOGDIR}/mismatch_${CO}.log"
  echo "  → ${CO}"
  if [ "$DRY" = "1" ]; then
    echo "     COHORT=$CO KALLISTO_ROOT=${REVISED_ROOT}/${CO} METADATA_DIR=$METADATA_DIR NSET=$NSET NMIS=$NMIS N_CORES=$N_CORES Rscript diag_mismatch_slope.R"
    continue
  fi
  T0=$(date +%s)
  if env COHORT="$CO" \
         KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
         METADATA_DIR="$METADATA_DIR" \
         NSET="$NSET" NMIS="$NMIS" N_CORES="$N_CORES" \
         Rscript diag_mismatch_slope.R > "$LOG" 2>&1; then
    echo "     ok  ($(( $(date +%s) - T0 )) 秒)  $LOG"
  else
    echo "     **失敗**  $LOG"
    FAILED+=("$CO")
  fi
done

if [ "$DRY" = "1" ]; then echo; echo "DRY=1 なので要約は作りません。"; exit 0; fi

echo
hr
echo "まとめ（mismatch_summary.R）"
hr
{
  echo "切片のペア崩し陰性対照   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "NSET=${NSET}  NMIS=${NMIS}  所要 $(( $(date +%s) - START )) 秒"
  echo
  echo "############################################################"
  echo "# 本物の値（方針書 §0 (D)、NSET=200）"
  echo "#   GSE244679  切片 +0.0250 ± 0.0178 (z=1.4)"
  echo "#   GSE127165  切片 +0.0553 ± 0.0073 (z=7.6)"
  echo "#   GSE144269  切片 +0.0659 ± 0.0122 (z=5.4)"
  echo "#   プール     +0.0544 ± 0.0059 (z=9.2, Q p=0.163)"
  echo "#   GSE40419   切片 +0.7256 ± 0.0114（交絡。除外済み）"
  echo "#"
  echo "# **下の表の true 列は NSET=${NSET} で計算し直した値なので、"
  echo "#   上の NSET=200 の値とは小数第 3 位あたりが動く。それは正常。**"
  echo "#   比較は同じ表の中（true 対 崩し中位）で行うこと。"
  echo "############################################################"
  echo
  COHORTS_CSV=$(echo "$COHORTS" | tr ' ' ',')
  COHORTS="$COHORTS_CSV" Rscript mismatch_summary.R 2>&1
  echo
  echo "############################################################"
  echo "# 失敗した実行"
  echo "############################################################"
  if [ "${#FAILED[@]}" -eq 0 ]; then echo "  なし"; else printf '  %s\n' "${FAILED[@]}"; fi
} > "$SUMMARY" 2>&1

cat "$SUMMARY"
echo
echo "報告書: ${SUMMARY}"

if [ "${#FAILED[@]}" -ne 0 ]; then
  echo "**${#FAILED[@]} 件失敗しています。** 個別のログ: ${LOGDIR}/"
  exit 1
fi

cat <<'EOS'

結果の使い方

  プールした切片（(2) の節）が判定の中心である。

  ** 比（崩し中位 / 本物）が 0 付近、経験的 p が下限に張り付く **
     → 中核の主張に対する陰性対照が通った。§0 (B)(C) に無かった穴が埋まる。
       Methods の妥当性確認として書き、Limitations にも 1 文置く。
       **ここまで来たら Results を書き始めてよい。**

  ** 比が 0.5 を超える **
     → 条件づけとマッチングが効果を作っている可能性が残る。
       **Results を書き始めない。** 主張の立て方から見直す。
       次に見るのは (1) のコホート別で、どのコホートが崩れているか。

  ** GSE40419 の比が 1 付近 **
     → 想定どおり。条件 = バッチなので崩しても残る。
       3 コホートとの対比としてそのまま使える。
  ** GSE40419 の比も 0 付近だった場合 **
     → 手続きそのものを疑うこと。バッチ交絡があるコホートで消えるなら、
       この陰性対照は交絡を検出できていない。N_CORES=1 で再現を確認し、
       置換が本当に disease 側だけを入れ替えているかをログで確かめる。

EOS
