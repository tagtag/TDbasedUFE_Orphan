#!/usr/bin/env bash
# =============================================================================
# run_verify.sh — 作業項目 8（rev.37 で修正した 3 本の再実行と印字の確認）
#
# 何をするか
#   1. diag_slope.R を 3 コホートで実行（NSET=200。論文 Table 2 の数値）
#   2. 浅いライブラリを除いた感度分析（該当 2 コホートのみ。§0 (L)）
#   3. de_logfc_<COHORT>_voom.csv が無ければ de_framework.R を実行
#   4. rho_compare.R を 3 コホートで実行（camera 同一集合の版を出すため 3 の後）
#   5. 全ログから確認すべき行だけ抜き出して verify_summary.txt にまとめる
#
# 何を確認するか（rev.36 で立てた宿題）
#   (a) 切片が d_zr_all = 0 への外挿かどうか   ← diag_slope.R が判定を印字
#   (b) 対照集合の抽出ばらつきの大きさ          ← 同上
#   (c) rho_compare.R の両側 p が方針書の表と一致するか（片側→両側の修正確認）
#   (d) 修正で Table 2 の数値が動いていないこと
#
# 使い方（リポジトリのルートで）
#   bash run_verify.sh                    # 全部走らせる
#   N_CORES=12 bash run_verify.sh         # 並列数を指定（強く推奨）
#   DRY=1 bash run_verify.sh              # コマンドを表示するだけ
#   COHORTS="GSE244679" bash run_verify.sh        # 1 コホートだけ
#   EXTRAS=1 bash run_verify.sh           # DESeq2 と MIN_COUNT 掃引も走らせる
#   FORCE_DE=1 bash run_verify.sh         # de_framework.R を必ず走らせ直す
#
# 所要時間の目安
#   diag_slope.R は NSET=200 なので患者 1 人あたり 200 回対照集合を引く。
#   GSE144269（70 人）は 14,000 回になるので N_CORES を上げること。
#   de_framework.R は abundance.tsv の読み込みが支配的だが、
#   counts_<COHORT>.rds にキャッシュされるので 2 回目以降は速い。
#
# 出力
#   logs/<日時>/  に各実行のログ
#   logs/<日時>/verify_summary.txt  に確認用の抜粋（これを読めば済む）
# =============================================================================

set -uo pipefail

# ------------------------------------------------------------------ 設定
REVISED_ROOT="${REVISED_ROOT:-Revised}"
METADATA_DIR="${METADATA_DIR:-${REVISED_ROOT}/metadata}"
COHORTS="${COHORTS:-GSE244679 GSE127165 GSE144269}"
N_CORES="${N_CORES:-1}"
NSET_SLOPE="${NSET_SLOPE:-200}"   # 論文 Table 2 は 200。既定の 20 ではない
NSET_RHO="${NSET_RHO:-20}"        # 論文 Table 5 は 20
DRY="${DRY:-0}"
EXTRAS="${EXTRAS:-0}"
FORCE_DE="${FORCE_DE:-0}"

# 浅いライブラリ（§0 (A)(L)）。GSE244679 は該当なし
declare -A SHALLOW=(
  [GSE127165]="SRR8631680"
  [GSE144269]="SRR10969201"
)

STAMP="$(date +%Y%m%d-%H%M%S)"
LOGDIR="logs/${STAMP}"
SUMMARY="${LOGDIR}/verify_summary.txt"

# ------------------------------------------------------------------ 前提確認
die() { echo "エラー: $*" >&2; exit 1; }

command -v Rscript >/dev/null 2>&1 || die "Rscript が見つかりません。"

for f in R/config.R R/00_functions.R R/07_interaction_permutation.R; do
  [ -f "$f" ] || die "$f がありません。リポジトリのルートで実行してください。
  （diag_slope.R / de_framework.R / rho_compare.R はこの 3 本を
    source() と readLines() で参照します。）"
done

for f in diag_slope.R de_framework.R rho_compare.R; do
  [ -f "$f" ] || die "$f がありません。rev.37 で修正した版を置いてください。"
done

[ -d "$REVISED_ROOT" ] || die "$REVISED_ROOT がありません。REVISED_ROOT で指定してください。"
[ -d "$METADATA_DIR" ] || die "$METADATA_DIR がありません。METADATA_DIR で指定してください。"

for CO in $COHORTS; do
  [ -d "${REVISED_ROOT}/${CO}" ] || die "${REVISED_ROOT}/${CO} がありません。"
done

mkdir -p "$LOGDIR"

# ------------------------------------------------------------------ 実行補助
FAILED=()

hr() { printf '%.0s=' {1..78}; echo; }

# run <ログ名> <環境変数の並び...> -- <スクリプト>
run() {
  local name="$1"; shift
  local log="${LOGDIR}/${name}.log"
  echo
  hr
  echo "▶ ${name}"
  echo "  $*"
  if [ "$DRY" = "1" ]; then
    echo "  （DRY=1 なので実行しません）"
    return 0
  fi
  local t0 t1
  t0=$(date +%s)
  if env "$@" > "$log" 2>&1; then
    t1=$(date +%s)
    echo "  完了（$((t1 - t0)) 秒）→ ${log}"
  else
    t1=$(date +%s)
    echo "  **失敗**（$((t1 - t0)) 秒）→ ${log}"
    echo "  ログの末尾:"
    tail -n 15 "$log" | sed 's/^/    /'
    FAILED+=("$name")
    return 1
  fi
}

# ------------------------------------------------------------------ 1. 傾き・切片
echo
hr
echo "1. diag_slope.R（論文 Table 2 の量。NSET=${NSET_SLOPE}）"
hr

for CO in $COHORTS; do
  run "slope_${CO}" \
    COHORT="$CO" \
    KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
    METADATA_DIR="$METADATA_DIR" \
    NSET="$NSET_SLOPE" \
    N_CORES="$N_CORES" \
    Rscript diag_slope.R
done

# ------------------------------------------------------------------ 2. 浅いライブラリ除外
echo
hr
echo "2. diag_slope.R（浅いライブラリ除外。§0 (L)）"
hr

for CO in $COHORTS; do
  LIBS="${SHALLOW[$CO]:-}"
  if [ -z "$LIBS" ]; then
    echo "  ${CO}: 該当ライブラリなし（中位の 1/5 未満のものがない）。飛ばします。"
    continue
  fi
  run "slope_${CO}_excl" \
    COHORT="$CO" \
    KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
    METADATA_DIR="$METADATA_DIR" \
    NSET="$NSET_SLOPE" \
    N_CORES="$N_CORES" \
    EXCLUDE_LIBS="$LIBS" \
    Rscript diag_slope.R
done

# ------------------------------------------------------------------ 3. DE 枠組み
echo
hr
echo "3. de_framework.R（rho_compare.R の camera 同一集合版に必要）"
hr

for CO in $COHORTS; do
  DEFILE="de_logfc_${CO}_voom.csv"
  if [ -f "$DEFILE" ] && [ "$FORCE_DE" != "1" ]; then
    echo "  ${CO}: ${DEFILE} が既にあります。飛ばします（FORCE_DE=1 で再実行）。"
    continue
  fi
  run "de_${CO}_voom" \
    COHORT="$CO" \
    KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
    METADATA_DIR="$METADATA_DIR" \
    ENGINE=voom \
    NSET="$NSET_RHO" \
    N_CORES="$N_CORES" \
    Rscript de_framework.R
done

# ------------------------------------------------------------------ 4. 患者単位検定
echo
hr
echo "4. rho_compare.R（両側化の確認。論文 Table 5 の量。NSET=${NSET_RHO}）"
hr

for CO in $COHORTS; do
  DEFILE="de_logfc_${CO}_voom.csv"
  if [ ! -f "$DEFILE" ] && [ "$DRY" != "1" ]; then
    echo "  ${CO}: **警告** ${DEFILE} が無いので camera 同一集合の版は出ません。"
    echo "         全集合の版だけになります（方針書の主表は camera 同一集合）。"
  fi
  run "rho_${CO}" \
    COHORT="$CO" \
    KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
    METADATA_DIR="$METADATA_DIR" \
    NSET="$NSET_RHO" \
    N_CORES="$N_CORES" \
    Rscript rho_compare.R
done

# ------------------------------------------------------------------ 5. 追加（任意）
if [ "$EXTRAS" = "1" ]; then
  echo
  hr
  echo "5. 追加の感度分析（GSE244679 のみ。§0 (J)(K)）"
  hr

  if echo "$COHORTS" | grep -q GSE244679; then
    run "de_GSE244679_deseq2" \
      COHORT=GSE244679 \
      KALLISTO_ROOT="${REVISED_ROOT}/GSE244679" \
      METADATA_DIR="$METADATA_DIR" \
      ENGINE=deseq2 \
      NSET="$NSET_RHO" \
      N_CORES="$N_CORES" \
      Rscript de_framework.R

    for MC in 5 1; do
      run "de_GSE244679_voom_mc${MC}" \
        COHORT=GSE244679 \
        KALLISTO_ROOT="${REVISED_ROOT}/GSE244679" \
        METADATA_DIR="$METADATA_DIR" \
        ENGINE=voom \
        MIN_COUNT="$MC" \
        NSET="$NSET_RHO" \
        N_CORES="$N_CORES" \
        Rscript de_framework.R
    done
  else
    echo "  COHORTS に GSE244679 が入っていないので飛ばします。"
  fi
fi

# ------------------------------------------------------------------ 6. 抜粋
if [ "$DRY" = "1" ]; then
  echo
  echo "DRY=1 なので要約は作りません。"
  exit 0
fi

pick() {  # pick <ログ> <パターン...>
  local log="$1"; shift
  [ -f "$log" ] || { echo "    （ログなし）"; return; }
  grep -aE "$1" "$log" | sed 's/^/    /' || echo "    （該当行なし）"
}

{
  echo "作業項目 8 の確認用抜粋   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "NSET: diag_slope=${NSET_SLOPE}  rho_compare=${NSET_RHO}  N_CORES=${N_CORES}"
  echo

  echo "############################################################"
  echo "# (d) Table 2 の数値が動いていないか"
  echo "#     方針書 rev.37 §0 (D) の表"
  echo "#       GSE244679  orphan 4.081±0.487  対照 1.155±0.027  増幅率 3.53"
  echo "#                  傾きの差 2.926±0.483  切片 +0.0250±0.0178"
  echo "#       GSE127165  orphan 2.385±0.145  対照 1.565±0.056  増幅率 1.52"
  echo "#                  傾きの差 0.820±0.133  切片 +0.0553±0.0073"
  echo "#       GSE144269  orphan 1.366±0.103  対照 1.076±0.064  増幅率 1.27"
  echo "#                  傾きの差 0.290±0.077  切片 +0.0659±0.0122"
  echo "############################################################"
  for CO in $COHORTS; do
    echo
    echo "## ${CO}"
    pick "${LOGDIR}/slope_${CO}.log" '^  (orphan   傾き|対照     傾き|傾き \(orphan|切片  |増幅率 =|解析対象|shortfall 0|注意: shortfall)'
  done

  echo
  echo "############################################################"
  echo "# (a) 切片は外挿か（d_zr_all = 0 がデータ範囲の内側か）"
  echo "#     外側なら Limitations に一文。方針書 §0 (D) の 1 と (F)"
  echo "############################################################"
  for CO in $COHORTS; do
    echo
    echo "## ${CO}"
    pick "${LOGDIR}/slope_${CO}.log" 'd_zr_all の範囲|切片を取る点|0 に最も近い患者|外挿である'
  done

  echo
  echo "############################################################"
  echo "# (b) 報告する ± に含まれない抽出ばらつきの大きさ"
  echo "#     Methods に一文。方針書 §0 (D) の定義小節と (F)"
  echo "############################################################"
  for CO in $COHORTS; do
    echo
    echo "## ${CO}"
    pick "${LOGDIR}/slope_${CO}.log" '対照集合の抽出による|切片の回帰 SE|無視できない大きさ'
  done

  echo
  echo "############################################################"
  echo "# (c) rho_compare.R の両側 p が方針書 §0 (I) の表と一致するか"
  echo "#     camera 同一集合（主表）"
  echo "#       GSE244679  Δ̄ -0.1169±0.0304  t p 8.2e-04  Wilcoxon 1.4e-03  Δ<0 18/24"
  echo "#       GSE127165  Δ̄ -0.0194±0.0078  t p 1.6e-02  Wilcoxon 2.5e-02  Δ<0 35/57"
  echo "#       GSE144269  Δ̄ -0.0455±0.0065  t p 1.4e-09  Wilcoxon 3.8e-08  Δ<0 54/70"
  echo "#     全集合（Supplementary）"
  echo "#       GSE244679  Δ̄ -0.0886±0.0240  t p 1.2e-03  Δ<0 17/24"
  echo "#       GSE127165  Δ̄ -0.0123±0.0027  t p 2.1e-05  Δ<0 43/57"
  echo "#       GSE144269  Δ̄ -0.0224±0.0024  t p 1.5e-13  Δ<0 61/70"
  echo "#"
  echo "#     片側のままなら t p がほぼ半分になる（8.2e-04 → 4.1e-04）。"
  echo "#     その場合は修正版が使われていない。"
  echo "############################################################"
  for CO in $COHORTS; do
    echo
    echo "## ${CO}"
    pick "${LOGDIR}/rho_${CO}.log" '^【|Delta = orphan|95% CI|1 標本 t p|Wilcoxon 符号付き順位 p|Delta < 0 の患者|rho_bar|orphan        [-+]|マッチ対照    中位'
  done

  echo
  echo "############################################################"
  echo "# (L) 浅いライブラリ除外（§0 (L)）"
  echo "#     GSE127165  増幅率 1.52→1.52  傾きの差 0.820→0.869  切片 +0.0553→+0.0549"
  echo "#     GSE144269  増幅率 1.27→1.44  傾きの差 0.290→0.472  切片 +0.0659→+0.0742"
  echo "############################################################"
  for CO in $COHORTS; do
    [ -n "${SHALLOW[$CO]:-}" ] || continue
    echo
    echo "## ${CO}（${SHALLOW[$CO]} 除外）"
    pick "${LOGDIR}/slope_${CO}_excl.log" '^  (傾き \(orphan|切片  |増幅率 =)|EXCLUDE_LIBS により'
  done

  if [ "$EXTRAS" = "1" ]; then
    echo
    echo "############################################################"
    echo "# 追加（§0 (J)(K)。GSE244679 のみ）"
    echo "#   DESeq2      log2FC 差 -0.389  camera p 0.185"
    echo "#   MIN_COUNT=5 通過 89.3%  log2FC 差 -0.297  camera p 0.140"
    echo "#   MIN_COUNT=1 通過 93.7%  log2FC 差 -0.265  camera p 0.150"
    echo "############################################################"
    for n in de_GSE244679_deseq2 de_GSE244679_voom_mc5 de_GSE244679_voom_mc1; do
      echo
      echo "## ${n}"
      pick "${LOGDIR}/${n}.log" 'orphan      +[0-9]+ →|log2FC 平均差|camera p  中位|camera Direction|推定した転写産物間相関'
    done
  fi

  echo
  echo "############################################################"
  echo "# 失敗した実行"
  echo "############################################################"
  if [ "${#FAILED[@]}" -eq 0 ]; then
    echo "  なし"
  else
    printf '  %s\n' "${FAILED[@]}"
  fi
} > "$SUMMARY" 2>&1

echo
hr
echo "要約: ${SUMMARY}"
hr
cat "$SUMMARY"

echo
if [ "${#FAILED[@]}" -ne 0 ]; then
  echo "**${#FAILED[@]} 件失敗しています。** 個別のログを見てください: ${LOGDIR}/"
  exit 1
fi

cat <<'EOS'

次にやること
  1. 上の (a) で「外側」と出ていたら、方針書 §0 (F) の Limitations に一文足す。
     「切片は観測範囲外への外挿である」と自分で書く。査読者に先に言われない。
  2. (b) の比が 20% を超えていたら、NSET を増やすか誤差に合成する。
     20% 未満なら Methods に「報告した標準誤差に含まれない」と一文書いて終わり。
  3. (c) で t p が表と一致すれば両側化の修正は成功。半分の値なら古い版が
     使われている（which Rscript と rho_compare.R の中身を確認）。
  4. (d) と (L) が表と一致すれば、修正で数値は動いていない。
     動いていたら原因を特定するまで本文を書き始めないこと。

EOS
