#!/usr/bin/env bash
## ---------------------------------------------------------------------------
## run_batch_audit.sh
##
## バッチ交絡まわりの宿題を一度に片付ける。
##
##   1. diag_batch.R      4 コホート全部（新 Table 1 の「バッチ分離」「技術指標」列）
##   2. batch_extra.R     4 コホート全部（浅いライブラリ全件 ＋ 陰性対照の両方向
##                        ＋ 本物のペアとの比較）
##   3. diag_slope.R      **GSE40419 でも回す**（除外したコホートで
##                        マッチ対照がどう動いたかを定量する）
##   4. 要約を batch_audit_summary.txt に抜粋
##
## 3 が要点である。GSE40419 は条件 = バッチなので、そこで出る切片は
## **純粋な技術交絡が作る切片**である。残す 3 コホートの +0.0544 と
## 並べれば、「この統計量が交絡でどこまで動くか」を数字で示せる。
## 査読者に先に指摘される前に自分で出す。
##
## 【使い方】リポジトリのルート（R/ の隣）で
##   ./run_batch_audit.sh
##   REVISED_ROOT=Revised N_CORES=12 NSET=200 ./run_batch_audit.sh
##   COHORTS="GSE40419" ./run_batch_audit.sh        # 1 コホートだけ
##   DRY=1 ./run_batch_audit.sh                      # コマンドを表示するだけ
##
## 【環境変数】
##   REVISED_ROOT  既定 Revised            kallisto 出力の親
##   METADATA_DIR  既定 Revised/metadata   サンプルシート
##   COHORTS       既定 4 コホート
##   N_CORES       既定 4
##   NSET          既定 200（diag_slope.R。論文の値と揃える）
##   NREP          既定 200（batch_extra.R の陰性対照の反復）
##   SKIP_SLOPE    1 で 3 を飛ばす（時間がかかるので）
## ---------------------------------------------------------------------------
set -u

REVISED_ROOT="${REVISED_ROOT:-Revised}"
METADATA_DIR="${METADATA_DIR:-${REVISED_ROOT}/metadata}"
COHORTS="${COHORTS:-GSE244679 GSE127165 GSE144269 GSE40419}"
N_CORES="${N_CORES:-4}"
NSET="${NSET:-200}"
NREP="${NREP:-200}"
DRY="${DRY:-0}"
SKIP_SLOPE="${SKIP_SLOPE:-0}"

LOGDIR="logs_batch_audit"
SUMMARY="batch_audit_summary.txt"
mkdir -p "$LOGDIR"

hr() { printf '%.0s-' {1..74}; echo; }
FAILED=()

run() {  # run <ログ名> <環境変数...> <コマンド...>
  local name="$1"; shift
  local log="${LOGDIR}/${name}.log"
  echo "  → ${name}"
  if [ "$DRY" = "1" ]; then printf '     %s\n' "$*"; return 0; fi
  if env "$@" > "$log" 2>&1; then
    echo "     ok  ($log)"
  else
    echo "     **失敗**  ($log)"
    FAILED+=("$name")
  fi
}

for CO in $COHORTS; do
  if [ ! -d "${REVISED_ROOT}/${CO}" ]; then
    echo "**${REVISED_ROOT}/${CO} がありません。REVISED_ROOT を確認してください。**"
    exit 1
  fi
done

echo
hr
echo "1. diag_batch.R（accession 分離・技術指標・浅いライブラリ・陰性対照）"
hr
for CO in $COHORTS; do
  run "batch_${CO}" \
    COHORT="$CO" \
    KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
    METADATA_DIR="$METADATA_DIR" \
    Rscript diag_batch.R
done

echo
hr
echo "2. batch_extra.R（浅いライブラリ全件・陰性対照の両方向・本物のペアとの比較）"
hr
for CO in $COHORTS; do
  run "extra_${CO}" \
    COHORT="$CO" \
    KALLISTO_ROOT="${REVISED_ROOT}/${CO}" \
    METADATA_DIR="$METADATA_DIR" \
    NREP="$NREP" \
    Rscript batch_extra.R
done

if [ "$SKIP_SLOPE" != "1" ]; then
  echo
  hr
  echo "3. diag_slope.R を GSE40419 で実行（交絡が作る切片の大きさを測る）"
  hr
  if echo "$COHORTS" | grep -q GSE40419; then
    run "slope_GSE40419" \
      COHORT=GSE40419 \
      KALLISTO_ROOT="${REVISED_ROOT}/GSE40419" \
      METADATA_DIR="$METADATA_DIR" \
      NSET="$NSET" \
      N_CORES="$N_CORES" \
      Rscript diag_slope.R
  else
    echo "  COHORTS に GSE40419 が無いので飛ばします。"
  fi
fi

if [ "$DRY" = "1" ]; then
  echo
  echo "DRY=1 なので要約は作りません。"
  exit 0
fi

pick() {  # pick <ログ> <パターン>
  local log="$1"; shift
  [ -f "$log" ] || { echo "    （ログなし）"; return; }
  grep -aE "$1" "$log" | sed 's/^/    /' || echo "    （該当行なし）"
}

{
  echo "バッチ交絡の監査   $(date '+%Y-%m-%d %H:%M:%S')"
  echo "NSET=${NSET}  NREP=${NREP}  N_CORES=${N_CORES}"
  echo
  echo "############################################################"
  echo "# (1) accession の分離と技術指標   → 新 Table 1 の 2 列"
  echo "#  既知（GSE40419）: normal 160121-160196 / disease 164550-164634"
  echo "#    深さ 2.398e7 / 4.452e7   マップ率 60.0% / 76.5%"
  echo "#    ゼロ率 全体 0.3237 / 0.3632   orphan 0.4626 / 0.8767"
  echo "############################################################"
  for CO in $COHORTS; do
    echo; echo "## ${CO}"
    pick "${LOGDIR}/batch_${CO}.log" '^  (normal   |disease  |disease の範囲|normal の範囲|→ |深さ|マップ率|ゼロ率)'
  done

  echo
  echo "############################################################"
  echo "# (2) 浅いライブラリの条件別内訳"
  echo "#  方針書 R3-1 の英文「seven libraries with fewer than 4e6"
  echo "#  reads are all on the normal side」は誤り。"
  echo "#  ERR164584 が disease で 3.86e6。正しい言い方は"
  echo "#  「最も浅い N 本がすべて normal 側」。下の → 行を使うこと。"
  echo "############################################################"
  for CO in $COHORTS; do
    echo; echo "## ${CO}"
    pick "${LOGDIR}/extra_${CO}.log" '【|条件別:|→ 最も浅い'
  done

  echo
  echo "############################################################"
  echo "# (3) 陰性対照 ― 本物のペア vs 同一条件のランダムな組"
  echo "#  本物の中位が帰無分布に埋もれるか、sd 比が 1 付近かを見る。"
  echo "#  **GSE40419 は条件 = バッチなので、本物が大きく出ても"
  echo "#    それは疾患効果ではない。** 他 3 コホートとの対比に使う。"
  echo "############################################################"
  for CO in $COHORTS; do
    echo; echo "## ${CO}"
    pick "${LOGDIR}/extra_${CO}.log" '本物のペア|どうし \(|経験的 p|sd 比'
  done

  if [ "$SKIP_SLOPE" != "1" ] && echo "$COHORTS" | grep -q GSE40419; then
    echo
    echo "############################################################"
    echo "# (4) GSE40419 の傾きと切片 ― 交絡が作る値の大きさ"
    echo "#  残す 3 コホートの実測（方針書 §0 (D)）"
    echo "#    GSE244679  傾きの差 2.926±0.483   切片 +0.0250±0.0178"
    echo "#    GSE127165  傾きの差 0.820±0.133   切片 +0.0553±0.0073"
    echo "#    GSE144269  傾きの差 0.290±0.077   切片 +0.0659±0.0122"
    echo "#    プール切片 +0.0544±0.0059 (z=9.2, Q p=0.163)"
    echo "#"
    echo "#  GSE40419 の切片がこれらと同程度かそれ以上なら、"
    echo "#  **この統計量は純粋な技術交絡だけでその値を作れる**ことになる。"
    echo "#  Limitations に必ず書く。隠すと査読で崩される。"
    echo "############################################################"
    echo
    pick "${LOGDIR}/slope_GSE40419.log" '^  (orphan   傾き|対照     傾き|傾き \(orphan|切片  |増幅率 =|解析対象|shortfall|注意: shortfall)|ゼロ率の変化|全転写産物  |orphan      |マッチ対照  |差 \(orphan'
  fi

  echo
  echo "############################################################"
  echo "# 失敗した実行"
  echo "############################################################"
  if [ "${#FAILED[@]}" -eq 0 ]; then echo "  なし"; else printf '  %s\n' "${FAILED[@]}"; fi
} > "$SUMMARY" 2>&1

echo
hr
echo "要約: ${SUMMARY}"
hr
cat "$SUMMARY"

echo
if [ "${#FAILED[@]}" -ne 0 ]; then
  echo "**${#FAILED[@]} 件失敗しています。** 個別のログ: ${LOGDIR}/"
  exit 1
fi

cat <<'EOS'

次にやること
  1. (1) の 4 コホート分を新 Table 1 の「バッチ分離」「技術指標の条件差」列に写す。
  2. (2) の「→ 最も浅い N 本がすべて ... 側」を response letter の該当文に使う。
     現行の「seven ... are all on the normal side」は誤りなので必ず差し替える。
  3. (3) で残す 3 コホートの sd 比が 1 付近なら、患者ごとのばらつきは個体差で
     説明できるということ。§0 (B)(C) の結論と整合する。
  4. (4) の GSE40419 の切片が +0.05 前後かそれ以上なら、**中核の主張と同じ
     大きさの値が技術交絡だけで出る**ことになる。その場合は
       - Limitations に明記する
       - 残す 3 コホートに accession の分離が無いこと（§0 (A)）を
         Methods で明示的に守りの論拠として書く
       - GSE40419 を「除外した」ではなく「交絡の実例」として提示する
     の 3 点を必ずやる。査読者に先に指摘されると主張の中核が崩れる。

EOS
