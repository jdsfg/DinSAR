#!/bin/bash
#=============================================================================
# PowerChina-1 DInSAR 一键自检脚本
# 功能：清理旧环境 -> 重装 -> 冒烟验证 -> 真实数据快速回归
#=============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INSTALLER_DIR="$SCRIPT_DIR"
VERIFY_DIR="$SCRIPT_DIR/verify_dinsar_run"

# 可选参数：指定原始数据目录
INPUT_DATA_DIR="${1:-}"

log() { echo "[SELF-CHECK] $*"; }
fail() { echo "❌ $*"; exit 1; }

log "开始一键自检"

# 保护条件
[ -d "$INSTALLER_DIR" ] || fail "安装目录不存在: $INSTALLER_DIR"

# 1) 环境清理（仅清理 GMTSAR 相关）
log "1/5 清理旧环境"
find /usr/local/bin -maxdepth 1 -type l \( -name 'make_slc_*' -o -name 'p2p_*.csh' -o -name 'xcorr' -o -name 'snaphu' \) -delete 2>/dev/null || true
rm -rf /usr/local/GMTSAR 2>/dev/null || true
rm -rf "$HOME/.local/gmtsar"

# 2) 重装
log "2/5 执行安装"
cd "$INSTALLER_DIR"
chmod +x ./install.sh
./install.sh --user
source ~/.bashrc >/dev/null 2>&1 || true
export PATH="$HOME/.local/gmtsar/bin:$PATH"

# 3) 冒烟验证
log "3/5 冒烟验证"
which make_slc_dj1 >/dev/null || fail "make_slc_dj1 不在 PATH"
which xcorr >/dev/null || fail "xcorr 不在 PATH"
snaphu >/tmp/self_check_snaphu.out 2>&1 || true
gmt --version >/dev/null || fail "gmt 不可用"

# 4) 选择测试数据并回归
log "4/5 回归测试"
if [ -n "$INPUT_DATA_DIR" ]; then
  CANDIDATES=("$INPUT_DATA_DIR")
else
  CANDIDATES=(
    "$INSTALLER_DIR/test_project_11days/raw"
    "$INSTALLER_DIR/test_project_004509/raw"
    "$INSTALLER_DIR/test_project/raw"
    "$INSTALLER_DIR/test_project_v22/raw"
  )
fi

DATA_DIR=""
for d in "${CANDIDATES[@]}"; do
  [ -d "$d" ] || continue
  xml_n=$(find "$d" -maxdepth 1 -type f -name '*.xml' | wc -l)
  tif_n=$(find "$d" -maxdepth 1 -type f \( -name '*.tiff' -o -name '*.tif' \) | wc -l)
  if [ "$xml_n" -ge 2 ] && [ "$tif_n" -ge 2 ]; then
    DATA_DIR="$d"
    break
  fi
done

[ -n "$DATA_DIR" ] || fail "未找到可用测试数据目录（需至少 2 个 XML + 2 个 TIFF/TIF）。"
log "使用测试数据: $DATA_DIR"

rm -rf "$VERIFY_DIR"
mkdir -p "$VERIFY_DIR/raw"
cp "$DATA_DIR"/*.xml "$VERIFY_DIR/raw/"
cp "$DATA_DIR"/*.tiff "$VERIFY_DIR/raw/" 2>/dev/null || cp "$DATA_DIR"/*.tif "$VERIFY_DIR/raw/"

chmod +x "$INSTALLER_DIR/run_dinsar.sh"
(time "$INSTALLER_DIR/run_dinsar.sh" --mode A --auto "$VERIFY_DIR")

# 5) 验证成果
log "5/5 校验结果"
if [ -f "$VERIFY_DIR/results/los_displacement.grd" ] || [ -f "$VERIFY_DIR/results/displacement.jpg" ] || [ -f "$VERIFY_DIR/results/displacement.png" ]; then
  echo "✅ 一键自检通过：部署与回归处理均正常。"
  ls -lh "$VERIFY_DIR/results"
else
  fail "未发现关键成果文件（los_displacement.grd / displacement.jpg / displacement.png）。"
fi
