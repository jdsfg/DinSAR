#!/bin/bash
#=============================================================================
# PowerChina-1 DInSAR 系统部署脚本 v2.1.2
# 
# 合并更新版: 包含 GCC 11+ 修复 + ALOS extern 修复 + run_dinsar.sh 安装
#
# 作者: PowerChina-1 DInSAR Team
# 日期: 2026-02-01
#=============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE="${SCRIPT_DIR}/deploy_$(date +%Y%m%d_%H%M%S).log"
GMTSAR_REPO="https://github.com/DONGYUSEN/gmtsar.git"
GMTSAR_BRANCH="master"
INSTALL_PREFIX="/usr/local/GMTSAR"
USER_INSTALL=false
OFFLINE_MODE=false
CHECK_ONLY=false
UNINSTALL=false
VERBOSE=false

# 颜色
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

log() { echo -e "${GREEN}[INFO]${NC} $1" | tee -a "$LOG_FILE"; }
error() { echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"; exit 1; }

check_dependencies() {
    log "检查系统依赖..."
    
    local missing=()
    local commands=("gcc" "g++" "make" "gmt" "autoconf")
    
    for cmd in "${commands[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        error "缺少必要命令: ${missing[*]}\n请安装: sudo apt-get install build-essential gmt gmt-gshhg autoconf libtiff-dev libhdf5-dev"
    fi
    
    log "依赖检查通过"
}

# 参数解析
for arg in "$@"; do
    case $arg in
        --user) USER_INSTALL=true; INSTALL_PREFIX="$HOME/.local/gmtsar" ;;
        --offline) OFFLINE_MODE=true ;;
    esac
done

fix_source_code() {
    log "应用 GCC 11+ 兼容性修补..."
    
    local src_dir="$GMTSAR_SRC"
    
    # --- Fix 1: config.mk -fcommon ---
    # 此步骤已移除，通过 CFLAGS 传递给 configure 脚本处理


    # --- Fix 2: gmtsar/sio_struct.c GMT_LONG typedef conflict ---
    local sio_c="$src_dir/gmtsar/sio_struct.c"
    if [ -f "$sio_c" ]; then
        # 移除可能冲突的 typedef long GMT_LONG;
        if grep -q "typedef long GMT_LONG;" "$sio_c"; then
             sed -i 's/typedef long GMT_LONG;/\/* typedef long GMT_LONG; *\//g' "$sio_c"
             log "已注释 gmtsar/sio_struct.c 中的 typedef GMT_LONG"
        fi
    fi

    # --- Fix 3: ALOS_preproc image_sio.h extern vars ---
    # 这是一个关键修复，因为多个 .c 文件包含了这个头文件，导致多重定义
    local alos_inc="$src_dir/preproc/ALOS_preproc/include/image_sio.h"
    local alos_lib="$src_dir/preproc/ALOS_preproc/lib_src"
    
    if [ -f "$alos_inc" ]; then
        if ! grep -q "extern int verbose" "$alos_inc"; then
            log "修补 image_sio.h 为 extern 声明..."
            sed -i 's/^int verbose;/extern int verbose;/g' "$alos_inc"
            sed -i 's/^int debug;/extern int debug;/g' "$alos_inc"
            sed -i 's/^int roi;/extern int roi;/g' "$alos_inc"
            sed -i 's/^int swap;/extern int swap;/g' "$alos_inc"
            sed -i 's/^int quad_pol;/extern int quad_pol;/g' "$alos_inc"
            sed -i 's/^int ALOS_format;/extern int ALOS_format;/g' "$alos_inc"
            sed -i 's/^int force_slope;/extern int force_slope;/g' "$alos_inc"
            sed -i 's/^int dopp;/extern int dopp;/g' "$alos_inc"
            sed -i 's/^int quiet_flag;/extern int quiet_flag;/g' "$alos_inc"
            sed -i 's/^int SAR_mode;/extern int SAR_mode;/g' "$alos_inc"
            sed -i 's/^int prefix_off;/extern int prefix_off;/g' "$alos_inc"
            
            sed -i 's/^double forced_slope;/extern double forced_slope;/g' "$alos_inc"
            sed -i 's/^double tbias;/extern double tbias;/g' "$alos_inc"
            sed -i 's/^double rbias;/extern double rbias;/g' "$alos_inc"
            sed -i 's/^double slc_fact;/extern double slc_fact;/g' "$alos_inc"
        else
            log "image_sio.h 已是 extern 声明"
        fi
    fi

    # --- Fix 4: Create globals_alos.c ---
    local alos_globals="$alos_lib/globals_alos.c"
    if [ ! -f "$alos_globals" ]; then
        log "创建 globals_alos.c..."
        cat > "$alos_globals" <<EOF
/*
 * globals_alos.c
 * Explicit definition of global variables declared extern in image_sio.h
 */
#include "image_sio.h"

int verbose = 0;
int debug = 0;
int roi = 0;
int swap = 0;
int quad_pol = 0;
int ALOS_format = 0;
int force_slope = 0;
int dopp = 0;
int quiet_flag = 0;
int SAR_mode = 0;
int prefix_off = 0;

double forced_slope = 0.0;
double tbias = 0.0;
double rbias = 0.0;
double slc_fact = 1.0;
EOF
    fi

    # --- Fix 5: Update ALOS Makefile to include globals_alos.o ---
    local alos_make="$alos_lib/Makefile"
    if [ -f "$alos_make" ]; then
        if ! grep -q "globals_alos" "$alos_make"; then
            log "更新 ALOS/lib_src/Makefile..."
            # 更健壮的替换：在 SRCS 行后添加，而不是替换
            if grep -q "^SRCS" "$alos_make"; then
                sed -i '/^SRCS/s/$/ globals_alos.c/' "$alos_make"
            else
                # 如果没有 SRCS 行，在 Makefile 开头添加
                sed -i '1i SRCS = globals_alos.c utils.c' "$alos_make"
            fi
        else
            log "ALOS Makefile 已包含 globals_alos.c"
        fi
        
        # 确保 Makefile 使用 Tab 缩进
        sed -i 's/^        /\t/g' "$alos_make"
    fi
}

main() {
    log "开始部署 v2.1.2 (合并更新版)..."
    
    # 0. 检查依赖
    check_dependencies
    
    # 1. 设置源码目录
    GMTSAR_SRC="$SCRIPT_DIR/gmtsar-master"
    if [ ! -d "$GMTSAR_SRC" ]; then
        error "找不到源码目录: $GMTSAR_SRC"
    fi

    # 2. 应用修复
    fix_source_code

    # 3. 编译
    log "开始编译..."
    cd "$GMTSAR_SRC"
    
    # 清理旧的编译结果，确保我们的修复生效
    make clean > /dev/null 2>&1
    
    # 配置 (通过 CFLAGS 添加 -fcommon)
    CFLAGS="-fcommon -O2" ./configure --prefix="$INSTALL_PREFIX" >> "$LOG_FILE" 2>&1
    
    # 编译
    make >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        error "编译失败，请查看日志: $LOG_FILE"
    fi
    
    # 4. 安装
    log "安装中..."
    make install >> "$LOG_FILE" 2>&1
    
    # 5. 安装 run_dinsar.sh
    log "安装数据处理脚本..."
    if [ -f "$SCRIPT_DIR/run_dinsar.sh" ]; then
        mkdir -p "$INSTALL_PREFIX/bin"
        cp "$SCRIPT_DIR/run_dinsar.sh" "$INSTALL_PREFIX/bin/"
        chmod +x "$INSTALL_PREFIX/bin/run_dinsar.sh"
        log "run_dinsar.sh 已安装到 $INSTALL_PREFIX/bin/"
    else
        log "警告: 未找到 run_dinsar.sh"
    fi
    
    # 6. 自动配置 PATH
    log "配置环境变量..."
    local shell_rc="$HOME/.bashrc"
    if ! grep -q "$INSTALL_PREFIX/bin" "$shell_rc" 2>/dev/null; then
        echo "" >> "$shell_rc"
        echo "# PowerChina-1 DInSAR System (added by install.sh)" >> "$shell_rc"
        echo "export PATH=\"$INSTALL_PREFIX/bin:\$PATH\"" >> "$shell_rc"
        log "PATH 已添加到 ~/.bashrc"
    else
        log "PATH 已配置"
    fi
    
    # 7. 验证安装
    log "验证安装..."
    export PATH="$INSTALL_PREFIX/bin:$PATH"
    
    local verify_ok=true
    for cmd in "make_slc_dj1" "xcorr" "phasediff" "snaphu"; do
        if ! command -v "$cmd" &>/dev/null; then
            log "警告: $cmd 未找到"
            verify_ok=false
        fi
    done
    
    if [ "$verify_ok" = true ]; then
        log "✅ 安装验证通过"
    else
        log "⚠️  部分命令未找到，请检查编译日志"
    fi
    
    log "部署成功!"
    log ""
    log "===========================================" 
    log "使用说明:"
    log "===========================================" 
    log "1. 加载环境变量:"
    log "   source ~/.bashrc"
    log ""
    log "2. 验证安装:"
    log "   make_slc_dj1 --help"
    log "   run_dinsar.sh --help"
    log ""
    log "3. 处理数据:"
    log "   run_dinsar.sh /path/to/project"
    log "   run_dinsar.sh --auto --mode A /path/to/project"
    log ""
    log "安装位置: $INSTALL_PREFIX"
    log "日志文件: $LOG_FILE"
}

main "$@"
