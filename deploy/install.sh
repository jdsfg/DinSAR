#!/bin/bash
#=============================================================================
# PowerChina-1 DInSAR 系统部署脚本 v2.1
# 
# 合并更新版: 包含 GCC 11+ 修复 + ALOS extern 修复
#
# 作者: PowerChina-1 DInSAR Team
# 日期: 2026-01-28
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
            sed -i 's/^[[:space:]]*int verbose;/extern int verbose;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int debug;/extern int debug;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int roi;/extern int roi;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int swap;/extern int swap;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int quad_pol;/extern int quad_pol;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int ALOS_format;/extern int ALOS_format;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int force_slope;/extern int force_slope;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int dopp;/extern int dopp;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int quiet_flag;/extern int quiet_flag;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int SAR_mode;/extern int SAR_mode;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*int prefix_off;/extern int prefix_off;/g' "$alos_inc"
            
            sed -i 's/^[[:space:]]*double forced_slope;/extern double forced_slope;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*double tbias;/extern double tbias;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*double rbias;/extern double rbias;/g' "$alos_inc"
            sed -i 's/^[[:space:]]*double slc_fact;/extern double slc_fact;/g' "$alos_inc"
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
        if ! grep -q "globals_alos.o" "$alos_make"; then
            log "更新 ALOS/lib_src/Makefile..."
            # 更通用：为 SRCS 变量前置 globals_alos.c，保留其余源码列表
            sed -i 's/^[[:space:]]*SRCS[[:space:]]*=[[:space:]]*/SRCS = globals_alos.c /' "$alos_make"
        fi
        
        # 确保 Makefile 使用 Tab 缩进 (之前遇到过 8空格的问题)
        sed -i 's/^        /\t/g' "$alos_make"
    fi
}

check_system_deps() {
    log "检查系统依赖..."
    local deps=("csh" "make" "gcc") 
    local install_cmd=""
    
    if command -v apt-get &>/dev/null; then
        install_cmd="apt-get install -y"
    elif command -v yum &>/dev/null; then
        install_cmd="yum install -y"
    fi

    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &>/dev/null; then
            log "缺少依赖: $dep"
            if [ -n "$install_cmd" ] && [ "$(id -u)" -eq 0 ]; then
                log "尝试自动安装 $dep..."
                $install_cmd "$dep" >> "$LOG_FILE" 2>&1
            else
                error "请手动安装 $dep (sudo apt install $dep / sudo yum install $dep) 后重试"
            fi
        fi
    done
}

main() {
    log "开始部署 v2.3 ..."
    
    check_system_deps

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

    # 4.1 显式安装 DJ1 预处理工具（上游默认流程可能未包含）
    if [ -d "$GMTSAR_SRC/preproc/DJ1_preproc" ]; then
        log "编译并安装 DJ1_preproc (make_slc_dj1)..."
        make -C "$GMTSAR_SRC/preproc/DJ1_preproc" install >> "$LOG_FILE" 2>&1 || error "DJ1_preproc 安装失败，请查看日志: $LOG_FILE"
    else
        log "未找到 DJ1_preproc 目录，跳过 make_slc_dj1 安装"
    fi
    
    # 强制 bin 目录权限
    if [ -d "$INSTALL_PREFIX/bin" ]; then
        chmod +x "$INSTALL_PREFIX/bin"/* 2>/dev/null
    fi

    # 4.2 可选创建系统软链接（便于未 source ~/.bashrc 的场景）
    if [ "$(id -u)" -eq 0 ] && [ "$INSTALL_PREFIX" = "/usr/local/GMTSAR" ] && [ -x "$INSTALL_PREFIX/bin/make_slc_dj1" ]; then
        ln -sf "$INSTALL_PREFIX/bin/make_slc_dj1" /usr/local/bin/make_slc_dj1 || true
    fi
    
    # 5. 自动写入 PATH（避免用户手动配置）
    local shell_rc="$HOME/.bashrc"
    local path_line="export PATH=\"\$PATH:$INSTALL_PREFIX/bin\""
    if [ -f "$shell_rc" ]; then
        if ! grep -Fq "$INSTALL_PREFIX/bin" "$shell_rc"; then
            echo "$path_line" >> "$shell_rc"
            log "已写入 PATH 到 $shell_rc"
        fi
    else
        echo "$path_line" > "$shell_rc"
        log "已创建 $shell_rc 并写入 PATH"
    fi

    log "部署成功!"
    log "已自动更新 PATH，重新登录或执行 'source ~/.bashrc' 后可直接使用 GMTSAR 命令"
    log ""
    log "使用数据处理脚本:"
    log "  ./run_dinsar.sh /path/to/project"
    log "  ./run_dinsar.sh --auto --mode A /path/to/project"
}

main "$@"
