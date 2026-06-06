#!/bin/bash
#=============================================================================
# PowerChina-1 DInSAR 数据处理脚本 v2.1
#
# 支持: DJ1 (电建一号), BC3 (涪城一号/天仪) X波段 SAR 数据
# 
# 用法:
#   ./run_dinsar.sh                           # 交互式向导
#   ./run_dinsar.sh /path/to/project          # 指定项目目录
#   ./run_dinsar.sh --auto --mode A /path     # 全自动处理
#
# 作者: PowerChina-1 DInSAR Team
# 日期: 2026-01-28
#=============================================================================

set -euo pipefail

#-----------------------------------------------------------------------------
# 全局变量
#-----------------------------------------------------------------------------
SCRIPT_VERSION="2.1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 处理参数
PROJECT_DIR=""
AUTO_MODE=false
INTERACTIVE=true

# 数据参数
MASTER_DATE=""
SLAVE_DATE=""
MASTER_XML=""
SLAVE_XML=""
SATELLITE_TYPE="BC3"
WAVELENGTH=0.0311  # X波段默认值 (m)

# 处理参数
PROCESS_MODE="FAST"  # FAST/ STANDARD/ QUALITY
UNWRAP_MODE="A"      # A=稳定模式, B=高精度模式（由 PROCESS_MODE 映射）
ATM_MODE="EMPIRICAL" # NONE/ EMPIRICAL (Phase-Elev)/ GACOS / ST_FILTER
GACOS_PATH=""        # GACOS 数据目录 (.ztd + .rsc)
REF_PX=""            # 参考点 Range 坐标 (可选，GACOS 模式使用)
REF_PY=""            # 参考点 Azimuth 坐标 (可选，GACOS 模式使用)

if [ -n "${DOWNSAMPLE_FACTOR+x}" ]; then
    DOWNSAMPLE_FACTOR_USER_SET=true
else
    DOWNSAMPLE_FACTOR_USER_SET=false
fi
DOWNSAMPLE_FACTOR="${DOWNSAMPLE_FACTOR:-4}"

SNAPHU_TILES="2 2"
SNAPHU_OVERLAP="200"
SNAPHU_NPROC=8
DETREND_ORDER=6

# 运行时变量
LOG_FILE=""
START_TIME=""

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

#-----------------------------------------------------------------------------
# 日志函数
#-----------------------------------------------------------------------------
log_info()  { echo -e "${GREEN}[✓]${NC} $*" | tee -a "$LOG_FILE"; }
log_warn()  { echo -e "${YELLOW}[⚠]${NC} $*" | tee -a "$LOG_FILE"; }
log_error() { echo -e "${RED}[✗]${NC} $*" | tee -a "$LOG_FILE"; }
log_step()  { echo -e "${BLUE}[→]${NC} $*" | tee -a "$LOG_FILE"; }

error_exit() {
    log_error "$1"
    log_error "详细日志: $LOG_FILE"
    exit "${2:-1}"
}

#-----------------------------------------------------------------------------
# 帮助信息
#-----------------------------------------------------------------------------
show_help() {
    cat << EOF
PowerChina-1 DInSAR 数据处理脚本 v${SCRIPT_VERSION}

用法: $0 [选项] [数据目录]

选项:
  --auto             全自动模式 (无交互)
    --mode fast|standard|quality
                                        全流程模式: Fast(默认) / Standard / Quality
                                        兼容旧写法: A=Fast, B=Quality
    --atm none|empirical|gacos
                                        大气校正模式: None / Empirical (默认) / GACOS
    --gacos-path <path>  指定 GACOS 数据目录
  --help             显示此帮助

全流程模式:
    Fast (默认):     优先时效性，缩小搜索窗与控制点密度
    Standard:        平衡时效与质量，推荐常规作业
    Quality:         复杂地形高精度，耗时更长

示例:
  $0                                    # 交互式向导
  $0 /data/project                      # 指定数据目录
    $0 --auto --mode fast /data/project   # 全自动 Fast 模式

数据目录结构:
  project/
  ├── raw/              # 放置 XML + TIFF 文件
  │   ├── image1.xml
  │   ├── image1.tiff
  │   ├── image2.xml
  │   └── image2.tiff
  └── topo/             # 放置 DEM (可选)
      └── dem.grd

EOF
    exit 0
}

#-----------------------------------------------------------------------------
# 参数解析
#-----------------------------------------------------------------------------
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --auto)
                AUTO_MODE=true
                INTERACTIVE=false
                shift
                ;;
            --mode)
                local mode_arg="$2"
                case "${mode_arg^^}" in
                    FAST|F)
                        PROCESS_MODE="FAST"
                        ;;
                    STANDARD|STD|S)
                        PROCESS_MODE="STANDARD"
                        ;;
                    A)
                        PROCESS_MODE="FAST"
                        ;;
                    QUALITY|HIGH|Q|B)
                        PROCESS_MODE="QUALITY"
                        ;;
                    *)
                        error_exit "无效 --mode: $2 (fast|standard|quality 或 A/B)"
                        ;;
                esac
                shift 2
                ;;
            --atm)
                case "${2^^}" in
                    NONE) ATM_MODE="NONE" ;;
                    EMPIRICAL|EMP) ATM_MODE="EMPIRICAL" ;;
                    GACOS) ATM_MODE="GACOS" ;;
                    ST_FILTER|ST|APS) ATM_MODE="ST_FILTER" ;;
                    *) error_exit "无效 --atm: $2 (none|empirical|gacos|st_filter)" ;;
                esac
                shift 2
                ;;
            --gacos-path)
                GACOS_PATH="$2"
                shift 2
                ;;
            --help|-h)
                show_help
                ;;
            -*)
                error_exit "未知选项: $1\n使用 --help 查看帮助"
                ;;
            *)
                PROJECT_DIR="$1"
                shift
                ;;
        esac
    done
    
    # 默认项目目录为当前目录
    if [ -z "$PROJECT_DIR" ]; then
        PROJECT_DIR="$(pwd)"
    fi
    
    # 转为绝对路径
    PROJECT_DIR="$(cd "$PROJECT_DIR" 2>/dev/null && pwd)" || error_exit "无效的项目目录: $PROJECT_DIR"
}

#----------------------------------------------------------------------------- 
# 分级策略：Fast / Standard / Quality
#-----------------------------------------------------------------------------
apply_tier_defaults() {
    local default_downsample=4
    case "$PROCESS_MODE" in
        FAST)
            UNWRAP_MODE="A"
            default_downsample=16
            SNAPHU_TILES="2 2"
            SNAPHU_NPROC=4
            ;;
        STANDARD)
            UNWRAP_MODE="A"
            default_downsample=2
            SNAPHU_TILES="2 2"
            SNAPHU_NPROC=4
            ;;
        QUALITY)
            UNWRAP_MODE="B"
            default_downsample=8
            SNAPHU_TILES="3 3"
            SNAPHU_OVERLAP="250"
            SNAPHU_NPROC=4
            ;;
        *)
            PROCESS_MODE="FAST"
            UNWRAP_MODE="A"
            default_downsample=4
            ;;
    esac

    if [ "$DOWNSAMPLE_FACTOR_USER_SET" = false ]; then
        DOWNSAMPLE_FACTOR="$default_downsample"
    fi
}

#-----------------------------------------------------------------------------
# 系统检测
#-----------------------------------------------------------------------------
check_system() {
    log_step "检测系统环境..."
    
    # 检查内存
    local total_mem_kb=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    local total_mem_gb=$((total_mem_kb / 1024 / 1024))
    
    log_info "系统内存: ${total_mem_gb}GB"
    
    # 根据内存自动调整参数
    if [ "$total_mem_gb" -lt 16 ]; then
        log_warn "内存较低，强制使用 Fast 模式"
        PROCESS_MODE="FAST"
        UNWRAP_MODE="A"
        DOWNSAMPLE_FACTOR=4
    fi
    
    # 检查必要命令
    local required_cmds=("gmt" "make_slc_dj1" "xcorr" "phasediff" "snaphu")
    local missing=()
    
    for cmd in "${required_cmds[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            missing+=("$cmd")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log_error "缺少必要命令: ${missing[*]}"
        log_error "请先运行 install.sh 安装系统"
        return 1
    fi
    
    log_info "系统环境检查通过"
    return 0
}

#-----------------------------------------------------------------------------
# 数据发现模块 (Filename Agnostic)
#-----------------------------------------------------------------------------
discover_data() {
    log_step "扫描数据目录..."
    
    local raw_dir="$PROJECT_DIR/raw"
    
    # 检查 raw 目录
    if [ ! -d "$raw_dir" ]; then
        log_error "数据目录不存在: $raw_dir"
        log_error "请创建 raw/ 目录并放入 XML + TIFF 文件"
        return 1
    fi
    
    # 查找所有 XML 文件 (不区分大小写)
    shopt -s nullglob nocaseglob
    local xml_files=("$raw_dir"/*.xml)
    shopt -u nullglob nocaseglob
    
    if [ ${#xml_files[@]} -eq 0 ]; then
        log_error "未找到 XML 文件"
        return 1
    fi
    
    log_info "发现 ${#xml_files[@]} 个 XML 文件"
    
    # 解析每个 XML 并配对 TIFF
    declare -A data_pairs
    
    for xml in "${xml_files[@]}"; do
        local basename=$(basename "$xml")
        local stem="${basename%.*}"
        
        # 尝试各种 TIFF 扩展名
        local tiff=""
        for ext in tiff TIFF tif TIF; do
            if [ -f "$raw_dir/${stem}.${ext}" ]; then
                tiff="$raw_dir/${stem}.${ext}"
                break
            fi
        done
        
        if [ -z "$tiff" ]; then
            log_warn "孤立 XML (无匹配 TIFF): $basename"
            continue
        fi
        
        # 从 XML 内部解析元数据 (不依赖文件名)
        local acq_date=""
        local sat_id=""
        
        # 尝试多种日期标签格式
        acq_date=$(grep -ioP '(?<=<startTime>|<AcquisitionDate>|<SceneDate>|<CenterTime>)[0-9TZ:\-]+' "$xml" 2>/dev/null | head -1)
        if [ -z "$acq_date" ]; then
            acq_date=$(grep -oP '\d{4}-\d{2}-\d{2}' "$xml" 2>/dev/null | head -1)
        fi
        
        # 标准化日期 YYYYMMDD
        acq_date=$(echo "$acq_date" | tr -d 'TZ:-' | cut -c1-8)
        
        # 卫星ID
        sat_id=$(grep -ioP '(?<=<SatelliteID>|<SensorID>|<Satellite>)[^<]+' "$xml" 2>/dev/null | head -1)
        [ -z "$sat_id" ] && sat_id="BC3"
        
        if [ -n "$acq_date" ]; then
            data_pairs["$acq_date"]="$xml|$tiff|$sat_id"
            log_info "[$sat_id] $acq_date: $basename"
        else
            log_warn "无法解析日期: $basename"
        fi
    done
    
    # 检查是否有足够的数据
    if [ ${#data_pairs[@]} -lt 2 ]; then
        log_error "至少需要 2 景影像进行干涉处理"
        return 1
    fi
    
    # 按日期排序选择主辅影像
    local sorted_dates=($(echo "${!data_pairs[@]}" | tr ' ' '\n' | sort))
    
    MASTER_DATE="${sorted_dates[0]}"
    SLAVE_DATE="${sorted_dates[1]}"
    
    local master_info="${data_pairs[$MASTER_DATE]}"
    local slave_info="${data_pairs[$SLAVE_DATE]}"
    
    MASTER_XML=$(echo "$master_info" | cut -d'|' -f1)
    MASTER_TIFF=$(echo "$master_info" | cut -d'|' -f2)
    SATELLITE_TYPE=$(echo "$master_info" | cut -d'|' -f3)
    
    SLAVE_XML=$(echo "$slave_info" | cut -d'|' -f1)
    SLAVE_TIFF=$(echo "$slave_info" | cut -d'|' -f2)
    
    log_info "主影像 (Master): $MASTER_DATE"
    log_info "辅影像 (Slave):  $SLAVE_DATE"
    log_info "卫星类型: $SATELLITE_TYPE"
    
    return 0
}

#-----------------------------------------------------------------------------
# 处理流程
#-----------------------------------------------------------------------------
setup_directories() {
    log_step "创建工作目录..."
    
    mkdir -p "$PROJECT_DIR"/{SLC,intf,topo,results,debug,logs}
    
    # 设置日志文件
    LOG_FILE="$PROJECT_DIR/logs/processing_$(date +%Y%m%d_%H%M%S).log"
    echo "处理开始: $(date)" > "$LOG_FILE"
    
    log_info "工作目录已创建"
}

generate_trans_dat() {
    cd "$INTF_DIR" || return 1
    [ -f trans.dat ] && return 0
    
    local dem_grd="$PROJECT_DIR/topo/dem.grd"
    if [ ! -f "$dem_grd" ] && [ -n "$DEM_PATH" ] && [ -f "$DEM_PATH" ]; then
        log_info "将 DEM 重投影到 EPSG:4326 (WGS84)..."
        gdalwarp -t_srs EPSG:4326 -r bilinear -of netCDF "$DEM_PATH" "${dem_grd}.tmp" >> "$LOG_FILE" 2>&1
        gmt grdconvert "${dem_grd}.tmp" "$dem_grd" >> "$LOG_FILE" 2>&1 || cp "${dem_grd}.tmp" "$dem_grd"
        rm -f "${dem_grd}.tmp"
    fi

    if [ -f "$dem_grd" ] && command -v SAT_llt2rat &>/dev/null && [ -f "${MASTER_DATE}.PRM" ]; then
        log_info "生成地理编码查找表 trans.dat (可能需数分钟)..."
        local rshift_orig=$(grep rshift "${MASTER_DATE}.PRM" 2>/dev/null | tail -1 | awk '{print $3}')
        command -v update_PRM &>/dev/null && update_PRM "${MASTER_DATE}.PRM" rshift 0 >> "$LOG_FILE" 2>&1
        gmt grd2xyz --FORMAT_FLOAT_OUT=%lf "$dem_grd" -s 2>/dev/null | SAT_llt2rat "${MASTER_DATE}.PRM" 1 -bod > trans.dat 2>> "$LOG_FILE"
        if [ -n "$rshift_orig" ] && command -v update_PRM &>/dev/null; then
            update_PRM "${MASTER_DATE}.PRM" rshift "$rshift_orig" >> "$LOG_FILE" 2>&1
        fi
        [ -s trans.dat ] && return 0
    fi
    return 1
}

generate_slc() {
    log_step "[1/8] 生成 SLC..."
    
    cd "$PROJECT_DIR/SLC" || error_exit "无法进入 SLC 目录"
    
    # 生成主影像 SLC
    log_info "生成主影像 SLC ($MASTER_DATE)..."
    if ! make_slc_dj1 "$MASTER_XML" "$MASTER_TIFF" "$MASTER_DATE" >> "$LOG_FILE" 2>&1; then
        log_warn "make_slc_dj1 失败，尝试 make_slc_s1a..."
        make_slc_s1a "$MASTER_XML" "$MASTER_TIFF" "$MASTER_DATE" >> "$LOG_FILE" 2>&1 || error_exit "主影像 SLC 生成失败"
    fi
    
    # 生成辅影像 SLC
    log_info "生成辅影像 SLC ($SLAVE_DATE)..."
    if ! make_slc_dj1 "$SLAVE_XML" "$SLAVE_TIFF" "$SLAVE_DATE" >> "$LOG_FILE" 2>&1; then
        log_warn "make_slc_dj1 失败，尝试 make_slc_s1a..."
        make_slc_s1a "$SLAVE_XML" "$SLAVE_TIFF" "$SLAVE_DATE" >> "$LOG_FILE" 2>&1 || error_exit "辅影像 SLC 生成失败"
    fi
    
    # 核心修复: 自动修复 PRM 文件参数 (earth_radius, SC_vel, SC_height, tie_point)
    if [ -f "$PROJECT_DIR/SLC/fix_prm.py" ]; then
        log_info "修复主影像 PRM 文件参数..."
        python3 "$PROJECT_DIR/SLC/fix_prm.py" "${MASTER_DATE}.PRM" "$MASTER_XML" "${MASTER_DATE}.LED" >> "$LOG_FILE" 2>&1
        log_info "修复辅影像 PRM 文件参数..."
        python3 "$PROJECT_DIR/SLC/fix_prm.py" "${SLAVE_DATE}.PRM" "$SLAVE_XML" "${SLAVE_DATE}.LED" >> "$LOG_FILE" 2>&1
    fi
    
    log_info "SLC 生成完成"
}

run_alignment() {
    log_step "[2/8] 影像配准..."
    
    cd "$PROJECT_DIR/SLC" || error_exit "无法进入 SLC 目录"
    
    local master_prm="${MASTER_DATE}.PRM"
    local slave_prm="${SLAVE_DATE}.PRM"
    
    if [ ! -f "$master_prm" ] || [ ! -f "$slave_prm" ]; then
        error_exit "PRM 文件不存在"
    fi
    
    local xsearch=128
    local ysearch=128
    local nx_corr=30
    local ny_corr=50
    local threshold_snr=10

    case "$PROCESS_MODE" in
        FAST)
            xsearch=64
            ysearch=64
            nx_corr=5
            ny_corr=10
            threshold_snr=8
            ;;
        STANDARD)
            xsearch=128
            ysearch=128
            nx_corr=100
            ny_corr=50
            threshold_snr=15
            ;;
        QUALITY)
            xsearch=256
            ysearch=256
            nx_corr=200
            ny_corr=50
            threshold_snr=20
            ;;
    esac

    log_info "运行 xcorr（${PROCESS_MODE}）..."
    xcorr "$master_prm" "$slave_prm" -xsearch "$xsearch" -ysearch "$ysearch" -nx "$nx_corr" -ny "$ny_corr" >> "$LOG_FILE" 2>&1
    
    if [ $? -ne 0 ]; then
        log_warn "xcorr 返回非零，检查日志..."
    fi

    # 若存在 fitoffset.csh 且生成了 freq_xcorr.dat，则用 threshold_snr 拟合并写回 PRM
    if command -v fitoffset.csh &>/dev/null && [ -f freq_xcorr.dat ]; then
        if fitoffset.csh 3 3 freq_xcorr.dat "$threshold_snr" >> "$slave_prm" 2>> "$LOG_FILE"; then
            log_info "fitoffset 已应用 (threshold_snr=$threshold_snr)"
        else
            log_warn "fitoffset 失败，使用 xcorr 原始偏移"
        fi
    else
        log_warn "fitoffset.csh 不可用或 freq_xcorr.dat 不存在，跳过偏移拟合"
    fi
    
    log_info "配准完成"
}

run_interferogram() {
    log_step "[3/8] 生成干涉图..."
    
    # Checkpoint
    local intf_dir="$PROJECT_DIR/intf/${MASTER_DATE}_${SLAVE_DATE}"
    if [ -f "$intf_dir/phase.grd" ] && [ -f "$intf_dir/amp.grd" ]; then
        log_info "已存在 phase.grd, 跳过干涉图生成 (Checkpoint)"
        INTF_DIR="$intf_dir"
        return 0
    fi
    
    # 创建干涉目录
    local intf_dir="$PROJECT_DIR/intf/${MASTER_DATE}_${SLAVE_DATE}"
    mkdir -p "$intf_dir"
    cd "$intf_dir" || error_exit "无法进入干涉目录"
    
    # 链接 PRM 和 SLC 文件
    ln -sf "$PROJECT_DIR/SLC/${MASTER_DATE}.PRM" .
    ln -sf "$PROJECT_DIR/SLC/${MASTER_DATE}.SLC" .
    ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}.PRM" .
    ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}.SLC" .
    ln -sf "$PROJECT_DIR/SLC/${MASTER_DATE}.LED" . 2>/dev/null
    ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}.LED" . 2>/dev/null

    # 兼容 GMTSAR 标准流程：在 phasediff 前补齐基线参数
    if command -v SAT_baseline &>/dev/null; then
        SAT_baseline "${MASTER_DATE}.PRM" "${SLAVE_DATE}.PRM" | tail -n 9 >> "${SLAVE_DATE}.PRM" 2>/dev/null || true
        SAT_baseline "${MASTER_DATE}.PRM" "${MASTER_DATE}.PRM" | grep height >> "${MASTER_DATE}.PRM" 2>/dev/null || true
    fi
    
    # 核心修复: 生成并去除地形相位
    local dem_grd="$PROJECT_DIR/topo/dem.grd"
    local topo_arg=""
    if [ -f "$dem_grd" ]; then
        if [ ! -f dem_ra.grd ]; then
            log_info "生成雷达坐标 DEM (dem_ra.grd) 以去除地形相位..."
            bash /home/mihu/gmtsar/DInSAR_v2.3_Delivery/fast_dem2topo_ra.sh "${MASTER_DATE}.PRM" "$dem_grd" "$LOG_FILE"
        fi
        if [ -f dem_ra.grd ]; then
            # 核心修复: 替换 dem_ra.grd 中的 NaN 值为 0，防止 phasediff.c 的 calc_average_topo 算得 NaN 导致全图变 NaN/0
            log_info "清洗雷达坐标 DEM (dem_ra.grd) 以免除 NaN 引起的计算链崩溃..."
            gmt grdmath dem_ra.grd 0 AND = dem_ra_clean.grd >> "$LOG_FILE" 2>&1
            mv dem_ra_clean.grd dem_ra.grd
            topo_arg="-topo dem_ra.grd"
            log_info "将去除地形相位 (Differential InSAR 模式)"
        fi
    else
        log_warn "未找到 $dem_grd，只能生成含地形相位的原始干涉图"
    fi
    
    log_info "运行 phasediff..."
    phasediff "${MASTER_DATE}.PRM" "${SLAVE_DATE}.PRM" -imag imag.grd -real real.grd $topo_arg >> "$LOG_FILE" 2>&1
    
    if [ ! -f real.grd ] || [ ! -f imag.grd ]; then
        error_exit "phasediff 失败"
    fi
    
    # 计算相位和振幅
    log_info "计算相位和振幅..."
    gmt grdmath real.grd imag.grd ATAN2 = phase.grd
    gmt grdmath real.grd imag.grd HYPOT = amp.grd
    
    # 真实的相干性应当由滤波步骤 (phasefilt) 计算
    log_info "干涉图生成完成"
    
    # 保存干涉目录路径
    INTF_DIR="$intf_dir"
}

run_filter() {
    log_step "[4/8] 相位滤波..."
    
    cd "$INTF_DIR" || error_exit "无法进入干涉目录"
    
    if [ -f "phasefilt.grd" ]; then
        log_info "已存在 phasefilt.grd, 跳过滤波 (Checkpoint)"
        FILT_PHASE="phasefilt.grd"
        return 0
    fi
    
    # 降采样 (如果需要)
    if [ "$DOWNSAMPLE_FACTOR" -gt 1 ]; then
        log_info "抗混叠降采样 ${DOWNSAMPLE_FACTOR}x (Gaussian Filter)..."
        local xinc yinc xinc_dec yinc_dec
        xinc=$(gmt grdinfo phase.grd -C 2>/dev/null | awk '{print $8}')
        yinc=$(gmt grdinfo phase.grd -C 2>/dev/null | awk '{print $9}')
        xinc_dec=$(awk -v inc="$xinc" -v f="$DOWNSAMPLE_FACTOR" 'BEGIN{printf "%.12g", inc*f}')
        yinc_dec=$(awk -v inc="$yinc" -v f="$DOWNSAMPLE_FACTOR" 'BEGIN{printf "%.12g", inc*f}')
        
        local fw=$((DOWNSAMPLE_FACTOR))
        # 核心优化1：使用高斯滤波进行纯复数域空间降采样，消除频域混叠
        gmt grdfilter real.grd -Fg${fw}/${fw} -D0 -I${xinc_dec}/${yinc_dec} -Greal_dec.grd >> "$LOG_FILE" 2>&1
        gmt grdfilter imag.grd -Fg${fw}/${fw} -D0 -I${xinc_dec}/${yinc_dec} -Gimag_dec.grd >> "$LOG_FILE" 2>&1
        
        gmt grdmath real_dec.grd imag_dec.grd ATAN2 = phase_dec.grd
        gmt grdmath real_dec.grd imag_dec.grd HYPOT = amp_dec.grd
        
        PHASE_FILE="phase_dec.grd"
        REAL_FILE="real_dec.grd"
        IMAG_FILE="imag_dec.grd"
        AMP_FILE="amp_dec.grd"
    else
        PHASE_FILE="phase.grd"
        REAL_FILE="real.grd"
        IMAG_FILE="imag.grd"
        AMP_FILE="amp.grd"
    fi
    
    # 核心优化2：增强型 Goldstein 滤波 (Alpha 强自适应)
    log_info "强力自适应 Goldstein 滤波 (Alpha=0.8)..."
    if command -v phasefilt &>/dev/null; then
        if [ ! -s "$REAL_FILE" ] || [ ! -s "$IMAG_FILE" ]; then
            error_exit "phasefilt 输入文件缺失或为空 (real/imag)"
        fi

        # 剥离无意义的 amp 输入，强制使用固定高强度的 alpha 以适应高山峡谷 (X-band)
        log_info "执行: phasefilt -imag $IMAG_FILE -real $REAL_FILE -alpha 0.8 -psize 32"
        phasefilt -imag "$IMAG_FILE" -real "$REAL_FILE" -alpha 0.8 -psize 32 >> "$LOG_FILE" 2>&1
        local phasefilt_rc=$?
        if [ $phasefilt_rc -ne 0 ]; then
            error_exit "phasefilt 失败 (退出码: $phasefilt_rc)，请检查日志"
        fi
        
        if [ -s filtphase.grd ]; then
            mv -f filtphase.grd phasefilt.grd 2>/dev/null
        fi
        if [ -s phasefilt.grd ]; then
            FILT_PHASE="phasefilt.grd"
        else
            error_exit "phasefilt 未生成有效输出 (filtphase.grd/phasefilt.grd)"
        fi
        
        # 核心优化3：重构真实物理相干性评估 (Pearson Correlation)
        log_info "计算基于空间包络滤波的真实干涉相干性 (Pearson Correlation)..."
        # 1. 对复数干涉图的实部和虚部进行空间平滑，计算分子强度
        gmt grdfilter "$REAL_FILE" -Fg20/20 -D0 -Greal_mean.grd >> "$LOG_FILE" 2>&1
        gmt grdfilter "$IMAG_FILE" -Fg20/20 -D0 -Gimag_mean.grd >> "$LOG_FILE" 2>&1
        gmt grdmath real_mean.grd imag_mean.grd HYPOT = num_amp.grd
        # 2. 对复数干涉图的绝对振幅进行空间平滑，计算分母
        gmt grdmath "$REAL_FILE" "$IMAG_FILE" HYPOT = raw_amp.grd
        gmt grdfilter raw_amp.grd -Fg20/20 -D0 -Gden_amp.grd >> "$LOG_FILE" 2>&1
        # 3. 计算相干系数 (防止除零，并将上限截断到 1.0，安全替换 NaN 为 0)
        gmt grdmath num_amp.grd den_amp.grd DIV 1.0 MIN 0 AND = corr.grd
    else
        FILT_PHASE="$PHASE_FILE"
        log_warn "phasefilt 不可用，跳过滤波"
    fi
    
    log_info "滤波完成"
}

run_unwrap() {
    log_step "[5/8] 相位解缠 (模式 $UNWRAP_MODE)..."
    
    cd "$INTF_DIR" || error_exit "无法进入干涉目录"

    # 解缠恢复：若已存在 unwrap.grd 则跳过解缠，便于中断后从去趋势步骤继续
    if [ -f unwrap.grd ]; then
        log_info "已存在解缠结果 unwrap.grd，跳过解缠步骤"
        return 0
    fi
    
    # 获取影像尺寸
    local ncols=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $10}')
    local nrows=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $11}')
    
    log_info "影像尺寸: ${ncols} x ${nrows}"
    
    # 相干性掩膜（低相干区域置零，降低解缠负担）
    local phase_for_unwrap="$FILT_PHASE"
    local corr_source="corr.grd"
    
    if [ -f "$corr_source" ]; then
        log_info "生成相干性硬掩膜 (阈值=0.20): $corr_source"
        gmt grdmath "$corr_source" 0.20 GE = coh_mask.grd >> "$LOG_FILE" 2>&1
        gmt grdmath "$FILT_PHASE" coh_mask.grd MUL 0 AND = phase_masked.grd >> "$LOG_FILE" 2>&1
        if [ -f phase_masked.grd ]; then
            phase_for_unwrap="phase_masked.grd"
        fi
    fi

    # 转换为二进制格式并安全替换 NaN 值为 0 (防止 snaphu 崩溃)
    log_info "准备 snaphu 相位与相干性双输入阵列..."
    gmt grdmath "$phase_for_unwrap" 0 AND = phase_no_nan.grd >> "$LOG_FILE" 2>&1
    gmt grd2xyz phase_no_nan.grd -ZTLf > phase.bin 2>> "$LOG_FILE"
    if [ -f "$corr_source" ]; then
        gmt grdmath "$corr_source" 0 AND = corr_no_nan.grd >> "$LOG_FILE" 2>&1
        gmt grd2xyz corr_no_nan.grd -ZTLf > corr.bin 2>> "$LOG_FILE"
    fi
    
    # 构建 snaphu 命令
    local snaphu_cmd="snaphu phase.bin $ncols -o unwrap.bin"
    local snaphu_conf=""
    local conf_candidates=(
        "/usr/local/GMTSAR/share/gmtsar/snaphu/config/snaphu.conf.brief"
        "$SCRIPT_DIR/gmtsar-master/share/gmtsar/snaphu/config/snaphu.conf.brief"
        "$SCRIPT_DIR/dinsar_deploy_v2.2_offline/gmtsar-master/share/gmtsar/snaphu/config/snaphu.conf.brief"
    )
    for c in "${conf_candidates[@]}"; do
        if [ -f "$c" ]; then
            snaphu_conf="$c"
            break
        fi
    done
    if [ -n "$snaphu_conf" ]; then
        snaphu_cmd="$snaphu_cmd -f $snaphu_conf"
    else
        log_warn "未找到 snaphu.conf.brief，使用 snaphu 默认配置"
    fi
    
    if [ "$UNWRAP_MODE" = "A" ]; then
        # 稳定模式: 分块 + 多线程
        snaphu_cmd="$snaphu_cmd --tile $SNAPHU_TILES $SNAPHU_OVERLAP $SNAPHU_OVERLAP --nproc $SNAPHU_NPROC"
        log_info "稳定模式: --tile $SNAPHU_TILES overlap=${SNAPHU_OVERLAP} --nproc $SNAPHU_NPROC"
    else
        # 高精度模式: 整图单线程
        log_info "高精度模式: 整图解缠 (需要较长时间)"
    fi
    
    snaphu_cmd="$snaphu_cmd -d"
    
    # 核心优化4：将相干性输入 SNAPHU 以构建高精度成本流网络 (Cost Graph)
    if [ -f corr.bin ]; then
        snaphu_cmd="$snaphu_cmd -c corr.bin"
        log_info "启用基于相关性的统计成本函数模型 (-c corr.bin)"
    fi
    
    log_info "运行 snaphu..."
    log_info "命令: $snaphu_cmd"
    
    # 记录开始时间
    local unwrap_start=$(date +%s)
    
    # 运行 snaphu
    $snaphu_cmd >> "$LOG_FILE" 2>&1
    
    if [ $? -ne 0 ]; then
        log_error "snaphu 失败，查看日志: $LOG_FILE"
        return 1
    fi
    
    # 记录耗时
    local unwrap_end=$(date +%s)
    local unwrap_time=$((unwrap_end - unwrap_start))
    log_info "snaphu 完成，耗时: $((unwrap_time / 60)) 分 $((unwrap_time % 60)) 秒"
    
    # 校验 unwrap.bin 尺寸（必须等于 ncols*nrows*4 bytes）
    if [ ! -s unwrap.bin ]; then
        error_exit "snaphu 未生成 unwrap.bin"
    fi
    local expected_bytes=$((ncols * nrows * 4))
    local actual_bytes
    actual_bytes=$(stat -c%s unwrap.bin 2>/dev/null || echo 0)
    if [ "$actual_bytes" -ne "$expected_bytes" ]; then
        error_exit "unwrap.bin 尺寸异常: 期望 ${expected_bytes} bytes (${ncols}x${nrows}x4)，实际 ${actual_bytes}"
    fi

    # 转换回 GRD 格式
    log_info "转换解缠结果..."
    local xmin=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $2}')
    local xmax=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $3}')
    local ymin=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $4}')
    local ymax=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $5}')
    local xinc=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $8}')
    local yinc=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $9}')
    
    gmt xyz2grd unwrap.bin -Gunwrap.grd -I$xinc/$yinc -R$xmin/$xmax/$ymin/$ymax -r -ZTLf >> "$LOG_FILE" 2>&1
    if [ $? -ne 0 ]; then
        error_exit "gmt xyz2grd 失败，无法生成 unwrap.grd"
    fi
    
    if [ ! -f unwrap.grd ]; then
        error_exit "解缠结果转换失败"
    fi
    
    # 备份原始解缠结果
    cp unwrap.grd "$PROJECT_DIR/debug/unwrap_raw.grd"
    
    log_info "相位解缠完成"
}

run_atm_correction() {
    log_step "[5.5/8] 大气校正 (${ATM_MODE})..."
    
    if [ "$ATM_MODE" = "NONE" ]; then
        log_info "大气校正被禁用 (ATM_MODE=NONE)"
        return 0
    fi

    cd "$INTF_DIR" || error_exit "无法进入干涉目录"
    
    if [ ! -f unwrap.grd ]; then
        log_warn "未找到解缠结果，跳过大气校正"
        return 0
    fi
    
    # 检查是否有 DEM
    local dem_grd="$PROJECT_DIR/topo/dem.grd"
    
    # [模式 A] 经验改正
    if [ "$ATM_MODE" = "EMPIRICAL" ]; then
        if [ ! -f "$dem_grd" ]; then
            log_warn "未找到 topo/dem.grd，无法进行地形相关大气校正"
            cp unwrap.grd unwrap_atm.grd
            return 0
        fi
        
        # 生成雷达坐标 DEM (dem_ra.grd) - 使用快速算法
        if [ ! -f dem_ra.grd ]; then
            log_info "正在生成雷达坐标 DEM (dem_ra.grd) [快速算法]..."
            bash /home/mihu/gmtsar/DInSAR_v2.3_Delivery/fast_dem2topo_ra.sh "${MASTER_DATE}.PRM" "$dem_grd" "$LOG_FILE"
        fi
        
        if [ ! -f dem_ra.grd ]; then
            log_warn "无法获取雷达坐标 DEM，跳过经验大气校正"
            return 0
        fi

        log_info "执行 Re_atm.csh 进行经验校正..."
        if [ -f "$SCRIPT_DIR/gmtsar-master/gmtsar/csh/Re_atm.csh" ]; then
            csh -f "$SCRIPT_DIR/gmtsar-master/gmtsar/csh/Re_atm.csh" unwrap.grd dem_ra.grd unwrap_atm.grd >> "$LOG_FILE" 2>&1
        else
            log_warn "未找到 Re_atm.csh，跳过大气校正"
            return 0
        fi

    # [模式 B] GACOS 校正
    elif [ "$ATM_MODE" = "GACOS" ]; then
        if [ -z "$GACOS_PATH" ]; then
            log_error "GACOS 模式需要指定 --gacos-path"
            return 1
        fi
        
        # GACOS 必须先有 trans.dat
        if ! generate_trans_dat; then
            log_error "GACOS 模式需要 trans.dat，但生成失败 (检查 DEM 是否存在)"
            return 1
        fi

        # 确定参考点 (默认为图像中心)
        if [ -z "$REF_PX" ] || [ -z "$REF_PY" ]; then
            local nx=$(gmt grdinfo unwrap.grd -C | awk '{print $10}')
            local ny=$(gmt grdinfo unwrap.grd -C | awk '{print $11}')
            REF_PX=$((nx / 2))
            REF_PY=$((ny / 2))
            log_info "使用默认参考点 (图像中心): $REF_PX, $REF_PY"
        fi

        log_info "执行 make_gacos_correction.csh..."
        local gacos_script="$SCRIPT_DIR/gmtsar-master/gmtsar/csh/make_gacos_correction.csh"
        if [ ! -f "$gacos_script" ]; then
            log_error "未找到 GACOS 脚本: $gacos_script"
            return 1
        fi
        
        # 参数: intf_dir GACOS_path ref_range ref_azimuth dem.grd
        # 注意: make_gacos_correction.csh 内部会 cd $1，所以传入 $PWD
        csh -f "$gacos_script" "$PWD" "$GACOS_PATH" "$REF_PX" "$REF_PY" "$dem_grd" >> "$LOG_FILE" 2>&1
        
        if [ -f unwrap_gacos_corrected.grd ]; then
            cp unwrap_gacos_corrected.grd unwrap_atm.grd
        else
            log_error "GACOS 校正未生成预期结果 (unwrap_gacos_corrected.grd)"
            return 1
        fi
    
    # [模式 C] 时空滤波 (ST_FILTER)
    elif [ "$ATM_MODE" = "ST_FILTER" ]; then
        log_info "执行 ST_Filter.csh 进行时空/空间自适应滤波..."
        local st_script="$SCRIPT_DIR/gmtsar-master/gmtsar/csh/ST_Filter.csh"
        if [ ! -f "$st_script" ]; then
            st_script="ST_Filter.csh"
        fi
        
        if ! command -v "$st_script" &>/dev/null && [ ! -f "$st_script" ]; then
            log_warn "未找到 ST_Filter.csh，跳过滤波"
            cp unwrap.grd unwrap_atm.grd
        else
            csh -f "$st_script" unwrap.grd unwrap_atm.grd "$PROJECT_DIR" >> "$LOG_FILE" 2>&1
        fi
    fi

    if [ -f unwrap_atm.grd ]; then
        log_info "大气校正完成 ($ATM_MODE)"
        mv unwrap.grd unwrap_pre_atm.grd
        cp unwrap_atm.grd unwrap.grd
    else
        log_warn "大气校正未产生输出"
    fi
}

run_detrend() {
    log_step "[6/8] 去轨道倾斜..."
    
    cd "$INTF_DIR" || error_exit "无法进入干涉目录"
    
    if [ -f "unwrap_detrended.grd" ]; then
        log_info "已存在 unwrap_detrended.grd, 跳过去倾斜 (Checkpoint)"
        return 0
    fi
    
    if [ ! -f unwrap.grd ]; then
        log_warn "未找到解缠结果，跳过去倾斜"
        return 0
    fi
    
    # 关键: 在相位域进行去倾斜 (物理正确)
    log_info "拟合轨道残余趋势面 (N=$DETREND_ORDER)..."
    
    gmt grdtrend unwrap.grd -N${DETREND_ORDER}r \
        -T"$PROJECT_DIR/debug/orbital_ramp.grd" \
        -Dunwrap_detrended.grd \
        -V >> "$LOG_FILE" 2>&1
    
    if [ $? -ne 0 ] || [ ! -f unwrap_detrended.grd ]; then
        log_warn "grdtrend 失败，使用原始解缠结果"
        cp unwrap.grd unwrap_detrended.grd
    fi
    
    log_info "轨道残余已保存: debug/orbital_ramp.grd"
    log_info "去倾斜完成"
}

run_displacement() {
    log_step "[7/8] 计算位移..."
    
    cd "$INTF_DIR" || error_exit "无法进入干涉目录"
    
    local input_phase="unwrap_detrended.grd"
    if [ ! -f "$input_phase" ]; then
        input_phase="unwrap.grd"
    fi
    
    if [ ! -f "$input_phase" ]; then
        log_warn "未找到解缠相位，跳过位移计算"
        return 0
    fi
    
    # 转换为位移 (mm)
    # 公式: d = φ × λ / (-4π) × 1000
    log_info "转换相位为 LOS 位移..."
    
    gmt grdmath "$input_phase" $WAVELENGTH MUL -4 PI MUL DIV 1000 MUL = los_displacement.grd
    
    if [ ! -f los_displacement.grd ]; then
        error_exit "位移计算失败"
    fi
    
    # 统计位移范围
    local disp_info=$(gmt grdinfo los_displacement.grd -C)
    local disp_min=$(echo "$disp_info" | awk '{printf "%.1f", $6}')
    local disp_max=$(echo "$disp_info" | awk '{printf "%.1f", $7}')
    
    log_info "LOS 位移范围: $disp_min ~ $disp_max mm"
    
    # 记录处理参数 (可追溯性)
    cat >> "$PROJECT_DIR/debug/processing_params.txt" << EOF
Processing Parameters (v2.1)
============================
Date: $(date)
Master: $MASTER_DATE
Slave: $SLAVE_DATE
Satellite: $SATELLITE_TYPE
Wavelength: $WAVELENGTH m
Detrend Order: $DETREND_ORDER
Process Mode: $PROCESS_MODE
Unwrap Mode: $UNWRAP_MODE
Downsample: ${DOWNSAMPLE_FACTOR}x
Displacement Range: $disp_min ~ $disp_max mm
EOF
    
    log_info "位移计算完成"
}

generate_outputs() {
    log_step "[8/8] 生成输出..."
    
    cd "$INTF_DIR" || error_exit "无法进入干涉目录"
    
    # 核心修复: 确保 FILT_PHASE 变量在跳过 Checkpoint 时依旧被正确注入
    if [ -z "$FILT_PHASE" ]; then
        if [ -f "phasefilt.grd" ]; then
            FILT_PHASE="phasefilt.grd"
        elif [ -f "phase_dec.grd" ]; then
            FILT_PHASE="phase_dec.grd"
        else
            FILT_PHASE="phase.grd"
        fi
    fi
    
    local results_dir="$PROJECT_DIR/results"
    mkdir -p "$results_dir"
    
    # 若缺少 dem.grd，尝试从 DEM_PATH 或 topo/*.tif 自动生成
    if [ ! -f "$PROJECT_DIR/topo/dem.grd" ]; then
        local src_dem=""
        if [ -n "$DEM_PATH" ] && [ -f "$DEM_PATH" ]; then
            src_dem="$DEM_PATH"
        else
            src_dem=$(ls "$PROJECT_DIR"/topo/*.{tif,TIF,tiff,TIFF} 2>/dev/null | head -n 1)
        fi
        if [ -n "$src_dem" ]; then
            log_info "检测到 DEM: $src_dem，强制重投影为 EPSG:4326 以保证地理编码不越界..."
            gdalwarp -t_srs EPSG:4326 -r bilinear -of netCDF "$src_dem" "${PROJECT_DIR}/topo/dem_tmp.nc" >> "$LOG_FILE" 2>&1
            gmt grdconvert "${PROJECT_DIR}/topo/dem_tmp.nc" "$PROJECT_DIR/topo/dem.grd" >> "$LOG_FILE" 2>&1 || cp "${PROJECT_DIR}/topo/dem_tmp.nc" "$PROJECT_DIR/topo/dem.grd"
            rm -f "${PROJECT_DIR}/topo/dem_tmp.nc"
        fi
    fi

    # GMTSAR bin（proj_ra2ll.csh 等）
    local gmtsar_bin="$SCRIPT_DIR/gmtsar-master/bin"
    
    # -------------------------------------------------------------------------
    # [Stage 6] 前端图层与数据导出 (GeoTIFF Export) — 强制地理编码
    # 需求：RA -> LL 地理编码，输出带 EPSG:4326 的 GeoTIFF，确保 ArcGIS 可读
    # -------------------------------------------------------------------------
    
    local have_trans=0
    if [ -f trans.dat ]; then
        have_trans=1
    elif generate_trans_dat; then
        have_trans=1
    fi
    
    if [ "$have_trans" -eq 1 ]; then
        log_info "[Stage 6] 执行地理编码 (RA -> LL)..."
        local proj_csh="$gmtsar_bin/proj_ra2ll.csh"
        if [ ! -x "$proj_csh" ]; then
            proj_csh="proj_ra2ll.csh"
        fi
        # 雷达坐标 -> 经纬度 GeoTIFF (使用高可用性 Python 算法)
        local py_geocode="/home/mihu/gmtsar/DInSAR_v2.3_Delivery/geocode_gdal.py"
        if [ -f unwrap_detrended.grd ]; then
            python3 "$py_geocode" trans.dat unwrap_detrended.grd "$results_dir/unwrap_detrended.tif" 0.0001 >> "$LOG_FILE" 2>&1 || true
        fi
        local corr_src=""
        [ -f filtcorr.grd ] && corr_src="filtcorr.grd" || corr_src="corr.grd"
        if [ -n "$corr_src" ] && [ -f "$corr_src" ]; then
            python3 "$py_geocode" trans.dat "$corr_src" "$results_dir/corr.tif" 0.0001 >> "$LOG_FILE" 2>&1 || true
        fi
        if [ -f los_displacement.grd ]; then
            python3 "$py_geocode" trans.dat los_displacement.grd "$results_dir/final_data_values.tif" 0.0001 >> "$LOG_FILE" 2>&1 || true
        fi
        
        # ----- 输出 A：前端图层 (Visual) - RGBA GeoTIFF，背景透明 -----
        # 位移：Jet 色表；相干：Gray 色表；grdimage -A -Q 生成带地理坐标的渲染图
        if [ -f los_displacement_ll.grd ]; then
            local dmin dmax
            dmin=$(gmt grdinfo los_displacement_ll.grd -C 2>/dev/null | awk '{print $6}')
            dmax=$(gmt grdinfo los_displacement_ll.grd -C 2>/dev/null | awk '{print $7}')
            [ -z "$dmin" ] && dmin=-100
            [ -z "$dmax" ] && dmax=100
            gmt makecpt -Cjet -T"$dmin"/"$dmax" -Z > def_layer.cpt 2>> "$LOG_FILE"
            gmt grdimage los_displacement_ll.grd -Rlos_displacement_ll.grd -JM15c -Cdef_layer.cpt -Q -A"$results_dir/final_layer_deformation.tif" >> "$LOG_FILE" 2>&1
            log_info "前端图层(位移): results/final_layer_deformation.tif"
        fi
        if [ -f corr_ll.grd ]; then
            local cmin cmax
            cmin=$(gmt grdinfo corr_ll.grd -C 2>/dev/null | awk '{print $6}')
            cmax=$(gmt grdinfo corr_ll.grd -C 2>/dev/null | awk '{print $7}')
            [ -z "$cmin" ] && cmin=0
            [ -z "$cmax" ] && cmax=1
            gmt makecpt -Cgray -T"$cmin"/"$cmax" -Z > coh_layer.cpt 2>> "$LOG_FILE"
            gmt grdimage corr_ll.grd -Rcorr_ll.grd -JM15c -Ccoh_layer.cpt -Q -A"$results_dir/final_layer_coherence.tif" >> "$LOG_FILE" 2>&1
            log_info "前端图层(相干): results/final_layer_coherence.tif"
        fi
        
        # ----- 输出 B：后端数据 (Data) - 纯数值 GeoTIFF，EPSG:4326 -----
        if [ -f los_displacement_ll.grd ]; then
            gmt grdconvert los_displacement_ll.grd "$results_dir/final_data_values.tif=gd:GTiff" >> "$LOG_FILE" 2>&1
            if command -v gdal_edit.py &>/dev/null; then
                gdal_edit.py -a_srs EPSG:4326 "$results_dir/final_data_values.tif" 2>> "$LOG_FILE" || true
            elif command -v gdal_translate &>/dev/null; then
                gdal_translate -a_srs EPSG:4326 "$results_dir/final_data_values.tif" "$results_dir/final_data_values_epsg.tif" 2>> "$LOG_FILE" && mv "$results_dir/final_data_values_epsg.tif" "$results_dir/final_data_values.tif" || true
            fi
            log_info "数据 GeoTIFF(数值): results/final_data_values.tif (EPSG:4326)"
        fi
    else
        log_warn "未找到 trans.dat 且无法从 DEM 生成，跳过地理编码；仅输出雷达坐标结果。提供 topo/dem.grd 可启用 GeoTIFF 地理编码。"
    fi

    # PNG 输出（使用动态范围）
    log_info "生成可视化图像..."
    
    # 干涉条纹图
    if [ -f "$FILT_PHASE" ]; then
        local pmin=$(gmt grdinfo "$FILT_PHASE" -C 2>/dev/null | awk '{print $6}')
        local pmax=$(gmt grdinfo "$FILT_PHASE" -C 2>/dev/null | awk '{print $7}')
        if [ -n "$pmin" ] && [ -n "$pmax" ] && [ "$pmin" != "NaN" ] && [ "$pmax" != "NaN" ] && [ "$pmin" != "$pmax" ]; then
            gmt begin "$results_dir/interferogram" png E100 >> "$LOG_FILE" 2>&1
            gmt makecpt -Ccyclic -T"$pmin"/"$pmax" > phase.cpt 2>> "$LOG_FILE"
            gmt grdimage "$FILT_PHASE" -JX15c -Cphase.cpt -Baf >> "$LOG_FILE" 2>&1
            gmt colorbar -Cphase.cpt -DJBC+w12c/0.4c -Baf+l"Phase (rad)" >> "$LOG_FILE" 2>&1
            gmt end >> "$LOG_FILE" 2>&1
            log_info "干涉条纹图: results/interferogram.png"
        else
            log_warn "相位数据范围不足或全为零，跳过干涉条纹图生成"
        fi
    fi
    
    # 位移图
    if [ -f los_displacement.grd ]; then
        local dmin2=$(gmt grdinfo los_displacement.grd -C 2>/dev/null | awk '{printf "%.1f", $6}')
        local dmax2=$(gmt grdinfo los_displacement.grd -C 2>/dev/null | awk '{printf "%.1f", $7}')
        if [ -n "$dmin2" ] && [ -n "$dmax2" ] && [ "$dmin2" != "NaN" ] && [ "$dmax2" != "NaN" ] && [ "$dmin2" != "$dmax2" ]; then
            gmt begin "$results_dir/displacement" png E100 >> "$LOG_FILE" 2>&1
            gmt makecpt -Cpolar -T"$dmin2"/"$dmax2" > disp.cpt 2>> "$LOG_FILE"
            gmt grdimage los_displacement.grd -JX15c -Cdisp.cpt -Baf >> "$LOG_FILE" 2>&1
            gmt colorbar -Cdisp.cpt -DJBC+w12c/0.4c -Baf+l"LOS Displacement (mm)" >> "$LOG_FILE" 2>&1
            gmt end >> "$LOG_FILE" 2>&1
            log_info "位移图: results/displacement.png"
        else
            log_warn "位移数据范围不足或全为零，跳过位移图生成"
        fi
        # 雷达坐标 GeoTIFF（无地理编码时保留）
        if [ "$have_trans" -eq 0 ]; then
            gmt grdconvert los_displacement.grd "$results_dir/displacement.tif=gd:GTiff" >> "$LOG_FILE" 2>&1
            log_info "GeoTIFF(雷达坐标): results/displacement.tif"
        fi
    fi
    
    # 复制关键文件到结果目录
    cp -f los_displacement.grd "$results_dir/" 2>/dev/null
    cp -f unwrap_detrended.grd "$results_dir/" 2>/dev/null
    [ -f los_displacement_ll.grd ] && cp -f los_displacement_ll.grd "$results_dir/" 2>/dev/null
    [ -f unwrap_detrended_ll.grd ] && cp -f unwrap_detrended_ll.grd "$results_dir/" 2>/dev/null
    [ -f corr_ll.grd ] && cp -f corr_ll.grd "$results_dir/" 2>/dev/null
    
    log_info "输出文件已保存到: $results_dir"
}

#-----------------------------------------------------------------------------
# 主函数
#-----------------------------------------------------------------------------
main() {
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  PowerChina-1 DInSAR 数据处理脚本 v${SCRIPT_VERSION}                ║"
    echo "║  支持: DJ1 (电建一号) / BC3 (涪城一号)                      ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    
    START_TIME=$(date +%s)
    
    # 解析参数
    parse_args "$@"
    apply_tier_defaults
    
    # 设置日志
    mkdir -p "$PROJECT_DIR/logs"
    LOG_FILE="$PROJECT_DIR/logs/processing_$(date +%Y%m%d_%H%M%S).log"
    
    log_info "项目目录: $PROJECT_DIR"
    log_info "处理模式: $PROCESS_MODE"
    log_info "解缠模式: $UNWRAP_MODE"
    log_info "降采样: ${DOWNSAMPLE_FACTOR}x"
    
    # 检查系统
    check_system || exit 1
    apply_tier_defaults
    
    # 设置目录
    setup_directories
    
    # 数据发现
    discover_data || exit 1
    
    # 确认开始
    if [ "$INTERACTIVE" = true ]; then
        echo ""
        read -p "开始处理? [Y/n] " confirm
        if [[ "$confirm" =~ ^[Nn]$ ]]; then
            echo "已取消"
            exit 0
        fi
    fi
    
    echo ""
    echo "========================================"
    echo "  开始处理"
    echo "========================================"
    echo ""
    
    # 执行处理流程
    generate_slc || exit 1
    run_alignment || exit 1
    run_interferogram || exit 1
    run_filter || exit 1
    run_unwrap || exit 1
    run_atm_correction || exit 1
    run_detrend || exit 1
    run_displacement || exit 1
    generate_outputs || exit 1
    
    # 完成
    local end_time=$(date +%s)
    local total_time=$((end_time - START_TIME))
    
    echo ""
    echo "========================================"
    echo "  处理完成"
    echo "========================================"
    echo ""
    log_info "总耗时: $((total_time / 60)) 分 $((total_time % 60)) 秒"
    log_info "结果目录: $PROJECT_DIR/results/"
    log_info "日志文件: $LOG_FILE"
    echo ""
}

# 运行主函数
main "$@"
