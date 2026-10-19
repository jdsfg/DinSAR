#!/bin/bash
#=============================================================================
# PowerChina-1 DInSAR 数据处理脚本 v2.5
#
# 支持: DJ1 (电建一号, X波段), BC3 (涪城一号/天仪, C波段), LT1 (陆探一号, L波段)
# 
# 用法:
#   ./run_dinsar.sh                           # 交互式向导
#   ./run_dinsar.sh /path/to/project          # 指定项目目录
#   ./run_dinsar.sh --auto --mode fast /path  # 全自动处理
#
# 作者: PowerChina-1 DInSAR Team
# 日期: 2026-01-28 (v2.5: 2026-10 并入 LT1 分支)
#=============================================================================

set -o pipefail

#-----------------------------------------------------------------------------
# 全局变量
#-----------------------------------------------------------------------------
SCRIPT_VERSION="2.5"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#-----------------------------------------------------------------------------
# 环境配置：把 GMTSAR 目录追加到 PATH 末尾，仅作兜底。
# 已在 PATH 中的命令（如 ~/.bashrc 配置的安装目录）优先，不会被源码目录里的同名程序覆盖。
#-----------------------------------------------------------------------------
_add_path_if_dir() {
    local d="$1"
    [ -d "$d" ] || return 0
    case ":$PATH:" in
        *":$d:"*) ;;
        *) PATH="$PATH:$d" ;;
    esac
}
for _root in "$HOME/gmtsar/gmtsar-master" "$SCRIPT_DIR/gmtsar-master"; do
    _add_path_if_dir "$_root/bin"
    _add_path_if_dir "$_root/gmtsar"
    _add_path_if_dir "$_root/preproc/DJ1_preproc/src_swath"
    _add_path_if_dir "$_root/preproc/LT1_preproc/src"
    _add_path_if_dir "$_root/preproc/LT1_preproc/src_orbit"
done
unset _root
export PATH

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
WAVELENGTH=""  # 将从 PRM 文件自动读取

# 处理参数
PROCESS_MODE="FAST"  # FAST/ STANDARD/ QUALITY
UNWRAP_MODE="A"      # A=稳定模式, B=高精度模式（由 PROCESS_MODE 映射）
ATM_MODE="NONE" # P1-6: 默认关闭全局线性大气校正，防止引入新误差 (支持 NONE, EMPIRICAL, GACOS, ST_FILTER, RAYTRACE)
WEATHER_DIR=""       # 物理射线追踪气象数据目录 (.nc 立方体, 可选)
COREG_MODE="legacy"  # legacy=xcorr+fitoffset (默认); multiscale=多尺度配准 ms_register.py (须显式 --coreg multiscale)
COREG_USER_SET=false # 命令行是否给了 --coreg (新版电建一号 XML+HDF5 在没给时自动改用 multiscale, 见 discover_data)
QA_MODE="on"
INPUT_FORMAT="tiff"  # tiff=旧版 XML+TIFF; h5=新版电建一号 XML+HDF5 (discover_data 里按同名 .h5 识别)
MS_REGISTER_ARGS="${MS_REGISTER_ARGS:-}"  # 额外传给 ms_register.py register 的参数 (高级用法, 如 "--scales 512,256,128")
GACOS_PATH=""        # GACOS 数据目录 (.ztd + .rsc)
REF_PX=""            # 参考点 Range 坐标 (可选，GACOS 模式使用)
REF_PY=""            # 参考点 Azimuth 坐标 (可选，GACOS 模式使用)


if [ -n "${DOWNSAMPLE_FACTOR+x}" ]; then
    DOWNSAMPLE_FACTOR_USER_SET=true
else
    DOWNSAMPLE_FACTOR_USER_SET=false
fi
DOWNSAMPLE_FACTOR="${DOWNSAMPLE_FACTOR:-4}"

# 新更加速模式标志
CROP_MODE=false


SNAPHU_TILES="2 2"
SNAPHU_OVERLAP="200"
SNAPHU_OVERLAP_USER_SET=false
SNAPHU_NPROC=8
SNAPHU_CONF_MODE="auto"
SNAPHU_CONF_PATH=""
SNAPHU_CONF_USED=""
CORR_FILTER_WIDTH=20

# grdtrend 的 -N: 1 = 只减去均值 (不去倾斜), 3 = 平面, 4 = 双线性, 6 = 二次曲面。
# 环境变量 DETREND_ORDER 明确指定的值优先于新版数据的自动选择 (见 discover_data)。
if [ -n "${DETREND_ORDER+x}" ]; then
    DETREND_ORDER_USER_SET=true
else
    DETREND_ORDER_USER_SET=false
fi
DETREND_ORDER="${DETREND_ORDER:-1}"

# 位移图 (displacement.png) 的色标范围: 取 |位移| 的第 N 百分位作为正负限, 100 = 用最大绝对值 (旧行为)。
# 环境变量 DISP_STRETCH_PCT 明确指定的值优先于新版数据的自动选择 (见 discover_data)。
if [ -n "${DISP_STRETCH_PCT+x}" ]; then
    DISP_STRETCH_PCT_USER_SET=true
else
    DISP_STRETCH_PCT_USER_SET=false
fi
DISP_STRETCH_PCT="${DISP_STRETCH_PCT:-100}"

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
    --atm none|empirical|gacos|st_filter|raytrace
                                        大气校正模式: None (默认) / Empirical / GACOS / ST_Filter / Raytrace
    --weather-dir <path> 指定气象场 NetCDF 目录 (--atm raytrace 可选)
    --gacos-path <path>  指定 GACOS 数据目录 (--atm gacos 时必需)
    --coreg legacy|multiscale
                                        配准方式: legacy (默认, xcorr + fitoffset) /
                                        multiscale (多尺度配准 ms_register.py: 不依赖 PRM 里的初始偏移,
                                        配准后用独立点位复核, 未通过则中止, 不会静默退回 legacy)
                                        额外参数用环境变量 MS_REGISTER_ARGS 传给 ms_register.py register
                                        新版电建一号 (XML+HDF5) 不加此项时自动用 multiscale, 其余数据仍默认 legacy
  --qa on|off|strict  产品质量检查 (默认 on；strict 时 QA FAIL 会中止流程)
  --snaphu-conf auto|base|<path>
                                        SNAPHU 配置: auto (默认, 按 PRM 派生) / base (旧配置) / 指定文件
  --snaphu-overlap N                   SNAPHU 分块重叠像元数 (正整数, 默认由处理模式决定)
  --crop             仅处理影像中心 50%x50% 区域 (加速验证用)
  --help             显示此帮助

支持卫星 (从 XML 自动识别): DJ1 / BC3 / LT1 (LT1 自动启用 PRF·DC 归一化与 L 波段配准参数)

环境变量:
  DETREND_ORDER      去轨道残余倾斜的 grdtrend -N: 1=只减均值, 3=平面, 6=二次曲面 (旧数据默认 1; 新版电建一号默认 3)
  DISP_STRETCH_PCT   位移图色标正负限取 |位移| 的第几百分位, 100=最大绝对值 (旧数据默认 100; 新版电建一号默认 98)
  TOPO_RA_DEC        雷达坐标 DEM 相对 SLC 的降采样倍数: 1=整景分辨率, N=粗网格, auto (旧数据默认 1; 新版电建一号默认 auto)
  DOWNSAMPLE_FACTOR  干涉图降采样倍数 (默认 4)
  MS_REGISTER_ARGS   额外传给 ms_register.py register 的参数

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
  ├── raw/              # 放置 XML + TIFF 文件 (新版电建一号: XML + 同名 .h5, 同名 .h5 优先于 .tiff)
  │   ├── image1.xml
  │   ├── image1.tiff
  │   ├── image2.xml
  │   └── image2.tiff
  └── topo/             # 放置 DEM (可选): dem.grd, 或 dem.tif 等 (自动转 EPSG:4326 的 dem.grd)
      └── dem.grd       # 高程须是相对 WGS84 椭球的高 (SRTM/Copernicus 等正高数据要先加大地水准面起伏)

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
            --crop)
                CROP_MODE=true
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
                    RAYTRACE|RAY) ATM_MODE="RAYTRACE" ;;
                    *) error_exit "无效 --atm: $2 (none|empirical|gacos|st_filter|raytrace)" ;;
                esac
                shift 2
                ;;
            --weather-dir)
                WEATHER_DIR="$2"
                shift 2
                ;;
            --gacos-path)
                GACOS_PATH="$2"
                shift 2
                ;;
            --coreg)
                case "${2,,}" in
                    legacy|xcorr) COREG_MODE="legacy" ;;
                    multiscale|ms) COREG_MODE="multiscale" ;;
                    *) error_exit "无效 --coreg: $2 (legacy|multiscale)" ;;
                esac
                COREG_USER_SET=true
                shift 2
                ;;
            --qa)
                case "${2,,}" in
                    on|off|strict) QA_MODE="${2,,}" ;;
                    *) error_exit "无效 --qa: $2 (on|off|strict)" ;;
                esac
                shift 2
                ;;
            --snaphu-conf)
                [ "$#" -ge 2 ] || error_exit "--snaphu-conf 需要 auto|base|<path>"
                case "${2,,}" in
                    auto|base)
                        SNAPHU_CONF_MODE="${2,,}"
                        SNAPHU_CONF_PATH=""
                        ;;
                    *)
                        SNAPHU_CONF_MODE="path"
                        if [[ "$2" = /* ]]; then
                            SNAPHU_CONF_PATH="$2"
                        else
                            local conf_dir
                            conf_dir=$(cd "$(dirname "$2")" 2>/dev/null && pwd) ||
                                error_exit "无效的 SNAPHU 配置路径: $2"
                            SNAPHU_CONF_PATH="$conf_dir/$(basename "$2")"
                        fi
                        ;;
                esac
                shift 2
                ;;
            --snaphu-overlap)
                [ "$#" -ge 2 ] && [[ "$2" =~ ^[1-9][0-9]*$ ]] ||
                    error_exit "无效 --snaphu-overlap: ${2:-} (必须为正整数)"
                SNAPHU_OVERLAP="$2"
                SNAPHU_OVERLAP_USER_SET=true
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
            [ "$SNAPHU_OVERLAP_USER_SET" = true ] || SNAPHU_OVERLAP="250"
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
        # 用户用环境变量 DOWNSAMPLE_FACTOR 明确指定的降采样倍数优先, 不被内存规则覆盖
        [ "$DOWNSAMPLE_FACTOR_USER_SET" = true ] || DOWNSAMPLE_FACTOR=4
    fi
    
    # 检查必要命令（SLC 预处理器 make_slc_* 因卫星而异，在 generate_single_slc 中按需检查）
    local required_cmds=("gmt" "xcorr" "phasediff" "snaphu" "resamp")
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
normalize_satellite_type() {
    local value
    value=$(echo "$1" | tr '[:lower:]' '[:upper:]')
    case "$value" in
        *LT1*|*LUTAN*|*LU_TAN*|*陆探*) echo "LT1" ;;
        *DJ1*|*POWERCHINA*|*DIANJIAN*) echo "DJ1" ;;
        *BC3*|*FC1*|*FUCHENG*|*SPACETY*|*TIANYI*) echo "BC3" ;;
        *) echo "BC3" ;;
    esac
}

# 兼容 LT1 的 xxx.meta.xml -> xxx.tiff 命名
find_matching_tiff() {
    local raw_dir="$1"
    local stem="$2"
    local candidate ext short_stem
    for ext in tiff TIFF tif TIF; do
        candidate="$raw_dir/${stem}.${ext}"
        [ -f "$candidate" ] && { echo "$candidate"; return 0; }
    done
    short_stem="${stem%.meta}"
    if [ "$short_stem" != "$stem" ]; then
        for ext in tiff TIFF tif TIF; do
            candidate="$raw_dir/${short_stem}.${ext}"
            [ -f "$candidate" ] && { echo "$candidate"; return 0; }
        done
    fi
    candidate=$(find "$raw_dir" -maxdepth 1 -type f \( -iname "${stem}*.tif" -o -iname "${stem}*.tiff" -o -iname "${short_stem}*.tif" -o -iname "${short_stem}*.tiff" \) | sort | head -1)
    [ -n "$candidate" ] && echo "$candidate"
    return 0
}

# 新版电建一号 (XML + HDF5): 同名 .h5 才含复数 SLC。同目录里的同名 .tiff 只是 8 bit 幅度预览图, 不能当 SLC 用,
# 所以有 .h5 时优先它。
find_matching_h5() {
    local raw_dir="$1"
    local stem="$2"
    local candidate ext
    for ext in h5 H5 hdf5 HDF5; do
        candidate="$raw_dir/${stem}.${ext}"
        [ -f "$candidate" ] && { echo "$candidate"; return 0; }
    done
    return 0
}

is_h5_input() {
    case "$1" in
        *.h5|*.H5|*.hdf5|*.HDF5) return 0 ;;
        *) return 1 ;;
    esac
}

extract_acquisition_date() {
    local xml="$1"
    local acq_date
    acq_date=$(grep -ioP '(?<=<startTime>|<AcquisitionDate>|<SceneDate>|<CenterTime>|<timeUTC>|<acquisition_start_utc>)[0-9TZ:.\-]+' "$xml" 2>/dev/null | head -1)
    [ -z "$acq_date" ] && acq_date=$(grep -oP '\d{4}-\d{2}-\d{2}' "$xml" 2>/dev/null | head -1)
    [ -z "$acq_date" ] && acq_date=$(basename "$xml" | grep -oP '(?<!\d)(19|20)\d{6}(?!\d)' | head -1)
    echo "$acq_date" | tr -d 'TZ:.-' | cut -c1-8
}

detect_satellite_type() {
    local xml="$1"
    local sat_id
    sat_id=$(grep -ioP '(?<=<SatelliteID>|<SensorID>|<Satellite>|<mission>|<missionInfo>|<satellite_name>)[^<]+' "$xml" 2>/dev/null | head -1)
    [ -z "$sat_id" ] && sat_id=$(basename "$xml")
    normalize_satellite_type "$sat_id"
}

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
        
        # 变量名沿用 tiff: 新版电建一号这里装的是同名 .h5 (见 find_matching_h5)
        local tiff
        tiff=$(find_matching_h5 "$raw_dir" "$stem")
        [ -z "$tiff" ] && tiff=$(find_matching_tiff "$raw_dir" "$stem")
        
        if [ -z "$tiff" ]; then
            log_warn "孤立 XML (无匹配 TIFF/h5): $basename"
            continue
        fi
        
        # 从 XML 内部解析元数据 (不依赖文件名)
        local acq_date=""
        local sat_id=""
        acq_date=$(extract_acquisition_date "$xml")
        sat_id=$(detect_satellite_type "$xml")
        
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

    INPUT_FORMAT="tiff"
    if is_h5_input "$MASTER_TIFF" || is_h5_input "$SLAVE_TIFF"; then
        if ! { is_h5_input "$MASTER_TIFF" && is_h5_input "$SLAVE_TIFF"; }; then
            log_error "主辅影像的数据格式不一致 (一景是 HDF5, 另一景是 TIFF): $MASTER_TIFF / $SLAVE_TIFF"
            return 1
        fi
        INPUT_FORMAT="h5"
    fi
    
    log_info "主影像 (Master): $MASTER_DATE"
    log_info "辅影像 (Slave):  $SLAVE_DATE"
    log_info "卫星类型: $SATELLITE_TYPE"
    log_info "数据格式: $INPUT_FORMAT"

    # 新版产品两景 SLC 的初始偏移常有数百~上千像素 (实测 SM 约 400 px, SP 约 1600x1900 px),
    # legacy 的 xcorr 搜索窗 (±64~256 px) 覆盖不了, 又不会报错, 只会悄悄给出噪声干涉图, 所以默认走 multiscale。
    if [ "$INPUT_FORMAT" = "h5" ] && [ "$COREG_USER_SET" = false ]; then
        COREG_MODE="multiscale"
        log_info "配准方式: multiscale (新版 XML+HDF5 默认; 需要旧流程请加 --coreg legacy)"
    fi

    # 新版产品的轨道常是实时轨道 (米~分米级误差), 在 X 波段 (半波长 1.56 cm) 的干涉图里留下几十 rad 的残余倾斜
    # (实测 SM 像对沿方位向约 11.7 周, 去掉后的残差标准差仅约 1 cm)。grdtrend -N1 只减均值, 倾斜原样保留, 位移图
    # 会被倾斜主导, 所以默认去平面 (-N3); 位移图色标取 98 分位, 否则极少数低相干离群值 (±数百 mm) 会把图冲成一片白。
    if [ "$INPUT_FORMAT" = "h5" ]; then
        if [ "$DETREND_ORDER_USER_SET" = false ]; then
            DETREND_ORDER=3
            log_info "去趋势: 平面 (grdtrend -N3; 新版 XML+HDF5 默认; 环境变量 DETREND_ORDER=1 可退回只减均值)"
        fi
        if [ "$DISP_STRETCH_PCT_USER_SET" = false ]; then
            DISP_STRETCH_PCT=98
            log_info "位移图色标: |位移| 的 98 分位 (新版 XML+HDF5 默认; 环境变量 DISP_STRETCH_PCT=100 可改用最大绝对值)"
        fi
    fi
    
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

find_dj1h5_converter() {
    local c
    for c in "$SCRIPT_DIR/gmtsar-master/preproc/DJ1_preproc/src_swath/make_slc_dj1h5.py" "$HOME/gmtsar/gmtsar-master/preproc/DJ1_preproc/src_swath/make_slc_dj1h5.py"; do
        [ -f "$c" ] && { echo "$c"; return 0; }
    done
    c=$(command -v make_slc_dj1h5.py 2>/dev/null) && { echo "$c"; return 0; }
    return 1
}

# 转换器需要 numpy 和 h5py; 系统 python 里没有 h5py 时, 它自己会去 ~/gmtsar/_tools/pylibs 找。
generate_single_slc_h5() {
    local role="$1"
    local xml="$2"
    local h5="$3"
    local out_stem="$4"
    local converter

    converter=$(find_dj1h5_converter) || error_exit "未找到 make_slc_dj1h5.py (应在 gmtsar-master/preproc/DJ1_preproc/src_swath/ 下)"
    command -v python3 &>/dev/null || error_exit "转换新版电建一号 HDF5 需要 python3 (含 numpy、h5py)"

    log_info "生成${role} SLC ($out_stem, DJ1 XML+HDF5, make_slc_dj1h5.py)..."
    python3 "$converter" "$xml" "$h5" "$out_stem" >> "$LOG_FILE" 2>&1 \
        || error_exit "${role} SLC 生成失败 (make_slc_dj1h5.py; 需要 python 模块 h5py, 详见日志)"
}

generate_single_slc() {
    local role="$1"
    local xml="$2"
    local tiff="$3"
    local out_stem="$4"
    local processor

    if is_h5_input "$tiff"; then
        generate_single_slc_h5 "$role" "$xml" "$tiff" "$out_stem"
        return
    fi

    case "$SATELLITE_TYPE" in
        LT1) processor="make_slc_lt1" ;;
        *)   processor="make_slc_dj1" ;;
    esac

    command -v "$processor" &>/dev/null || error_exit "未找到 $processor，请先编译/安装对应预处理器"

    log_info "生成${role} SLC ($out_stem, $SATELLITE_TYPE, $processor)..."
    if [ "$processor" = "make_slc_lt1" ]; then
        # 第 4 个参数为 SLC 幅度缩放因子
        "$processor" "$xml" "$tiff" "$out_stem" 10.0 >> "$LOG_FILE" 2>&1 || error_exit "${role} SLC 生成失败"
    elif ! "$processor" "$xml" "$tiff" "$out_stem" >> "$LOG_FILE" 2>&1; then
        log_warn "$processor 失败，尝试 make_slc_s1a..."
        make_slc_s1a "$xml" "$tiff" "$out_stem" >> "$LOG_FILE" 2>&1 || error_exit "${role} SLC 生成失败"
    fi
}

# 新版电建一号: 给两景 PRM 追加 SC_vel / earth_radius / SC_height*。要成对做 (从影像沿用主影像的 earth_radius),
# 所以不放在单景转换里 (见 make_slc_dj1h5.py 文档)。calc_dop_orb 末尾的 0 把 fd1 强制成 0 (产品是零多普勒时间),
# 转换器已经写了 fd1=0, 这里滤掉 calc_dop_orb 输出里重复的 fd1 行。须在 SLC 目录里调用。
add_orbit_params_to_prm() {
    command -v calc_dop_orb &>/dev/null || error_exit "未找到 calc_dop_orb (GMTSAR)"
    local m_prm="${MASTER_DATE}.PRM"
    local s_prm="${SLAVE_DATE}.PRM"
    local m_tmp="${MASTER_DATE}.orb.tmp"
    local s_tmp="${SLAVE_DATE}.orb.tmp"
    local re

    log_info "用 calc_dop_orb 补齐 PRM 的轨道参数 (SC_vel / earth_radius / SC_height*)..."
    calc_dop_orb "$m_prm" "$m_tmp" 0 0 >> "$LOG_FILE" 2>&1 || error_exit "calc_dop_orb 失败 (主影像)"
    re=$(awk '/^earth_radius/ {print $3; exit}' "$m_tmp")
    [ -n "$re" ] || error_exit "calc_dop_orb 没有给出主影像的 earth_radius"
    calc_dop_orb "$s_prm" "$s_tmp" "$re" 0 >> "$LOG_FILE" 2>&1 || error_exit "calc_dop_orb 失败 (辅影像)"
    grep -v '^fd1' "$m_tmp" >> "$m_prm"
    grep -v '^fd1' "$s_tmp" >> "$s_prm"
    rm -f "$m_tmp" "$s_tmp"
}

generate_slc() {
    log_step "[1/8] 生成 SLC..."
    
    cd "$PROJECT_DIR/SLC" || error_exit "无法进入 SLC 目录"
    
    generate_single_slc "主影像" "$MASTER_XML" "$MASTER_TIFF" "$MASTER_DATE"
    generate_single_slc "辅影像" "$SLAVE_XML" "$SLAVE_TIFF" "$SLAVE_DATE"
    
    # 核心修复: 自动修复 PRM 文件参数 (earth_radius, SC_vel, SC_height, tie_point)
    # fix_prm.py 依赖天仪/电建 XML 结构，LT1 不适用。优先用项目 SLC 目录里的副本，其次用脚本同目录的。
    # 新版电建一号 (XML+HDF5) 的 XML 结构不同, 也不用 fix_prm.py: 由 calc_dop_orb 按 GMTSAR 惯例从轨道算这些参数。
    local fix_prm_py=""
    if [ "$INPUT_FORMAT" = "h5" ]; then
        add_orbit_params_to_prm
    elif [ "$SATELLITE_TYPE" != "LT1" ]; then
        for c in "$PROJECT_DIR/SLC/fix_prm.py" "$SCRIPT_DIR/fix_prm.py"; do
            [ -f "$c" ] && { fix_prm_py="$c"; break; }
        done
    fi
    if [ -n "$fix_prm_py" ]; then
        log_info "修复主影像 PRM 文件参数 ($fix_prm_py)..."
        python3 "$fix_prm_py" "${MASTER_DATE}.PRM" "$MASTER_XML" "${MASTER_DATE}.LED" >> "$LOG_FILE" 2>&1
        log_info "修复辅影像 PRM 文件参数..."
        python3 "$fix_prm_py" "${SLAVE_DATE}.PRM" "$SLAVE_XML" "${SLAVE_DATE}.LED" >> "$LOG_FILE" 2>&1
    fi
    
    # 动态读取雷达波长 (P0: 波长硬编码修复)
    local prm_wave
    prm_wave=$(grep "radar_wavelength" "${MASTER_DATE}.PRM" 2>/dev/null | awk '{print $3}')
    if [ -n "$prm_wave" ]; then
        WAVELENGTH="$prm_wave"
        log_info "从 PRM 动态读取雷达波长: $WAVELENGTH m"
    else
        WAVELENGTH=0.0311
        log_warn "未找到雷达波长，使用默认值: $WAVELENGTH m"
    fi
    
    log_info "SLC 生成完成"
}

#-----------------------------------------------------------------------------
# LT1 专用预处理：PRF 归一化 + 多普勒中心频率(DC)归一化
# LT1 两景 PRF 可相差 1% 以上，不先归一化会使方位向配准系统性失配。
#-----------------------------------------------------------------------------
run_lt1_preprocess() {
    [ "$SATELLITE_TYPE" = "LT1" ] || return 0

    log_step "[1.5/8] LT1 专用预处理 (PRF/DC 归一化)..."

    cd "$PROJECT_DIR/SLC" || error_exit "无法进入 SLC 目录"

    local master_prm="${MASTER_DATE}.PRM"
    local slave_prm="${SLAVE_DATE}.PRM"
    local prf_m prf_s prf_ratio

    prf_m=$(awk -F= '/^PRF/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$master_prm")
    prf_s=$(awk -F= '/^PRF/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$slave_prm")

    if [ -n "$prf_m" ] && [ -n "$prf_s" ]; then
        prf_ratio=$(awk -v m="$prf_m" -v s="$prf_s" 'BEGIN {r = (m-s)/s; if (r<0) r=-r; printf("%.6f", r)}')
        log_info "PRF: master=$prf_m slave=$prf_s (相对差异=${prf_ratio})"

        if awk -v r="$prf_ratio" 'BEGIN {exit !(r > 0.001)}'; then
            log_info "PRF 差异 >0.1%，重采样 slave SLC 到 master PRF..."
            cp -n "$slave_prm" "${slave_prm}.before_prf_norm"

            python3 - "$prf_m" "$prf_s" "$SLAVE_DATE" >> "$LOG_FILE" 2>&1 << 'PYPRF' || error_exit "LT1 PRF 归一化失败"
import os, re, sys
import numpy as np

prf_m, prf_s, stem = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
prm_text = open(f"{stem}.PRM").read()
n_rng = int(re.search(r"num_rng_bins\s*=\s*(\d+)", prm_text).group(1))
n_samples = n_rng * 2  # int16 实/虚交织

slc = np.memmap(f"{stem}.SLC", dtype=np.int16, mode="r")
n_lines = len(slc) // n_samples
slc = slc[: n_lines * n_samples].reshape(n_lines, n_samples)

ratio = prf_s / prf_m
new_n = int(n_lines / ratio)
src = np.arange(new_n) * ratio
i0 = np.clip(src.astype(np.int64), 0, n_lines - 2)
frac = (src - i0).astype(np.float32)

tmp = f"{stem}.SLC.prfnorm"
chunk = 2048
with open(tmp, "wb") as fo:
    for s in range(0, new_n, chunk):
        e = min(s + chunk, new_n)
        a = slc[i0[s:e]].astype(np.float32)
        b = slc[i0[s:e] + 1].astype(np.float32)
        f = frac[s:e, None]
        out = np.clip(np.round(a * (1 - f) + b * f), -32768, 32767).astype(np.int16)
        fo.write(out.tobytes())
del slc
os.replace(tmp, f"{stem}.SLC")

prm_text = re.sub(r"(PRF\s*=\s*)\S+", lambda m: m.group(1) + repr(prf_m), prm_text)
prm_text = re.sub(r"(num_lines\s*=\s*)\S+", lambda m: m.group(1) + str(new_n), prm_text)
open(f"{stem}.PRM", "w").write(prm_text)
print(f"PRF normalized: {n_lines} -> {new_n} lines, PRF {prf_s} -> {prf_m}")
PYPRF
            log_info "PRF 归一化完成"
        else
            log_info "PRF 差异 <0.1%，无需归一化"
        fi
    fi

    local fd1_m fd1_s
    fd1_m=$(awk -F= '/^fd1/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$master_prm")
    fd1_s=$(awk -F= '/^fd1/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$slave_prm")
    if [ -n "$fd1_m" ] && [ -n "$fd1_s" ]; then
        log_info "DC 归一化: slave fd1=$fd1_s -> master fd1=$fd1_m"
        if command -v update_PRM &>/dev/null; then
            update_PRM "$slave_prm" fd1 "$fd1_m" >> "$LOG_FILE" 2>&1
        else
            sed -i "s/^fd1.*/fd1 = $fd1_m/" "$slave_prm"
        fi
    fi

    log_info "LT1 预处理完成"
}

#-----------------------------------------------------------------------------
# 多尺度配准 (--coreg multiscale): 估计 -> resamp 重采样 -> 独立点位复核 一步完成。
# 输出文件名与 legacy 流程相同 (${SLAVE_DATE}RS.PRM / RS.SLC), 下游步骤不必区分; 未通过复核时直接中止,
# 不退回 legacy, 避免把没通过质量门限的配准结果悄悄交给后面的步骤。
#-----------------------------------------------------------------------------
ms_register_python_engine_requested() {
    local args=" ${MS_REGISTER_ARGS:-} "
    [[ "$args" == *" --engine python "* || "$args" == *" --engine=python "* ]]
}

preflight_ms_register() {
    [ "$COREG_MODE" = "multiscale" ] || return 0
    ms_register_python_engine_requested && return 0

    local candidate
    if [ -n "${MS_REGISTER_BIN:-}" ] && [ -f "$MS_REGISTER_BIN" ] && [ -x "$MS_REGISTER_BIN" ]; then
        return 0
    fi
    candidate=$(type -P ms_register 2>/dev/null) || candidate=""
    if [ -n "$candidate" ] && [ -f "$candidate" ] && [ -x "$candidate" ]; then
        return 0
    fi
    candidate="$SCRIPT_DIR/gmtsar-master/gmtsar/ms_register"
    if [ -f "$candidate" ] && [ -x "$candidate" ]; then
        return 0
    fi

    error_exit "多尺度配准需要原生 ms_register 二进制。请运行 make -C gmtsar-master/gmtsar ms_register，或设置 MS_REGISTER_ARGS=\"--engine python\""
}

find_ms_register() {
    local c
    if c=$(command -v ms_register.py 2>/dev/null); then
        echo "$c"
        return 0
    fi
    for c in "$SCRIPT_DIR/gmtsar-master/gmtsar/csh/ms_register.py" "$HOME/gmtsar/gmtsar-master/gmtsar/csh/ms_register.py"; do
        [ -f "$c" ] && { echo "$c"; return 0; }
    done
    return 1
}

run_alignment_multiscale() {
    local master_prm="$1"
    local slave_prm="$2"
    local slave_rs_prm="${SLAVE_DATE}RS.PRM"
    local slave_rs_slc="${SLAVE_DATE}RS.SLC"
    local ms_register rc

    ms_register=$(find_ms_register) || error_exit "未找到 ms_register.py (应在 gmtsar-master/gmtsar/csh/ 下)"
    command -v python3 &>/dev/null || error_exit "多尺度配准需要 python3 (含 numpy)"

    log_info "运行多尺度配准 (ms_register.py register)..."
    # --prm-slc-file: resamp 输出的 PRM 沿用从影像原来的 SLC_file, run_interferogram 会把 RS.SLC 链接成这个名字。
    # MS_REGISTER_ARGS 故意不加引号, 以便传入多个参数。
    # shellcheck disable=SC2086
    python3 "$ms_register" register "$master_prm" "$slave_prm" \
        --out-slc "$slave_rs_slc" --out-prm "$slave_rs_prm" --prm-slc-file "${SLAVE_DATE}.SLC" \
        --est-prm "${SLAVE_DATE}_est.PRM" \
        --report coreg_report.json --log coreg_offsets.log --points coreg_points.dat \
        $MS_REGISTER_ARGS >> "$LOG_FILE" 2>&1
    rc=$?
    if [ $rc -ne 0 ] || [ ! -f "$slave_rs_prm" ] || [ ! -f "$slave_rs_slc" ]; then
        error_exit "多尺度配准未通过 (退出码 $rc), 详见 $PROJECT_DIR/SLC/coreg_report.json 和日志; 不会自动退回 legacy 配准"
    fi
    log_info "多尺度配准完成: $slave_rs_prm (报告: SLC/coreg_report.json)"
}

run_alignment() {
    log_step "[2/8] 影像配准..."
    
    cd "$PROJECT_DIR/SLC" || error_exit "无法进入 SLC 目录"
    
    local master_prm="${MASTER_DATE}.PRM"
    local slave_prm="${SLAVE_DATE}.PRM"
    
    if [ ! -f "$master_prm" ] || [ ! -f "$slave_prm" ]; then
        error_exit "PRM 文件不存在"
    fi

    if [ "$COREG_MODE" = "multiscale" ]; then
        run_alignment_multiscale "$master_prm" "$slave_prm"
        log_info "配准完成"
        return 0
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

    # LT1 专用：L 波段偏移量大，先按 near_range 差估计粗距离向偏移作为搜索中心，再用大窗口 xcorr。
    # 符号：rshift = (nr_master - nr_slave) / 距离向像元间距（与 FFT 互相关输出方向相反，见 迭代.md 10b）
    if [ "$SATELLITE_TYPE" = "LT1" ]; then
        local nr_m nr_s rsr coarse_rshift
        nr_m=$(awk -F= '/near_range/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$master_prm")
        nr_s=$(awk -F= '/near_range/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$slave_prm")
        rsr=$(awk -F= '/rng_samp_rate/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$master_prm")
        if [ -n "$nr_m" ] && [ -n "$nr_s" ] && [ -n "$rsr" ]; then
            coarse_rshift=$(awk -v nm="$nr_m" -v ns="$nr_s" -v rs="$rsr" \
                'BEGIN {c=299792458.0; pix=c/(2*rs); printf("%d", (nm-ns)/pix)}')
            log_info "LT1 粗距离向偏移 rshift=$coarse_rshift (nr_m=$nr_m, nr_s=$nr_s)"
            command -v update_PRM &>/dev/null && update_PRM "$slave_prm" rshift "$coarse_rshift" >> "$LOG_FILE" 2>&1
        else
            log_warn "LT1: 无法读取 near_range/rng_samp_rate，跳过粗偏移估计"
        fi
        xsearch=512
        ysearch=512
        nx_corr=30
        ny_corr=30
        threshold_snr=20
        log_info "LT1: 使用 L 波段配准参数 (xsearch=$xsearch ysearch=$ysearch nx=$nx_corr ny=$ny_corr)"
    fi

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
    
    # 核心修复 [P0-1]: 执行 resamp 共配准重采样 (没有重采样会极大降低相干性)
    log_info "执行 resamp 共配准重采样..."
    local slave_rs_prm="${SLAVE_DATE}RS.PRM"
    local slave_rs_slc="${SLAVE_DATE}RS.SLC"
    resamp "$master_prm" "$slave_prm" "$slave_rs_prm" "$slave_rs_slc" 4 >> "$LOG_FILE" 2>&1
    if [ ! -f "$slave_rs_prm" ] || [ ! -f "$slave_rs_slc" ]; then
        error_exit "resamp 失败 (无法生成 $slave_rs_prm)"
    fi
    log_info "重采样完成: $slave_rs_prm"
    
    log_info "配准完成"
}

# 把 $DEM_PATH 或 topo/ 下的 DEM (GeoTIFF 等) 重投影成 EPSG:4326 的 topo/dem.grd。
# 没有 DEM 返回 0 (之后只得到含地形相位的干涉图); 有 DEM 却没能转成 dem.grd 返回 1。
# 新版数据要赶在生成干涉图之前做: run_interferogram 靠 dem.grd 去地形相位, 此时还没有的话只能得到含地形相位的干涉图。
prepare_dem() {
    local dem_grd="$PROJECT_DIR/topo/dem.grd"
    [ -f "$dem_grd" ] && return 0

    local src_dem=""
    if [ -n "$DEM_PATH" ] && [ -f "$DEM_PATH" ]; then
        src_dem="$DEM_PATH"
    else
        src_dem=$(ls "$PROJECT_DIR"/topo/*.{tif,TIF,tiff,TIFF} 2>/dev/null | head -n 1)
    fi
    [ -n "$src_dem" ] || return 0

    log_info "检测到 DEM: $src_dem，重投影为 EPSG:4326 (topo/dem.grd)..."
    mkdir -p "$PROJECT_DIR/topo"
    local dem_tmp="$PROJECT_DIR/topo/dem_tmp.nc"
    rm -f "$dem_tmp"
    gdalwarp -t_srs EPSG:4326 -r bilinear -of netCDF "$src_dem" "$dem_tmp" >> "$LOG_FILE" 2>&1
    [ -f "$dem_tmp" ] || return 1
    gmt grdconvert "$dem_tmp" "$dem_grd" >> "$LOG_FILE" 2>&1 || cp "$dem_tmp" "$dem_grd"
    rm -f "$dem_tmp"
    [ -f "$dem_grd" ]
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
    if [ -f "$PROJECT_DIR/SLC/${SLAVE_DATE}RS.PRM" ]; then
        ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}RS.PRM" "${SLAVE_DATE}.PRM"
        ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}RS.SLC" "${SLAVE_DATE}.SLC"
        log_info "使用已重采样的从影像进行干涉"
    else
        ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}.PRM" .
        ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}.SLC" .
    fi
    ln -sf "$PROJECT_DIR/SLC/${MASTER_DATE}.LED" . 2>/dev/null
    ln -sf "$PROJECT_DIR/SLC/${SLAVE_DATE}.LED" . 2>/dev/null

    # 兼容 GMTSAR 标准流程：在 phasediff 前补齐基线参数
    if command -v SAT_baseline &>/dev/null; then
        SAT_baseline "${MASTER_DATE}.PRM" "${SLAVE_DATE}.PRM" | tail -n 9 >> "${SLAVE_DATE}.PRM" 2>/dev/null || true
        SAT_baseline "${MASTER_DATE}.PRM" "${MASTER_DATE}.PRM" | grep height >> "${MASTER_DATE}.PRM" 2>/dev/null || true
    fi
    
    # 核心修复: 生成并去除地形相位
    # 新版数据在这里先把 GeoTIFF 等 DEM 转成 dem.grd (旧流程的 DEM 转换仍留在地理编码阶段, 这里不改)
    if [ "$INPUT_FORMAT" = "h5" ]; then
        prepare_dem || error_exit "DEM 转换失败 (gdalwarp / gmt grdconvert, 详见日志), 无法去除地形相位"
    fi
    local dem_grd="$PROJECT_DIR/topo/dem.grd"
    local topo_arg=""
    if [ -f "$dem_grd" ]; then
        if [ ! -f dem_ra.grd ]; then
            log_info "生成雷达坐标 DEM (dem_ra.grd) 以去除地形相位..."
            # 整景分辨率的 dem_ra.grd 与 phasediff 的 real/imag 同样大 (新版 SM 整景 1.45 GB), 8~16 GB 内存的机器上
            # 会把内存挤爆。phasediff 允许 topo 网格比 SLC 粗 (SLC 尺寸须是它的整数倍), 新版产品用 4x 粗网格 (auto: 4/2/1 里能整除的最大值)。
            # 其余数据保持原行为 (1 = 整景分辨率)。
            local topo_dec=1
            [ "$INPUT_FORMAT" = "h5" ] && topo_dec="auto"
            bash "$SCRIPT_DIR/fast_dem2topo_ra.sh" "${MASTER_DATE}.PRM" "$dem_grd" "$LOG_FILE" "${TOPO_RA_DEC:-$topo_dec}"
        fi
        if [ -f dem_ra.grd ]; then
            # 核心修复: 替换 dem_ra.grd 中的 NaN 值为 0，防止 phasediff.c 的 calc_average_topo 算得 NaN 导致全图变 NaN/0
            log_info "清洗雷达坐标 DEM (dem_ra.grd) 以免除 NaN 引起的计算链崩溃..."
            gmt grdmath dem_ra.grd 0 AND = dem_ra_clean.grd >> "$LOG_FILE" 2>&1
            mv dem_ra_clean.grd dem_ra.grd
            topo_arg="-topo dem_ra.grd"
            log_info "将去除地形相位 (Differential InSAR 模式)"
        elif [ "$INPUT_FORMAT" = "h5" ]; then
            # 给了 DEM 却没能去地形相位, 悄悄往下走只会得到带地形条纹的干涉图, 所以新版数据在这里直接中止
            error_exit "雷达坐标 DEM (dem_ra.grd) 生成失败, 无法去除地形相位 (详见日志)"
        else
            log_warn "雷达坐标 DEM (dem_ra.grd) 生成失败, 干涉图将含地形相位 (详见日志)"
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
        local phasefilt_cmd=("phasefilt" "-imag" "$IMAG_FILE" "-real" "$REAL_FILE" "-alpha" "0.8" "-psize" "32")
        log_info "执行: ${phasefilt_cmd[*]}"
        "${phasefilt_cmd[@]}" >> "$LOG_FILE" 2>&1
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
        gmt grdfilter "$REAL_FILE" -Fg${CORR_FILTER_WIDTH}/${CORR_FILTER_WIDTH} -D0 -Greal_mean.grd >> "$LOG_FILE" 2>&1
        gmt grdfilter "$IMAG_FILE" -Fg${CORR_FILTER_WIDTH}/${CORR_FILTER_WIDTH} -D0 -Gimag_mean.grd >> "$LOG_FILE" 2>&1
        gmt grdmath real_mean.grd imag_mean.grd HYPOT = num_amp.grd
        # 2. 对复数干涉图的绝对振幅进行空间平滑，计算分母
        gmt grdmath "$REAL_FILE" "$IMAG_FILE" HYPOT = raw_amp.grd
        gmt grdfilter raw_amp.grd -Fg${CORR_FILTER_WIDTH}/${CORR_FILTER_WIDTH} -D0 -Gden_amp.grd >> "$LOG_FILE" 2>&1
        # 3. 计算相干系数 (防止除零，并将上限截断到 1.0，安全替换 NaN 为 0)
        gmt grdmath num_amp.grd den_amp.grd DIV 1.0 MIN 0 AND = corr.grd
    else
        FILT_PHASE="$PHASE_FILE"
        log_warn "phasefilt 不可用，跳过滤波"
    fi
    
    log_info "滤波完成"
}

find_snaphu_conf() {
    local candidate
    local conf_candidates=(
        "/usr/local/GMTSAR/share/gmtsar/snaphu/config/snaphu.conf.brief"
        "$SCRIPT_DIR/gmtsar-master/share/gmtsar/snaphu/config/snaphu.conf.brief"
        "$SCRIPT_DIR/dinsar_deploy_v2.2_offline/gmtsar-master/share/gmtsar/snaphu/config/snaphu.conf.brief"
    )
    for candidate in "${conf_candidates[@]}"; do
        if [ -f "$candidate" ]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

run_unwrap() {
    log_step "[5/8] 相位解缠 (模式 $UNWRAP_MODE)..."
    
    cd "$INTF_DIR" || error_exit "无法进入干涉目录"

    case "$SNAPHU_CONF_MODE" in
        auto)
            local base_conf master_prm conf_summary
            base_conf=$(find_snaphu_conf) || error_exit "未找到 snaphu.conf.brief，无法生成自动配置"
            master_prm="$PROJECT_DIR/SLC/${MASTER_DATE}.PRM"
            [ -f "$master_prm" ] || error_exit "自动 SNAPHU 配置缺少主影像 PRM: $master_prm"
            SNAPHU_CONF_USED="$INTF_DIR/snaphu_auto.conf"
            conf_summary=$(python3 "$SCRIPT_DIR/snaphu_conf_gen.py" \
                --base "$base_conf" \
                --prm "$master_prm" \
                --downsample "$DOWNSAMPLE_FACTOR" \
                --corr-filter-width "$CORR_FILTER_WIDTH" \
                --out "$SNAPHU_CONF_USED" 2>&1) ||
                error_exit "SNAPHU 自动配置生成失败: $conf_summary"
            log_info "SNAPHU 配置(auto): $SNAPHU_CONF_USED"
            log_info "SNAPHU 派生参数: $conf_summary"
            ;;
        base)
            SNAPHU_CONF_USED=$(find_snaphu_conf || true)
            if [ -n "$SNAPHU_CONF_USED" ]; then
                log_info "SNAPHU 配置(base): $SNAPHU_CONF_USED"
            else
                log_warn "未找到 snaphu.conf.brief，使用 snaphu 默认配置"
            fi
            ;;
        path)
            SNAPHU_CONF_USED="$SNAPHU_CONF_PATH"
            [ -n "$SNAPHU_CONF_USED" ] && [ -r "$SNAPHU_CONF_USED" ] ||
                error_exit "SNAPHU 配置不可读: ${SNAPHU_CONF_USED:-<empty>}"
            log_info "SNAPHU 配置(path): $SNAPHU_CONF_USED"
            ;;
        *)
            error_exit "未知 SNAPHU 配置模式: $SNAPHU_CONF_MODE"
            ;;
    esac

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
    
    if [ "$CROP_MODE" = true ]; then
        log_info "[加速模式] 裁剪影像中心 50%x50% 区域以加速测试验证..."
        local cx=$((ncols / 2))
        local cy=$((nrows / 2))
        local dx=$((ncols / 4))
        local dy=$((nrows / 4))
        local xmin=$((cx - dx)); [ $xmin -lt 0 ] && xmin=0
        local xmax=$((cx + dx)); [ $xmax -gt $ncols ] && xmax=$ncols
        local ymin=$((cy - dy)); [ $ymin -lt 0 ] && ymin=0
        local ymax=$((cy + dy)); [ $ymax -gt $nrows ] && ymax=$nrows
        
        gmt grdcut "$FILT_PHASE" -R$xmin/$xmax/$ymin/$ymax -Gphasefilt_crop.grd >> "$LOG_FILE" 2>&1
        FILT_PHASE="phasefilt_crop.grd"
        phase_for_unwrap="$FILT_PHASE"
        ncols=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $10}')
        nrows=$(gmt grdinfo "$FILT_PHASE" -C | awk '{print $11}')
        
        if [ -f "$corr_source" ]; then
            gmt grdcut "$corr_source" -R$xmin/$xmax/$ymin/$ymax -Gcorr_crop.grd >> "$LOG_FILE" 2>&1
            corr_source="corr_crop.grd"
        fi
    fi

    if [ -f "$corr_source" ]; then
        log_info "跳过生成相干性硬掩膜(P1-4)，交由 snaphu 统计代价模型自适应处理低相干区"
    fi

    # 转换为二进制格式并安全替换 NaN 值为 0 (防止 snaphu 崩溃)
    log_info "准备 snaphu 相位与相干性双输入阵列..."
    gmt grdmath "$phase_for_unwrap" 0 AND = phase_no_nan.grd >> "$LOG_FILE" 2>&1
    gmt grd2xyz phase_no_nan.grd -ZTLf > phase.bin 2>> "$LOG_FILE"
    if [ -f "$corr_source" ]; then
        gmt grdmath "$corr_source" 0 AND = corr_no_nan.grd >> "$LOG_FILE" 2>&1
        gmt grd2xyz corr_no_nan.grd -ZTLf > corr.bin 2>> "$LOG_FILE"
    fi
    
    # 构建 snaphu 命令 (P0-2: 添加 -s 强制开启 SMOOTH 形变代价模式)
    local snaphu_cmd=("snaphu" "phase.bin" "$ncols" "-o" "unwrap.bin" "-s")
    if [ -n "$SNAPHU_CONF_USED" ]; then
        snaphu_cmd+=("-f" "$SNAPHU_CONF_USED")
    else
        log_warn "SNAPHU 配置未设置，使用 snaphu 默认配置"
    fi
    
    if [ "$UNWRAP_MODE" = "A" ]; then
        # 稳定模式: 分块 + 多线程
        local tiles=($SNAPHU_TILES)
        local overlap=$SNAPHU_OVERLAP
        if [ "$CROP_MODE" = true ]; then
            tiles=(1 1)
            overlap=0
        fi
        snaphu_cmd+=("--tile" "${tiles[@]}" "$overlap" "$overlap" "--nproc" "$SNAPHU_NPROC")
        log_info "稳定模式: --tile ${tiles[*]} overlap=${overlap} --nproc $SNAPHU_NPROC"
    else
        # 高精度模式: 整图单线程
        log_info "高精度模式: 整图解缠 (需要较长时间)"
    fi
    
    snaphu_cmd+=("-d")
    
    # 核心优化4：将相干性输入 SNAPHU 以构建高精度成本流网络 (Cost Graph)
    if [ -f corr.bin ]; then
        snaphu_cmd+=("-c" "corr.bin")
        log_info "启用基于相关性的统计成本函数模型 (-c corr.bin)"
    fi
    
    log_info "运行 snaphu..."
    log_info "命令: ${snaphu_cmd[*]}"
    
    # 记录开始时间
    local unwrap_start=$(date +%s)
    
    # 运行 snaphu
    "${snaphu_cmd[@]}" >> "$LOG_FILE" 2>&1
    
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

find_atm_script() {
    local name="$1"
    local p="$SCRIPT_DIR/gmtsar-master/gmtsar/csh/$name"
    if [ -f "$p" ]; then
        echo "$p"
        return 0
    fi
    p="$(command -v "$name" 2>/dev/null)"
    if [ -n "$p" ] && [ -f "$p" ]; then
        echo "$p"
        return 0
    fi
    return 1
}

run_atm_correction() {
    log_step "[5.5/8] 大气校正 (${ATM_MODE})..."

    if [ "$ATM_MODE" = "NONE" ]; then
        if [ -e "$INTF_DIR/unwrap_pre_atm.grd" ] || [ -L "$INTF_DIR/unwrap_pre_atm.grd" ] ||
           [ -e "$INTF_DIR/atm_correction.state" ] || [ -L "$INTF_DIR/atm_correction.state" ]; then
            log_error "已有大气校正检查点，NONE 不会自动还原原始相位；请在独立工程中重新处理"
            return 1
        fi
        log_info "大气校正被禁用 (ATM_MODE=NONE)"
        return 0
    fi

    cd "$INTF_DIR" || error_exit "无法进入干涉目录"
    if { [ -e unwrap_pre_atm.grd ] || [ -L unwrap_pre_atm.grd ] || [ -e atm_correction.state ] || [ -L atm_correction.state ]; } &&
       { [ ! -f unwrap_pre_atm.grd ] || [ ! -f atm_correction.state ] || [ ! -f unwrap.grd ]; }; then
        log_error "已有大气校正记录但检查点不完整，保留全部文件；请在独立工程中重新处理"
        return 1
    fi
    if [ ! -f unwrap.grd ]; then
        log_warn "未找到解缠结果，跳过大气校正"
        return 0
    fi

    local dem_grd="$PROJECT_DIR/topo/dem.grd"
    local atm_script script_name
    case "$ATM_MODE" in
        EMPIRICAL)
            script_name=Re_atm.csh
            if [ ! -f "$dem_grd" ]; then
                log_warn "未找到 topo/dem.grd，无法进行经验大气校正"
                [ ! -f unwrap_pre_atm.grd ] || return 1
                return 0
            fi
            if [ ! -f dem_ra.grd ]; then
                bash "$SCRIPT_DIR/fast_dem2topo_ra.sh" "${MASTER_DATE}.PRM" "$dem_grd" "$LOG_FILE" >> "$LOG_FILE" 2>&1 || return 1
            fi
            [ -s dem_ra.grd ] || { log_error "无法获取雷达坐标 DEM"; return 1; }
            ;;
        GACOS)
            script_name=make_gacos_correction.csh
            if [ -z "$GACOS_PATH" ] || [ ! -d "$GACOS_PATH" ] || [ ! -f "$dem_grd" ]; then
                log_error "GACOS 模式需要有效的 --gacos-path 和 topo/dem.grd"
                return 1
            fi
            GACOS_PATH=$(realpath -- "$GACOS_PATH") || return 1
            generate_trans_dat && [ -s trans.dat ] || { log_error "GACOS 模式需要有效的 trans.dat"; return 1; }
            if [ -z "$REF_PX" ] || [ -z "$REF_PY" ]; then
                local center
                center=$(gmt grdinfo unwrap.grd -C | awk '{print ($2+$3)/2, ($4+$5)/2}') || return 1
                read -r REF_PX REF_PY <<< "$center"
                [ -n "$REF_PX" ] && [ -n "$REF_PY" ] || return 1
                log_info "使用网格坐标中心作为参考点: $REF_PX, $REF_PY (未验证稳定性)"
            fi
            ;;
        RAYTRACE)
            script_name=Raytrace_atm.csh
            if [ ! -f "$dem_grd" ]; then
                log_error "RAYTRACE 模式需要有效的 topo/dem.grd"
                return 1
            fi
            generate_trans_dat && [ -s trans.dat ] || { log_error "RAYTRACE 模式需要有效的 trans.dat"; return 1; }
            ;;
        ST_FILTER) script_name=ST_Filter.csh ;;
        *) log_error "未知大气模式: $ATM_MODE"; return 1 ;;
    esac
    atm_script=$(find_atm_script "$script_name") || { log_error "未找到 $script_name"; return 1; }
    atm_script=$(realpath -- "$atm_script") || return 1

    local source_phase="$INTF_DIR/unwrap.grd" f signature output_hash
    [ -f unwrap_pre_atm.grd ] && source_phase="$INTF_DIR/unwrap_pre_atm.grd"
    local inputs=("$source_phase" "$atm_script") prms=("$INTF_DIR"/*.PRM)
    local gacos_prm gacos_led st_pixel_m
    if [ "$ATM_MODE" = "EMPIRICAL" ]; then
        inputs+=("$INTF_DIR/dem_ra.grd")
    elif [ "$ATM_MODE" = "GACOS" ]; then
        inputs+=("$dem_grd" "$INTF_DIR/trans.dat")
        for f in "$MASTER_DATE" "$SLAVE_DATE"; do
            inputs+=("$GACOS_PATH/$f.ztd" "$GACOS_PATH/$f.ztd.rsc")
        done
        gacos_prm=${prms[0]}
        [ -f "$gacos_prm" ] || { log_error "GACOS 缺少 PRM"; return 1; }
        local led candidate
        led=$(awk '/led_file/ {print $3; exit}' "$gacos_prm") || return 1
        [ -n "$led" ] || { log_error "GACOS PRM 缺少 led_file"; return 1; }
        if [[ "$led" = /* ]]; then
            [ ! -f "$led" ] || gacos_led="$led"
        else
            for candidate in "$INTF_DIR/$led" "$PROJECT_DIR/SLC/$led" "$PROJECT_DIR/raw/$led" "$PROJECT_DIR"/*/raw/"$led"; do
                if [ -f "$candidate" ]; then
                    gacos_led=$(realpath -- "$candidate") || return 1
                    break
                fi
            done
        fi
        [ -n "$gacos_led" ] || { log_error "GACOS 找不到 PRM 引用的轨道: $led"; return 1; }
        inputs+=("$gacos_prm" "$gacos_led")
    elif [ "$ATM_MODE" = "RAYTRACE" ]; then
        inputs+=("$dem_grd" "$INTF_DIR/trans.dat")
        for f in "${prms[@]}"; do
            [ ! -f "$f" ] || inputs+=("$f")
        done
        if [ -n "$WEATHER_DIR" ] && [ -d "$WEATHER_DIR" ]; then
            while IFS= read -r -d '' f; do
                inputs+=("$f")
            done < <(find "$WEATHER_DIR" -maxdepth 2 -type f \( -name "*.nc" -o -name "*.grib2" \) -print0 | sort -z)
        elif [ -d "$PROJECT_DIR/weather" ]; then
            while IFS= read -r -d '' f; do
                inputs+=("$f")
            done < <(find "$PROJECT_DIR/weather" -maxdepth 2 -type f \( -name "*.nc" -o -name "*.grib2" \) -print0 | sort -z)
        fi
    else
        # ST_FILTER 只读取首个 PRM 的像元尺度；地理编码会重写其余字段。
        st_pixel_m=${ST_PIXEL_M-3.0}
        if [ -f "${prms[0]}" ]; then
            local fs prf vel
            fs=$(awk '/rng_samp_rate/ {print $3; exit}' "${prms[0]}") || return 1
            prf=$(awk '/^PRF/ {print $3; exit}' "${prms[0]}") || return 1
            vel=$(awk '/SC_vel/ {print $3; exit}' "${prms[0]}") || return 1
            if [ -n "$fs" ] && [ -n "$prf" ] && [ -n "$vel" ]; then
                st_pixel_m=$(awk -v fs="$fs" -v prf="$prf" -v vel="$vel" 'BEGIN {printf "%.3f", (299792458/(2*fs)+vel/prf)/2}') || return 1
            fi
        fi
        if [ -d "$PROJECT_DIR/intf" ]; then
            while IFS= read -r -d '' f; do
                [ "$f" -ef "$source_phase" ] || [ "$f" -ef "$INTF_DIR/unwrap.grd" ] || inputs+=("$f")
            done < <(find "$PROJECT_DIR/intf" -name unwrap.grd -print0 | sort -z)
        fi
    fi
    signature=$(
        {
            printf '%s\0' "$ATM_MODE" "$REF_PX" "$REF_PY" "$st_pixel_m" "$GACOS_PATH" "$WEATHER_DIR"
            for f in "${inputs[@]}"; do
                [ -f "$f" ] || exit 1
                sha256sum < "$f" || exit 1
            done
        } | sha256sum
    ) || { log_error "大气校正输入缺失或无法读取"; return 1; }
    signature=${signature%% *}

    if [ -f unwrap_pre_atm.grd ]; then
        local state=()
        [ ! -f atm_correction.state ] || mapfile -t state < atm_correction.state
        output_hash=$(sha256sum < unwrap.grd) || return 1
        if [ "${state[0]:-}" != "$ATM_MODE" ] || [ "${state[1]:-}" != "$signature" ] || [ "${state[2]:-}" != "${output_hash%% *}" ]; then
            log_error "大气校正检查点缺少来源记录或输入/模式/结果已变；保留全部文件，请在独立工程中重新处理"
            return 1
        fi
        log_info "大气校正输入与输出校验通过，跳过重复校正 (Checkpoint)"
        return 0
    fi
    if [ -f unwrap_detrended.grd ]; then
        log_error "已存在去趋势结果，不能在其上游追加校正后继续使用旧缓存；请在独立工程中重新处理"
        return 1
    fi

    # 各校正脚本会覆写固定临时文件名，必须隔离旧产物和原始输入。
    local work
    work=$(mktemp -d "$PROJECT_DIR/debug/.atm-work.XXXXXX") || return 1
    if ! (
        trap 'rm -rf -- "$work"' EXIT
        local stage="$work/$(basename "$INTF_DIR")"
        mkdir "$stage" || exit 1
        cp -- "$source_phase" "$stage/unwrap.grd" || exit 1
        for f in "$INTF_DIR"/*.PRM "$INTF_DIR"/*.LED; do
            [ ! -f "$f" ] || cp -L -- "$f" "$stage/" || exit 1
        done
        if [ "$ATM_MODE" = "GACOS" ]; then
            cp -- "$gacos_led" "$stage/.atm-orbit.LED" || exit 1
            awk '$1 == "led_file" {$3 = ".atm-orbit.LED"} {print}' "$gacos_prm" > "$stage/${gacos_prm##*/}" || exit 1
        fi
        for f in dem_ra.grd trans.dat; do
            [ ! -f "$INTF_DIR/$f" ] || ln -s -- "$INTF_DIR/$f" "$stage/$f" || exit 1
        done
        cd "$stage" || exit 1
        case "$ATM_MODE" in
            EMPIRICAL) csh -f "$atm_script" unwrap.grd dem_ra.grd unwrap_atm.grd || exit 1 ;;
            GACOS)
                csh -f "$atm_script" "$stage" "$GACOS_PATH" "$REF_PX" "$REF_PY" "$dem_grd" || exit 1
                cp -- unwrap_gacos_corrected.grd unwrap_atm.grd || exit 1
                ;;
            ST_FILTER) csh -f "$atm_script" "$source_phase" unwrap_atm.grd "$PROJECT_DIR" || exit 1 ;;
            RAYTRACE) csh -f "$atm_script" "$source_phase" unwrap_atm.grd "$PROJECT_DIR" "$stage" "$WEATHER_DIR" || exit 1 ;;
        esac
        python3 - "$source_phase" unwrap_atm.grd <<'PY'
import sys
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
a, b = [gdal.Open(p) for p in sys.argv[1:]]
if (a.RasterXSize, a.RasterYSize, a.RasterCount) != (b.RasterXSize, b.RasterYSize, b.RasterCount) or a.RasterCount != 1:
    raise ValueError("Atmospheric result dimensions differ from input")
if not np.allclose(a.GetGeoTransform(), b.GetGeoTransform(), rtol=0, atol=1e-9):
    raise ValueError("Atmospheric result grid is shifted or resampled")
pa, pb = a.GetProjection(), b.GetProjection()
if bool(pa) != bool(pb) or (pa and not osr.SpatialReference(wkt=pa).IsSame(osr.SpatialReference(wkt=pb))):
    raise ValueError("Atmospheric result coordinate reference differs from input")
count, sum_squared = 0, 0.0
for row in range(0, a.RasterYSize, 256):
    rows = min(256, a.RasterYSize - row)
    arrays = [d.GetRasterBand(1).ReadAsArray(0, row, d.RasterXSize, rows).astype(np.float64) for d in (a, b)]
    masks = []
    for d, v in zip((a, b), arrays):
        mask = np.isfinite(v)
        nodata = d.GetRasterBand(1).GetNoDataValue()
        if nodata is not None:
            mask &= v != nodata
        masks.append(mask)
    if not np.array_equal(*masks):
        raise ValueError("Atmospheric result changes the valid-data footprint")
    delta = arrays[1][masks[0]] - arrays[0][masks[0]]
    count += delta.size
    sum_squared += float(np.dot(delta, delta))
if not count:
    raise ValueError("Atmospheric result has no valid pixels")
print("Validated atmospheric result: %d pixels; correction RMS %.9g rad (not an accuracy score)" % (count, np.sqrt(sum_squared / count)))
PY
        [ $? -eq 0 ] || exit 1
        cp -- unwrap_atm.grd "$INTF_DIR/unwrap_atm.grd" || exit 1
        if [ "$ATM_MODE" = "GACOS" ]; then
            cp -- unwrap_gacos_corrected.grd "$INTF_DIR/unwrap_gacos_corrected.grd" || exit 1
        fi
        cp -- "$source_phase" "$INTF_DIR/unwrap_pre_atm.grd" || exit 1
        mv -- unwrap_atm.grd "$INTF_DIR/unwrap.grd" || exit 1
        output_hash=$(sha256sum < "$INTF_DIR/unwrap.grd") || exit 1
        printf '%s\n' "$ATM_MODE" "$signature" "${output_hash%% *}" > "$INTF_DIR/atm_correction.state" || exit 1
    ) >> "$LOG_FILE" 2>&1; then
        log_error "大气校正执行或结果校验失败；未接受未经校验的输出，详见日志"
        return 1
    fi
    log_info "大气校正输出已校验并保存 ($ATM_MODE)"
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
    # LT1 (L 波段, 导师反馈): 使用 N5 曲面，其余卫星沿用 DETREND_ORDER
    if [ "$SATELLITE_TYPE" = "LT1" ]; then
        DETREND_ORDER=5
        log_info "LT1: 使用 N5 曲面去趋势"
    fi
    log_info "拟合轨道残余趋势面 (grdtrend -N$DETREND_ORDER; 1=只减均值, 3=平面, 6=二次曲面)..."
    
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
Processing Parameters (v${SCRIPT_VERSION})
============================
Date: $(date)
Master: $MASTER_DATE
Slave: $SLAVE_DATE
Satellite: $SATELLITE_TYPE
Wavelength: $WAVELENGTH m
Detrend Order: $DETREND_ORDER
Process Mode: $PROCESS_MODE
Coregistration: $COREG_MODE
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
        local py_geocode="$SCRIPT_DIR/geocode_gdal.py"
        # 整景 SLC 的尺寸 (主影像 PRM): 输入网格被降采样过时, geocode_gdal.py 靠它把 trans.dat 的坐标缩放到输入网格上
        local full_ncols full_nrows
        full_ncols=$(awk '$1 == "num_rng_bins" {print $3; exit}' "${MASTER_DATE}.PRM" 2>/dev/null)
        full_nrows=$(awk '$1 == "num_lines" {print $3; exit}' "${MASTER_DATE}.PRM" 2>/dev/null)
        local full_dims=""
        [ -n "$full_ncols" ] && [ -n "$full_nrows" ] && full_dims="$full_ncols $full_nrows"
        if [ -f unwrap_detrended.grd ]; then
            python3 "$py_geocode" trans.dat unwrap_detrended.grd "$results_dir/unwrap_detrended.tif" 0.0001 $full_dims >> "$LOG_FILE" 2>&1 \
                || log_warn "unwrap_detrended 地理编码失败 (详见日志)"
        fi
        local corr_src=""
        [ -f filtcorr.grd ] && corr_src="filtcorr.grd" || corr_src="corr.grd"
        if [ -n "$corr_src" ] && [ -f "$corr_src" ]; then
            python3 "$py_geocode" trans.dat "$corr_src" "$results_dir/corr.tif" 0.0001 $full_dims >> "$LOG_FILE" 2>&1 \
                || log_warn "相干性地理编码失败 (详见日志)"
        fi
        if [ -f los_displacement.grd ]; then
            python3 "$py_geocode" trans.dat los_displacement.grd "$results_dir/final_data_values.tif" 0.0001 $full_dims >> "$LOG_FILE" 2>&1 \
                || log_warn "LOS 位移地理编码失败 (详见日志)"
        fi
        
    else
        log_warn "未找到 trans.dat 且无法从 DEM 生成，跳过地理编码；仅输出雷达坐标结果。提供 topo/dem.grd 可启用 GeoTIFF 地理编码。"
    fi

    # PNG 输出（使用动态范围）
    log_info "生成可视化图像 (Safe Python Renderer)..."
    
    cat > make_safe_thumbnails.py << 'EOF'
import sys
import numpy as np
from osgeo import gdal
from PIL import Image

def arr_to_img(data, cmap_name, vmin, vmax, width=800):
    h, w = data.shape
    if w > width:
        step = max(1, w // width)
        data = data[::step, ::step]
    
    new_h, new_w = data.shape
    
    # Replace NaNs with 0 for color mapping
    nan_mask = np.isnan(data)
    data = np.nan_to_num(data, nan=0.0)
    
    norm = (data - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0, 1)
    
    rgb = np.zeros((*data.shape, 3), dtype=np.uint8)
    
    if cmap_name == 'hsv':
        # Vectorized HSV to RGB
        h_ang = norm * 6.0
        i = np.floor(h_ang)
        f = h_ang - i
        # Since S=1, V=1: p=0, q=1-f, t=f
        p = np.zeros_like(h_ang)
        q = 1.0 - f
        t = f
        i = (i % 6).astype(int)
        
        c1 = (i == 0); rgb[c1, 0], rgb[c1, 1], rgb[c1, 2] = 255, (t[c1]*255).astype(np.uint8), 0
        c2 = (i == 1); rgb[c2, 0], rgb[c2, 1], rgb[c2, 2] = (q[c2]*255).astype(np.uint8), 255, 0
        c3 = (i == 2); rgb[c3, 0], rgb[c3, 1], rgb[c3, 2] = 0, 255, (t[c3]*255).astype(np.uint8)
        c4 = (i == 3); rgb[c4, 0], rgb[c4, 1], rgb[c4, 2] = 0, (q[c4]*255).astype(np.uint8), 255
        c5 = (i == 4); rgb[c5, 0], rgb[c5, 1], rgb[c5, 2] = (t[c5]*255).astype(np.uint8), 0, 255
        c6 = (i == 5); rgb[c6, 0], rgb[c6, 1], rgb[c6, 2] = 255, 0, (q[c6]*255).astype(np.uint8)
    elif cmap_name == 'RdBu_r':
        rgb[:,:,0] = np.clip((norm * 2 * 255), 0, 255).astype(np.uint8)
        rgb[:,:,1] = np.clip((1 - np.abs(norm*2 - 1)) * 255, 0, 255).astype(np.uint8)
        rgb[:,:,2] = np.clip(((1-norm) * 2 * 255), 0, 255).astype(np.uint8)
    else:
        val = (norm * 255).astype(np.uint8)
        rgb = np.stack([val, val, val], axis=-1)
        
    # Paint NaNs as dark grey
    rgb[nan_mask] = [50, 50, 50]
        
    img = Image.fromarray(rgb)
    return img

try:
    ds1 = gdal.Open(sys.argv[1])
    if ds1:
        pf = ds1.GetRasterBand(1).ReadAsArray()
        img1 = arr_to_img(pf, 'hsv', -np.pi, np.pi)
        img1.save(sys.argv[2])
except Exception as e:
    print("Phase render error:", e)

try:
    ds2 = gdal.Open(sys.argv[3])
    if ds2:
        disp = ds2.GetRasterBand(1).ReadAsArray()
        # 第 5 个参数: 色标正负限取 |位移| 的第几百分位; 100 (默认) = 最大绝对值
        pct = float(sys.argv[5]) if len(sys.argv) > 5 else 100.0
        finite = np.abs(disp[np.isfinite(disp)])
        if finite.size == 0:
            limit = 50
        elif pct >= 100:
            limit = float(finite.max())
        else:
            limit = float(np.percentile(finite, pct))
        if not limit > 0:
            limit = float(finite.max()) if finite.size and finite.max() > 0 else 1.0
        print("Displacement colour scale: +/-%.1f mm (percentile %g of |displacement|)" % (limit, pct))
        img2 = arr_to_img(disp, 'RdBu_r', -limit, limit)
        img2.save(sys.argv[4])
except Exception as e:
    print("Disp render error:", e)
EOF
    python3 make_safe_thumbnails.py "$FILT_PHASE" "$results_dir/interferogram.png" "los_displacement.grd" "$results_dir/displacement.png" "$DISP_STRETCH_PCT" >> "$LOG_FILE" 2>&1
    log_info "干涉条纹图: results/interferogram.png"
    log_info "位移图: results/displacement.png"
        # 雷达坐标 GeoTIFF（无地理编码时保留）
        if [ "$have_trans" -eq 0 ] && [ -f los_displacement.grd ]; then
            gmt grdconvert los_displacement.grd "$results_dir/displacement.tif=gd:GTiff" >> "$LOG_FILE" 2>&1
            log_info "GeoTIFF(雷达坐标): results/displacement.tif"
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
run_qa() {
    if [ "$QA_MODE" = "off" ]; then
        log_info "产品 QA 已关闭 (--qa off)"
        return 0
    fi

    local report_dir="$PROJECT_DIR/results/qa"
    local report_path="$report_dir/qa_report.json"
    local master_prm="$PROJECT_DIR/SLC/${MASTER_DATE}.PRM"
    local qa_rc
    local overall=""
    local qa_tiles="$SNAPHU_TILES"
    local qa_overlap="$SNAPHU_OVERLAP"
    if [ "$CROP_MODE" = true ]; then
        qa_tiles="1 1"
        qa_overlap=0
    fi
    local qa_cmd=(
        python3 "$SCRIPT_DIR/dinsar_qa.py"
        --intf-dir "$INTF_DIR"
        --project-dir "$PROJECT_DIR"
        --master-prm "$master_prm"
        --wavelength "$WAVELENGTH"
        --downsample "$DOWNSAMPLE_FACTOR"
        --snaphu-tiles "$qa_tiles"
        --snaphu-overlap "$qa_overlap"
        --out-dir "$report_dir"
    )

    mkdir -p "$report_dir" || {
        if [ "$QA_MODE" = "strict" ]; then
            error_exit "无法创建 QA 输出目录: $report_dir"
        fi
        log_warn "无法创建 QA 输出目录: $report_dir"
        return 0
    }
    if [ -n "$SNAPHU_CONF_USED" ]; then
        qa_cmd+=(--snaphu-conf "$SNAPHU_CONF_USED")
    fi
    if [ "$QA_MODE" = "strict" ]; then
        qa_cmd+=(--strict)
    fi

    log_step "运行产品质量检查 (--qa $QA_MODE)..."
    if "${qa_cmd[@]}" >> "$LOG_FILE" 2>&1; then
        qa_rc=0
    else
        qa_rc=$?
    fi
    if [ -f "$report_path" ]; then
        overall=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["overall"])' "$report_path" 2>/dev/null)
    fi
    if [ -n "$overall" ]; then
        log_info "QA 总体结论: $overall; 报告: $report_path"
    else
        log_warn "QA 未生成有效报告 (rc=$qa_rc); 预期报告: $report_path"
    fi

    if [ "$QA_MODE" = "strict" ] && { [ "$qa_rc" -ne 0 ] || [ "$overall" = "FAIL" ]; }; then
        error_exit "严格 QA 未通过; 报告: $report_path" 2
    fi
    if [ "$QA_MODE" = "on" ] && { [ "$qa_rc" -ne 0 ] || [ "$overall" = "FAIL" ]; }; then
        log_warn "产品 QA 未通过，流程继续; 报告: $report_path"
    fi
    return 0
}

main() {
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║  PowerChina-1 DInSAR 数据处理脚本 v${SCRIPT_VERSION}                ║"
    echo "║  支持: DJ1 (电建一号) / BC3 (涪城一号) / LT1 (陆探一号)     ║"
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
    log_info "配准方式: $COREG_MODE"
    log_info "解缠模式: $UNWRAP_MODE"
    log_info "降采样: ${DOWNSAMPLE_FACTOR}x"
    log_info "产品 QA: $QA_MODE"
    
    # 检查系统
    check_system || exit 1
    
    # 设置目录
    setup_directories
    
    # 数据发现
    discover_data || exit 1
    preflight_ms_register || exit 1
    
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
    run_lt1_preprocess || exit 1
    run_alignment || exit 1
    run_interferogram || exit 1
    run_filter || exit 1
    run_unwrap || exit 1
    run_atm_correction || exit 1
    run_detrend || exit 1
    run_displacement || exit 1
    generate_outputs || exit 1
    run_qa || exit 1
    
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

# 运行主函数 (被 source 时只加载函数, 便于单元测试)
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
