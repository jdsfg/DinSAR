#!/bin/bash
# 用法: fast_dem2topo_ra.sh master.PRM dem.grd logfile [dec]
#   dem.grd 的高程须是相对 WGS84 椭球的高度; 输出 dem_ra.grd (雷达坐标下的地形, GMTSAR 的 topo_ra 约定:
#   高程是相对 PRM 的 earth_radius 球面的, 与 phasediff 的球面几何一致)。
#   dec (可选): 输出网格相对 SLC 的降采样倍数。
#     1 (默认)  整景分辨率;
#     N > 1     每 N 个像元一个节点的粗网格。phasediff 允许 topo 比 SLC 粗, 前提是 SLC 尺寸是它的整数倍, 读的时候
#               第 k 个像元直接取第 k/N (整除) 个节点。整景分辨率的 dem_ra.grd 与 phasediff 的 real/imag 同样大,
#               大场景 (如 1.45 GB) 会把内存挤爆;
#     auto      4/2/1 里能同时整除 SLC 两个维度的最大值。
#   DEM 点投到雷达坐标后由 topo_ra_grid.py 做 Delaunay 线性内插 (与 gmt triangulate 的结果一致, 但快得多)。
set -o pipefail

PRM=$1
DEM_GRD=$2
LOG=$3
DEC=${4:-1}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

echo "--- FAST Radar Mapping ---" >> "$LOG"

# Get dimensions from PRM
num_rng=$(awk '$1 == "num_rng_bins" {print $3; exit}' "$PRM")
num_az=$(awk '$1 == "num_lines" {print $3; exit}' "$PRM")
if [ -z "$num_rng" ] || [ -z "$num_az" ]; then
    echo "PRM 里没有 num_rng_bins / num_lines: $PRM" >> "$LOG"
    exit 1
fi

if [ "$DEC" = "auto" ]; then
    DEC=1
    for d in 4 2; do
        if [ $((num_rng % d)) -eq 0 ] && [ $((num_az % d)) -eq 0 ]; then
            DEC=$d
            break
        fi
    done
fi
case "$DEC" in
    ''|*[!0-9]*|0) echo "无效的 dec=$DEC, 改用 1" >> "$LOG"; DEC=1 ;;
esac
if [ $((num_rng % DEC)) -ne 0 ] || [ $((num_az % DEC)) -ne 0 ]; then
    echo "SLC ${num_rng}x${num_az} 不能被 $DEC 整除, 改用整景分辨率 (dec=1)" >> "$LOG"
    DEC=1
fi

echo "SLC: ${num_rng}x${num_az}  dec: $DEC" >> "$LOG"

# Generate mapping points from DEM
# SAT_llt2rat -bod 每条记录输出 5 个 double: 距离向像元、方位向像元、相对 earth_radius 球面的高程、经度、纬度。
if [ ! -s dem.xyz ]; then
    echo "Generating mapping points from DEM..." >> "$LOG"
    rm -f dem.xyz dem.xyz.tmp
    if ! gmt grd2xyz --FORMAT_FLOAT_OUT=%lf "$DEM_GRD" -s 2>> "$LOG" | SAT_llt2rat "$PRM" 1 -bod > dem.xyz.tmp 2>> "$LOG" \
            || [ ! -s dem.xyz.tmp ]; then
        echo "DEM 投到雷达坐标失败 (没有点落在影像范围内, 或 SAT_llt2rat 出错)" >> "$LOG"
        rm -f dem.xyz.tmp
        exit 1
    fi
    mv dem.xyz.tmp dem.xyz
fi

# 网格第 j 行 = 第 j 段方位线, 用 -ZTLf 读进来正好落在 GMT 内存的第 j 行 (phasediff 按这个顺序读 topo), 所以这里不需要再翻转。
rm -f dem_ra.grd topo_ra.f32
if ! python3 "$HERE/topo_ra_grid.py" dem.xyz "$num_rng" "$num_az" "$DEC" topo_ra.f32 >> "$LOG" 2>&1 \
        || ! gmt xyz2grd topo_ra.f32 -R0/$num_rng/0/$num_az -I$DEC/$DEC -r -ZTLf -Gdem_ra.grd >> "$LOG" 2>&1; then
    echo "雷达坐标 DEM 网格生成失败" >> "$LOG"
    rm -f dem_ra.grd topo_ra.f32 dem.xyz
    exit 1
fi

rm -f topo_ra.f32 dem.xyz
echo "Fast radar mapping completed." >> "$LOG"
