#!/bin/csh -f
#
#  Re_atm.csh - Atmospheric Phase Correction (Empirical Phase-Elevation Model)
#  
#  Part of PowerChina-1 (Spacety) DInSAR System
#  Module: Atmospheric Correction
#
#  Usage: Re_atm.csh unwrap.grd dem_ra.grd [output.grd]
#
#  Input:  unwrap.grd  - Unwrapped phase in radar coordinates (rad)
#          dem_ra.grd  - DEM in radar coordinates (m)
#  Output: output.grd  - Atmospheric corrected phase (rad)
#

if ( $#argv < 2 ) then
    echo ""
    echo "Usage: Re_atm.csh unwrap.grd dem_ra.grd [output.grd]"
    echo ""
    exit 1
endif

set unw_file = $1
set dem_file = $2
set out_file = "unwrap_atm.grd"

if ( $#argv >= 3 ) then
    set out_file = $3
endif

# Check input files
if ( ! -f $unw_file ) then
    echo "  ERROR: Unwrapped phase file $unw_file not found!"
    exit 1
endif

if ( ! -f $dem_file ) then
    echo "  ERROR: Radar-coded DEM file $dem_file not found!"
    exit 1
endif

echo "  --- Atmospheric Phase Screen (APS) Correction ---"
echo "  Method: Empirical Phase-Elevation Regression"
echo "  Input Phase: $unw_file"
echo "  Input DEM:   $dem_file"

# 1. 采样与准备数据 (排除 NaN 值)
# 我们需要将相位和高程对齐
echo "  [1/3] Sampling and correlating phase with elevation..."

# 获取 DEM 的分辨率并重采样相位图以对齐 (防止网格不匹配)
set R = `gmt grdinfo $unw_file -I-`
gmt grdsample $dem_file $R -Gtmp_dem_matched.grd

# 将两者转为 XYZ 并合并，过滤 NaN
gmt grd2xyz $unw_file -s > tmp_phase.xyz
gmt grd2xyz tmp_dem_matched.grd -s > tmp_dem.xyz

# 使用 blockmean 减少数据量并提高回归稳定性 (降采样到 10x10 窗口)
# 获取范围
set REG = `gmt grdinfo $unw_file -I- -C | awk '{print "-R"$2"/"$3"/"$4"/"$5}'`
set INC = `gmt grdinfo $unw_file -I- -C | awk '{print "-I"$8*10"/"$9*10}'`

paste tmp_phase.xyz tmp_dem.xyz | awk '{if($3!="NaN" && $6!="NaN") print $1,$2,$3,$6}' > tmp_combined.xyz

# 2. 线性回归计算系数
# 只有在样本量足够时才运行
set n_pts = `wc -l < tmp_combined.xyz`
if ( $n_pts < 100 ) then
    echo "  WARNING: Too few valid points ($n_pts) for regression. Skipping."
    cp $unw_file $out_file
    goto cleanup
endif

echo "  [2/3] Performing linear regression on $n_pts samples..."
# 格式: x y phase elev -> 我们对 elev(index 3) 和 phase(index 2) 回归
# -Ey: Vertical misfit (standard OLS)
# -i3,2: X=col 3(elev), Y=col 2(phase)
# -Fp: Output only model parameters in a single line
gmt regress tmp_combined.xyz -i3,2 -Ey -Fp > tmp_reg_result.txt

# -Fp 输出格式: N meanX meanY angle misfit slope icept sigma_slope sigma_icept r R N_effective
# Slope 是第 6 个字段，Intercept 是第 7 个字段
set slope = `awk '{print $6}' tmp_reg_result.txt`
set intercept = `awk '{print $7}' tmp_reg_result.txt`

if ( "$slope" == "" ) then
    echo "  ERROR: Failed to estimate regression coefficients."
    cp $unw_file $out_file
    goto cleanup
endif

echo "  Regression result: Phase = ($slope) * Elevation + ($intercept)"

# 3. 应用修正
echo "  [3/3] Applying correction to grid..."
# Formula: phase_corr = phase - (slope * elevation + intercept)
# 有时我们只扣除相关项，保留常数项(归零由后期处理)
gmt grdmath $unw_file $dem_file $slope MUL SUB = $out_file

echo "  Atmospheric correction applied: $out_file"

cleanup:
rm -f tmp_phase.xyz tmp_dem.xyz tmp_combined.xyz tmp_dem_matched.grd
echo "  Done."
echo ""

