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
#  Pixels where the phase or the DEM is 0/NaN are treated as invalid: they are
#  excluded from the regression and passed through unchanged.
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

# 1. 采样与准备数据
echo "  [1/3] Sampling and correlating phase with elevation..."

# 将 DEM 重采样到相位网格 (防止网格不匹配)
gmt grdsample $dem_file -R$unw_file -Gtmp_dem_matched.grd
if ( $status != 0 || ! -f tmp_dem_matched.grd ) then
    echo "  ERROR: Failed to resample DEM onto the phase grid."
    cp $unw_file $out_file
    goto cleanup
endif

# 不使用 grd2xyz -s: 两个网格各自丢弃 NaN 行会使 paste 后的行错位。
# 无效像元 (相位=0 为掩膜/填充值, DEM=0 为 NaN 清洗值) 必须剔除, 否则它们主导回归。
set nxy = `gmt grdinfo $unw_file -C | awk '{print $10*$11}'`
set stride = `echo $nxy | awk '{print int($1/1000000)+1}'`

gmt grd2xyz $unw_file > tmp_phase.xyz
gmt grd2xyz tmp_dem_matched.grd > tmp_dem.xyz

paste tmp_phase.xyz tmp_dem.xyz | awk -v s=$stride '{if($3!="NaN" && $6!="NaN" && $3!=0 && $6!=0 && (++k % s)==0) print $1,$2,$3,$6}' > tmp_combined.xyz

# 2. 线性回归计算系数
# 只有在样本量足够时才运行
set n_pts = `wc -l < tmp_combined.xyz`
if ( $n_pts < 100 ) then
    echo "  WARNING: Too few valid points ($n_pts) for regression. Skipping."
    cp $unw_file $out_file
    goto cleanup
endif

echo "  [2/3] Performing linear regression on $n_pts samples (stride $stride)..."
# 格式: x y phase elev -> 对 elev(index 3) 和 phase(index 2) 回归
# -Ey: Vertical misfit (standard OLS); 稳健回归 (-Nr/-Nw) 在 10^5 点以上耗时数分钟，不采用
# -Fp: Output only model parameters in a single line
gmt regress tmp_combined.xyz -i3,2 -Ey -Fp > tmp_reg_result.txt

# -Fp 输出格式: N meanX meanY angle misfit slope icept sigma_slope sigma_icept r R N_effective
set mean_elev = `awk '{print $2}' tmp_reg_result.txt`
set slope = `awk '{print $6}' tmp_reg_result.txt`
set intercept = `awk '{print $7}' tmp_reg_result.txt`
set corr = `awk '{print $10}' tmp_reg_result.txt`

if ( "$slope" == "" || "$mean_elev" == "" ) then
    echo "  ERROR: Failed to estimate regression coefficients."
    cp $unw_file $out_file
    goto cleanup
endif

echo "  Regression result: Phase = ($slope) * Elevation + ($intercept), r = $corr"

# 3. 应用修正
echo "  [3/3] Applying correction to grid..."
# phase_corr = phase - slope * (elevation - mean_elevation)
# 减去样本平均高程后, 修正量不依赖 DEM 的高程基准 (radar-coded DEM 可能相对地球半径),
# 相位整体均值不被平移。
gmt grdmath $unw_file tmp_dem_matched.grd $mean_elev SUB $slope MUL SUB = tmp_corrected.grd
# 无效像元保持原值
gmt grdmath $unw_file 0 NEQ tmp_dem_matched.grd 0 NEQ MUL tmp_corrected.grd $unw_file IFELSE = $out_file

echo "  Atmospheric correction applied: $out_file"

cleanup:
rm -f tmp_phase.xyz tmp_dem.xyz tmp_combined.xyz tmp_dem_matched.grd tmp_corrected.grd tmp_reg_result.txt
echo "  Done."
echo ""

