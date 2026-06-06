#!/bin/csh -f
#
#  ST_Filter.csh - Spatio-Temporal Atmospheric Filtering (Phase C)
#  
#  Part of PowerChina-1 (Spacety) DInSAR System
#  Module: Advanced Atmospheric Correction
#
#  Usage: ST_Filter.csh input.grd output.grd [project_dir]
#
#  Method: 
#    - If multiple interferograms exist, uses a common-point stacking approach 
#      to estimate and remove transient atmospheric phase screen (APS).
#    - If single pair, applies a spatial high-pass/low-pass Gaussian filter 
#      to separate turbulent APS from deformation.
#

if ( $#argv < 2 ) then
    echo ""
    echo "Usage: ST_Filter.csh input.grd output.grd [project_dir]"
    echo ""
    exit 1
endif

set infile = $1
set outfile = $2
set proj_dir = "."
if ( $#argv >= 3 ) then
    set proj_dir = $3
endif

if ( ! -f $infile ) then
    echo "  ERROR: Input file $infile not found!"
    exit 1
endif

echo "  --- Spatio-Temporal APS Filtering (Phase C) ---"

# Check for stacking possibility
set intf_list = `find $proj_dir/intf -name "unwrap.grd" | grep -v "$infile"`
set n_stack = `echo $intf_list | wc -w`

if ( $n_stack >= 3 ) then
    echo "  [Mode] Multi-temporal Stacking Filter ($n_stack peers found)"
    # Simple Temporal APS estimation: 
    # APS_i = unw_i - median(all unw mapped to same time span)
    # Here we use a simpler common-mode approach: 
    # High-frequency spatial components that are inconsistent over time are APS.
    
    # 1. Calculate stack mean (common deformation signal)
    # We use a localized stack for efficiency
    gmt grdmath $intf_list[1] $intf_list[2] ADD $intf_list[3] ADD 3 DIV = tmp_stack_mean.grd
    
    # 2. Difference current with stack mean
    # The residual contains transient APS
    gmt grdmath $infile tmp_stack_mean.grd SUB = tmp_residual.grd
    
    # 3. Spatial filter the residual to get APS
    # APS is spatially smooth (~1-5km)
    gmt grdfilter tmp_residual.grd -Gtmp_aps.grd -Fg2000 -D4
    
    # 4. Subtract APS from original
    gmt grdmath $infile tmp_aps.grd SUB = $outfile
    
    rm -f tmp_stack_mean.grd tmp_residual.grd tmp_aps.grd
    echo "  Multi-temporal filter applied."
else
    echo "  [Mode] Spatial Adaptive Filter (Single-pair mode)"
    # Atmospheric turbulence follows a power law. 
    # We use a Gaussian high-pass filter to separate small-scale turbulence 
    # from larger-scale deformation trends (if detrending was already done).
    
    # Separating APS (high-frequency spatial) from long-term signal
    # We assume deformation is smoother than turbulent APS
    gmt grdfilter $infile -Gtmp_lowpass.grd -Fg3000 -D4
    gmt grdmath $infile tmp_lowpass.grd SUB = tmp_highpass.grd
    
    # In p2p, we often want to keep the lowpass (deformation) and remove highpass (noise)
    # But APS can be mid-frequency. 
    # Standard practice: Adaptive Goldstein-like but on unwrapped phase.
    
    # For DJ1/BC3 (X-band), we'll do a simple low-pass smoothing to suppress high-freq APS
    gmt grdfilter $infile -G$outfile -Fg500 -D4
    echo "  Spatial adaptive smoothing applied (Gaussian 500m)."
endif

echo "  ST_Filter.csh completed."
echo ""
