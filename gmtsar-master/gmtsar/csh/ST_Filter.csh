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
#    - If >=3 other interferograms of identical grid size exist under
#      <project_dir>/intf, uses a common-point stacking approach to estimate
#      and remove transient atmospheric phase screen (APS).
#    - Otherwise applies a spatial Gaussian low-pass to suppress turbulent APS.
#
#  Note: input.grd is in radar coordinates (x = range, y = azimuth, in pixels),
#        so all filter widths are converted from metres to grid units and
#        grdfilter runs in Cartesian mode (-D0). The pixel size is read from
#        a *.PRM in the current directory when available, else ST_PIXEL_M
#        (env, metres) or 3 m is used.
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

# Ground-equivalent pixel size (m), used to turn metre widths into grid units
set pix_m = 3.0
if ( $?ST_PIXEL_M ) then
    set pix_m = $ST_PIXEL_M
endif
set nonomatch
set prm_list = ( *.PRM )
unset nonomatch
if ( -f "$prm_list[1]" ) then
    set fs  = `grep rng_samp_rate $prm_list[1] | head -1 | awk '{print $3}'`
    set prf = `grep '^PRF' $prm_list[1] | head -1 | awk '{print $3}'`
    set vel = `grep SC_vel $prm_list[1] | head -1 | awk '{print $3}'`
    if ( "$fs" != "" && "$prf" != "" && "$vel" != "" ) then
        set pix_m = `echo $fs $prf $vel | awk '{dr=299792458/(2*$1); da=$3/$2; printf("%.3f", (dr+da)/2)}'`
    endif
endif
echo "  Pixel size used for filter widths: $pix_m m"

# Peers: other interferograms of the project with the same grid as the input.
# The input itself must be excluded by resolved path, not by name, because
# every pair directory contains a file called unwrap.grd.
set self_abs = `readlink -f $infile`
set ref_dims = `gmt grdinfo $infile -C | awk '{print $10, $11}'`
set intf_list = ()
if ( -d $proj_dir/intf ) then
    foreach f ( `find $proj_dir/intf -name "unwrap.grd" | sort` )
        if ( "`readlink -f $f`" == "$self_abs" ) continue
        set dims = `gmt grdinfo $f -C | awk '{print $10, $11}'`
        if ( "$dims" != "$ref_dims" ) continue
        set intf_list = ( $intf_list $f )
    end
endif
set n_stack = $#intf_list

if ( $n_stack >= 3 ) then
    echo "  [Mode] Multi-temporal Stacking Filter ($n_stack peers found)"
    # Common-mode approach: the mean of the peers approximates the signal that
    # is stable over time (only meaningful for peers with similar time spans);
    # what the current pair has on top of it, once smoothed, is taken as APS.
    
    # 1. Calculate stack mean (common deformation signal)
    gmt grdmath $intf_list[1] $intf_list[2] ADD $intf_list[3] ADD 3 DIV = tmp_stack_mean.grd
    
    # 2. Difference current with stack mean
    # The residual contains transient APS
    gmt grdmath $infile tmp_stack_mean.grd SUB = tmp_residual.grd
    
    # 3. Spatial filter the residual to get APS
    # APS is spatially smooth (~1-5km)
    set w_px = `echo $pix_m | awk '{printf("%.1f", 2000/$1)}'`
    gmt grdfilter tmp_residual.grd -Gtmp_aps.grd -Fg$w_px -D0
    
    # 4. Subtract APS from original
    gmt grdmath $infile tmp_aps.grd SUB = $outfile
    
    rm -f tmp_stack_mean.grd tmp_residual.grd tmp_aps.grd
    echo "  Multi-temporal filter applied."
else
    echo "  [Mode] Spatial Adaptive Filter (Single-pair mode, $n_stack comparable peers)"
    # APS is spatially smooth relative to decorrelation noise; for X-band
    # (DJ1/BC3) a Gaussian low-pass suppresses the high-frequency part.
    set w_px = `echo $pix_m | awk '{printf("%.1f", 500/$1)}'`
    gmt grdfilter $infile -G$outfile -Fg$w_px -D0
    echo "  Spatial adaptive smoothing applied (Gaussian 500 m = $w_px px)."
endif

if ( ! -f $outfile ) then
    echo "  ERROR: ST_Filter failed to produce $outfile"
    exit 1
endif

echo "  ST_Filter.csh completed."
echo ""
