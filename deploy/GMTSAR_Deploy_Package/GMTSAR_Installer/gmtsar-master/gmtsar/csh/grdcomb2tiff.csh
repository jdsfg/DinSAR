#!/bin/csh -f
#
#  ysdong@cug, 2026
#  Fixed by AI Agent - 2026-01-22
#
unset noclobber
#
# script to convert a grd file to a geotiff file
# Supports both geographic and radar coordinates
#
if ($#argv < 5 || $#argv > 5) then
 echo " "
 echo "Usage: grdcomb2tiff.csh grd_file1 cptfile1 grd_file2 cptfile2  output"
 echo "grd_file2 以40%的透明度叠加到grd_file1之上。"
 echo "Example: grdcomb2tiff.csh final-amp_ll final-amp.cpt phasefilt_ll phase.cpt phase_amp_ll"
 echo " "
 exit 1
endif 
#
# Get grid info
set INFO = `gmt grdinfo $1 -C`
set XMIN = `echo $INFO | cut -d' ' -f2`
set XMAX = `echo $INFO | cut -d' ' -f3`
set DX = `echo $INFO | cut -d' ' -f8`

# Check if geographic or radar coordinates
set XRANGE = `gmt math -Q $XMAX $XMIN SUB = `
echo "X range: $XRANGE"

# Use Cartesian projection for radar coords (range > 360)
if (`gmt math -Q $XRANGE 360 GT = ` == 1) then
  echo "Using Cartesian projection (radar coordinates)"
  set PROJ = "-JX20c/15c"
  set DPI = 150
else
  echo "Using geographic projection"
  set PROJ = "-Jx1id"
  set DPI = `gmt math -Q $DX INV RINT = `
endif

echo "DPI: $DPI"
gmt set COLOR_MODEL = hsv
gmt set PS_MEDIA = tabloid
#
  gmt grdimage $1 -C$2 $PROJ -P -Y2i -X2i -Q -K -V >  $5.ps
  gmt grdimage $3 -C$4 $PROJ -t40 -Q -O -V >> $5.ps
#
#   now make the geotiff 
#
echo "Make $5.tiff"
gmt psconvert $5.ps -W+g+t"$1" -E$DPI -P -A -Tt
rm -f $5.ps gmt.conf gmt.history
#
