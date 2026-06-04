#!/bin/csh -f
#
#  Geo_all.csh - Geocoding for Final Deformation Map
#  
#  Part of PowerChina-1 (Spacety) DInSAR System
#  Module: Deformation Inversion
#
#  Usage: Geo_all.csh los_def.grd trans.dat [resolution]
#
#  Input:  los_def.grd  - LOS deformation in radar coordinates (meters)
#          trans.dat    - Radar to geographic coordinate transformation
#          resolution   - Output resolution in arc-seconds (default: 4)
#
#  Output: los_def_ll.grd - Geocoded LOS deformation (WGS84)
#          los_def_ll.kml - Google Earth overlay
#

if ( $#argv < 2 ) then
    echo ""
    echo "Usage: Geo_all.csh los_def.grd trans.dat [resolution]"
    echo ""
    echo "  Input:  los_def.grd  - LOS deformation (radar coords)"
    echo "          trans.dat    - Coordinate transformation table"
    echo "          resolution   - Output resolution in arc-sec (default: 4)"
    echo ""
    echo "  Output: los_def_ll.grd - Geocoded deformation (WGS84)"
    echo "          los_def_ll.kml - Google Earth overlay"
    echo ""
    exit 1
endif

set deffile = $1
set transfile = $2
set res = 4

if ( $#argv >= 3 ) then
    set res = $3
endif

echo "Geo_all.csh - Geocoding Deformation Map"
echo "  Input:      $deffile"
echo "  Trans:      $transfile"
echo "  Resolution: ${res} arc-seconds"

# Check input files
if ( ! -f $deffile ) then
    echo "ERROR: Input file $deffile not found!"
    exit 1
endif

if ( ! -f $transfile ) then
    echo "ERROR: Transformation file $transfile not found!"
    exit 1
endif

# Get geographic bounds from trans.dat
echo "  Determining geographic bounds..."
set bounds = `gmt gmtinfo $transfile -bi5d -C`
set west = `echo $bounds | awk '{printf "%.4f", $7}'`
set east = `echo $bounds | awk '{printf "%.4f", $8}'`
set south = `echo $bounds | awk '{printf "%.4f", $9}'`
set north = `echo $bounds | awk '{printf "%.4f", $10}'`

echo "  Bounds: $west/$east/$south/$north"

# Project to geographic coordinates
echo "  Projecting to geographic coordinates..."

# Sample deformation at trans.dat points
gmt grdtrack $transfile -G$deffile -bi5d -bo6d > def_track.dat

# Extract lon, lat, deformation
gmt gmtconvert def_track.dat -bi6d -bo3d -o3,4,5 > lldef.dat

# Grid the data
set inc = `echo $res | awk '{print $1/3600}'`
gmt blockmedian lldef.dat -R$west/$east/$south/$north -I$inc -bi3f -bo3f > lldef_med.dat
gmt surface lldef_med.dat -R$west/$east/$south/$north -I$inc -bi3f -Glos_def_ll.grd -T0.25 -N1000

# Add metadata
gmt grdedit los_def_ll.grd -D"LOS Deformation"/"meters"/"PowerChina-1 DInSAR"

echo "  Creating visualization..."

# Create color palette (blue = subsidence, red = uplift)
gmt grd2cpt los_def_ll.grd -Cpolar -Z -D > def.cpt

# Create KML
gmt grdimage los_def_ll.grd -JX5d -Cdef.cpt -Q > los_def_ll.ps
gmt psconvert los_def_ll.ps -TG -A -W+k

# Cleanup
rm -f def_track.dat lldef.dat lldef_med.dat def.cpt los_def_ll.ps

# Statistics
set stats = `gmt grdinfo los_def_ll.grd -L1 -C`
set def_min = `echo $stats | awk '{printf "%.4f", $6}'`
set def_max = `echo $stats | awk '{printf "%.4f", $7}'`

echo "  Geocoded deformation range: $def_min to $def_max meters"
echo "  Output: los_def_ll.grd"
echo "  Output: los_def_ll.kml"
echo "Geo_all.csh completed."

