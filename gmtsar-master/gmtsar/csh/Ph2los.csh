#!/bin/csh -f
#
#  Ph2los.csh - Phase to LOS Deformation Conversion
#  
#  Part of PowerChina-1 (Spacety) DInSAR System
#  Module: Deformation Inversion
#
#  Usage: Ph2los.csh unwrap.grd wavelength
#
#  Input:  unwrap.grd  - Unwrapped phase (radians)
#          wavelength  - Radar wavelength in meters (default: 0.0555 for C-band)
#
#  Output: los_def.grd - LOS deformation in meters
#
#  Formula: d = (lambda / 4*pi) * phi
#           where phi is unwrapped phase in radians
#

if ( $#argv < 1 ) then
    echo ""
    echo "Usage: Ph2los.csh unwrap.grd [wavelength]"
    echo ""
    echo "  Input:  unwrap.grd  - Unwrapped phase (radians)"
    echo "          wavelength  - Radar wavelength in meters"
    echo "                        Default: 0.0555 (C-band, Tianyi BC3)"
    echo ""
    echo "  Output: los_def.grd - LOS deformation (meters)"
    echo ""
    echo "  Formula: d = (lambda / 4*pi) * phi"
    echo ""
    exit 1
endif

set unwrapfile = $1
set lambda = 0.0555

if ( $#argv >= 2 ) then
    set lambda = $2
endif

# Calculate conversion factor: lambda / (4 * pi)
set factor = `echo $lambda | awk '{printf "%.10f", $1 / (4 * 3.14159265359)}'`

echo "Ph2los.csh - Phase to LOS Deformation"
echo "  Input:      $unwrapfile"
echo "  Wavelength: $lambda m"
echo "  Factor:     $factor m/rad"

# Check input file exists
if ( ! -f $unwrapfile ) then
    echo "ERROR: Input file $unwrapfile not found!"
    exit 1
endif

# Convert phase to deformation
# d = (lambda / 4*pi) * phi
echo "  Converting phase to LOS deformation..."
gmt grdmath $unwrapfile $factor MUL = los_def.grd

# Get statistics
set stats = `gmt grdinfo los_def.grd -L1 -C`
set def_min = `echo $stats | awk '{printf "%.4f", $6}'`
set def_max = `echo $stats | awk '{printf "%.4f", $7}'`

echo "  Deformation range: $def_min to $def_max meters"
echo "  Output: los_def.grd"
echo "Ph2los.csh completed."

