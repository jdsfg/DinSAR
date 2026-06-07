#!/bin/bash
PRM=$1
DEM_GRD=$2
LOG=$3

echo "--- FAST Radar Mapping ---" >> "$LOG"

# Get dimensions from PRM
num_rng=$(grep num_rng_bins "$PRM" | awk '{print $3}')
num_az=$(grep num_lines "$PRM" | awk '{print $3}')

REGION="-R0/$num_rng/0/$num_az"
# Coarse spacing for interpolation speed
INC_COARSE="-I4/4"
INC_FINE="-I1/1"

echo "Region: $REGION" >> "$LOG"

# Generate mapping points from DEM
if [ ! -f dem.xyz ]; then
    echo "Generating mapping points from DEM..." >> "$LOG"
    gmt grd2xyz --FORMAT_FLOAT_OUT=%lf "$DEM_GRD" -s | SAT_llt2rat "$PRM" 1 -bod > dem.xyz
fi

# Blockmedian to reduce points
gmt blockmedian dem.xyz $REGION $INC_COARSE -bi3d -bo3d -r > temp.rat

# Triangulate at coarse resolution
gmt triangulate temp.rat $REGION $INC_COARSE -bi3d -Gtopo_ra_coarse.grd -r >> "$LOG" 2>&1

# Sample back to fine resolution
gmt grdsample topo_ra_coarse.grd $REGION $INC_FINE -Gtopo_ra_fine.grd >> "$LOG" 2>&1

# Move directly without FLIPUD
mv topo_ra_fine.grd dem_ra.grd

rm -f temp.rat topo_ra_coarse.grd topo_ra_fine.grd dem.xyz
echo "Fast radar mapping completed." >> "$LOG"
