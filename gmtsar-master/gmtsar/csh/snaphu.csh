#!/bin/csh -f
#
#  snaphu.csh - Phase Unwrapping Wrapper Script
#  
#  Part of PowerChina-1 (Spacety) DInSAR System
#  Module: Phase Processing
#
#  Usage: snaphu.csh phasefilt.grd corr.grd [threshold]
#
#  Input:  phasefilt.grd  - Filtered wrapped phase (radians)
#          corr.grd       - Coherence for cost calculation
#          threshold      - Coherence threshold (default: 0.1)
#
#  Output: unwrap.grd     - Unwrapped phase (radians)
#

if ( $#argv < 2 ) then
    echo ""
    echo "Usage: snaphu.csh phasefilt.grd corr.grd [corr_threshold]"
    echo ""
    echo "  Input:  phasefilt.grd  - Filtered wrapped phase"
    echo "          corr.grd       - Coherence map"
    echo "          corr_threshold - Mask threshold (default: 0.1)"
    echo ""
    echo "  Output: unwrap.grd     - Unwrapped phase (radians)"
    echo ""
    exit 1
endif

set phasefile = $1
set corrfile = $2
set threshold = 0.1

if ( $#argv >= 3 ) then
    set threshold = $3
endif

echo "snaphu.csh - Phase Unwrapping"
echo "  Input phase: $phasefile"
echo "  Input corr:  $corrfile"
echo "  Threshold:   $threshold"

# Get dimensions (ncols = $10, nrows = $11)
set ncols = `gmt grdinfo $phasefile -C | awk '{print $10}'`
set nrows = `gmt grdinfo $phasefile -C | awk '{print $11}'`
echo "  Dimensions:  $ncols cols x $nrows rows"

# Convert to SNAPHU binary format using Python (handles NaN properly)
echo "  Converting to SNAPHU format..."

python3 << EOF
import subprocess
import numpy as np

# Get phase data
result = subprocess.run(['gmt', 'grd2xyz', '$phasefile', '-ZBLf'], capture_output=True)
phase = np.frombuffer(result.stdout, dtype=np.float32).copy()

# Get corr data
result = subprocess.run(['gmt', 'grd2xyz', '$corrfile', '-ZBLf'], capture_output=True)
corr = np.frombuffer(result.stdout, dtype=np.float32).copy()

# Replace NaN with 0
phase[np.isnan(phase)] = 0.0
corr[np.isnan(corr)] = 0.0

# Apply coherence threshold mask
phase[corr < $threshold] = 0.0
corr[corr < $threshold] = 0.0

phase.tofile('phase.in')
corr.tofile('corr.in')
print(f"  Phase range: {phase.min():.3f} to {phase.max():.3f}")
print(f"  Corr range: {corr.min():.3f} to {corr.max():.3f}")
EOF

if ( $status != 0 ) then
    echo "ERROR: Failed to convert to SNAPHU format"
    exit 1
endif

# Run SNAPHU (linelength = nrows for row-major data)
echo "  Running SNAPHU unwrapping..."
snaphu phase.in $nrows -d -c corr.in -o unwrap.out

if ( ! -f unwrap.out ) then
    echo "ERROR: SNAPHU failed to produce output"
    exit 1
endif

# Convert back to GMT grid using Python
echo "  Converting to GMT grid..."

python3 << EOF
import numpy as np
import subprocess

# Read unwrapped phase
unwrap = np.fromfile('unwrap.out', dtype=np.float32)
unwrap = unwrap.reshape($nrows, $ncols)

# Get grid info
result = subprocess.run(['gmt', 'grdinfo', '$phasefile', '-C'], capture_output=True, text=True)
info = result.stdout.strip().split()
xmin, xmax = float(info[1]), float(info[2])
ymin, ymax = float(info[3]), float(info[4])
xinc = (xmax - xmin) / ($ncols - 1) if $ncols > 1 else 1
yinc = (ymax - ymin) / ($nrows - 1) if $nrows > 1 else 1

# Write as xyz for GMT
with open('unwrap.xyz', 'w') as f:
    for i in range($nrows):
        for j in range($ncols):
            x = xmin + j * xinc
            y = ymax - i * yinc
            f.write(f"{x} {y} {unwrap[i,j]}\n")
EOF

# Create grid
set region = `gmt grdinfo $phasefile -I-`
set inc = `gmt grdinfo $phasefile -I`
gmt xyz2grd unwrap.xyz -Gunwrap.grd $region $inc

# Cleanup
rm -f phase.in corr.in unwrap.out unwrap.xyz

echo "  Output: unwrap.grd"
echo "snaphu.csh completed."
