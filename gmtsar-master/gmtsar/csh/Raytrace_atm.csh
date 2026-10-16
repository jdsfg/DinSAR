#!/bin/csh -f
# Raytrace_atm.csh
# Wrapper for GMTSAR ray-tracing physical atmospheric correction workflow

if ( $#argv < 4 ) then
    echo "Usage: Raytrace_atm.csh <source_phase.grd> <output_atm.grd> <project_dir> <stage_dir> [weather_dir]"
    exit 1
endif

set source_phase = "$1"
set output_phase = "$2"
set project_dir  = "$3"
set stage_dir    = "$4"
set weather_dir  = ""
if ( $#argv >= 5 ) then
    set weather_dir = "$5"
endif

set script_dir = "$0:h"
set py_workflow = "$script_dir/raytrace_atm_workflow.py"

if ( ! -f "$py_workflow" ) then
    echo "Error: Cannot find $py_workflow"
    exit 1
endif

python3 "$py_workflow" \
    --source-phase "$source_phase" \
    --output-phase "$output_phase" \
    --project-dir "$project_dir" \
    --stage-dir "$stage_dir" \
    --weather-dir "$weather_dir"

exit $status
