#!/bin/csh -f

set script_dir = "$0:h"
if ( -x "$script_dir/atm_raytrace" ) then
    exec "$script_dir/atm_raytrace" $argv:q
else if ( -x "$script_dir/../atm_raytrace" ) then
    exec "$script_dir/../atm_raytrace" $argv:q
else
    exec atm_raytrace $argv:q
endif
