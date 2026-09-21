#!/bin/bash

rm -f *.nc
rm -f *.exe
rm -f nanuq.output* output.namelist* nanuq.stat* nanuq.step layout_nanuq.dat communication_nanuq.txt timing_nanuq.output layout.dat timing_nanuq_gnuplot.sh communication_report.txt
rm -f namelist_dom_cfg namelist_ice_cfg ; # they should be s-links!
rm -f axis_def_nanuq.xml context_nanuq.xml domain_def_nanuq.xml field_def_nanuq.xml grid_def_nanuq.xml namelist_dom_ref namelist_ice_ref
