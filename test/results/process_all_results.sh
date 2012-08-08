#!/bin/bash
#
# This script wipes the processed results directory and runs each script.
#
# NB they shouldn't be run individually as they append results to existing files if the directory isn't wiped.
#
# Last updated 14/2/12

if cd processed; then
    rm -f processed_*.dat
fi
cd ..

./process_overall_results_WT.sh
./process_overall_results_files.sh
./process_proliferation_only_results_files.sh

