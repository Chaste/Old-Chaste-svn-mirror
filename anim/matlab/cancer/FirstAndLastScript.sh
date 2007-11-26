#!/bin/sh
for (( i = 30; i <= 80; i++ ))
do
	cd /local/pmxgm/Simulation_Results/16_stem_cell_Meineke_recreate/sunter3/2007-11-22-19-19/MeinekeLabellingExperiment/results_from_time_"$i"0.667/vis_results/
	head -1 results.viznodes > first_line.txt | tail -1 results.viznodes > last_line.txt
done
