#!/bin/sh
for (( i = 30; i <= 81; i++ ))
do
	cd /local/pmxaw/MeinekeLabellingExperiment/results_from_time_"$i"0.667/vis_results/
	head -1 results.viznodes > first_line.txt | tail -1 results.viznodes > last_line.txt
done
