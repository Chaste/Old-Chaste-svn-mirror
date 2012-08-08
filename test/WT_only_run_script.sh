#!/bin/bash
#
# Script to illustrate running batch jobs on cluster and passing in arguments.
#
# NB At the moment you have to manually make sure you only submit the correct number of jobs for number of processors
# This can be replaced with a call to a queueing software/grid engine thing when one is installed.
#
# This script assumes that the following has been run successfully:
# scons compile_only=1 build=GccOpt test_suite=projects/CryptInvasion/test/TestProliferationOnly.hpp
#

SIM_RUN[0]="1"
SIM_RUN[1]="2"
SIM_RUN[2]="3"
SIM_RUN[3]="4"

for (( i=0 ; i<${#SIM_RUN[*]} ; i++))
do
	echo "Beginning run = ${SIM_RUN[$i]}."
    # NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
	# ">" directs std::cout to the file.
	# "2>&1" directs std::cerr to the same place.
	# "&" on the end lets the script carry on and not wait until this has finished.
    nice -20 ../build/optimised_ndebug/TestWildTypeMonoclonalityRunner ${SIM_RUN[$i]} > WT_${SIM_RUN[$i]}_output.txt 2>&1 &
done
echo "Jobs submitted"
