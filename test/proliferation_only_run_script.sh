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

#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[0]="0.3"
MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[0]="0.4"
#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[2]="0.5"
#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[3]="0.6"
#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[0]="0.7"
#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[1]="0.8"
#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[2]="0.9"
#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[3]="1.0"

for (( i=0 ; i<${#MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[*]} ; i++))
do
	echo "Beginning run for mutant proliferative height fraction = ${MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[$i]}."
    # NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
	# ">" directs std::cout to the file.
	# "2>&1" directs std::cerr to the same place.
	# "&" on the end lets the script carry on and not wait until this has finished.
    nice -25 ../build/optimised_native_ndebug/TestProliferationOnlyRunner ${MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[$i]} > proliferation_only_${MUTANT_PROLIFERATIVE_HEIGHT_FRACTION[$i]}_output.txt 2>&1 &
done
echo "Jobs submitted"
