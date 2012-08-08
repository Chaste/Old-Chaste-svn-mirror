#!/bin/bash
#
# Script to illustrate running batch jobs on cluster and passing in arguments.
#
# NB At the moment you have to manually make sure you only submit the correct number of jobs for number of processors
# This can be replaced with a call to a queueing software/grid engine thing when one is installed.
#
# This script assumes that the following has been run successfully:
# scons compile_only=1 build=GccOpt test_suite=projects/CryptInvasion/test/TestAdhesionOnly.hpp
#

MUTANT_DAMPING_MULTIPLIER[0]="0.4"
MUTANT_DAMPING_MULTIPLIER[1]="0.3"

for (( j=0 ; j<${#MUTANT_DAMPING_MULTIPLIER[*]} ; j++))
  do
  	echo "Beginning run for mutant damping multiplier = ${MUTANT_DAMPING_MULTIPLIER[$j]}."
	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
	# ">" directs std::cout to the file.
	# "2>&1" directs std::cerr to the same place.
	# "&" on the end lets the script carry on and not wait until this has finished.
    nice -20 ../build/optimised_native_ndebug/TestAdhesionOnlyRunner ${MUTANT_DAMPING_MULTIPLIER[$j]} > adhesion_only_${MUTANT_DAMPING_MULTIPLIER[$j]}_output.txt 2>&1 &
  done
echo "Jobs submitted"
