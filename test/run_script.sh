#!/bin/bash
#
# Script to illustrate running batch jobs on cluster and passing in arguments.
#
# NB At the moment you have to manually make sure you only submit the correct number of jobs for number of processors
# This can be replaced with a call to a queueing software/grid engine thing when one is installed.
#
# This script assumes that the following has been run successfully:
# scons compile_only=1 build=GccOpt test_suite=projects/CryptInvasion/test/TestCryptInvasion.hpp
#
# This parameter set started at 18:30 on Fri 30 Oct by AlexF

MUTANT_STATE[0]="1"
MUTANT_STATE[1]="2"
MUTANT_STATE[2]="3"
MUTANT_HEIGHT_FRACTION[0]="0"
MUTANT_HEIGHT_FRACTION[1]="1"
MUTANT_DAMPING_MULTIPLIER[0]="2.0"

for (( k=0 ; k<${#MUTANT_HEIGHT_FRACTION[*]} ; k++))
do
  for (( j=0 ; j<${#MUTANT_DAMPING_MULTIPLIER[*]} ; j++))
  do
  	for (( i=0 ; i<${#MUTANT_STATE[*]} ; i++))
  	do
  		echo "Beginning run for mutation state = ${MUTANT_STATE[$i]}, mutant height fraction = ${MUTANT_HEIGHT_FRACTION[$k]}, mutant damping multiplier = ${MUTANT_DAMPING_MULTIPLIER[$j]}."
	    # NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
	    # ">" directs std::cout to the file.
	    # "2>&1" directs std::cerr to the same place.
	    # "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../build/optimised_native_ndebug/TestCryptInvasionRunner ${MUTANT_STATE[$i]} ${MUTANT_HEIGHT_FRACTION[$k]} ${MUTANT_DAMPING_MULTIPLIER[$j]} > crypt_invasion_${MUTANT_STATE[$i]}_${MUTANT_HEIGHT_FRACTION[$k]}_${MUTANT_DAMPING_MULTIPLIER[$j]}_output_native.txt 2>&1 &
	done
  done
done
echo "Jobs submitted"
