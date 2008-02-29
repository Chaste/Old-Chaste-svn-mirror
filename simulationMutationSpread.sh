#$ -j y
#$ -cwd
#$ -m be 
#$ -M gary.mirams@gmail.com
export LD_LIBRARY_PATH=/opt/petsc-2.2.1-with-mpi/lib/libO_c++/linux-gnu:/usr/local/gcc-3.4.5/lib/:/home/pmxgm/lib:/opt/boost/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/ufs/pmxgm/lib
#Make sure all the work goes on the local disk in a unique folder
export CHASTE_TEST_OUTPUT=/scratch/pmxgm/$JOB_ID
mkdir $CHASTE_TEST_OUTPUT
date
echo "Memory before run is: "
df -h
echo "Writing to local disk at " $CHASTE_TEST_OUTPUT
#Copy archive across
cp -r MutationSpread $CHASTE_TEST_OUTPUT
./TestMutationSpreadRunner

date
echo "Memory after run is: "
df -h
#Copy results back and clean up
cp -r $CHASTE_TEST_OUTPUT/MutationSpread MutationSpreadAfter 
rm -rf $CHASTE_TEST_OUTPUT
