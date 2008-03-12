#$ -j y
#$ -cwd
#$ -m be 
#$ -M pmxaw@nottingham.ac.uk
export LD_LIBRARY_PATH=/opt/petsc-2.2.1-with-mpi/lib/libO_c++/linux-gnu:/usr/local/gcc-3.4.5/lib/:/home/pmxaw/lib:/opt/boost/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/ufs/pmxaw/lib
#Make sure all the work goes on the local disk in a unique folder
export CHASTE_TEST_OUTPUT=/scratch/pmxaw/$JOB_ID
mkdir $CHASTE_TEST_OUTPUT
date
echo "Memory before run is: "
df -h
echo "Writing to local disk at " $CHASTE_TEST_OUTPUT
#Copy archive across
cp -r Noddy_No_No $CHASTE_TEST_OUTPUT
./TestMakeNoNoRunner
date
echo "Memory after run is: "
df -h
#Copy results back and clean up
cp -r $CHASTE_TEST_OUTPUT . 
rm -rf $CHASTE_TEST_OUTPUT
