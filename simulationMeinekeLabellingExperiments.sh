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

#echo "Disk usage before run is: "
#df -h
#du -sh /scratch/pmxgm/$JOB_ID

echo "Writing to local disk at " $CHASTE_TEST_OUTPUT
#Copy archive across

cp -r MeinekeLabellingExperiment $CHASTE_TEST_OUTPUT
./TestMeinekeLabellingExperimentsRunner

date

#echo "Disk usage after run is: "
#df -h
#du -sh /scratch/pmxgm/$JOB_ID

#Copy results back and clean up
cp -r $CHASTE_TEST_OUTPUT . 
rm -rf $CHASTE_TEST_OUTPUT

# Get rid of the $JOB_ID folder and copy things back to the original place.
rm -rf MeinekeLabellingExperiment/
mv $JOB_ID/* .
rm -rf $JOB_ID/


