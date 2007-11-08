#$ -j y
#$ -cwd
#$ -m be 
#$ -M gary.mirams@gmail.com
export LD_LIBRARY_PATH=/opt/petsc-2.2.1-with-mpi/lib/libO_c++/linux-gnu:/usr/local/gcc-3.4.5/lib/:/home/pmxgm/lib
./TestMeinekeLabellingExperimentsRunner

