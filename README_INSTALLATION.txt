CHASTE DEVELOPER INSTALLATION GUIDE

Go to a terminal and try typing
$HOME
if it doesn't return a directory then you need to define the variable e.g.
HOME=~

The following packages and libraries are all compulsory for chaste to run 
(the version numbers of packages and libraries might need updating). 
Please install them as described below.

==========SCONS:=============
(Python is a prerequisite for this)
Use your package manager to install scons or

wget http://mesh.dl.sourceforge.net/sourceforge/scons/scons-1.2.0.tar.gz
gunzip scons-1.2.0.tar.gz
tar -xf scons-1.2.0.tar
python setup.py install --prefix=$HOME

========MPI and PETSC:=======

cd $HOME
wget ftp://ftp.mcs.anl.gov/pub/mpi/mpich.tar.gz
tar -xzvf mpich.tar.gz
cd mpich-1.2.7p1
./configure --prefix=${HOME}/mpi --with-comm=shared --with-device=ch_shmem --enable-sharedlib --disable-f77
make
cd examples/test/
make testing
cd ../..
make install
cd
rm -rf mpich-1.2.7p1/  mpich.tar.gz 

cd $HOME
wget ftp://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-2.3.3-p15.tar.gz
tar -xzvf petsc-2.3.3-p15.tar.gz 
rm petsc-2.3.3-p15.tar.gz 
cd petsc-2.3.3-p15/
export PETSC_DIR=`pwd`
./config/configure.py  --download-c-blas-lapack=1 --with-mpi-dir=${HOME}/mpi --with-x=false  -PETSC_ARCH=linux-gnu --with-clanguage=cxx
make all
./config/configure.py  --download-c-blas-lapack=1 --with-mpi-dir=${HOME}/mpi --with-x=false   --with-debugging=0 -PETSC_ARCH=linux-gnu-opt --with-clanguage=cxx
make all
./config/configure.py  --download-c-blas-lapack=1 --with-mpi-dir=${HOME}/mpi --with-x=false  -PETSC_ARCH=linux-gnu-profile --with-clanguage=cxx --CFLAGS=-pg -CXXFLAGS=-pg -LDFLAGS=-pg
make all

=========BOOSTJAM and BOOST:===========

Download BoostJam from http://sourceforge.net/project/showfiles.php?group_id=7586&package_id=72941 according to your platform into your ${HOME} directory.
Download Boost from http://sourceforge.net/project/showfiles.php?group_id=7586&package_id=8041 according to your platform into your ${HOME} directory.
Unpak and install them both:

tar -xzvf boost-jam-3.1.12-1-linuxx86.tgz
tar -zxvf boost_1_33_1.tar.gz
cd boost_1_33_1
../boost-jam-3.1.12-1-linuxx86/bjam "-sTOOLS=gcc" --prefix=$HOME/boost install

==============HDF5:====================

cd $HOME
wget ftp://ftp.hdfgroup.org/HDF5/prev-releases/hdf5-1.6.6/src/hdf5-1.6.6.tar.gz
tar -zxf hdf5-1.6.6.tar.gz 
cd hdf5-1.6.6
CC=mpicc ./configure --enable-parallel --prefix=${HOME}/hdf5
make
cd test
make check
cd ../testpar
make check
cd ..
make install


================METIS:==================

cd $HOME
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.0pre2.tar.gz
tar -zxf metis-5.0pre2.tar.gz
cd metis-5.0pre2

Edit the file Makefile.in, changing the line CC = cc to CC = gcc.
Also edit the file include/metis.h and specify the width (32 or 64 bits) of the elementary data type used in METIS. This is controled by the IDXTYPEWIDTH constant. For now, on a 32 bit architecture you can only specify a width of 32, whereas for a 64 bit architecture you can specify a width of either 32 or 64 bits.
Then run:

make all

============= XSD and XML =============
Installing Apache Xerces XML libraries:

wget http://apache.mirror.infiniteconflict.com/xerces/c/3/sources/xerces-c-3.0.1.tar.gz
tar -zxf xerces-c-3.0.1.tar.gz
cd xerces-c-3.0.1
./configure
make
cd ..

If you experience any problems then you can modify the configuration options
for your system according to http://xerces.apache.org/xerces-c/build-3.html

Installing XSD:

cd ~
wget http://codesynthesis.com/download/xsd/2.3/linux-gnu/i686/xsd-2.3.1-i686-linux-gnu.tar.bz2
tar -xjf xsd-2.3.1-i686-linux-gnu.tar.bz2
ln -s ~/xsd-2.3.1-i686-linux-gnu/bin/xsd ~/bin/xsd


===========SET ENVIRONMENTAL VARIABLES AND PATHS:===================

Edit the .bashrc file in the home directory and append the following:

export LD_LIBRARY_PATH=$HOME/petsc-2.2.1/lib/libg_c++/linux-gnu/

NB. change petsc-2.2.1 to the vesions you installed in earlier steps NB. You will need to logout and in again for the above to take effect. 

Also, alter ~/.bash_profile to set the PATH as
PATH=$HOME/bin:$PATH

/// \todo Check that the below instructions are valid once release is finalised.

You will need to edit the python/hostconfig/default.py to give the correct paths to each of the libraries.
 

========== Compiling chaste ==========
Now we should have all the necessary libraries to compile chaste. The following commands 
should generate static and dynamic chaste libraries.

scons compile_only=1 chaste_libs=1 static=0 build=GccOpt exe=1 apps
scons compile_only=1 chaste_libs=1 static=1 build=GccOpt exe=1 apps

The executable Chaste or Chaste.o can be found in the apps/src folder.
