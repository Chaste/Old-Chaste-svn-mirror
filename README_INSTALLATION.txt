############################ ISSUES TO BE RESOLVED ########################
##
## Notes
##  1. We changed "X=Y" to "export X=Y" everywhere
##  2. Changed bash_profile to bashrc
##  3. Changed the xerces instructions to the website xerces instructions
##
## Errors:
## 1. Error in parallel checks of HDF5 (the second 'make check'). See #763 for error messages
##      - TestHdf5DataReader doesn't compile
## 2. Triangle is not mentioned at all - or on the website??? 
##      - Tests which require triangle fail "Remeshing (by calling triangle) failed"
##      - Tests which require tetgen fail 
##          - "tetgen: tetgen.cxx:4918: void tetgenmesh::findedge(tetgenmesh::triface*, double*, double*): Assertion `i < 3' failed."
## 3. The readme said output would be written in <chaste_dir>/testoutput if CHASTE_TEST_OUTPUT is not set, but after running scons
##    the output was written to /tmp/<username>/testoutout. Updated readme to say the latter. Bug?
##
## Issues:
## **1**. Not at all clear how to edit default.py once you get there - easy to edit the petsc path
##        but then had to just start running scons and adding the unfound libraries to other_libpaths
## 2. Maybe get all libraries put in the same place, $CHASTE_LIBS/libs, rather than
##    $CHASTE_LIBS/libs and $CHASTE_LIBS/boost/libs and $CHASTE_LIBS/blah/libs etc.
## 3. Are the two edits in the Metis section needed. We purposefully didn't do the first
##    (CC=cc -> CC=gcc) and forgot the second (edit include/metis.h) but had no problems.
## 4. Intel production not mentioned - suggest we create a public webpage for this 
##     (installing MKL, optimised petsc, etc) and give a link here.
## 5. Installation takes too long to be done at the workshop??? 
##    
## 
##
###########################################################################




CHASTE DEVELOPER INSTALLATION GUIDE

First define the folder where you want to install the libraries which chaste depends on, e.g.
export CHASTE_LIBS=~/chaste-libs

If it doesn't already exist, create this directory:
mkdir $CHASTE_LIBS



===================================================================================
=
=   Installing prerequisites
=
===================================================================================

The following packages and libraries are all compulsory for chaste to run 
(the version numbers of packages and libraries might need updating). 
Please install them as described below.


==========SCONS:=============
(Python is a prerequisite for this)
Use your package manager to install scons or

cd $CHASTE_LIBS
wget http://mesh.dl.sourceforge.net/sourceforge/scons/scons-1.2.0.tar.gz
tar zxf scons-1.2.0.tar.gz
cd scons-1.2.0
python setup.py install --prefix=$CHASTE_LIBS
cd ..
rm -rf scons-1.2.0.tar.gz scons-1.2.0

========MPI and PETSC:=======

cd $CHASTE_LIBS
wget ftp://ftp.mcs.anl.gov/pub/mpi/mpich.tar.gz
tar -xzvf mpich.tar.gz
cd mpich-1.2.7p1
./configure --prefix=$CHASTE_LIBS --with-comm=shared --with-device=ch_shmem --enable-sharedlib --disable-f77
make
cd examples/test/
make testing
cd ../..
make install
cd ..
rm -rf mpich-1.2.7p1/  mpich.tar.gz 

cd $CHASTE_LIBS
wget ftp://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-2.3.3-p15.tar.gz
tar -xzvf petsc-2.3.3-p15.tar.gz
rm -rf petsc-2.3.3-p15.tar.gz
cd petsc-2.3.3-p15/
export PETSC_DIR=`pwd`
./config/configure.py  --download-c-blas-lapack=1 --with-mpi-dir=$CHASTE_LIBS/mpi --with-x=false  -PETSC_ARCH=linux-gnu --with-clanguage=cxx
make all
./config/configure.py  --download-c-blas-lapack=1 --with-mpi-dir=$CHASTE_LIBS/mpi --with-x=false   --with-debugging=0 -PETSC_ARCH=linux-gnu-opt --with-clanguage=cxx
make all
./config/configure.py  --download-c-blas-lapack=1 --with-mpi-dir=$CHASTE_LIBS/mpi --with-x=false  -PETSC_ARCH=linux-gnu-profile --with-clanguage=cxx --CFLAGS=-pg -CXXFLAGS=-pg -LDFLAGS=-pg
make all
cd ..

=========BOOSTJAM and BOOST:===========

cd $CHASTE_LIBS
wget http://fastbull.dl.sourceforge.net/sourceforge/boost/boost-jam-3.1.17-1-linuxx86.tgz
wget http://garr.dl.sourceforge.net/sourceforge/boost/boost_1_33_1.tar.gz
tar -xzvf boost-jam-3.1.17-1-linuxx86.tgz
tar -zxvf boost_1_33_1.tar.gz
cd boost_1_33_1
../boost-jam-3.1.17-1-linuxx86/bjam "-sTOOLS=gcc" --prefix=$CHASTE_LIBS/boost install
cd ..
rm -rf boost_1_33_1.tar.gz boost-jam-3.1.17-1-linuxx86.tgz

==============HDF5:====================

cd $CHASTE_LIBS
wget ftp://ftp.hdfgroup.org/HDF5/prev-releases/hdf5-1.6.6/src/hdf5-1.6.6.tar.gz
tar -zxf hdf5-1.6.6.tar.gz
cd hdf5-1.6.6
CC=mpicc ./configure --enable-parallel --prefix=$CHASTE_LIBS/hdf5
make
cd test
make check
cd ../testpar
make check                                                 ### FAILS HERE
cd ..
make install


================METIS:==================

cd $CHASTE_LIBS
wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.0pre2.tar.gz
tar -zxf metis-5.0pre2.tar.gz
cd metis-5.0pre2

Edit the file Makefile.in, changing the line CC = cc to CC = gcc.
Also edit the file include/metis.h and specify the width (32 or 64 bits) of the elementary data type used in METIS. This is controled by the IDXTYPEWIDTH constant. For now, on a 32 bit architecture you can only specify a width of 32, whereas for a 64 bit architecture you can specify a width of either 32 or 64 bits.
Then run:

make all
cd ..
rm -f metis-5.0pre2.tar.gz

============= XSD and XML =============
Installing Apache Xerces XML libraries:

cd $CHASTE_LIBS
wget http://archive.apache.org/dist/xml/xerces-c/Xerces-C_2_7_0/source/xerces-c-src_2_7_0.tar.gz
tar -zxf xerces-c-src_2_7_0.tar.gz
cd xerces-c-src_2_7_0/
export XERCESCROOT=`pwd`
cd $XERCESCROOT/src/xercesc
./runConfigure -plinux -cgcc -xg++ -P$CHASTE_LIB/xerces
make
make install
cd $CHASTE_LIBS

If you experience any problems then you can modify the configuration options
for your system according to http://xerces.apache.org/xerces-c/build-3.html


Installing XSD:

cd $CHASTE_LIBS
wget http://codesynthesis.com/download/xsd/2.3/linux-gnu/i686/xsd-2.3.1-i686-linux-gnu.tar.bz2
tar -xjf xsd-2.3.1-i686-linux-gnu.tar.bz2
ln -s $CHASTE_LIBS/xsd-2.3.1-i686-linux-gnu/bin/xsd $CHASTE_LIBS/bin/xsd
rm -f xsd-2.3.1-i686-linux-gnu.tar.bz2


=============== TETGEN =====================

cd $CHASTE_LIBS
wget http://www.wias-berlin.de/people/si/tetgen1.3.4.tar.gz
tar -xfvz tetgen1.3.4.tar.gz
cd tetgen1.3.4
make
mv tetgen $CHASTE_LIBS/bin/
cd ..
rm -rf tetgen1.3.4*


===========SET ENVIRONMENTAL VARIABLES AND PATHS:===================

Edit the ~/.bashrc file in the home directory and append the following:

export CHASTE_LIBS=<THE FOLDER YOU CALLED $CHASTE_LIBS ABOVE>
export PATH=$CHASTE_LIBS/bin:$PATH
export LD_LIBRARY_PATH=$CHASTE_LIBS/petsc-2.3.3/lib/libg_c++/linux-gnu/

Note that test output will be put in the folder given by the environment variable CHASTE_TEST_OUTPUT. If this has not been
assigned to anything, it will be (temporarily) set to be the directory "/tmp/<USER_NAME>/testoutput/" whenever scons is run.
If you would prefer output in a different directory also add the following to the .bashrc file

export CHASTE_TEST_OUTPUT=<YOUR CHOICE OF DIRECTORY>


NB. You will need to logout and in again for the above to take effect (or instead do "source ~/.bashrc" to implement the above changes for that particular terminal).



///\todo Check that the below instructions are valid once release is finalised.




===================================================================================
=
=   Running Chaste
=
===================================================================================



========== Compiling and running chaste ==========

First, return to the main Chaste directory (the folder with the code and this README file)

You will need to edit the python/hostconfig/default.py to give the correct paths to each of the libraries.

Now we should have all the necessary libraries to compile chaste. To run the Chaste tests to see if everything compiles and the tests pass
type one of the following (note that scons should always be run from the main Chaste directory):
 
scons          # runs all the tests - will take a while
scons global   # just runs the tests in the global folder
scons io       # just runs the tests in the io folder
scons mesh     # just runs the tests in the mesh folder
scons ode      # just runs the tests in the ode folder
scons pde      # just runs the tests in the pde folder
scons heart    # just runs the tests in the heart folder
scons cancer   # just runs the tests in the cancer folder

If the tests compile and run, to check the results, open a web browser and go to the following (local) page:

file:///<CHASTE_CODE_DIRECTORY>/testoutput/<NAME_OF_YOUR_COMPUTER>.default/index.html

for example: file:///home/scratch/chaste-code/testoutput/msc30.ecs.ox.ac.uk.default/index.html
Note that this is just a file in the Chaste directory.


To run a single test, do (for example):
scons test_suite=heart/test/bidomain/TestBidomainProblem.hpp

To run a single test in parallel with 2 processors, do:
scons build=Parallel test_suite=heart/test/bidomain/TestBidomainProblem.hpp   # 2 processors
scons build=Parallel3 test_suite=heart/test/bidomain/TestBidomainProblem.hpp  # 3 processors

To run a single test with an optimised gcc build, do
scons build=GccOpt test_suite=heart/test/bidomain/TestBidomainProblem.hpp 

To just compile a single test, do (for example):
scons test_suite=heart/test/bidomain/TestBidomainProblem.hpp compile_only=1

///\todo mention the existence of an intel production build and give a link to installation instructions??

Note that if you scons build=<SOMETHING>, the webpage to go to is
file:///<CHASTE_CODE_DIRECTORY>/testoutput/<NAME_OF_YOUR_COMPUTER>.<SOMETHING>/index.html
for example: file:///home/scratch/chaste-code/testoutput/msc30.ecs.ox.ac.uk.GccOpt/index.html




========== Creating an executable ==========

*IF* you want to create the standalone cardiac chaste executable:

The following commands should generate static and dynamic chaste libraries.

scons compile_only=1 chaste_libs=1 static=0 build=GccOpt exe=1 apps
scons compile_only=1 chaste_libs=1 static=1 build=GccOpt exe=1 apps

The executable Chaste or Chaste.o can be found in the apps/src folder.

In order to run a simulation edit the ChasteParameters.xml file according to your needs and, 
from the chaste directory, type

apps/src/Chaste ChasteParameters.xml

