INCS = -I. -Icxxtest \
-Icancer/src -Icancer/src/common -Icancer/src/mesh -Icancer/src/odes -Icancer/src/tissue -Icancer/src/tissue/cell -Icancer/src/tissue/killers -Icancer/src/tissue/cell/cycle \
-Iglobal/src \
-Ilinalg/src -Ilinalg/src/common \
-Ipde/src -Ipde/src/problem -Ipde/src/common -Ipde/src/solver -Ipde/src/problem/elasticity -Ipde/src/problem/common -Ipde/src/solver/elasticity -Ipde/src/solver/common \
-Imesh/src -Imesh/src/common -Imesh/src/voronoi -Imesh/src/writer -Imesh/src/reader  \
-Iio/src -Iio/src/writer -Iio/src/reader \
-Iode/src -Iode/src/problem -Iode/src/common -Iode/src/solver -Iode/src/problem/cardiac \

INCS += -I/opt/boost/include/boost-1_33_1

LIBS=cancer/src/common/CancerParameters.o \
cancer/src/common/CancerEventHandler.o \
cancer/src/tissue/CryptStatistics.o \
cancer/src/tissue/cell/cycle/AbstractCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractSimpleCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractSimpleMeinekeCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractOdeBasedCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractWntOdeBasedCellCycleModel.o \
cancer/src/tissue/cell/cycle/TysonNovakCellCycleModel.o \
cancer/src/tissue/cell/cycle/SimpleWntCellCycleModel.o \
cancer/src/tissue/cell/cycle/StochasticCellCycleModel.o \
cancer/src/tissue/cell/cycle/FixedCellCycleModel.o \
cancer/src/tissue/cell/cycle/WntCellCycleModel.o \
cancer/src/tissue/cell/cycle/IngeWntSwatCellCycleModel.o \
cancer/src/tissue/cell/TissueCell.o \
cancer/src/odes/WntCellCycleOdeSystem.o \
cancer/src/odes/IngeWntSwatCellCycleOdeSystem.o \
cancer/src/odes/TysonNovak2001OdeSystem.o \
cancer/src/tissue/cell/cycle/WntGradient.o \
global/src/LogFile.o \
global/src/Exception.o  \
global/src/OutputFileHandler.o \
global/src/RandomNumberGenerator.o \
global/src/SimulationTime.o \
global/src/TimeStepper.o \
io/src/writer/ColumnDataWriter.o \
ode/src/solver/AbstractOneStepIvpOdeSolver.o \
ode/src/solver/RungeKutta4IvpOdeSolver.o \
ode/src/common/AbstractOdeSystem.o \

CXXFLAGS = -DSPECIAL_SERIAL -O3 -march=opteron ${INCS} 
#On userpc44
#LDFLAGS =   -lboost_serialization

#On engels in Nottingham
LDFLAGS =   -L/opt/boost/lib -lboost_serialization-gcc

default:	TestGenerateSteadyStateCryptRunner TestCryptSimulation2dRunner

FRESH_DIR=`date +%F-%H-%M`

# This test generates the archives which are used in the profiling test Test2DCryptRepresentativeSimulation.hpp

TestGenerateSteadyStateCryptRunner.cpp:	cancer/test/TestGenerateSteadyStateCrypt.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestGenerateSteadyStateCryptRunner.cpp cancer/test/TestGenerateSteadyStateCrypt.hpp

TestGenerateSteadyStateCryptRunner: TestGenerateSteadyStateCryptRunner.o ${LIBS}
	g++ TestGenerateSteadyStateCryptRunner.o ${LIBS} -o TestGenerateSteadyStateCryptRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
	cp TestGenerateSteadyStateCryptRunner ${FRESH_DIR} ;\
	cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
	cp ../simulationGenerateSteadyStateCrypt.sh .  ;\
	mv simulationGenerateSteadyStateCrypt.sh simulation.sh
	
TestCryptSimulation2dRunner.cpp:	cancer/test/TestCryptSimulation2d.hpp
	cxxtest/cxxtestgen.py  --error-printer  -o TestCryptSimulation2dRunner.cpp cancer/test/TestCryptSimulation2d.hpp

TestCryptSimulation2dRunner: TestCryptSimulation2dRunner.o ${LIBS}
	
# This runs the test which generates MeinekeLabellingExperiment data.

TestMeinekeLabellingExperimentsRunner.cpp:	projects/GaryM/test/TestMeinekeLabellingExperimentsSunterData.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMeinekeLabellingExperimentsRunner.cpp projects/GaryM/test/TestMeinekeLabellingExperimentsSunterData.hpp

TestMeinekeLabellingExperimentsRunner: TestMeinekeLabellingExperimentsRunner.o ${LIBS}
	g++ TestMeinekeLabellingExperimentsRunner.o ${LIBS} -o TestMeinekeLabellingExperimentsRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
	# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/MeinekeLabellingExperiment; mkdir ${FRESH_DIR}/MeinekeLabellingExperiment/archive ;\
	cd ${FRESH_DIR}/MeinekeLabellingExperiment/archive ;\
	cp ../../../projects/GaryM/test/data/SteadyStateSimpleWnt/sunter1_archive/* . ;\
	cd ../.. ;\
	# Finished copying archives across.
	cp TestMeinekeLabellingExperimentsRunner ${FRESH_DIR} ;\
	cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
	cp ../simulationMeinekeLabellingExperiments.sh .  ;\
	mv simulationMeinekeLabellingExperiments.sh simulation.sh

# A more useful test to label a cell near the bottom at random and follow mutation's progress.

TestMutationSpreadRunner.cpp:	projects/GaryM/test/TestMutationSpread.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMutationSpreadRunner.cpp projects/GaryM/test/TestMutationSpread.hpp

TestMutationSpreadRunner: TestMutationSpreadRunner.o ${LIBS}
	g++ TestMutationSpreadRunner.o ${LIBS} -o TestMutationSpreadRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/MutationSpread; mkdir ${FRESH_DIR}/MutationSpread/archive ;\
	cd ${FRESH_DIR}/MutationSpread/archive ;\
	cp ../../../projects/GaryM/test/data/SteadyStateSimpleWnt/sunter3_archive/* . ;\
	cd ../.. ;\
# Finished copying archives across.
	cp TestMutationSpreadRunner ${FRESH_DIR} ;\
	cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
	cp ../simulationMutationSpread.sh . ;\
	mv simulationMutationSpread.sh simulation.sh 
# End of different test.


FULL_INCS = -isystem /home/chaste/petsc-2.3.2-p4/bmake/linux-intel-opt-mkl \
-isystem /home/chaste/petsc-2.3.2-p4/include \
-isystem ../../../xsd-2.3.1-i686-linux-gnu/libxsd \
${INCS} \
-I heart/src/problem \
-I heart/src/solver \
-I heart/src/pdes \
-I heart/src/odes \
-I heart/src/stimulus \
-I heart/src/io \
-I global/src

FULL_LIBS= -Llinklib -Lheart/build/intel -Lheart \
-L/home/chaste/petsc-2.3.2-p4/lib/linux-intel-opt-mkl \
-L/opt/intel/cc/9.1.039/lib -L/home/chaste/petsc-2.3.2-p4/externalpackages/f2cblaslapack/linux-gnu \
-L/opt/intel/mkl/9.1.023/lib/32 \
-Llib -ltestheart -lheart -lode -lmesh -llinalg -lio -lglobal \
-lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec\
-lpetsc -lmkl_lapack -lmkl -lboost_serialization -lxerces-c

Chaste:	Chaste.o heart/build/intel/src/io/ChasteParameters.o
	mpicxx -CC=icpc -o Chaste Chaste.o heart/build/intel/src/io/ChasteParameters.o ${FULL_LIBS}

Chaste.o: Chaste.cpp #Need more dependencies here! 
	mpicxx -CC=icpc ${FULL_INCS} -c Chaste.cpp

heart/build/intel/src/io/ChasteParameters.o:	heart/src/io/ChasteParameters.cpp
	mpicxx -CC=icpc ${FULL_INCS} -c  heart/src/io/ChasteParameters.cpp -o heart/build/intel/src/io/ChasteParameters.o

clean:
	rm *.o  */src/*/*.o */src/*/*/*.o */src/*/*/*/*.o Test*.cpp
