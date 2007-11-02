INCS = -I. -Icxxtest \
-Icancer/src -Icancer/src/common -Icancer/src/mesh -Icancer/src/odes -Icancer/src/tissue -Icancer/src/tissue/cell -Icancer/src/tissue/killers -Icancer/src/tissue/cell/cycle \
-Iglobal/src \
-Ilinalg/src -Ilinalg/src/common \
-Ipde/src -Ipde/src/problem -Ipde/src/common -Ipde/src/solver -Ipde/src/problem/elasticity -Ipde/src/problem/common -Ipde/src/solver/elasticity -Ipde/src/solver/common \
-Imesh/src -Imesh/src/common -Imesh/src/voronoi -Imesh/src/writer -Imesh/src/reader -Imesh/src/decimator  \
-Iio/src -Iio/src/writer -Iio/src/reader \
-Iode/src -Iode/src/problem -Iode/src/common -Iode/src/solver -Iode/src/problem/cardiac \

INCS += -I/opt/boost/include/boost-1_33_1

LIBS=cancer/src/common/CancerParameters.o \
cancer/src/tissue/cell/cycle/AbstractCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractSimpleCellCycleModel.o \
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

CXXFLAGS = -DSPECIAL_SERIAL -O3 ${INCS}

#On userpc44
#LDFLAGS =   -lboost_serialization

#On engels in Nottingham
LDFLAGS =   -L/opt/boost/lib -lboost_serialization-gcc

default:	TestMakeNiceCryptSimsRunner TestCryptSimulation2dRunner

FRESH_DIR=`date +%F-%H-%M`

# This test generates the archives which are used in the profiling test Test2DCryptRepresentativeSimulation.hpp

TestMakeNiceCryptSimsRunner.cpp:	cancer/test/TestMakeNiceCryptSims.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMakeNiceCryptSimsRunner.cpp cancer/test/TestMakeNiceCryptSims.hpp

TestMakeNiceCryptSimsRunner: TestMakeNiceCryptSimsRunner.o ${LIBS}
	g++ TestMakeNiceCryptSimsRunner.o ${LIBS} -o TestMakeNiceCryptSimsRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
	cp TestMakeNiceCryptSimsRunner ${FRESH_DIR} ;\
	cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
	cp ../simulationNiceCryptSims.sh .  ;\
	mv simulationNiceCryptSims.sh simulation.sh

TestCryptSimulation2dRunner.cpp:	cancer/test/TestCryptSimulation2d.hpp
	cxxtest/cxxtestgen.py  --error-printer  -o TestCryptSimulation2dRunner.cpp cancer/test/TestCryptSimulation2d.hpp

TestCryptSimulation2dRunner: TestCryptSimulation2dRunner.o ${LIBS}

# A more useful test to label a cell near the bottom at random and follow mutation's progress.

TestMutationSpreadRunner.cpp:	projects/GaryM/test/TestMutationSpread.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMutationSpreadRunner.cpp cancer/test/TestMutationSpread.hpp

TestMutationSpreadRunner: TestMutationSpreadRunner.o ${LIBS}
	g++ TestMutationSpreadRunner.o ${LIBS} -o TestMutationSpreadRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/NiceCryptSim; mkdir ${FRESH_DIR}/NiceCryptSim/archive ;\
	cd ${FRESH_DIR}/NiceCryptSim/archive ;\
	cp ../../../cancer/test/data/NiceCryptSim/archive/* . ;\
	cd ../.. ;\
# Finished copying archives across.
	cp TestMutationSpreadRunner ${FRESH_DIR} ;\
	cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
	cp ../simulationMutationSpread.sh . ;\
	mv simulationMutationSpread.sh simulation.sh 

# End of different test.

clean:
	rm *.o  */src/*/*.o */src/*/*/*.o */src/*/*/*/*.o Test*.cpp
