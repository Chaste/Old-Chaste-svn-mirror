INCS = -I. -Icxxtest -Ipde/src -Ipde/src/problem -Ipde/src/common -Ipde/src/solver -Ipde/src/problem/elasticity -Ipde/src/problem/common \
-Ipde/src/solver/elasticity -Ipde/src/solver/common -Imesh/src -Imesh/src/voronoi -Imesh/src/writer -Imesh/src/common -Imesh/src/reader \
-Imesh/src/decimator  \
-Iio/src -Iio/src/writer -Iio/src/reader \
-Iode/src -Iode/src/problem -Iode/src/common -Iode/src/solver -Iode/src/problem/cancer -Iode/src/problem/cardiac -Iglobal/src -Ilinalg/src \
-Ilinalg/src/common -Icoupled/src -Icoupled/src/problem -Icoupled/src/common -Icoupled/src/solver -Icoupled/src/problem/cancer \
-Icoupled/src/problem/elasticity -Icoupled/src/problem/cardiac -Icoupled/src/solver/cancer -Icoupled/src/solver/elasticity \
-Icoupled/src/solver/cardiac -Imodels/src -Imodels/src/crypt -Imodels/src/crypt/cell -Imodels/src/crypt/killers -Imodels/src/crypt/cell/cycle

INCS += -I/opt/boost/include/boost-1_33_1

LIBS=global/src/Exception.o  \
global/src/LogFile.o \
global/src/CancerParameters.o \
models/src/crypt/cell/cycle/AbstractCellCycleModel.o \
models/src/crypt/cell/cycle/TysonNovakCellCycleModel.o \
models/src/crypt/cell/cycle/StochasticCellCycleModel.o \
models/src/crypt/cell/cycle/FixedCellCycleModel.o \
models/src/crypt/cell/cycle/WntCellCycleModel.o \
models/src/crypt/cell/MeinekeCryptCell.o \
ode/src/problem/cancer/WntCellCycleOdeSystem.o \
models/src/crypt/WntGradient.o \
global/src/OutputFileHandler.o \
global/src/RandomNumberGenerator.o \
global/src/SimulationTime.o \
global/src/TimeStepper.o \
io/src/writer/ColumnDataWriter.o \
ode/src/problem/cancer/TysonNovak2001OdeSystem.o \
ode/src/solver/AbstractOneStepIvpOdeSolver.o \
ode/src/solver/RungeKutta4IvpOdeSolver.o \
ode/src/common/AbstractOdeSystem.o \

CXXFLAGS = -DSPECIAL_SERIAL -O3 ${INCS}

#On userpc44
#LDFLAGS =   -lboost_serialization

#On engels in Nottingham
LDFLAGS =   -L/opt/boost/lib -lboost_serialization-gcc

default:	TestMakeNiceCryptSimsRunner TestCryptSimulation2DPeriodicRunner

FRESH_DIR=`date +%F-%H-%M`

# This test generates the archives which are used in the profiling test Test2DCryptRepresentativeSimulation.hpp

TestMakeNiceCryptSimsRunner.cpp:	models/test/TestMakeNiceCryptSims.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMakeNiceCryptSimsRunner.cpp models/test/TestMakeNiceCryptSims.hpp

TestMakeNiceCryptSimsRunner: TestMakeNiceCryptSimsRunner.o ${LIBS}
	g++ TestMakeNiceCryptSimsRunner.o ${LIBS} -o TestMakeNiceCryptSimsRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
	cp TestMakeNiceCryptSimsRunner ${FRESH_DIR} ;\
	cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
	cp ../simulationNiceCryptSims.sh .  ;\
	mv simulationNiceCryptSims.sh simulation.sh

TestCryptSimulation2DPeriodicRunner.cpp:	models/test/TestCryptSimulation2DPeriodic.hpp
	cxxtest/cxxtestgen.py  --error-printer  -o TestCryptSimulation2DPeriodicRunner.cpp models/test/TestCryptSimulation2DPeriodic.hpp

TestCryptSimulation2DPeriodicRunner: TestCryptSimulation2DPeriodicRunner.o ${LIBS}

# A more useful test to label a cell near the bottom at random and follow mutation's progress.

TestMutationSpreadRunner.cpp:	models/test/TestMutationSpread.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMutationSpreadRunner.cpp models/test/TestMutationSpread.hpp

TestMutationSpreadRunner: TestMutationSpreadRunner.o ${LIBS}
	g++ TestMutationSpreadRunner.o ${LIBS} -o TestMutationSpreadRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	mkdir ${FRESH_DIR} ; mkdir ${FRESH_DIR}/bin ;\
# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/NiceCryptSim; mkdir ${FRESH_DIR}/NiceCryptSim/archive ;\
	cd ${FRESH_DIR}/NiceCryptSim/archive ;\
	cp ../../../models/test/data/NiceCryptSim/archive/* . ;\
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
