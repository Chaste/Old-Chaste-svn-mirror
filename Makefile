LOCALINCS = -I. -Icxxtest -Ipde/src -Ipde/src/problem -Ipde/src/common -Ipde/src/solver -Ipde/src/problem/elasticity -Ipde/src/problem/common \
-Ipde/src/solver/elasticity -Ipde/src/solver/common -Imesh/src -Imesh/src/voronoi -Imesh/src/writer -Imesh/src/common -Imesh/src/reader \
-Imesh/src/decimator  \
-Iio/src -Iio/src/writer -Iio/src/reader \
-Iode/src -Iode/src/problem -Iode/src/common -Iode/src/solver -Iode/src/problem/cancer -Iode/src/problem/cardiac -Iglobal/src -Ilinalg/src \
-Ilinalg/src/common -Icoupled/src -Icoupled/src/problem -Icoupled/src/common -Icoupled/src/solver -Icoupled/src/problem/cancer \
-Icoupled/src/problem/elasticity -Icoupled/src/problem/cardiac -Icoupled/src/solver/cancer -Icoupled/src/solver/elasticity \
-Icoupled/src/solver/cardiac -Imodels/src -Imodels/src/crypt -Imodels/src/crypt/cell -Imodels/src/crypt/killers -Imodels/src/crypt/cell/cycle

INCS=${LOCALINCS} -I/home/chaste/petsc-2.2.1/include/  -I/home/chaste/petsc-2.2.1/bmake/linux-gnu/ -I/home/chaste/petsc-2.2.1/include/mpiuni/

LIBS=global/src/Exception.o  \
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
LDFLAGS =  -L/home/chaste/petsc-2.2.1/lib/libg/linux-gnu/ -lboost_serialization -lpetsc -lmpiuni


default:	TestMakeNiceCryptSimsRunner TestCryptSimulation2DPeriodicRunner



TestMakeNiceCryptSimsRunner.cpp:	models/test/TestMakeNiceCryptSims.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMakeNiceCryptSimsRunner.cpp models/test/TestMakeNiceCryptSims.hpp

TestMakeNiceCryptSimsRunner: TestMakeNiceCryptSimsRunner.o ${LIBS}


TestCryptSimulation2DPeriodicRunner.cpp:	models/test/TestCryptSimulation2DPeriodic.hpp
	cxxtest/cxxtestgen.py  --error-printer  -o TestCryptSimulation2DPeriodicRunner.cpp models/test/TestCryptSimulation2DPeriodic.hpp

TestCryptSimulation2DPeriodicRunner: TestCryptSimulation2DPeriodicRunner.o ${LIBS}

clean:
	rm *.o  */src/*/*.o */src/*/*/*.o */src/*/*/*/*.o Test*.cpp
