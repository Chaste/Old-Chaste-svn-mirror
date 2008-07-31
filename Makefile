INCS = -I. -Icxxtest \
-Icancer/src -Icancer/src/common -Icancer/src/mesh -Icancer/src/odes -Icancer/src/tissue -Icancer/src/tissue/statistics -Icancer/src/tissue/mechanics -Icancer/src/tissue/cell -Icancer/src/tissue/killers -Icancer/src/tissue/cell/cycle \
-Iglobal/src \
-Ilinalg/src -Ilinalg/src/common \
-Ipde/src -Ipde/src/problem -Ipde/src/common -Ipde/src/solver -Ipde/src/problem/elasticity -Ipde/src/problem/common -Ipde/src/solver/elasticity -Ipde/src/solver/common \
-Imesh/src -Imesh/src/common -Imesh/src/voronoi -Imesh/src/writer -Imesh/src/reader -Imesh/src/triangle \
-Iio/src -Iio/src/writer -Iio/src/reader \
-Iode/src -Iode/src/problem -Iode/src/common -Iode/src/solver -Iode/src/problem/cardiac \

INCS += -I/opt/boost/include/boost-1_33_1

LIBS=cancer/src/common/CancerParameters.o \
cancer/src/common/CancerEventHandler.o \
cancer/src/mesh/Cylindrical2dMesh.o \
cancer/src/mesh/HoneycombMeshGenerator.o \
cancer/src/tissue/statistics/AbstractCryptStatistics.o \
cancer/src/tissue/statistics/CryptStatistics.o \
cancer/src/tissue/cell/cycle/AbstractCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractSimpleCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractSimpleMeinekeCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractOdeBasedCellCycleModel.o \
cancer/src/tissue/cell/cycle/AbstractWntOdeBasedCellCycleModel.o \
cancer/src/tissue/cell/cycle/TysonNovakCellCycleModel.o \
cancer/src/tissue/cell/cycle/SimpleWntCellCycleModel.o \
cancer/src/tissue/cell/cycle/StochasticCellCycleModel.o \
cancer/src/tissue/cell/cycle/StochasticWntCellCycleModel.o \
cancer/src/tissue/cell/cycle/FixedCellCycleModel.o \
cancer/src/tissue/cell/cycle/WntCellCycleModel.o \
cancer/src/tissue/cell/cycle/WntConcentration.o \
cancer/src/tissue/cell/cycle/IngeWntSwatCellCycleModel.o \
cancer/src/tissue/cell/TissueCell.o \
cancer/src/tissue/killers/RadialSloughingCellKiller.o \
cancer/src/tissue/killers/SloughingCellKiller.o \
cancer/src/tissue/mechanics/CryptProjectionSpringSystem.o \
cancer/src/odes/WntCellCycleOdeSystem.o \
cancer/src/odes/IngeWntSwatCellCycleOdeSystem.o \
cancer/src/odes/TysonNovak2001OdeSystem.o \
global/src/LogFile.o \
global/src/Exception.o  \
global/src/OutputFileHandler.o \
global/src/RandomNumberGenerator.o \
global/src/SimulationTime.o \
global/src/TimeStepper.o \
io/src/writer/ColumnDataWriter.o \
mesh/src/triangle/triangle.o \
ode/src/solver/AbstractOneStepIvpOdeSolver.o \
ode/src/solver/RungeKutta4IvpOdeSolver.o \
ode/src/common/AbstractOdeSystem.o \

CXXFLAGS = -DTRILIBRARY -DANSI_DECLARATORS -DSPECIAL_SERIAL -O3 -march=opteron ${INCS} 
#On userpc44
LDFLAGS =   -lboost_serialization

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
	mkdir ${FRESH_DIR} ;\
	cp TestGenerateSteadyStateCryptRunner ${FRESH_DIR} ;\
	cp simulationGenerateSteadyStateCrypt.sh ${FRESH_DIR}  ;\
	mv ${FRESH_DIR}/simulationGenerateSteadyStateCrypt.sh ${FRESH_DIR}/simulation.sh

# This test generates Inge cells crypt

TestMakeNiceCryptSimsAlexWRunner.cpp:	projects/AlexW_New/test/TestMakeNiceCryptSimsAlexW_proj.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMakeNiceCryptSimsAlexWRunner.cpp projects/AlexW_New/test/TestMakeNiceCryptSimsAlexW_proj.hpp

TestMakeNiceCryptSimsAlexWRunner: TestMakeNiceCryptSimsAlexWRunner.o ${LIBS}
	g++ TestMakeNiceCryptSimsAlexWRunner.o ${LIBS} -o TestMakeNiceCryptSimsAlexWRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxaw@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ;\
	cp TestMakeNiceCryptSimsAlexWRunner ${FRESH_DIR} ;\
	cp simulationMakeNiceCryptSimsAlexW.sh ${FRESH_DIR}  ;\
	mv ${FRESH_DIR}/simulationMakeNiceCryptSimsAlexW.sh ${FRESH_DIR}/simulation.sh
	
# This test generates is same as above but with copying for yes area-drag, yes edge-based springs

TestMakeYesYesRunner.cpp:	projects/AlexW_New/test/TestMakeNiceCryptSimsAlexW_proj.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMakeYesYesRunner.cpp projects/AlexW_New/test/TestMakeNiceCryptSimsAlexW_proj.hpp

TestMakeYesYesRunner: TestMakeYesYesRunner.o ${LIBS}
	g++ TestMakeYesYesRunner.o ${LIBS} -o TestMakeYesYesRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxaw@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ;
	# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/Noddy_Yes_Yes; mkdir ${FRESH_DIR}/Noddy_Yes_Yes/archive ;\
	cd ${FRESH_DIR}/Noddy_Yes_Yes/archive ;\
	cp ../../../projects/AlexW_New/test/data/SteadyStateSimpleWnt/Noddy_Yes_Yes/archive/mesh_300.* . ;\
	cp ../../../projects/AlexW_New/test/data/SteadyStateSimpleWnt/Noddy_Yes_Yes/archive/tissue_sim_at_time_300.arch . ;\
	cd ../.. ;\
	# Finished copying archives across.
	cp TestMakeYesYesRunner ${FRESH_DIR} ;\
	cp simulationMakeYesYes.sh ${FRESH_DIR}  ;\
	mv ${FRESH_DIR}/simulationMakeYesYes.sh ${FRESH_DIR}/simulation.sh
	
# This test generates is same as above but with copying for No area-drag, No edge-based springs

TestMakeNoNoRunner.cpp:	projects/AlexW_New/test/TestMakeNiceCryptSimsAlexW_proj.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMakeNoNoRunner.cpp projects/AlexW_New/test/TestMakeNiceCryptSimsAlexW_proj.hpp

TestMakeNoNoRunner: TestMakeNoNoRunner.o ${LIBS}
	g++ TestMakeNoNoRunner.o ${LIBS} -o TestMakeNoNoRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxaw@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ;\
	# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/Noddy_No_No; mkdir ${FRESH_DIR}/Noddy_No_No/archive ;\
	cd ${FRESH_DIR}/Noddy_No_No/archive ;\
	cp ../../../projects/AlexW_New/test/data/SteadyStateSimpleWnt/Noddy_No_No/archive/mesh_300.* . ;\
	cp ../../../projects/AlexW_New/test/data/SteadyStateSimpleWnt/Noddy_No_No/archive/tissue_sim_at_time_300.arch . ;\
	cd ../.. ;\
	# Finished copying archives across.
	cp TestMakeNoNoRunner ${FRESH_DIR} ;\
	cp simulationMakeNoNo.sh  ${FRESH_DIR} ;\
	mv ${FRESH_DIR}/simulationMakeNoNo.sh ${FRESH_DIR}/simulation.sh	
	
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
	mkdir ${FRESH_DIR} ;\
	# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/MeinekeLabellingExperiment; mkdir ${FRESH_DIR}/MeinekeLabellingExperiment/archive ;\
	cd ${FRESH_DIR}/MeinekeLabellingExperiment/archive ;\
	cp ../../../projects/GaryM/test/data/SteadyStateMeinekeStochastic/sunter1_archive/mesh_300.* . ;\
	cp ../../../projects/GaryM/test/data/SteadyStateMeinekeStochastic/sunter1_archive/tissue_sim_at_time_300.arch . ;\
	cd ../.. ;\
	# Finished copying archives across.
	cp TestMeinekeLabellingExperimentsRunner ${FRESH_DIR} ;\
	cp simulationMeinekeLabellingExperiments.sh ${FRESH_DIR} ;\
	mv ${FRESH_DIR}/simulationMeinekeLabellingExperiments.sh ${FRESH_DIR}/simulation.sh

# This runs the test which generates MeinekeLabellingExperiment for AlexWdata.

TestMeinekeLabellingAlexWRunner.cpp:	projects/AlexW_New/test/TestSimpleWntLabellingExperimentsSunterDataAlexW.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMeinekeLabellingAlexWRunner.cpp projects/AlexW_New/test/TestSimpleWntLabellingExperimentsSunterDataAlexW.hpp

TestMeinekeLabellingAlexWRunner: TestMeinekeLabellingAlexWRunner.o ${LIBS}
	g++ TestMeinekeLabellingAlexWRunner.o ${LIBS} -o TestMeinekeLabellingAlexWRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxaw@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ;\
	# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/MeinekeLabellingAlexW; mkdir ${FRESH_DIR}/MeinekeLabellingAlexW/archive ;\
	cd ${FRESH_DIR}/MeinekeLabellingAlexW/archive ;\
	cp ../../../projects/AlexW_New/test/data/SteadyStateSimpleWnt/sunter1_archive/mesh_300.* . ;\
	cp ../../../projects/AlexW_New/test/data/SteadyStateSimpleWnt/sunter1_archive/tissue_sim_at_time_300.arch . ;\
	cd ../.. ;\
	# Finished copying archives across.
	cp TestMeinekeLabellingAlexWRunner ${FRESH_DIR} ;\
	cp simulationMeinekeLabellingAlexW.sh ${FRESH_DIR}  ;\
	mv ${FRESH_DIR}/simulationMeinekeLabellingAlexW.sh ${FRESH_DIR}/simulation.sh

# A more useful test to label a cell near the bottom at random and follow mutation's progress.

TestMutationSpreadRunner.cpp:	projects/GaryM/test/TestMutationSpread.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestMutationSpreadRunner.cpp projects/GaryM/test/TestMutationSpread.hpp

TestMutationSpreadRunner: TestMutationSpreadRunner.o ${LIBS}
	g++ TestMutationSpreadRunner.o ${LIBS} -o TestMutationSpreadRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ;
# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/MutationSpread; mkdir ${FRESH_DIR}/MutationSpread/archive ;\
	cd ${FRESH_DIR}/MutationSpread/archive ;\
	cp ../../../projects/GaryM/test/data/SteadyStateMeinekeStochastic/sunter3_archive/* . ;\
	cd ../.. ;\
# Finished copying archives across.
	cp TestMutationSpreadRunner ${FRESH_DIR} ;\
	cp simulationMutationSpread.sh ${FRESH_DIR} ;\
	mv ${FRESH_DIR}/simulationMutationSpread.sh ${FRESH_DIR}/simulation.sh 
# End of different test.

# A test to generate graphs of the time it takes for monoclonality
TestNicheSuccessionTimeDistributionsRunner.cpp:	projects/AlexF/test/TestNicheSuccessionTimeDistributions.hpp
	cxxtest/cxxtestgen.py  --error-printer -o TestNicheSuccessionTimeDistributionsRunner.cpp projects/AlexF/test/TestNicheSuccessionTimeDistributions.hpp

TestNicheSuccessionTimeDistributionsRunner: TestNicheSuccessionTimeDistributionsRunner.o ${LIBS}
	g++ TestNicheSuccessionTimeDistributionsRunner.o ${LIBS} -o TestNicheSuccessionTimeDistributionsRunner ${LDFLAGS};\
	echo "Making new experiment in ${FRESH_DIR} " ;\
	echo "Do scp -r -C ${FRESH_DIR} pmxgm@deimos.nottingham.ac.uk:" ;\
	echo "Then qsub simulation.sh on deimos";\
	echo "If 'owt funny happens when this is compiling type 'make clean' to do this from fresh" ;\
	mkdir ${FRESH_DIR} ;\
# Need to copy across the starting state of the simulation
	mkdir ${FRESH_DIR}/NicheSuccessionTime; mkdir ${FRESH_DIR}/NicheSuccessionTime/archive ;\
	cd ${FRESH_DIR}/NicheSuccessionTime/archive ;\
	cp ../../../projects/AlexF/test/data/SteadyStateCryptMeinekeGeometry/archive/* . ;\
	cd ../.. ;\
# Need the triangle binary for crypt projection simulations
	mkdir ${FRESH_DIR}/bin ;\
    cd ${FRESH_DIR}/bin ;\
	cp ../../bin/triangle triangle ;\
	cd .. ;\
# Finished copying archives across.
	cp TestNicheSuccessionTimeDistributionsRunner ${FRESH_DIR} ;\
	cp simulationNicheSuccession.sh ${FRESH_DIR} ;\
	mv ${FRESH_DIR}/simulationNicheSuccession.sh ${FRESH_DIR}/simulation.sh 
# End of monoclonality time test.

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
