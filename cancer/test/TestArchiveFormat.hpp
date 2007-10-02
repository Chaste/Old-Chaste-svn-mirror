#ifndef TESTARCHIVEFORMAT_HPP_
#define TESTARCHIVEFORMAT_HPP_


#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "MeinekeCryptCell.hpp"
#include "WntCellCycleModel.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "AbstractCellKiller.hpp"
#include "SloughingCellKiller.hpp"

class TestArchiveFormat : public CxxTest::TestSuite
{
public:
/**
 * This test is required because Test2DCryptRepresentativeSimulationn loads an archive
 * stored in cancer/test/data. When the archiving of TissueSimulation and associate classes
 * the stored archive needs to be update. This test checks that the archive can be loaded,
 * and will seg fault if not.
 * It does nothing more, so it runs quickly and can be in the continuous test pack.
 * 
 * 
 * When the archive file changes, re-run TestMakeNiceCryptSims, having commented out the
 * for loop which loads, runs and saves the simulation. One should be left with
 * the archive after the first save which can be copied to cancer/test/data.
 */
    void TestLoadArchive() throw (Exception)
    {        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        std::string test_to_profile = "NiceCryptSim";
        double t = 350;   // this is the folder and time that the stored results were archived (needed to know foldernames)
        
        // The archive needs to be copied from cancer/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation.     
        OutputFileHandler any_old_handler("",false);
        std::string test_output_directory = any_old_handler.GetTestOutputDirectory();
        std::string test_data_directory = "cancer/test/data/" + test_to_profile +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +" "+ test_output_directory +"/";     
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);
        
        TissueSimulation<2>* p_simulator = TissueSimulation<2>::Load(test_to_profile,t);
        delete p_simulator;
                
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTARCHIVEFORMAT_HPP_*/
