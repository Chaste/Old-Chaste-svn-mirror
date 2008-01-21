#ifndef TESTREPRESENTATIVESIMULATION_HPP_
#define TESTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "CryptSimulation2d.hpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
//#include "TissueCell.hpp"
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

class TestRepresentativeSimulation : public CxxTest::TestSuite
{
public:
void TestRepresentativeSimulationForProfiling() throw (Exception)
    {        
        SimulationTime::Instance()->SetStartTime(0.0);

        std::string test_to_load = "SteadyStateCrypt";
        std::string test_to_profile = "CryptProfiling";
        double t = 150;   // this is the folder and time that the stored results were archived (needed to know foldernames)
        double run_for = 10; // run for 10 hours.
        
        // create a new clean directory...
        OutputFileHandler file_handler(test_to_profile,true);   
        
        // The archive needs to be copied from cancer/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation.     
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string test_data_directory = "cancer/test/data/" + test_to_load +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +"* "+ test_output_directory +"/" + test_to_profile + "/";     
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);
        
        CryptSimulation2d* p_simulator = CryptSimulation2d::Load(test_to_profile,t);
        p_simulator->SetEndTime(t+run_for); // start time + duration
        p_simulator->Solve();
        delete p_simulator;
                
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTREPRESENTATIVESIMULATION_HPP_*/
