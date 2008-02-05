#ifndef TESTARCHIVEFORMAT_HPP_
#define TESTARCHIVEFORMAT_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>

#include "CryptSimulation2d.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "TissueCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "TrianglesMeshReader.cpp"
#include "ColumnDataReader.hpp"
#include "OutputFileHandler.hpp"


class TestArchiveFormat : public CxxTest::TestSuite
{
public:

    /**
     * This test is required because Test2DCryptRepresentativeSimulation loads 
     * an archive stored in cancer/test/data. When the archiving of 
     * TissueSimulation and associate classes is updated the stored archive 
     * needs to be updated. This test checks that the archive can be loaded,
     * and will seg fault if not. It does nothing more, so it runs quickly 
     * and can be in the continuous test pack.
     * 
     * IF THIS TEST FAILS:
     * - You have probably changed an archiving function somewhere
     * - You need to remake cancer/test/data/<test below>/archive/
     * - To do this re-run TestGenerateSteadyStateCrypt.hpp
     * - Archives produced can be copied to :
     *   cancer/test/data/<test below>/archive/
     * 
     * (it is a long test, currently just < 5hours, and could be 
     * run overnight - please do this rather than just moving it 
     * to the failing test pack(!) because these files are now 
     * the basis of some proper simulations for
     * the papers that are on the way...) 
     */
    void TestLoadArchive() throw (Exception)
    {
        SimulationTime::Instance()->SetStartTime(0.0);

        std::string test_to_profile = "SteadyStateCrypt";
        double t = 150;   // this is the folder and time that the stored results were archived (needed to know foldernames)
        
        // Open a new directory...
        OutputFileHandler file_handler(test_to_profile,true);   
        
        // The archive needs to be copied from cancer/test/data/<test_to_profile>
        // to the testoutput directory to continue running the simulation.     
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string test_data_directory = "cancer/test/data/" + test_to_profile +"/";
        std::string command = "cp -Rf --remove-destination " + test_data_directory +" "+ test_output_directory +"/";     
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value, 0);
        
        CryptSimulation2d* p_simulator = CryptSimulation2d::Load(test_to_profile,t);
        p_simulator->SetEndTime(t + 1);
        delete p_simulator;
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTARCHIVEFORMAT_HPP_*/
