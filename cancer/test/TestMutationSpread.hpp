#ifndef TESTMUTATIONSPREAD_HPP_
#define TESTMUTATIONSPREAD_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "CryptSimulation2d.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "OutputFileHandler.hpp"
#include "TissueCell.hpp"
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

class TestMutationSpread : public CxxTest::TestSuite
{
    
public:

    /*
     * The idea of this set of simulations is to label one of the cells at random
     * (as 'LABELLED' or later as a mutant). We then track the progress of its 
     * progeny throughout the crypt. When the labelled cells are swept away, 
     * or take over a crypt we end the simulation and log the result and time. 
     */
    void TestWhetherMutationsSpread() throw (Exception)
    {        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);

        /*
         * We load the steady state that the profiled test uses
         */
        std::string test_to_profile = "NiceCryptSim";
        double load_time = 350;   // this is the folder and time that the stored results were archived (needed to know foldernames)
        double time_of_each_run = 10; // run for 10 hours.
        double end_of_simulation = 1000;
        
        
        // Call a function to label a cell
        CryptSimulation2d* p_simulator = CryptSimulation2d::Load(test_to_profile,load_time);
        unsigned label_this = Label();
        p_simulator->rGetTissue().rGetCellAtNodeIndex(label_this).SetMutationState(LABELLED);
        p_simulator->Save();
        
        // write out to file which cell it was
        OutputFileHandler results_handler("NiceCryptSim",false);
        out_stream file=results_handler.OpenOutputFile("overall_results.dat");
        std::vector<double> position = p_simulator->GetNodeLocation(label_this);
        (*file) << "Node = " << label_this << " at x = " << position[0] << "\ty = " << position[1] << "\n" << std::flush;
  
        delete p_simulator;
        
        // run the simulator
        for (double t=load_time; t<end_of_simulation+0.5; t += time_of_each_run)
        {
            CryptSimulation2d* p_simulator = CryptSimulation2d::Load(test_to_profile,t);
            p_simulator->SetEndTime(t+time_of_each_run);
            p_simulator->Solve();
            p_simulator->Save();
            
            // Check the spread of the mutation (whether to stop?)
            c_vector<unsigned,5> cell_type_count = p_simulator->GetCellTypeCount();
            bool exit_now = false;
            
            (*file) << "Time = " << SimulationTime::Instance()->GetDimensionalisedTime() <<
                    "\t Labelled cells = " << cell_type_count[1] <<"\n"<<std::flush;
            
            // if labelled population has been swept away then stop.
            if (cell_type_count[1] == 0u
              ||cell_type_count[1] == p_simulator->rGetTissue().GetNumRealCells())
            {
                (*file) << "Finished simulation\n"<<std::flush;
                exit_now = true;   
            }
            
            delete p_simulator;
            if (exit_now) break;
        }
        (*file) << "Closing file\n"<<std::flush;        
        file->close();   
                
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

    /*
     * The cell labelling list. Here we make a list of all the cells saved 
     * in the simulation at t=350. We pick one at random to label and follow its
     * cell line
     */
    unsigned Label()
    {
        std::vector<unsigned> label_these;
        label_these.push_back(129);
        label_these.push_back(195);
        //label_these.push_back(82);
        label_these.push_back(106);
        label_these.push_back(98);
        label_these.push_back(77);
        label_these.push_back(89);
        label_these.push_back(227);
        label_these.push_back(155);
        label_these.push_back(72);
        //label_these.push_back(79);
        label_these.push_back(90);
        label_these.push_back(102);
        label_these.push_back(85);
        label_these.push_back(95);
        label_these.push_back(275);
        label_these.push_back(93);
        label_these.push_back(78);
        label_these.push_back(83);
        label_these.push_back(70);
        //label_these.push_back(206);
        label_these.push_back(40);
        label_these.push_back(41);
        label_these.push_back(92);
        label_these.push_back(71);
        label_these.push_back(289);
        label_these.push_back(86);
        label_these.push_back(84);
        //label_these.push_back(149);
        label_these.push_back(87);
        label_these.push_back(241);
        label_these.push_back(76);
        //label_these.push_back(81);
        label_these.push_back(99);
        //label_these.push_back(75);
        label_these.push_back(483);
        label_these.push_back(487);
        label_these.push_back(405);
        label_these.push_back(458);
        label_these.push_back(494);
        label_these.push_back(497);
        label_these.push_back(440);
        //label_these.push_back(448);
        label_these.push_back(505);
        label_these.push_back(504);
        label_these.push_back(465);
        //label_these.push_back(453);
        //label_these.push_back(451);
        label_these.push_back(454);
        //label_these.push_back(429);
        label_these.push_back(363);
        label_these.push_back(469);
        label_these.push_back(503);
        label_these.push_back(306);
        label_these.push_back(456);
        label_these.push_back(452);
        label_these.push_back(273);
        label_these.push_back(438);
        label_these.push_back(364);
        label_these.push_back(359);
        label_these.push_back(407);
        label_these.push_back(272);
        label_these.push_back(233);
        label_these.push_back(175);
        label_these.push_back(415);
        label_these.push_back(493);
        label_these.push_back(431);
        label_these.push_back(334);
        label_these.push_back(500);
        label_these.push_back(361);
        label_these.push_back(376);
        //label_these.push_back(263);
        label_these.push_back(378);
        label_these.push_back(468);
        
        unsigned this_index = RandomNumberGenerator::Instance()->randMod(label_these.size());
        
        return label_these[this_index];   
    }


};

#endif /*TESTMUTATIONSPREAD_HPP_*/
