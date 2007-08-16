#ifndef TESTTUMOURSPHEROID_HPP_
#define TESTTUMOURSPHEROID_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulation.cpp"

#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "RandomCellKiller.hpp"



class RadiusBasedCellKiller : public AbstractCellKiller<2>
{
private :
    c_vector<double,2> mCentre;
    double mTimeStep;

public :
    RadiusBasedCellKiller(Crypt<2>* pCrypt, c_vector<double,2> centre, double timeStep)
        : AbstractCellKiller<2>(pCrypt),
          mCentre(centre),
          mTimeStep(timeStep)
    {
    }
    
    virtual void TestAndLabelCellsForApoptosisOrDeath()
    {
        for(Crypt<2>::Iterator cell_iter = mpCrypt->Begin();
            cell_iter != mpCrypt->End();
            ++cell_iter)
        {
            const c_vector<double,2>& location = cell_iter.GetNode()->rGetLocation();
            double dist_to_centre = norm_2(location - mCentre);
            
            
            double prob_of_death = 2*mTimeStep - 1*mTimeStep*dist_to_centre;
            if (prob_of_death<=0.0)
            {
                prob_of_death=0.0;
            }
            
            if (!cell_iter->HasApoptosisBegun() &&
                RandomNumberGenerator::Instance()->ranf() < prob_of_death)
            {
                cell_iter->StartApoptosis();
            }    
        }
    }
};


class TestTumourSpheroid : public CxxTest::TestSuite
{  
public :
    void dontTest2dSpheroid() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
       
        int num_cells_depth = 10;
        int num_cells_width = 10;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(TRANSIT, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = -p_gen->ranf()*p_params->GetTransitCellCycleTime();
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
                
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("2dSpheroid");
        simulator.SetEndTime(5.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        //simulator.UseCutoffPoint(1.5);
        
        c_vector<double,2> centre(2);
        centre(0) = (double)num_cells_width/2.0;
        centre(1) = (double)num_cells_depth/2.0;
        
        AbstractCellKiller<2>* p_killer = new RadiusBasedCellKiller(&crypt, centre, simulator.GetDt());
        simulator.AddCellKiller(p_killer);
        
        simulator.Solve();
        
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
    void TestApoptosisSpringLengths() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
       
        int num_cells_depth = 1;
        int num_cells_width = 2;
        double crypt_length = num_cells_depth-0.0;
        double crypt_width = num_cells_width-0.0;
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<MeinekeCryptCell> cells;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            MeinekeCryptCell cell(TRANSIT, HEALTHY, 0, new FixedCellCycleModel());
            double birth_time = -p_gen->ranf()*p_params->GetTransitCellCycleTime();
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
                
        Crypt<2> crypt(*p_mesh, cells);
        crypt.SetGhostNodes(ghost_node_indices);

        TissueSimulation<2> simulator(crypt);

        simulator.SetOutputDirectory("2dSpheroidApoptosis");
        simulator.SetEndTime(1.0);
        simulator.SetMaxCells(50);
        simulator.SetMaxElements(100);
        
        CancerParameters::Instance()->SetApoptosisTime(2.0);
        crypt.rGetCellAtNodeIndex(14).StartApoptosis();
        crypt.rGetCellAtNodeIndex(15).StartApoptosis();
        simulator.SetNoBirth(true);
                        
        simulator.Solve();
        
        /* We track the locations of two dying cells (a and b) and two
         * live cells adjacent to them (c and d)
         * 
         * All cells begin distance 1 apart.
         * 
         * a and b move together to leave a gap of 0.
         * a and c (and b and d) move to a distance of 0.5 apart.
         */
        
        c_vector<double, 2> a_location = crypt.rGetMesh().GetNode(14)->rGetLocation();
        c_vector<double, 2> b_location = crypt.rGetMesh().GetNode(15)->rGetLocation();
        c_vector<double, 2> c_location = crypt.rGetMesh().GetNode(20)->rGetLocation();
        c_vector<double, 2> d_location = crypt.rGetMesh().GetNode(21)->rGetLocation();
        
        double a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                                (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        double a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                                (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        double c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                                (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));
        
        TS_ASSERT_DELTA(a_b_separation , 0.5, 1e-1);
        TS_ASSERT_DELTA(a_c_separation , 0.75, 1e-1);
        TS_ASSERT_DELTA(c_d_separation , 1.0, 1e-1);
        
        simulator.SetEndTime(1.99);
        simulator.Solve();
        
        a_location = crypt.rGetMesh().GetNode(14)->rGetLocation();
        b_location = crypt.rGetMesh().GetNode(15)->rGetLocation();
        c_location = crypt.rGetMesh().GetNode(20)->rGetLocation();
        d_location = crypt.rGetMesh().GetNode(21)->rGetLocation();
        
        a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                         (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                         (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                         (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));
        
        TS_ASSERT_DELTA(a_b_separation , 0.01, 1e-1);
        TS_ASSERT_DELTA(a_c_separation , 0.5, 1e-1);
        //TS_ASSERT_DELTA(c_d_separation , 1.0, 1e-1); no longer attached by a spring
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 4u);
        
        simulator.SetEndTime(2.01);
        simulator.Solve();
        
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 2u);
        
        
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
   
};

#endif /*TESTTUMOURSPHEROID_HPP_*/
