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
            
            double prob_of_death = 0.2*mTimeStep - 0.01*mTimeStep*dist_to_centre;
            
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
    void Test2dSpheroid() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
       
        int num_cells_depth = 10;
        int num_cells_width = 10;
        double crypt_length = num_cells_depth-1.0;
        double crypt_width = num_cells_width-1.0;
        
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

        simulator.SetOutputDirectory("2dSpheroid");
        simulator.SetEndTime(12.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
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
};

#endif /*TESTTUMOURSPHEROID_HPP_*/
