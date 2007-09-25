#ifndef TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulationWithNutrients.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "WntGradient.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "ColumnDataReader.hpp"
#include "SimulationTime.hpp"
#include "SloughingCellKiller.hpp"
#include "PetscTools.hpp"

// Possible types of Cell Cycle Model (just for CreateVectorOfCells method)
typedef enum CellCycleType_
{
    FIXED,
    STOCHASTIC,
    WNT,
    TYSONNOVAK
} CellCycleType;


class SimpleEllipticPde : public AbstractLinearEllipticPde<2>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return -1;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> , double )
    {
        return 0.0;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(ChastePoint<2> )
    {
        return identity_matrix<double>(2);
    }
        
};


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

class TestTissueSimulationWithNutrients : public CxxTest::TestSuite
{
    void CreateVectorOfCells(std::vector<MeinekeCryptCell>& rCells, 
                             ConformingTetrahedralMesh<2,2>& rMesh, 
                             CellCycleType cycleType, 
                             bool randomBirthTimes,
                             double y0 = 0.3,
                             double y1 = 2.0,
                             double y2 = 3.0,
                             double y3 = 4.0)
    {
        RandomNumberGenerator *p_random_num_gen=RandomNumberGenerator::Instance();
        unsigned num_cells = rMesh.GetNumNodes();

        AbstractCellCycleModel* p_cell_cycle_model = NULL;
        double typical_transit_cycle_time;
        double typical_stem_cycle_time;
        
        CancerParameters* p_params = CancerParameters::Instance();
        
        rCells.reserve(num_cells);
        
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;

            double y = rMesh.GetNode(i)->GetPoint().rGetLocation()[1];
            
            if (cycleType==FIXED)
            {
                p_cell_cycle_model = new FixedCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==STOCHASTIC)
            {
                p_cell_cycle_model = new StochasticCellCycleModel();
                typical_transit_cycle_time = p_params->GetTransitCellCycleTime();
                typical_stem_cycle_time = p_params->GetStemCellCycleTime();
            }
            else if (cycleType==WNT)
            {
                WntGradient wnt_gradient(LINEAR);
                double wnt = wnt_gradient.GetWntLevel(y);
                p_cell_cycle_model = new WntCellCycleModel(wnt);
                typical_transit_cycle_time = 16.0;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else if (cycleType==TYSONNOVAK)
            {
                p_cell_cycle_model = new TysonNovakCellCycleModel();
                typical_transit_cycle_time = 1.25;
                typical_stem_cycle_time = typical_transit_cycle_time;
            }
            else
            {
                EXCEPTION("Cell Cycle Type is not recognised");   
            }
            
            
            double birth_time = 0.0;
            
            if (y <= y0)
            {
                cell_type = STEM;
                generation = 0;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_stem_cycle_time; // hours
                }
            }
            else if (y < y1)
            {
                cell_type = TRANSIT;
                generation = 1;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y2)
            {
                cell_type = TRANSIT;
                generation = 2;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else if (y < y3)
            {
                cell_type = TRANSIT;
                generation = 3;
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
            }
            else
            {
                if(randomBirthTimes)
                {
                    birth_time = -p_random_num_gen->ranf()*typical_transit_cycle_time; // hours 
                }
                if(cycleType==WNT || cycleType==TYSONNOVAK)
                {
                    // There are no fully differentiated cells!
                    cell_type = TRANSIT;
                    
                }
                else
                {
                    cell_type = DIFFERENTIATED;
                }                
                generation = 4;
            }

            MeinekeCryptCell cell(cell_type, HEALTHY, generation, p_cell_cycle_model);
            
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            rCells.push_back(cell);
        }
    }
    
public:

    void TestWithOxygen() throw(Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }
        
        
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
        
        SimpleEllipticPde pde;

        TissueSimulationWithNutrients<2> simulator(crypt, &pde);

        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.5);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        //simulator.UseCutoffPoint(1.5);
        
        c_vector<double,2> centre(2);
        centre(0) = (double)num_cells_width/2.0;
        centre(1) = (double)num_cells_depth/2.0;
        
        AbstractCellKiller<2>* p_killer = new RadiusBasedCellKiller(&crypt, centre, simulator.GetDt());
        simulator.AddCellKiller(p_killer);
        
        simulator.Solve();
        
        // add get methods etc and test
        
        
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
