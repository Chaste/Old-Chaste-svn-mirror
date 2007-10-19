#ifndef TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <cxxtest/TestSuite.h>
#include "TissueSimulationWithNutrients.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <vector>
#include "FixedCellCycleModel.hpp"
#include "ColumnDataReader.hpp"
#include "SimulationTime.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "PetscTools.hpp"

class SimpleOxygenBasedCellCycleModel : public FixedCellCycleModel
{
private:
    double mTimeProgressingThroughCellCycle;
        
public:
    SimpleOxygenBasedCellCycleModel() 
        : FixedCellCycleModel(),
          mTimeProgressingThroughCellCycle(0.0)
    {
    }
    
    bool ReadyToDivide()
    {
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell);
        
        // the std::max is a hack, due to the choice of PDE
        // (a simple nonlinear PDE would ensure positivity)   
//        if (oxygen_concentration < 0.0)
//        {
//            EXCEPTION("Oxygen concentration has gone negative - check the oxygen PDE");
//        }
        mTimeProgressingThroughCellCycle = mTimeProgressingThroughCellCycle + std::max(oxygen_concentration,0.0)*SimulationTime::Instance()->GetTimeStep(); 
        
        bool result = false;        
        if ( mTimeProgressingThroughCellCycle > CancerParameters::Instance()->GetHepaOneCellG1Duration() + 
                                                        CancerParameters::Instance()->GetSG2MDuration() )
        {
            result = true;
        }
        
        return result;
    }
    
    AbstractCellCycleModel* CreateCellCycleModel()
    {
        return new SimpleOxygenBasedCellCycleModel();
    }
};

class SimpleLinearEllipticPde : public AbstractLinearEllipticPde<2>
{

public:
    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return -0.1; //-1.0;
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

class TestTissueSimulationWithNutrients : public CxxTest::TestSuite
{
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
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
                
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells
        std::vector<TissueCell> cells;
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(HEPA_ONE, HEALTHY, 0, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -p_gen->ranf()*(p_params->GetTransitCellG1Duration()
                                               +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
                
        Tissue<2> tissue(*p_mesh, cells);
        tissue.SetGhostNodes(ghost_node_indices);
        
        SimpleLinearEllipticPde pde;

        TissueSimulationWithNutrients<2> simulator(tissue, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.2);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        //simulator.UseCutoffPoint(1.5);
        
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
        
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);
        
        double start_time = std::clock();
        
        simulator.Solve();
        
        double end_time = std::clock();
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "Time to perform simulation was = " << elapsed_time << "\n";
                
        // add get methods etc and test
        
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
