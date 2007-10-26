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
#include "SimpleOxygenBasedCellCycleModel.hpp"

class OxygenPde : public AbstractNonlinearEllipticPde<2>
{
public:

    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> , double u)
    {
        return -u;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(ChastePoint<2> , double u)
    {
        return identity_matrix<double>(2);
    }
    
    c_matrix<double,2,2> ComputeDiffusionTermPrime(ChastePoint<2> , double u)
    {
        return zero_matrix<double>(2);
    }
    
    double ComputeNonlinearSourceTermPrime(ChastePoint<2> , double u)
    {
        return -1.0;
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
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
                
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
        
        OxygenPde pde;

        TissueSimulationWithNutrients<2> simulator(tissue, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.1);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.UseCutoffPoint(1.5);
        
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
