#ifndef TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <stdio.h>
#include <time.h>
#include <cxxtest/TestSuite.h>
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include <iostream>
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

class SimplePdeForTesting : public AbstractNonlinearEllipticPde<2>
{
public:

    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return -0.1;
    }
    
    double ComputeNonlinearSourceTerm(ChastePoint<2> , double u)
    {
        return 0.0;
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
        return 0.0;
    }
};

class SimpleOxygenPde : public AbstractNonlinearEllipticPde<2>
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

    /* 
     * Test the PostSolve method using the problem del squared C = 0.1
     * on the unit disc, with boundary condition C=1 on r=1, which has 
     * analytic solution C = 1-0.025*(1-r^2).
     */
    void TestPostSolveMethod() throw(Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }
        
        // instantiate singleton objects
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
       
       
        // set up mesh
        ConformingTetrahedralMesh<2,2>* p_mesh = new ConformingTetrahedralMesh<2,2>;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        p_mesh->ConstructFromMeshReader(mesh_reader);
            
        // set up cells
        std::vector<TissueCell> cells;        
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(HEPA_ONE, HEALTHY, 0, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -p_gen->ranf()*(p_params->GetHepaOneCellG1Duration()
                                               +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // set up tissue        
        Tissue<2> tissue(*p_mesh, cells);
        
        // set up cellwisedata and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);
        
        // set up PDE
        SimplePdeForTesting pde;
        
        // set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, &pde);
        simulator.SetOutputDirectory("TestPostSolve");
        simulator.SetEndTime(0.1);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.rGetMeinekeSystem().UseCutoffPoint(1.5);
        
        // set up cell killer and pass into simulation
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
        
        // solve PDE             
        simulator.PostSolve();        
        
        // check the correct solution was obtained           
        for (Tissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {   
            double radius = norm_2(tissue.GetLocationOfCell(*cell_iter));
            double analytic_solution = 1 - 0.025*(1 - pow(radius,2.0));
            
            TS_ASSERT_DELTA(p_data->GetValue(&(*cell_iter)), analytic_solution, 1e-2);
        }     
        
        // tidy up
        delete p_mesh;
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }

    void TestWithOxygen() throw(Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }
        
        // instantiate singleton objects
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
              
        // set up mesh
        int num_cells_depth = 10;
        int num_cells_width = 10;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
                    
        // set up cells
        std::vector<TissueCell> cells;        
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(HEPA_ONE, HEALTHY, 0, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -p_gen->ranf()*(p_params->GetHepaOneCellG1Duration()
                                               +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // set up tissue        
        Tissue<2> tissue(*p_mesh, cells);
        
        // set up cellwisedata and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);
        
        // set up PDE
        SimpleOxygenPde pde;
        
        // set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, &pde);
        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.01);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        simulator.rGetMeinekeSystem().UseCutoffPoint(1.5);
        
        // set up cell killer and pass into simulation
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
               
        // run tissue simulation and print out time taken    
        double start_time = std::clock();
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        double end_time = std::clock();
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "Time taken to perform simulation was = " << elapsed_time << "\n";
                        
        // tidy up
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
