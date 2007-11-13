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
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "Meineke2001SpringSystem.hpp" 

class SimplePdeForTesting : public AbstractNonlinearEllipticPde<2>
{
public:
    double ComputeLinearSourceTerm(ChastePoint<2> )
    {
        return -1.0;
    }

    double ComputeNonlinearSourceTerm(ChastePoint<2> , double )
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
        return -0.1*u;
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
     * A two-part test for the PostSolve() method.
     * 
     * Firstly, test the PDE solver using the problem del squared C = 0.1 
     * on the unit disc, with boundary condition C=1 on r=1, which has 
     * analytic solution C = 1-0.125*(1-r^2).
     * 
     * Secondly, test that cells' hypoxic durations are correctly updated when a 
     * nutrient distribution is prescribed.
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
        
        // change the hypoxic concentration, just for this test
        p_params->SetHepaOneCellHypoxicConcentration(0.9);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        // since we are not calling Solve(), we need to set up 
        // the simulation time manually
	    SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0/120.0, 1);   
       
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
        
        Meineke2001SpringSystem<2>* p_spring_system = new Meineke2001SpringSystem<2>(tissue); 
        p_spring_system->UseCutoffPoint(1.5); 
              
        // set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, p_spring_system, &pde); 
        simulator.SetOutputDirectory("TestPostSolveMethod");
        simulator.SetEndTime(1.0/120.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
        // set up cell killer and pass into simulation
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
        
        simulator.PostSolve();                
        
        // check the correct solution was obtained           
        for (Tissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {   
            double radius = norm_2(tissue.GetLocationOfCell(*cell_iter));
            double analytic_solution = 1 - 0.25*(1 - pow(radius,2.0));
            
            // First part of test - check that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(&(*cell_iter)), analytic_solution, 1e-2);
            
            // Second part of test - check that each cell's hypoxic duration is correctly updated
            if ( p_data->GetValue(&(*cell_iter)) >= p_params->GetHepaOneCellHypoxicConcentration() )
            {
            	TS_ASSERT_DELTA(cell_iter->GetHypoxicDuration(), 0.0, 1e-5);
            }
            else
            {
            	TS_ASSERT_DELTA(cell_iter->GetHypoxicDuration(), 1.0/120.0, 1e-5);	
            } 
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
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
              
        // set up mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
                    
        // set up cells
        std::vector<TissueCell> cells;        
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(HEPA_ONE, HEALTHY, 0, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*(p_params->GetHepaOneCellG1Duration()
                                                                       +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // set up tissue        
        Tissue<2> tissue(*p_mesh, cells);
        
        // set up Cellwiseata and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        p_data->SetTissue(tissue);
        
        // since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        // (note: it would really make more sense to put the PDE solver stuff in a PreSolve method)  
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
			p_data->SetValue(1.0, p_mesh->GetNode(i));
		}
        
        // set up PDE
        SimpleOxygenPde pde;
        
        Meineke2001SpringSystem<2>* p_spring_system = new Meineke2001SpringSystem<2>(tissue); 
        p_spring_system->UseCutoffPoint(1.5); 
                  
        // set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, p_spring_system, &pde); 
        simulator.SetOutputDirectory("TissueSimulationWithOxygen");
        simulator.SetEndTime(0.5);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
        // set up cell killer and pass into simulation
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
               
        double start_time = std::clock();       
        
        // run tissue simulation 
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        
        // record final mesh size for visualizer
        TS_ASSERT_THROWS_NOTHING(simulator.WriteFinalMeshSizeForVisualizer());
                        
        double end_time = std::clock();
        
        // print out time taken to run tissue simulation    
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "Time taken to perform simulation was = " << elapsed_time << "\n";
                        
        // tidy up
        delete p_killer;
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellwiseData<2>::Destroy();
    }
    
    /*
     * This test compares the visualizer output from the previous test with a known file.
     * 
     * Note - if the previous test is changed we need to update the file this test refers to. 
     */
    void TestWriteFinalMeshSizeForVisualizerMethod() throw (Exception)
    {
        // work out where the previous test wrote its files
        OutputFileHandler handler("TissueSimulationWithOxygen",false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0/vis_results";
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.vizelements cancer/test/data/TissueSimulationWithOxygen_vis/previous_results.vizelements").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.viznodes cancer/test/data/TissueSimulationWithOxygen_vis/previous_results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/results.vizsetup cancer/test/data/TissueSimulationWithOxygen_vis/previous_results.vizsetup").c_str()), 0);
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
