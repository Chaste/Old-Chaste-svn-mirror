#ifndef TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_
#define TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_

#include <stdio.h>
#include <time.h>
#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
    double ComputeLinearSourceTerm(const ChastePoint<2>& )
    {
        return -1.0;
    }

    double ComputeNonlinearSourceTerm(const ChastePoint<2>& , double )
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& , double u)
    {
        return identity_matrix<double>(2);
    }
    
    c_matrix<double,2,2> ComputeDiffusionTermPrime(const ChastePoint<2>& , double u)
    {
        return zero_matrix<double>(2);
    }
    
    double ComputeNonlinearSourceTermPrime(const ChastePoint<2>& , double u)
    {
        return 0.0;
    }
};

class SimpleOxygenPde : public AbstractNonlinearEllipticPde<2>
{
public:

    double ComputeLinearSourceTerm(const ChastePoint<2>& )
    {
        return 0.0;
    }
    
    double ComputeNonlinearSourceTerm(const ChastePoint<2>& , double u)
    {
        return -0.1*u;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& , double u)
    {
        return identity_matrix<double>(2);
    }
    
    c_matrix<double,2,2> ComputeDiffusionTermPrime(const ChastePoint<2>& , double u)
    {
        return zero_matrix<double>(2);
    }
    
    double ComputeNonlinearSourceTermPrime(const ChastePoint<2>& , double u)
    {
        return -1.0;
    }
};

class TestTissueSimulationWithNutrients : public CxxTest::TestSuite
{
    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);
        CancerParameters::Instance()->Reset();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        // Clear up singleton classes
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
    
public:

    /* 
     * A two-part test for the PostSolve() method.
     * 
     * Firstly, test the PDE solver using the problem del squared C = 0.1 
     * on the unit disc, with boundary condition C=1 on r=1, which has 
     * analytic solution C = 1-0.25*(1-r^2).
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
        
        // change the hypoxic concentration, just for this test
        CancerParameters::Instance()->SetHepaOneCellHypoxicConcentration(0.9);
        CancerParameters::Instance()->SetHepaOneParameters();

        // set up mesh
        ConformingTetrahedralMesh<2,2>* p_mesh = new ConformingTetrahedralMesh<2,2>;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        p_mesh->ConstructFromMeshReader(mesh_reader);
            
        // set up cells
        std::vector<TissueCell> cells;        
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                    (CancerParameters::Instance()->GetHepaOneCellG1Duration()
                                    +CancerParameters::Instance()->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cell.SetSymmetricDivision();
            cells.push_back(cell);
        }
        
        // set up tissue        
        Tissue<2> tissue(*p_mesh, cells);
        
        // set up cellwisedata and associate it with the tissue
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumNodesAndVars(p_mesh->GetNumNodes(), 1);
        p_data->SetTissue(tissue);
        
        // since values are first passed in to CellwiseData before it is updated in PostSolve(),
        // we need to pass it some initial conditions to avoid memory errors
        // (note: it would really make more sense to put the PDE solver stuff in a PreSolve method)  
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i));
        }
        
        // set up PDE
        SimplePdeForTesting pde;
        
        Meineke2001SpringSystem<2>* p_spring_system = new Meineke2001SpringSystem<2>(tissue);
        
        // use an extremely small cutoff so that no cells interact 
        // - this is to ensure that in the Solve method, the cells don't move
        // (we need to call Solve to set up the .viznutrient file)
        p_spring_system->UseCutoffPoint(0.0001); 
              
        // set up tissue simulation
        TissueSimulationWithNutrients<2> simulator(tissue, p_spring_system, &pde); 
        simulator.SetOutputDirectory("TestPostSolveMethod");
        simulator.SetEndTime(2.0/120.0);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
        // set up cell killer and pass into simulation
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
        
        simulator.Solve();                
        
        // check the correct solution was obtained           
        for (Tissue<2>::Iterator cell_iter = tissue.Begin();
             cell_iter != tissue.End();
             ++cell_iter)
        {   
            double radius = norm_2(tissue.GetLocationOfCell(*cell_iter));
            double analytic_solution = 1 - 0.25*(1 - pow(radius,2.0));
            
            // Get cell model
            AbstractCellCycleModel* p_abstract_model=cell_iter->GetCellCycleModel();
            SimpleOxygenBasedCellCycleModel* p_oxygen_model = static_cast <SimpleOxygenBasedCellCycleModel*>(p_abstract_model);
            
            // First part of test - check that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(&(*cell_iter)), analytic_solution, 1e-2);
            
            // Second part of test - check that each cell's hypoxic duration is correctly updated
            if ( p_data->GetValue(&(*cell_iter)) >= CancerParameters::Instance()->GetHepaOneCellHypoxicConcentration() )
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 0.0, 1e-5);
            }
            else
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 1/120.0, 1e-5);
            } 
        }     
        
        // tidy up
        delete p_mesh;
        delete p_killer;
        CellwiseData<2>::Destroy();
    }
        
    void TestWithOxygen() throw(Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }
        
        CancerParameters::Instance()->SetHepaOneParameters();
        
        // set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
                    
        // set up cells
        std::vector<TissueCell> cells;        
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*
                                    (CancerParameters::Instance()->GetHepaOneCellG1Duration()
                                    +CancerParameters::Instance()->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cell.SetSymmetricDivision();
            cells.push_back(cell);
        }
        
        // set up tissue        
        Tissue<2> tissue(*p_mesh, cells);
        
        // set up CellwiseData and associate it with the tissue
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
               
        
        // run tissue simulation 
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        
        // record final mesh size for visualizer
        TS_ASSERT_THROWS_NOTHING(simulator.WriteFinalMeshSizeForVisualizer());
                        
        // test positions        
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.4968, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.8635, 1e-4);
        
        std::vector<double> node_15_location = simulator.GetNodeLocation(15);
        TS_ASSERT_DELTA(node_15_location[0], 0.4976, 1e-4);
        TS_ASSERT_DELTA(node_15_location[1], 2.5977, 1e-4);
                
        // test the CellwiseData result
        TissueCell* p_cell = &(simulator.rGetTissue().rGetCellAtNodeIndex(5));
        TS_ASSERT_DELTA(CellwiseData<2>::Instance()->GetValue(p_cell), 0.9604, 1e-4);
        
        p_cell = &(simulator.rGetTissue().rGetCellAtNodeIndex(15));
        TS_ASSERT_DELTA(CellwiseData<2>::Instance()->GetValue(p_cell), 0.9584, 1e-4);
        
                        
        // tidy up
        delete p_killer;
        CellwiseData<2>::Destroy();
    }
    
    /*
     * This test compares the visualizer output from the previous test with a known file.
     * 
     * Note - if the previous test is changed we need to update the file this test refers to. 
     */
    void TestWriteNutrient() throw (Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }

        // work out where the previous test wrote its files
        OutputFileHandler handler("TissueSimulationWithOxygen",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/vis_results/results.viznutrient";         
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TissueSimulationWithOxygen_vis/results.viznutrient").c_str()), 0);     
    }
    
    void TestArchiving() throw (Exception)
    {
        if (!PetscTools::IsSequential())
        {
            TS_TRACE("This test does not pass in parallel yet.");
            return;
        }
        
        CancerParameters::Instance()->SetHepaOneParameters();
        
        // set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
                    
        // set up cells
        std::vector<TissueCell> cells;        
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel());
            double birth_time = -1.0 - ( (double) i/p_mesh->GetNumNodes() )*
                                            (CancerParameters::Instance()->GetHepaOneCellG1Duration()
                                             +CancerParameters::Instance()->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cell.SetSymmetricDivision();
            cells.push_back(cell);
        }
        
        // set up tissue        
        Tissue<2> tissue(*p_mesh, cells);
        
        // set up CellwiseData and associate it with the tissue
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
        simulator.SetOutputDirectory("TissueSimulationWithNutrientsSaveAndLoad");
        simulator.SetEndTime(0.2);
        simulator.SetMaxCells(400);
        simulator.SetMaxElements(800);
        
        // set up cell killer and pass into simulation
        AbstractCellKiller<2>* p_killer = new OxygenBasedCellKiller<2>(&tissue);
        simulator.AddCellKiller(p_killer);
        
        simulator.Solve();
        
        simulator.Save();
        
        TissueSimulationWithNutrients<2>* p_simulator = TissueSimulationWithNutrients<2>::Load("TissueSimulationWithNutrientsSaveAndLoad", 0.2);
        p_simulator->SetPde(&pde);
        
        p_simulator->SetEndTime(0.5);
        p_simulator->Solve();        
        
        // These results are from time 0.5 in TestWithOxygen.
        std::vector<double> node_5_location = p_simulator->GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.4968, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.8635, 1e-4);
        
        std::vector<double> node_15_location = p_simulator->GetNodeLocation(15);
        TS_ASSERT_DELTA(node_15_location[0], 0.4976, 1e-4);
        TS_ASSERT_DELTA(node_15_location[1], 2.5977, 1e-4);
        
        // test CellwiseData was set up correctly
        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(),true);
        
        // test the CellwiseData result
        TissueCell* p_cell = &(p_simulator->rGetTissue().rGetCellAtNodeIndex(5));
        TS_ASSERT_DELTA(CellwiseData<2>::Instance()->GetValue(p_cell), 0.9604, 1e-4);
        
        p_cell = &(p_simulator->rGetTissue().rGetCellAtNodeIndex(15));
        TS_ASSERT_DELTA(CellwiseData<2>::Instance()->GetValue(p_cell), 0.9584, 1e-4);
        
        delete p_killer;
        delete p_simulator;       
        CellwiseData<2>::Destroy();
        
    }

};
#endif /*TESTTISSUESIMULATIONWITHNUTRIENTS_HPP_*/
