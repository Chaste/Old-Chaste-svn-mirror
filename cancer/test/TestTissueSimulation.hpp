#ifndef TESTTISSUESIMULATION_HPP_
#define TESTTISSUESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <stdio.h>
#include <time.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <vector>
#include <iostream>

#include "TissueSimulation.cpp"
#include "FixedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "AbstractCancerTestSuite.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CryptProjectionSpringSystem.hpp"
#include "RandomCellKiller.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "SimpleWntCellCycleModel.hpp"

/** 
 *  Note: Most tests of TissueSimulation are in TestCryptSimulation2d
 */

class TestTissueSimulation : public AbstractCancerTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCancerTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCancerTestSuite::tearDown();
    }
    
public:

    void TestOutputStatistics() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools
        RandomNumberGenerator *p_gen=RandomNumberGenerator::Instance();
        CancerParameters::Instance()->SetHepaOneParameters();
        
        // Set up mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0u, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
                    
        // Set up cells
        std::vector<TissueCell> cells;        
        
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(STEM, HEALTHY, new TysonNovakCellCycleModel());
            double birth_time = -1.0*p_gen->ranf();
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // Set up tissue        
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        tissue.SetWriteTissueAreas(true); // record the spheroid radius and necrotic radius
              
        Meineke2001SpringSystem<2> spring_system(tissue);
        spring_system.UseCutoffPoint(1.5);
                  
        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, &spring_system);
        simulator.SetOutputDirectory("TissueSimulationWritingProteins");
        simulator.SetEndTime(0.5);
        simulator.SetOutputCellVariables(true);
        simulator.SetOutputCellCyclePhases(true);

        // Run tissue simulation 
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        OutputFileHandler handler("TissueSimulationWritingProteins",false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellvariables.dat";
        TS_ASSERT_EQUALS(system(("diff " + results_file + " cancer/test/data/TissueSimulationWritingProteins/cellvariables.dat").c_str()), 0);
    }
    
    /**
     *  Test a tissue simulation with a non-Meineke spring system. 
     * 
     *  This test consists of a standard crypt projection model simulation with a 
     *  radial sloughing cell killer, a crypt projection cell cycle model that 
     *  depends on a radial Wnt gradient, and the crypt projection model spring 
     *  system, and store the results for use in later archiving tests.
     */
    void TestTissueSimulationWithCryptProjectionSpringSystem() throw (Exception)
    {        
        CancerParameters *p_params = CancerParameters::Instance();        
        p_params->SetWntStemThreshold(0.95);
        
        double a = 0.2;
        double b = 2.0;      
        p_params->SetCryptProjectionParameterA(a);
        p_params->SetCryptProjectionParameterB(b);   
                
        // Set up mesh
        int num_cells_depth = 20;
        int num_cells_width = 20;     
        unsigned thickness_of_ghost_layer = 3;   
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);
        
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);
        
        p_mesh->Translate(-width_of_mesh/2,-height_of_mesh/2);                   
        
        // To start off with, set up all cells to be of type TRANSIT
        std::vector<TissueCell> cells;
        
        std::cout << "num nodes = " << p_mesh->GetNumNodes() << "\n" << std::flush;
        
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, new SimpleWntCellCycleModel());
            cell.InitialiseCellCycleModel();
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                  ( p_params->GetTransitCellG1Duration()
                                   +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
                                                          
        // Make a tissue
        MeshBasedTissueWithGhostNodes<2> crypt(*p_mesh, cells, ghost_node_indices);         

        // Set up the Wnt gradient 
        WntConcentration::Instance()->SetType(RADIAL); 
        WntConcentration::Instance()->SetTissue(crypt);   
        
        // Create the spring system
        CryptProjectionSpringSystem spring_system(crypt);
                
        // Make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt, &spring_system, false, false);
        
        // Create a radial cell killer and pass it in to the tissue simulation
        c_vector<double,2> centre = zero_vector<double>(2);
        double crypt_radius = pow(CancerParameters::Instance()->GetCryptLength()/a, 1.0/b);

        RadialSloughingCellKiller killer(&crypt, centre, crypt_radius);
        crypt_projection_simulator.AddCellKiller(&killer);
        
        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionSimulation");
        crypt_projection_simulator.SetEndTime(0.25);
        
        // Run the simulation
        TS_ASSERT_THROWS_NOTHING(crypt_projection_simulator.Solve());
        
        // These cells just divided and have been gradually moving apart.
        // These results are from time 0.25.
        std::vector<double> node_302_location = crypt_projection_simulator.GetNodeLocation(302);
        TS_ASSERT_DELTA(node_302_location[0], -0.0954, 1e-4);
        TS_ASSERT_DELTA(node_302_location[1], 0.2475, 1e-4);
        
        std::vector<double> node_506_location = crypt_projection_simulator.GetNodeLocation(506);
        TS_ASSERT_DELTA(node_506_location[0], -0.7046, 1e-4);
        TS_ASSERT_DELTA(node_506_location[1], 0.6178, 1e-4);

        // Test the Wnt gradient result
        TissueCell* p_cell = &(crypt.rGetCellAtNodeIndex(302));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 0.9991, 1e-4);
        p_cell = &(crypt.rGetCellAtNodeIndex(506));
        TS_ASSERT_DELTA(WntConcentration::Instance()->GetWntLevel(p_cell), 0.9898, 1e-4);
        
        // Tidy up
        WntConcentration::Destroy();
    }
    
    /**
     *  Test a tissue simulation with a cell killer.
     * 
     *  In this test, we solve a tissue simulation without ghost nodes and 
     *  check that the numbers of nodes and cells match at the end of the 
     *  simulation. 
     * 
     *  This test currently fails so is commented out. 
     */
    void xTestTissueSimulationWithCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0, false);
        ConformingTetrahedralMesh<2,2>* p_mesh = generator.GetMesh();
        
        // Set up cells, one for each node. Give each cell a random birth time.
        std::vector<TissueCell> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (CancerParameters::Instance()->GetStemCellG1Duration()
                                    + CancerParameters::Instance()->GetSG2MDuration() );        
            TissueCell cell(STEM, HEALTHY, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        // Create a tissue
        MeshBasedTissue<2> tissue(*p_mesh, cells);
        
        // Create a mechanics system
        Meineke2001SpringSystem<2> spring_system(tissue);
        spring_system.UseCutoffPoint(1.5);
                  
        // Set up tissue simulation
        TissueSimulation<2> simulator(tissue, &spring_system);
        simulator.SetOutputDirectory("TestTissueSimulationWithCellDeath");
        simulator.SetEndTime(0.5);
        
        // Add cell killer
        RandomCellKiller<2> random_cell_killer(&tissue, 0.05);
        simulator.AddCellKiller(&random_cell_killer);

        // Run tissue simulation 
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
        
        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetTissue().GetNumNodes(), simulator.rGetTissue().GetNumRealCells());
    }
    
};
#endif /*TESTTISSUESIMULATION_HPP_*/
