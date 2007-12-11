#ifndef TESTCRYPTPROJECTIONSIMULATION_HPP_
#define TESTCRYPTPROJECTIONSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "TissueSimulation.cpp"
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>
#include <ctime>
#include <vector>
#include "OutputFileHandler.hpp"
#include "ColumnDataReader.hpp"
#include "RadialSloughingCellKiller.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CryptProjectionSpringSystem.hpp"
#include "CompareMeshes.hpp"
    
class TestCryptProjectionSimulation : public CxxTest::TestSuite
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
    
    /**
     *  Test a standard crypt projection model simulation with a radial sloughing 
     *  cell killer, a crypt projection cell cycle model that depends on a radial 
     *  Wnt gradient, and the crypt projection model spring system, and store the 
     *  results for use in later archiving tests.
     */
    void TestTissueSimulationWithCryptProjectionSpringSystem() throw (Exception)
    {        
        CancerParameters *p_params = CancerParameters::Instance();        
        p_params->SetRadialWntThreshold(0.95);
        
        double a = 0.2;
        double b = 2.0;      
        p_params->SetCryptProjectionParameterA(a);
        p_params->SetCryptProjectionParameterB(b);   
                
        // Set up mesh
        int num_cells_depth = 20;
        int num_cells_width = 20;     
        unsigned thickness_of_ghost_layer = 3;   
        
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, thickness_of_ghost_layer, false);
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        std::set<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        c_vector<double,2> width_extremes = p_mesh->GetWidthExtremes(0u);
        c_vector<double,2> height_extremes = p_mesh->GetWidthExtremes(1u);
        
        double width_of_mesh = (num_cells_width/(num_cells_width+2.0*thickness_of_ghost_layer))*(width_extremes[1] - width_extremes[0]);
        double height_of_mesh = (num_cells_depth/(num_cells_depth+2.0*thickness_of_ghost_layer))*(height_extremes[1] - height_extremes[0]);
        
        p_mesh->Translate(-width_of_mesh/2,-height_of_mesh/2);                   
        
        // To start off with, set up all cells to be of type TRANSIT
        std::vector<TissueCell> cells;
        
        for(unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TissueCell cell(TRANSIT, HEALTHY, new SimpleWntCellCycleModel());            
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*(p_params->GetTransitCellG1Duration()
                                                +p_params->GetSG2MDuration());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cell.SetSymmetricDivision();
            cells.push_back(cell);
        }
                                                          
        // Make a tissue and pass in the Wnt gradient
        Tissue<2> crypt(*p_mesh, cells);         
        crypt.SetGhostNodes(ghost_node_indices);          
        
        WntGradient::Instance()->SetType(RADIAL); 
        WntGradient::Instance()->SetTissue(crypt);   
        
        // Create the spring system
        CryptProjectionSpringSystem* p_spring_system = new CryptProjectionSpringSystem(crypt);
                
        // Make a tissue simulation
        TissueSimulation<2> crypt_projection_simulator(crypt, p_spring_system);
        
        // Create a radial cell killer and pass it in to the tissue simulation        
        c_vector<double,2> centre = zero_vector<double>(2);
        double crypt_radius = pow(CancerParameters::Instance()->GetCryptLength()/a, 1.0/b);

        RadialSloughingCellKiller* p_killer = new RadialSloughingCellKiller(&crypt, centre, crypt_radius);
        crypt_projection_simulator.AddCellKiller(p_killer);
        
        // Set up the simulation
        crypt_projection_simulator.SetOutputDirectory("CryptProjectionSimulation");
        crypt_projection_simulator.SetEndTime(0.25);
        crypt_projection_simulator.SetMaxCells(1000);
        crypt_projection_simulator.SetMaxElements(2000);
        
        // Run the simulation
        double start_time = std::clock();        
        TS_ASSERT_THROWS_NOTHING(crypt_projection_simulator.Solve());                
        double end_time = std::clock();
        
        // Print out time taken to run crypt simulation    
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "Time to perform crypt projection simulation was = " << elapsed_time << "\n";
        
        // These cells just divided and have been gradually moving apart. These results 
        // are from time 0.25, which is also tested below after a save and a load.
        std::vector<double> node_302_location = crypt_projection_simulator.GetNodeLocation(302);
        TS_ASSERT_DELTA(node_302_location[0], -0.0954, 1e-4);
        TS_ASSERT_DELTA(node_302_location[1], 0.2475, 1e-4);
        std::vector<double> node_506_location = crypt_projection_simulator.GetNodeLocation(506);
        TS_ASSERT_DELTA(node_506_location[0], -0.7046, 1e-4);
        TS_ASSERT_DELTA(node_506_location[1], 0.6178, 1e-4);

        // Test the Wnt gradient result
        TissueCell* p_cell = &(crypt.rGetCellAtNodeIndex(302));
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(p_cell), 0.9991, 1e-4);
        p_cell = &(crypt.rGetCellAtNodeIndex(506));
        TS_ASSERT_DELTA(WntGradient::Instance()->GetWntLevel(p_cell), 0.9898, 1e-4);
        
        // Tidy up
        delete p_killer;
        WntGradient::Destroy();
    }
    
};


#endif /*TESTCRYPTPROJECTIONSIMULATION_HPP_*/
