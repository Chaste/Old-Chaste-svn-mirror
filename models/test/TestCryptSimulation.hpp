#ifndef TESTCRYPTSIMULATION_HPP_
#define TESTCRYPTSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation.hpp"

#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"

#include "CancerParameters.hpp"
#include "ColumnDataReader.hpp"


class TestCryptSimulation : public CxxTest::TestSuite
{
    void Make1dCryptMesh(std::string meshFilename, unsigned numNodesInEachDimension, double length)
    {
        OutputFileHandler output_file_handler("CryptMesh");
        out_stream p_node_file = output_file_handler.OpenOutputFile(meshFilename+".node");
        (*p_node_file) << std::scientific;
        
        out_stream p_elem_file = output_file_handler.OpenOutputFile(meshFilename+".ele");
        (*p_elem_file) << std::scientific;
        
        unsigned num_nodes    = numNodesInEachDimension;
        unsigned num_elements = num_nodes - 1;
        
        (*p_node_file) << num_nodes << "\t1\t0\t1" << std::endl;
        for (unsigned i = 0; i < num_nodes; i++)
        {
            int b = 0;
            if ((i == 0) || (i == num_nodes-1))
            {
                b = 1;
            }
            double x = length*i/(num_nodes-1);
            (*p_node_file) << i << "\t" << x << "\t" << b << std::endl;
        }
        p_node_file->close();
        
        (*p_elem_file) << num_elements << "\t2\t0" << std::endl;
        for (unsigned i = 0; i < num_elements; i++)
        {
            (*p_elem_file) << i << "\t" << i << "\t" << i+1 << std::endl;
        }
        p_elem_file->close();
    }
    
    void CheckAgainstPreviousRun(std::string resultDirectory, unsigned maxCells)
    {
        std::cout << "Comparing " << resultDirectory << std::endl << std::flush;
        
        ColumnDataReader computed_results = ColumnDataReader(resultDirectory,
                                                             "tabulated_results",
                                                             true);
                                                             
        ColumnDataReader expected_results = ColumnDataReader("models/test/data/" + resultDirectory,
                                                             "tabulated_results",
                                                             false);
                                                             
        for (unsigned cell=0; cell<maxCells; cell++)
        {
            std::stringstream cell_type_var_name;
            std::stringstream cell_position_var_name;
            cell_type_var_name << "cell_type_" << cell;
            cell_position_var_name << "cell_position_" << cell;
            
            // Vector of Cell Types
            std::vector<double> expected_cell_types = expected_results.GetValues(cell_type_var_name.str());
            std::vector<double> computed_cell_types = computed_results.GetValues(cell_type_var_name.str());
            
            //Vector of Cell Positions
            std::vector<double> expected_cell_positions = expected_results.GetValues(cell_position_var_name.str());
            std::vector<double> computed_cell_positions = computed_results.GetValues(cell_position_var_name.str());
            
            //Comparing expected and computed vector length
            TS_ASSERT_EQUALS(expected_cell_types.size(), computed_cell_types.size());
            TS_ASSERT_EQUALS(expected_cell_positions.size(), computed_cell_positions.size());
            
            //Walkthrough of the expected and computed vectors
            for (unsigned time_step = 0; time_step < expected_cell_types.size(); time_step++)
            {
                TS_ASSERT_EQUALS(expected_cell_types[time_step], computed_cell_types[time_step]);
                TS_ASSERT_DELTA(expected_cell_positions[time_step], computed_cell_positions[time_step],1e-6);
            }
        }
    }
    
public:

    void TestExceptions(void)
    {
        Make1dCryptMesh("1D_crypt_mesh", 22, 21);
        
        OutputFileHandler output_file_handler("");
        std::string testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // throws because start time not set on simulation time
        TS_ASSERT_THROWS_ANYTHING(CryptSimulation bad_simulator(mesh));
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptSimulation simulator(mesh);
        simulator.SetDt(1.0/120.0); // i.e. 30s: this is the default - just for coverage
        simulator.SetEndTime(0.25);
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());
        
    }
    
    // No birth because SetIncludeRandomBirth() was not called.
    // No death because the initial length is 21 and only cells
    // with position greater than 22 (=default crypt length)
    // are sloughed off.
    void Test1dChainWithNoBirth(void) throw(Exception)
    {
        CancerParameters::Instance()->Reset();
        
        Make1dCryptMesh("1D_crypt_mesh", 22, 21);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0]=10.5;
        mesh.SetNode(10, shifted_point);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptSimulation simulator(mesh);
        simulator.SetOutputDirectory("CryptWithNoBirthAndNoDeath");
        simulator.SetMaxCells(22);
        simulator.SetEndTime(0.25);
        simulator.Solve();
        
        // Note: converges very slowly so large tolerance of 0.1
        for (unsigned index = 0; index<mesh.GetNumNodes(); index++)
        {
            Node<1> *p_node = mesh.GetNode(index);
            TS_ASSERT(!p_node->IsDeleted());
            const c_vector<double,1>& r_node_loc = p_node->rGetLocation();
            TS_ASSERT_DELTA(r_node_loc[0], index, 1e-1);
        }
        
        CheckAgainstPreviousRun("CryptWithNoBirthAndNoDeath", 22);
    }
    
    // Death because this test has a cell starting at the end of
    // the crypt.
    void Test1dChainWithDeathAndNoBirth(void) throw(Exception)
    {
        CancerParameters::Instance()->Reset();

        Make1dCryptMesh("1D_crypt_mesh", 23, 22);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0]=10.5;
        mesh.SetNode(10, shifted_point);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        CryptSimulation simulator(mesh);
        simulator.SetOutputDirectory("CryptWithDeathButNoBirth");
        simulator.SetMaxCells(23);
        simulator.SetEndTime(0.25);
        simulator.Solve();
        
        unsigned dead_cells = 0;
        // Note: converges very slowly so small tolerance of 0.1
        for (unsigned index = 0; index<mesh.GetNumAllNodes(); index++)
        {
            Node<1> *p_node = mesh.GetNode(index);
            if (!p_node->IsDeleted())
            {
                const c_vector<double,1>& r_node_loc = p_node->rGetLocation();
                TS_ASSERT_DELTA(r_node_loc[0], index, 1e-1);
                // Note: since there is no birth, the position of each node should still settle
                // down to its index
            }
            else
            {
                dead_cells++;
            }
        }
        
        TS_ASSERT_LESS_THAN(0u, dead_cells);
        
        CheckAgainstPreviousRun("CryptWithDeathButNoBirth", 23);
    }
    
        
    
    // Create a chain of meineke cells (1 stem, 14 transit cells and 8 differentiated)
    // and pass into the simulation class
    void Test1DChainWithMeinekeCells() throw (Exception)
    {
        CancerParameters::Instance()->Reset();
        RandomNumberGenerator *p_rand_gen=RandomNumberGenerator::Instance();
        
        CancerParameters *p_params = CancerParameters::Instance();
        
        double crypt_length = 22.0;
        p_params->SetCryptLength(crypt_length);
        
        Make1dCryptMesh("1D_crypt_mesh", 23, crypt_length);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_rand_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new StochasticCellCycleModel);
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells);
        
        simulator.SetOutputDirectory("CryptWithCells");
        simulator.SetMaxCells(50);
        simulator.SetEndTime(10.0);
        
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithCells",50);
        RandomNumberGenerator::Destroy();
    }
    
    
    // same as Test1DChainWithBirthVariableRestLength but with Meineke cells.
    // (see comment for Test1DChainWithBirthVariableRestLength).
    void Test1DChainWithMeinekeCellsAndGrowth() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        RandomNumberGenerator *p_rand_gen = RandomNumberGenerator::Instance();
        
        double crypt_length = 22.0;
        p_params->SetCryptLength(crypt_length);
        
        Make1dCryptMesh("1D_crypt_mesh", 23, crypt_length);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_rand_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation simulator(mesh, cells);
        
        simulator.SetOutputDirectory("CryptWithCellsAndGrowth");
        simulator.SetMaxCells(50);
        simulator.SetEndTime(10.0);
        
        simulator.SetIncludeVariableRestLength();
        
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithCellsAndGrowth",50);
        RandomNumberGenerator::Destroy();
    }
    
    
      
    void Test1DChainWithTysonNovakCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        RandomNumberGenerator *p_rand_gen = RandomNumberGenerator::Instance();
        
        double crypt_length = 22.0;
        p_params->SetCryptLength(crypt_length);
        
        Make1dCryptMesh("1D_crypt_mesh", 23, crypt_length);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // For Tyson-Novak Cells
        double temp_stem = p_params->GetStemCellCycleTime();
        double temp_transit = p_params->GetTransitCellCycleTime();
        p_params->SetStemCellCycleTime(0.15);
        p_params->SetTransitCellCycleTime(0.15);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -p_rand_gen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -p_rand_gen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new TysonNovakCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation simulator(mesh, cells);
        
        simulator.SetOutputDirectory("CryptWithTysonNovakCells");
        simulator.SetMaxCells(500);
        simulator.SetEndTime(1.35);
        
        simulator.SetIncludeVariableRestLength();
        
        simulator.Solve();
        
        cells = simulator.GetCells();
        // check we have had loads of birth
        // N.B. that if this test is run for longer it will fail 
        // because they get too squashed in (T&N is too quick).
        TS_ASSERT_EQUALS(cells.size() , num_cells+65u);
        
        p_params->SetStemCellCycleTime(temp_stem);
        p_params->SetTransitCellCycleTime(temp_transit);
    }
    
    
    
    /////////////////////////////////////////////////////////////////////////////
    // create a chain of 1 stem cell and the test differentiated and check
    // that there is correct number of cells and they are in the correct order
    /////////////////////////////////////////////////////////////////////////////
    void Test1dChainCorrectCellNumbers()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        p_params->Reset();
        RandomNumberGenerator::Instance();
        
        // check the stem cell cycle time is still 24 hrs, otherwise
        // this test might not pass
        TS_ASSERT_DELTA(p_params->GetStemCellCycleTime(), 24, 1e-12);
        TS_ASSERT_DELTA(p_params->GetTransitCellCycleTime(), 12, 1e-12);
        
        
        double crypt_length = 5.0;
        p_params->SetCryptLength(crypt_length);
        
        Make1dCryptMesh("1D_crypt_mesh", 6, crypt_length);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            if (i == 0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = 0; //hours - doesn't matter for stem cell;
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, HEALTHY, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells);
        simulator.SetOutputDirectory("Crypt1dTestCorrectCellNumbers");
        simulator.SetMaxCells(50);
        simulator.SetEndTime(40);
        
        simulator.SetIncludeVariableRestLength();
        
        simulator.Solve();
        
        std::vector<MeinekeCryptCell> cells_after_simulation = simulator.GetCells();
        //Warning - there are 6 live cells and one dead one sloughed off still in existance.
        //n.b. throughout the simulation 2 cells are sloughed off but one place is reused
        TS_ASSERT_EQUALS((int) cells_after_simulation.size(),7);
        for (unsigned index = 0; index<mesh.GetNumAllNodes(); index++)
        {
            if (!mesh.GetNode(index)->IsDeleted())
            {
                MeinekeCryptCell cell = cells_after_simulation[index];
                double x = mesh.GetNode(index)->GetPoint()[0];
                if (fabs(x)<1e-2)
                {
                    CryptCellType type  = cell.GetCellType();
                    TS_ASSERT(type==STEM);
                }
                else if ((fabs(x-1)<1e-2)||(fabs(x-2)<1e-2))
                {
                    CryptCellType type  = cell.GetCellType();
                    TS_ASSERT(type==TRANSIT);
                }
                else if ((fabs(x-3)<1e-2)||(fabs(x-4)<1e-2)||(fabs(x-5)<1e-2))
                {
                    CryptCellType type  = cell.GetCellType();
                    TS_ASSERT(type==DIFFERENTIATED);
                }
                else
                {
                    //There shouldn't be any cells at non-integer positions provided resting length =1
                    TS_ASSERT(false);
                }
            }
        }
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCRYPTSIMULATION_HPP_*/

