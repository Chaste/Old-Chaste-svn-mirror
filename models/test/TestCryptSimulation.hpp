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
        
        CryptSimulation simulator(mesh);
        simulator.SetEndTime(0.25);
        TS_ASSERT_THROWS_ANYTHING(simulator.Solve());
    }
    
    // No birth because SetIncludeRandomBirth() was not called.
    // No death because the initial length is 21 and only cells
    // with position greater than 22 (=default crypt length)
    // are sloughed off.
    void test1dChainWithNoBirth(void) throw(Exception)
    {
        Make1dCryptMesh("1D_crypt_mesh", 22, 21);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> shifted_point;
        shifted_point.rGetLocation()[0]=10.5;
        mesh.SetNode(10, shifted_point);
        
        CryptSimulation simulator(mesh);
        simulator.SetOutputDirectory("CryptWithNoBirthAndNoDeath");
        simulator.SetMaxCells(22);
        simulator.SetEndTime(0.25);
        simulator.Solve();
        
        // Note: converges very slowly so large tolerance of 0.1
        for (int index = 0; index<mesh.GetNumNodes(); index++)
        {
            Node<1> *p_node = mesh.GetNode(index);
            TS_ASSERT(!p_node->IsDeleted());
            Point<1> point = p_node->rGetPoint();
            double position = point.rGetLocation()[0];
            TS_ASSERT_DELTA(position, index, 1e-1);
        }
        
        //CheckAgainstPreviousRun("CryptWithNoBirthAndNoDeath", 22);
    }
    
    // Death because this test has a cell starting at the end of
    // the crypt.
    void Test1dChainWithDeathAndNoBirth(void) throw(Exception)
    {
        Make1dCryptMesh("1D_crypt_mesh", 23, 22);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        Point<1> shifted_point;
        shifted_point.rGetLocation()[0]=10.5;
        mesh.SetNode(10, shifted_point);
        
        CryptSimulation simulator(mesh);
        simulator.SetOutputDirectory("CryptWithDeathButNoBirth");
        simulator.SetMaxCells(23);
        simulator.SetEndTime(0.25);
        simulator.Solve();
        
        unsigned dead_cells = 0;
        // Note: converges very slowly so small tolerance of 0.1
        for (int index = 0; index<mesh.GetNumAllNodes(); index++)
        {
            Node<1> *p_node = mesh.GetNode(index);
            if (!p_node->IsDeleted())
            {
                Point<1> point = p_node->rGetPoint();
                double position = point.rGetLocation()[0];
                TS_ASSERT_DELTA(position, index, 1e-1);
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
    
    // Includes random birth (random position in a random element), death and constant
    // rest length by default
    void Test1DChainWithBirthConstantRestLength() throw (Exception)
    {
        // Note that random numbers are reseeded with srandom(0) by the following constructor.
        RandomNumberGenerator rand_gen;
        CancerParameters *p_params = CancerParameters::Instance();
        
        Make1dCryptMesh("1D_crypt_mesh", 23, 22);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type=STEM;
            unsigned generation=0;
            double birth_time=-p_params->GetStemCellCycleTime()+1.0*i; //hours
            MeinekeCryptCell cell(cell_type, birth_time, generation, new StochasticCellCycleModel(&rand_gen));
            cell.SetNodeIndex(i);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells);
        simulator.SetOutputDirectory("CryptWithBirthConstantRestLength");
        simulator.SetMaxCells(33);
        simulator.SetIncludeRandomBirth();
        simulator.SetEndTime(10.0);
        simulator.Solve();
        
        CheckAgainstPreviousRun("CryptWithBirthConstantRestLength",33);
    }
    
    // Includes variable rest length but not fully working in the class (see
    // CryptSimulation.hpp and look for the
    // "(age1<1.0/time_scale && age2<1.0/time_scale && fabs(age1-age2)<1e-6)" if.
    void Test1DChainWithBirthVariableRestLength() throw (Exception)
    {
        RandomNumberGenerator rand_gen;
        Make1dCryptMesh("1D_crypt_mesh", 23, 22);
        
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type=STEM;
            unsigned generation=0;
            double birth_time= -1.0; //hours
            MeinekeCryptCell cell(cell_type, birth_time, generation, new StochasticCellCycleModel(&rand_gen));
            cell.SetNodeIndex(i);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells);
        simulator.SetOutputDirectory("CryptWithBirthVariableRestLength");
        simulator.SetMaxCells(25);
        simulator.SetIncludeRandomBirth();
        simulator.SetIncludeVariableRestLength();
        simulator.SetEndTime(10.0);
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithBirthVariableRestLength",25);
    }
    
    
    // Create a chain of meineke cells (1 stem, 14 transit cells and 8 differentiated)
    // and pass into the simulation class
    void Test1DChainWithMeinekeCells() throw (Exception)
    {
        RandomNumberGenerator rand_gen;
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
                birth_time = -rand_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -rand_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, 0.0, generation, new StochasticCellCycleModel(&rand_gen));
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
    }
    
    
    // same as Test1DChainWithBirthVariableRestLength but with Meineke cells.
    // (see comment for Test1DChainWithBirthVariableRestLength).
    void Test1DChainWithMeinekeCellsAndGrowth() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator rand_gen;
        
        double crypt_length = 22.0;
        p_params->SetCryptLength(crypt_length);
        
        Make1dCryptMesh("1D_crypt_mesh", 23, crypt_length);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<1,1> mesh_reader(testoutput_dir+"/CryptMesh/1D_crypt_mesh");
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
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
                birth_time = -rand_gen.ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -rand_gen.ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
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
    }
    
    

    /////////////////////////////////////////////////////////////////////////////
    // create a chain of 1 stem cell and the test differentiated and check
    // that there is correct number of cells and they are in the correct order 
    /////////////////////////////////////////////////////////////////////////////
    void Test1dChainCorrectCellNumbers()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerator rand_gen;
        
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
            MeinekeCryptCell cell(cell_type, 0.0, generation, new FixedCellCycleModel());
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
        
        CryptSimulation simulator(mesh, cells);
        
        simulator.SetOutputDirectory("Crypt1dTestCorrectCellNumbers");
        simulator.SetMaxCells(50);
        simulator.SetEndTime(40);
        
        simulator.SetIncludeVariableRestLength();
        
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        std::vector<MeinekeCryptCell> cells_after_simulation = simulator.GetCells();
        //Warning - there are 6 live cells and one dead one sloughed off still in existance.
        //n.b. throughout the simulation 2 cells are sloughed off but one place is reused
        TS_ASSERT_EQUALS((int) cells_after_simulation.size(),7);
        for (int index = 0; index<mesh.GetNumAllNodes(); index++)
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
    }
};

#endif /*TESTCRYPTSIMULATION_HPP_*/

