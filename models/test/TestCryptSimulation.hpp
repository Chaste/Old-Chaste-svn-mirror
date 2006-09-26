#ifndef TESTCRYPTSIMULATION_HPP_
#define TESTCRYPTSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include <iostream>
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
    
    void CheckAgainstPreviousRun(std::string resultDirectory)
    {
    	OutputFileHandler output_file_handler("");
    	std::string testoutput_dir = output_file_handler.GetTestOutputDirectory();
    	std::string compare_command = "cmp " + 
            testoutput_dir + "/" + resultDirectory + "/results " +
            "models/test/data/" + resultDirectory + "/results";
        TS_ASSERT_EQUALS(system(compare_command.c_str()), 0);
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
        simulator.SetEndTime(6.0);
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
        simulator.SetEndTime(6.0);
        simulator.Solve();
        
        
        // Note: converges very slowly so large tolerance of 0.1
        for (int index = 0; index<mesh.GetNumNodes(); index++)
        {
            Node<1> *p_node = mesh.GetNodeAt(index);
            TS_ASSERT(!p_node->IsDeleted());
            Point<1> point = p_node->rGetPoint();
            double position = point.rGetLocation()[0];
            TS_ASSERT_DELTA(position, index, 1e-1);
        }
        
        CheckAgainstPreviousRun("CryptWithNoBirthAndNoDeath");
    }
    
    // Death because this test has a cell starting at the end of
    // the crypt.
    void test1dChainWithDeathAndNoBirth(void) throw(Exception)
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
        simulator.SetEndTime(6.0);
        simulator.Solve();
        
        unsigned dead_cells = 0;
        // Note: converges very slowly so small tolerance of 0.1
        for (int index = 0; index<mesh.GetNumAllNodes(); index++)
        {
            Node<1> *p_node = mesh.GetNodeAt(index);
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
        
        CheckAgainstPreviousRun("CryptWithDeathButNoBirth");
    }
    
    // Includes random birth (random position in a random element), death and constant
    // rest length by default
    void Test1DChainWithBirthConstantRestLength()
    {
    	// Note that random numbers are reseeded with srandom(0) by the following constructor.
     	RandomNumberGenerators *pGen=new RandomNumberGenerators;
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
            MeinekeCryptCell cell(cell_type, birth_time, generation, new StochasticCellCycleModel(pGen));
            cell.SetNodeIndex(i);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells);
        simulator.SetOutputDirectory("CryptWithBirthConstantRestLength");
        simulator.SetIncludeRandomBirth();
        simulator.SetEndTime(240.0);
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithBirthConstantRestLength");
        delete pGen;
    }
    
    // Includes variable rest length but not fully working in the class (see
    // CryptSimulation.hpp and look for the
    // "(age1<1.0/time_scale && age2<1.0/time_scale && fabs(age1-age2)<1e-6)" if.
    void Test1DChainWithBirthVariableRestLength() throw (Exception)
    {
        RandomNumberGenerators *pGen=new RandomNumberGenerators;
               
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
            double birth_time= -1.0*i; //hours
            MeinekeCryptCell cell(cell_type, birth_time, generation, new StochasticCellCycleModel(pGen));
            cell.SetNodeIndex(i);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells); 
        
        simulator.SetOutputDirectory("CryptWithBirthVariableRestLength");
        simulator.SetIncludeRandomBirth();
        simulator.SetIncludeVariableRestLength();
        simulator.SetEndTime(240.0);
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithBirthVariableRestLength");
    }
    
    // Create a chain of meineke cells (1 stem, 14 transit cells and 8 differentiated)
    // and pass into the simulation class
    void Test1DChainWithMeinekeCells() throw (Exception)
    {
        RandomNumberGenerators *pGen=new RandomNumberGenerators;
  		CancerParameters *p_params = CancerParameters::Instance();

        double crypt_length = 22.0;
        
        
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
                birth_time = -pGen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -pGen->ranf()*p_params->GetTransitCellCycleTime(); //hours
            }
            else
            {
                cell_type = DIFFERENTIATED;
                generation = 4;
                birth_time = 0; //hours
            }
            MeinekeCryptCell cell(cell_type, 0.0, generation, new StochasticCellCycleModel(pGen));
            cell.SetNodeIndex(i);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }
              
        CryptSimulation simulator(mesh, cells);
        simulator.SetOutputDirectory("CryptWithCells");
        simulator.SetEndTime(240.0);
        simulator.SetCryptLength(crypt_length);
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithCells");
    }
    
    
    // same as Test1DChainWithBirthVariableRestLength but with Meineke cells.
    // (see comment for Test1DChainWithBirthVariableRestLength).
    void Test1DChainWithMeinekeCellsAndGrowth() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        RandomNumberGenerators *pGen = new RandomNumberGenerators;
        double crypt_length = 22.0;
        
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
                birth_time = -pGen->ranf()*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (i < 15)
            {
                cell_type = TRANSIT;
                generation = 1 + (i - 1) / 5;
                birth_time = -pGen->ranf()*p_params->GetTransitCellCycleTime(); //hours
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
        simulator.SetEndTime(240.0);
        simulator.SetCryptLength(crypt_length);
        simulator.SetIncludeVariableRestLength();
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        
        CheckAgainstPreviousRun("CryptWithCellsAndGrowth");
    }
    
    void DoNotTestForSingleStemCell() throw (Exception)
    {
    	RandomNumberGenerators *pGen = new RandomNumberGenerators;
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
            CryptCellType cell_type;
            unsigned generation;
            double birth_time=0; //hours
            MeinekeCryptCell cell(cell_type, birth_time, generation, new StochasticCellCycleModel(pGen));
            cell.SetNodeIndex(i);
            cells.push_back(cell);
        }
        
        
        CryptSimulation simulator(mesh, cells);
        simulator.SetOutputDirectory("CryptSingleStemCellCheck");
        simulator.SetIncludeRandomBirth();
        simulator.SetEndTime(48.0);
        TS_ASSERT_THROWS_NOTHING( simulator.Solve() );
        CheckAgainstPreviousRun("CryptSingleStemCellCheck");
        //Simulation Run ended
        
        // Testing for only one stem cell
    	ColumnDataReader *mpResultReader;
    	mpResultReader = new ColumnDataReader("/tmp/chaste/testoutput/CryptSingleStemCellCheck", "results",
                                                 false);

		delete pGen;
    }
    
};

#endif /*TESTCRYPTSIMULATION_HPP_*/

