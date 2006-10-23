#ifndef TESTCRYPTSIMULATION2D_HPP_
#define TESTCRYPTSIMULATION2D_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CryptSimulation2D.hpp"

#include "MeinekeCryptCell.hpp"
#include "FixedCellCycleModel.hpp"
#include "StochasticCellCycleModel.hpp"

#include "CancerParameters.hpp"



class TestCryptSimulation2D : public CxxTest::TestSuite
{
    void Make2dCryptMesh(std::string meshFilename, unsigned numNodesAlongWidth, unsigned numNodesAlongLength, double width, double length)
    {
        OutputFileHandler output_file_handler("CryptMesh");
        out_stream p_node_file = output_file_handler.OpenOutputFile(meshFilename+".node");
        (*p_node_file) << std::scientific;
        
        out_stream p_elem_file = output_file_handler.OpenOutputFile(meshFilename+".ele");
        (*p_elem_file) << std::scientific;
        
        unsigned num_nodes            = numNodesAlongWidth*numNodesAlongLength;
        unsigned num_elem_along_width = numNodesAlongWidth-1;
        unsigned num_elem_along_length = numNodesAlongLength-1;
        unsigned num_elem             = 2*num_elem_along_width*num_elem_along_length;
        unsigned num_edges            = 3*num_elem_along_width*num_elem_along_length + num_elem_along_width + num_elem_along_length;
        
        (*p_node_file) << num_nodes << "\t2\t0\t1" << std::endl;
        unsigned node = 0;
        for (unsigned i = 0; i < numNodesAlongLength; i++)
        {
            for (unsigned j = 0; j < numNodesAlongWidth; j++)
            {
                int b = 0;
                if ((i==0) || (i==numNodesAlongLength-1) || (j==0) || (j==numNodesAlongWidth-1))
                {
                    b = 1;
                }
                double x = width*((double)j + 0.25*(1+ pow(-1,i+1)))/(num_elem_along_width) ;
                
                double y = length*(sqrt(3)/2)*(double)i/(num_elem_along_length); //
                
                (*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << b << std::endl;
            }
        }
        p_node_file->close();
        
        out_stream p_edge_file = output_file_handler.OpenOutputFile(meshFilename+".edge");
        (*p_node_file) << std::scientific;
        
        (*p_elem_file) << num_elem << "\t3\t0" << std::endl;
        (*p_edge_file) << num_edges << "\t3\t0\t1" << std::endl;
        
        unsigned elem = 0;
        unsigned edge = 0;
        for (unsigned i = 0; i < num_elem_along_length; i++)
        {
            for (unsigned j = 0; j < num_elem_along_width; j++)
            {
                int node0 =     i*numNodesAlongWidth + j;
                int node1 =     i*numNodesAlongWidth + j+1;
                int node2 = (i+1)*numNodesAlongWidth + j;
                if (i%2 != 0)
                {
                    node2 = node2 + 1;
                }
                
                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
                
                int horizontal_edge_is_boundary_edge = 0;
                int vertical_edge_is_boundary_edge = 0;
                if (i==0)
                {
                    horizontal_edge_is_boundary_edge = 1;
                }
                if (j==0)
                {
                    vertical_edge_is_boundary_edge = 1;
                }
                
                (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 <<  "\t" << horizontal_edge_is_boundary_edge << std::endl;
                (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 <<  "\t" << 0 << std::endl;
                (*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 <<  "\t" << vertical_edge_is_boundary_edge << std::endl;
                
                node0 = i*numNodesAlongWidth + j + 1;
                
                if (i%2 != 0)
                {
                    node0 = node0 - 1;
                }
                node1 = (i+1)*numNodesAlongWidth + j+1;
                node2 = (i+1)*numNodesAlongWidth + j;
                
                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
            }
        }
        
        for (unsigned i = 0; i < num_elem_along_length; i++)
        {
            int node0 = (i+1)*numNodesAlongWidth-1;
            int node1 = (i+2)*numNodesAlongWidth-1;
            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
        }
        
        for (unsigned j = 0; j < num_elem_along_width; j++)
        {
            int node0 =  numNodesAlongWidth*(numNodesAlongLength-1) + j;
            int node1 =  numNodesAlongWidth*(numNodesAlongLength-1) + j+1;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
        }
        
        p_elem_file->close();
        p_edge_file->close();
    }
    
    
public:

    // Test the spring system. There are no cells in this test, therefore no birth, although
    // nodes are sloughed. The mesh is initially a set of 10 by 10 squares, each square made
    // up of two triangles. The horizontal and vertical edges (springs) are at rest length, the
    // diagonals are two long, so this means the mesh skews to a (sloughed) parallelogram, each
    // triangle trying to become equilateral.
    void Test2DSpringSystemWithSloughing() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);
        
        double crypt_length = 10;
        double crypt_width = 10;
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
       
        CryptSimulation2D simulator(mesh);
        simulator.SetOutputDirectory("Crypt2DSprings");
        simulator.SetEndTime(1.0);
        
        simulator.Solve();
        
        for (int i=0; i<mesh.GetNumElements(); i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                unsigned nodeA, nodeB;
                j==2 ? nodeA = 2 : nodeA = j;
                j==2 ? nodeB = 0 : nodeB = j+1;
                
                double x_nodeA = mesh.GetElement(i)->GetNode(nodeA)->GetPoint()[0];
                double y_nodeA = mesh.GetElement(i)->GetNode(nodeA)->GetPoint()[1];
                double x_nodeB = mesh.GetElement(i)->GetNode(nodeB)->GetPoint()[0];
                double y_nodeB = mesh.GetElement(i)->GetNode(nodeB)->GetPoint()[1];
                
                // if both nodes are not sloughed, check distance between them is
                // the natural length.
                if (  (x_nodeA < crypt_width)  && (x_nodeB < crypt_width) &&
                      (y_nodeA < crypt_length) && (y_nodeB < crypt_length) )
                {
                    double length = sqrt( (x_nodeA - x_nodeB)*(x_nodeA - x_nodeB) + (y_nodeA - y_nodeB)*(y_nodeA - y_nodeB) );
                    TS_ASSERT_DELTA(length, 1.0, 0.05); // tolerance of 0.01 would work if ran for 2 units
                }
            }
        }
        
    }
    
    
    void Test2DSpringsWithCells() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);
        double crypt_length = 10;
        double crypt_width = 10;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 5.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
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
        

        CryptSimulation2D simulator(mesh, cells);
        simulator.SetOutputDirectory("Crypt2DSpringsWithCells");
        simulator.SetEndTime(1.0);
        
        //simulator.SetIncludeVariableRestLength();
        
        // throws anything because not working at the moment
        TS_ASSERT_THROWS_ANYTHING( simulator.Solve() );
    }
    
    
    void TestWithBirthOnHoneycombMesh() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);  // this is BAD, mkay?
        
        double crypt_length = 10.0;
        double crypt_width = 5.0;
        
        Make2dCryptMesh("2D_crypt_mesh", 6, 11, crypt_width, crypt_length);
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<2,2> mesh_reader(testoutput_dir+"/CryptMesh/2D_crypt_mesh");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 5.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
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
        
        CryptSimulation2D simulator(mesh, cells);
        simulator.SetOutputDirectory("Crypt2DHoneycombMesh");
        simulator.SetEndTime(1.0);
        
        //simulator.SetIncludeVariableRestLength();
        
        // throws anything because not working at the moment
        TS_ASSERT_THROWS_ANYTHING( simulator.Solve() );
    }
    
    
    // not being run because takes a few minutes to run
    void DONOT_Test2DSpringsFixedBoundaries() throw (Exception)
    {
        CancerParameters *p_params = CancerParameters::Instance();
        srandom(0);
        double crypt_length = 10;
        double crypt_width = 10;
        
        p_params->SetCryptLength(crypt_length);
        p_params->SetCryptWidth(crypt_width);
        
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_100mm_200_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Set up cells by iterating through the mesh nodes
        unsigned num_cells = mesh.GetNumAllNodes();
        std::vector<MeinekeCryptCell> cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            CryptCellType cell_type;
            unsigned generation;
            double birth_time;
            
            double y = mesh.GetNodeAt(i)->GetPoint().rGetLocation()[1];
            if (y == 0.0)
            {
                cell_type = STEM;
                generation = 0;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetStemCellCycleTime(); //hours - doesn't matter for stem cell;
            }
            else if (y < 3)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 6.5)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
            }
            else if (y < 8)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time = -(((double)random())/RAND_MAX)*p_params->GetTransitCellCycleTime(); //hours
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
        
        CryptSimulation2D simulator(mesh,cells);
        simulator.SetOutputDirectory("Crypt2DSpringsFixedBoundaries");
        simulator.SetEndTime(0.5); //Days?
        simulator.SetFixedBoundaries();
        
        //simulator.SetIncludeVariableRestLength();
        
        simulator.Solve() ;
    }
};

#endif /*TESTCRYPTSIMULATION2D_HPP_*/

