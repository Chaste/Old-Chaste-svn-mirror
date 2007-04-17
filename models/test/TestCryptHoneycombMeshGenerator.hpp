#ifndef TESTCRYPTHONECOMBMESHGENERATOR_HPP_
#define TESTCRYPTHONECOMBMESHGENERATOR_HPP_
#include <cxxtest/TestSuite.h>
#include "CryptHoneycombMeshGenerator.hpp"

class TestCryptHoneycombMeshGenerator : public CxxTest::TestSuite
{
private:
    void Output2DNodesToFile(ConformingTetrahedralMesh<2,2>* p_mesh, std::string fileName)
    {
        OutputFileHandler handler("");
        out_stream file=handler.OpenOutputFile(fileName);
        
        unsigned num_nodes=p_mesh->GetNumNodes();
    
        for (unsigned i=0; i<num_nodes; i++)
        {
            c_vector<double, 2> location = p_mesh->GetNode(i)->rGetLocation();
            (*file) << location[0] << "\t" << location[1] << "\n" << std::flush;
        }
        
        file->close();
    }
    
    void Output2DNodesToFileCylindrical(Cylindrical2dMesh* p_mesh, std::string fileName)
    {
        OutputFileHandler handler("");
        out_stream file=handler.OpenOutputFile(fileName);
        
        unsigned num_nodes=p_mesh->GetNumNodes();
    
        for (unsigned i=0; i<num_nodes; i++)
        {
            c_vector<double, 2> location = p_mesh->GetNode(i)->rGetLocation();
            (*file) << location[0] << "\t" << location[1] << "\n" << std::flush;
        }
        
        file->close();
    }
    
public:

    void TestCryptHoneycombMeshGeneratorCylindricalRelaxed() throw(Exception)
    {
        unsigned num_cells_width = 8;
        unsigned num_cells_depth = 22;
        unsigned ghosts = 2;
        
        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth, ghosts, true);
         
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();             
        // check the mesh
        
        ConformingTetrahedralMesh<2,2>* p_mesh2;
        
        TS_ASSERT_THROWS_ANYTHING(p_mesh2=generator.GetMesh());
                
        
        Output2DNodesToFileCylindrical(p_mesh, "cylindrical_node_positions.dat");
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),(num_cells_width)*(num_cells_depth+2*ghosts));
                
        // zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0],0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1],-(double)ghosts*sqrt(3.0)/2.0 , 1e-5);
        
        // first real node
        int index = (num_cells_width)*ghosts; // 4 here is the number of ghost nodes in a row
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], 0.0,1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 0.0,1e-12);
        
        // last real node
        index = (ghosts+num_cells_depth)*(num_cells_width)-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], 7.5,1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 21.0*sqrt(3)/2.0,1e-4);
        
        // last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], 7.5,1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], (ghosts+num_cells_depth-1)*sqrt(3.0)/2.0,1e-4);
        
        // check the ghost nodes
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        for (unsigned i=0; i< num_cells_width*ghosts ; i++)
        {
            TS_ASSERT_EQUALS(ghost_node_indices[i],i);
        }
        TS_ASSERT_EQUALS(ghost_node_indices[num_cells_width*ghosts+1],(ghosts+num_cells_depth)*num_cells_width+1);
        
        CancerParameters* p_params = CancerParameters::Instance();
        TS_ASSERT_DELTA(p_params->GetCryptWidth(), (double)num_cells_width, 1e-7);
        TS_ASSERT_DELTA(p_params->GetCryptLength(), sqrt(3)*num_cells_depth/2.0, 1e-7);
    }
    
    void TestCryptHoneycombMeshGeneratorCylindricalCompressed() throw(Exception)
    {
        unsigned num_cells_width = 8;
        unsigned num_cells_depth = 22;
        double width = 7.0;
        unsigned ghosts = 2;
        
        double x_factor = width/(double)num_cells_width;
        
        CryptHoneycombMeshGenerator generator(num_cells_width, num_cells_depth, width, ghosts, true);
         
        Cylindrical2dMesh* p_mesh=generator.GetCylindricalMesh();             
        // check the mesh
        
        ConformingTetrahedralMesh<2,2>* p_mesh2;
        
        TS_ASSERT_THROWS_ANYTHING(p_mesh2=generator.GetMesh());
                
        
        Output2DNodesToFileCylindrical(p_mesh, "cylindrical_node_positions.dat");
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(),(num_cells_width)*(num_cells_depth+2*ghosts));
                
        // zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0],0.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1],-x_factor*(double)ghosts*sqrt(3.0)/2.0 , 1e-5);
        
        // first real node
        int index = (num_cells_width)*ghosts; // 4 here is the number of ghost nodes in a row
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], 0.0,1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], 0.0,1e-12);
        
        // last real node
        index = (ghosts+num_cells_depth)*(num_cells_width)-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], x_factor*7.5,1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], x_factor*21.0*sqrt(3)/2.0,1e-4);
        
        // last node
        int last_node = p_mesh->GetNumNodes()-1;
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], x_factor*7.5,1e-12);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], x_factor*(ghosts+num_cells_depth-1)*sqrt(3.0)/2.0,1e-4);
        
        // check the ghost nodes
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        for (unsigned i=0; i< num_cells_width*ghosts ; i++)
        {
            TS_ASSERT_EQUALS(ghost_node_indices[i],i);
        }
        TS_ASSERT_EQUALS(ghost_node_indices[num_cells_width*ghosts+1],(ghosts+num_cells_depth)*num_cells_width+1);
        
        CancerParameters* p_params = CancerParameters::Instance();
        TS_ASSERT_DELTA(p_params->GetCryptWidth(), x_factor*(double)num_cells_width, 1e-7);
        TS_ASSERT_DELTA(p_params->GetCryptLength(), x_factor*sqrt(3)*num_cells_depth/2.0, 1e-7);
    }
    
    void TestCryptOldPeriodicHoneycombMeshGeneratorRelaxed() throw(Exception)
    {
        int num_cells_width = 8;
        int num_cells_depth = 12;
        double width = 8.0;
        unsigned ghosts = 4;
        
        CryptHoneycombMeshGenerator generator(num_cells_width,num_cells_depth,width,ghosts,false);
        
        double length = (double)num_cells_depth*(sqrt(3)/2)*width/(double)num_cells_width;
        
        // check the mesh
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        TS_ASSERT_THROWS_ANYTHING(p_mesh=generator.GetCylindricalMesh());
        
        TS_ASSERT_EQUALS((unsigned)p_mesh->GetNumNodes(),(num_cells_width+1+2*ghosts)*(num_cells_depth+2*ghosts));
        
        // Scaling Factor
        double F = (width/(double)num_cells_width);
        double spooky = (double)ghosts;
        
        // zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0],-spooky*F, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1],-spooky*F*sqrt(3)/2,1e-6);
        
        unsigned this_many_ghosts_at_start = ((2*ghosts+num_cells_width+1)*ghosts+ghosts);
        
        // first real node
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[0], 0.0,1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[1], 0.0,1e-6);
        
        // last real node
        int index = (2*ghosts+num_cells_width+1)*(ghosts+num_cells_depth)+ghosts+num_cells_width;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], width,  1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], length, 1e-4);
        
        // last node
        int last_node = p_mesh->GetNumNodes()-1;
        double last_node_y = length+(spooky-1)*F*(sqrt(3)/2);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], width+(spooky+0.5)*F,1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], last_node_y,1e-5);
        
        // check the ghost nodes
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        for (unsigned i=0; i<this_many_ghosts_at_start; i++)
        {
            TS_ASSERT(ghost_node_indices[i]==i);
        }
        // Check that the next ghost node is the other side of the stem cells...
        TS_ASSERT_EQUALS(ghost_node_indices[this_many_ghosts_at_start],this_many_ghosts_at_start+num_cells_width+1)
                
        CancerParameters* p_params = CancerParameters::Instance();
        TS_ASSERT_DELTA(p_params->GetCryptWidth(), width, 1e-7);
        TS_ASSERT_DELTA(p_params->GetCryptLength(), length, 1e-7);
    }
    
    void TestCryptOldPeriodicHoneycombMeshGeneratorCompressed() throw(Exception)
    {
        int num_cells_width = 8;
        int num_cells_depth = 12;
        double width = 6.0;
        unsigned ghosts = 4;
        
        CryptHoneycombMeshGenerator generator(num_cells_width,num_cells_depth,width,ghosts,false);
        
        double length = (double)num_cells_depth*(sqrt(3)/2)*width/(double)num_cells_width;
        
        // check the mesh
        ConformingTetrahedralMesh<2,2>* p_mesh=generator.GetMesh();
        
        TS_ASSERT_THROWS_ANYTHING(p_mesh=generator.GetCylindricalMesh());
        
        TS_ASSERT_EQUALS((unsigned)p_mesh->GetNumNodes(),(num_cells_width+1+2*ghosts)*(num_cells_depth+2*ghosts));
        
        // Scaling Factor
        double F = (width/(double)num_cells_width);
        double spooky = (double)ghosts;
        
        // zeroth node
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[0],-spooky*F, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->GetPoint()[1],-spooky*F*sqrt(3)/2,1e-6);
        
        unsigned this_many_ghosts_at_start = ((2*ghosts+num_cells_width+1)*ghosts+ghosts);
        
        // first real node
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[0], 0.0,1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(this_many_ghosts_at_start)->GetPoint()[1], 0.0,1e-6);
        
        // last real node
        int index = (2*ghosts+num_cells_width+1)*(ghosts+num_cells_depth)+ghosts+num_cells_width;
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[0], width,  1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(index)->GetPoint()[1], length, 1e-4);
        
        // last node
        int last_node = p_mesh->GetNumNodes()-1;
        double last_node_y = length+(spooky-1)*F*(sqrt(3)/2);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[0], width+(spooky+0.5)*F,1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(last_node)->GetPoint()[1], last_node_y,1e-6);
        
        // check the ghost nodes
        std::vector<unsigned> ghost_node_indices = generator.GetGhostNodeIndices();
        
        for (unsigned i=0; i<this_many_ghosts_at_start; i++)
        {
            TS_ASSERT(ghost_node_indices[i]==i);
        }
        // Check that the next ghost node is the other side of the stem cells...
        TS_ASSERT_EQUALS(ghost_node_indices[this_many_ghosts_at_start],this_many_ghosts_at_start+num_cells_width+1)
                
        CancerParameters* p_params = CancerParameters::Instance();
        TS_ASSERT_DELTA(p_params->GetCryptWidth(), width, 1e-7);
        TS_ASSERT_DELTA(p_params->GetCryptLength(), length, 1e-7);
    }
    

};


#endif /*TESTCRYPTHONECOMBMESHGENERATOR_HPP_*/
