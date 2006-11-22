#ifndef CRYPTHONEYCOMBMESHGENERATOR_HPP_
#define CRYPTHONEYCOMBMESHGENERATOR_HPP_
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"

class CryptHoneycombMeshGenerator 
{
private:
    ConformingTetrahedralMesh<2,2>* mpMesh;
    std::vector<int> mGhostNodeIndices;
    std::string mMeshFilename; 
    double mCryptWidth;
    double mCryptDepth;
    
    void MakeHoneycombMeshFiles(unsigned numCellWidth, unsigned numCellDepth)
    {
        double x0=-2.0;
        double y0=-sqrt(3)/2;
        OutputFileHandler output_file_handler("");
        out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
        (*p_node_file) << std::scientific;
        
        out_stream p_elem_file = output_file_handler.OpenOutputFile(mMeshFilename+".ele");
        (*p_elem_file) << std::scientific;
        unsigned num_nodes_along_width = numCellWidth+4;
        unsigned num_nodes_along_depth = numCellDepth+4;
        unsigned total_num_nodes       = num_nodes_along_width*num_nodes_along_depth;
        unsigned num_elem_along_width = num_nodes_along_width-1;
        unsigned num_elem_along_depth = num_nodes_along_depth-1;
        unsigned num_elem             = 2*num_elem_along_width*num_elem_along_depth;
        unsigned num_edges            = 3*num_elem_along_width*num_elem_along_depth + num_elem_along_width + num_elem_along_depth;
        
        (*p_node_file) << total_num_nodes << "\t2\t0\t1" << std::endl;
        unsigned node = 0;
        for (unsigned i = 0; i < num_nodes_along_depth; i++)
        {
            for (unsigned j = 0; j < num_nodes_along_width; j++)
            {
                int b = 0;
                if ((i==0) || (i==num_nodes_along_depth-1) || (j==0) || (j==num_nodes_along_width-1))
                {
                    b = 1;
                }
                double x = x0 + ((double)j + 0.25*(1+ pow(-1,i+1))) ;
                
                double y = y0 + (sqrt(3)/2)*(double)i;
                
                (*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << b << std::endl;
            }
        }
        p_node_file->close();
        
        out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
        (*p_node_file) << std::scientific;
        
        (*p_elem_file) << num_elem << "\t3\t0" << std::endl;
        (*p_edge_file) << num_edges << "\t3\t0\t1" << std::endl;
        
        unsigned elem = 0;
        unsigned edge = 0;
        for (unsigned i = 0; i < num_elem_along_depth; i++)
        {
            for (unsigned j = 0; j < num_elem_along_width; j++)
            {
                int node0 =     i*num_nodes_along_width + j;
                int node1 =     i*num_nodes_along_width + j+1;
                int node2 = (i+1)*num_nodes_along_width + j;
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
                
                node0 = i*num_nodes_along_width + j + 1;
                
                if (i%2 != 0)
                {
                    node0 = node0 - 1;
                }
                node1 = (i+1)*num_nodes_along_width + j+1;
                node2 = (i+1)*num_nodes_along_width + j;
                
                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
            }
        }
        
        for (unsigned i = 0; i < num_elem_along_depth; i++)
        {
            int node0 = (i+1)*num_nodes_along_width-1;
            int node1 = (i+2)*num_nodes_along_width-1;
            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
        }
        
        for (unsigned j = 0; j < num_elem_along_width; j++)
        {
            int node0 =  num_nodes_along_width*(num_nodes_along_depth-1) + j;
            int node1 =  num_nodes_along_width*(num_nodes_along_depth-1) + j+1;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
        }
        
        p_elem_file->close();
        p_edge_file->close();
    }


    void ComputeGhostNodes()
    {
        assert(mpMesh!=NULL);
        for (int i=0; i<mpMesh->GetNumNodes(); i++)
        {
            double x = mpMesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = mpMesh->GetNode(i)->GetPoint().rGetLocation()[1];
            if ((x<0)||(x>mCryptWidth)||(y>mCryptDepth)||(y<0))
            {
                mGhostNodeIndices.push_back(i);
            }
        }
    }


public:
    CryptHoneycombMeshGenerator(int numCellWidth, int numCellDepth)
    {
        mCryptWidth = numCellWidth*1; //*1 because cells are considered to be size one
        mCryptDepth = sqrt(3)*numCellDepth/2;
        
        mMeshFilename = "2D_temporary_crypt_mesh";
        MakeHoneycombMeshFiles(numCellWidth,numCellDepth);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<2,2> mesh_reader(testoutput_dir+ mMeshFilename);
        
        mpMesh = new ConformingTetrahedralMesh<2,2>;
        mpMesh->ConstructFromMeshReader(mesh_reader);
        
        ComputeGhostNodes();
    }
    
    ConformingTetrahedralMesh<2,2>* GetMesh()
    {
        return mpMesh;   
    }
    
    std::vector<int> GetGhostNodeIndices()
    {
        return mGhostNodeIndices;
    }



};
#endif /*CRYPTHONEYCOMBMESHGENERATOR_HPP_*/
