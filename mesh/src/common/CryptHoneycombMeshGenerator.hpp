#ifndef CRYPTHONEYCOMBMESHGENERATOR_HPP_
#define CRYPTHONEYCOMBMESHGENERATOR_HPP_

#include "Cylindrical2dMesh.cpp"
#include "NodeMap.hpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>
#include "OutputFileHandler.hpp"
#include "CancerParameters.hpp"

/**
 *  Generator of honeycomb mesh
 *
 *  This class takes in options such as width, height, number of ghost nodes
 *  and generates a mesh and ghost node info. NOTE: the user should delete the
 *  mesh after use
 */
class CryptHoneycombMeshGenerator
{
private:
    ConformingTetrahedralMesh<2,2>* mpMesh;
    Cylindrical2dMesh* mpCylindricalMesh;
    std::vector<unsigned> mGhostNodeIndices;
    std::string mMeshFilename;
    double mCryptWidth;
    double mCryptDepth;
    double mBottom;
    double mTop;
    std::vector<unsigned> mTopNodes;
    std::vector<unsigned> mBottomNodes;
    unsigned mNumCellWidth;
    bool mCylindrical;
       
    //////////////////////////////////////////////////////////////
    // Periodic Honeycomb mesh maker
    /////////////////////////////////////////////////////////////
    void Make2dPeriodicCryptMesh(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, double width, unsigned ghosts)
    {
        OutputFileHandler output_file_handler("");
        
        out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
        (*p_node_file) << std::scientific;
        
        out_stream p_elem_file = output_file_handler.OpenOutputFile(mMeshFilename+".ele");
        (*p_elem_file) << std::scientific;
        
        double horizontal_spacing = width / (double)numNodesAlongWidth;
        double vertical_spacing = (sqrt(3)/2)*horizontal_spacing;
//		double vertical_spacing = (sqrt(3)/2)*length / (double)numNodesAlongLength;

        // This line needed to define ghost nodes later...
        mCryptDepth = (double)numNodesAlongLength * (sqrt(3)/2)* width /(double)numNodesAlongWidth;
        
        if(!mCylindrical)
        {
            // Add an extra node to the width (this is the repeated periodic node)
            numNodesAlongWidth++;
            numNodesAlongWidth = numNodesAlongWidth + 2*ghosts;
        }
        numNodesAlongLength = numNodesAlongLength + 2*ghosts;
        
        unsigned num_nodes            = numNodesAlongWidth*numNodesAlongLength;
        unsigned num_elem_along_width = numNodesAlongWidth-1;
        unsigned num_elem_along_length = numNodesAlongLength-1;
        unsigned num_elem             = 2*num_elem_along_width*num_elem_along_length;
        unsigned num_edges            = 3*num_elem_along_width*num_elem_along_length + num_elem_along_width + num_elem_along_length;
        
        double x0 = -horizontal_spacing*ghosts;
        if(mCylindrical)
        {
            x0 = 0;
        }
        double y0 = -vertical_spacing*ghosts;
        mBottom = -vertical_spacing*ghosts;
        mTop = mBottom + vertical_spacing*(numNodesAlongLength-1);
        
        (*p_node_file) << num_nodes << "\t2\t0\t1" << std::endl;
        unsigned node = 0;
        for (unsigned i = 0; i < numNodesAlongLength; i++)
        {
            for (unsigned j = 0; j < numNodesAlongWidth; j++)
            {
                unsigned boundary = 0;
                if ((i==0) || (i==numNodesAlongLength-1))
                {
                    boundary = 1;
                }
                if (!mCylindrical)
                {
                    if ((j==0) || (j==numNodesAlongWidth-1))
                    {
                        boundary = 1;
                    }
                }
                double x = x0 + horizontal_spacing*((double)j + 0.25*(1.0+ pow(-1,i+1)));
                
                double y = y0 + vertical_spacing*(double)i;
                
                (*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << boundary << std::endl;
                if(i==0)
                {
                    mBottomNodes.push_back(node-1);   
                }
                if(i==numNodesAlongLength-1)
                {
                    mTopNodes.push_back(node-1);   
                }
            }
        }
        p_node_file->close();
        
        out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
        (*p_node_file) << std::scientific;
        
        (*p_elem_file) << num_elem << "\t3\t0" << std::endl;
        (*p_edge_file) << num_edges << "\t3\t0\t1" << std::endl;
        
        unsigned elem = 0;
        unsigned edge = 0;
        for (unsigned i = 0; i < num_elem_along_length; i++)
        {
            for (unsigned j = 0; j < num_elem_along_width; j++)
            {
                unsigned node0 =     i*numNodesAlongWidth + j;
                unsigned node1 =     i*numNodesAlongWidth + j+1;
                unsigned node2 = (i+1)*numNodesAlongWidth + j;
                if (i%2 != 0)
                {
                    node2 = node2 + 1;
                }
                
                (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
                
                unsigned horizontal_edge_is_boundary_edge = 0;
                unsigned vertical_edge_is_boundary_edge = 0;
                if (i==0)
                {
                    horizontal_edge_is_boundary_edge = 1;
                }
                if (j==0 && !mCylindrical)
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
            unsigned node0 = (i+1)*numNodesAlongWidth-1;
            unsigned node1 = (i+2)*numNodesAlongWidth-1;
            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
        }
        
        
        for (unsigned j = 0; j < num_elem_along_width; j++)
        {
            unsigned node0 =  numNodesAlongWidth*(numNodesAlongLength-1) + j;
            unsigned node1 =  numNodesAlongWidth*(numNodesAlongLength-1) + j+1;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
        }
        
        p_elem_file->close();
        p_edge_file->close();
    }
    
    
    void ComputeGhostNodes1()
    {
        assert(mpMesh!=NULL);
        
        for (unsigned i=0; i<mpMesh->GetNumNodes(); i++)
        {
            double x = mpMesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = mpMesh->GetNode(i)->GetPoint().rGetLocation()[1];
            if ((x<0)||(x>mCryptWidth*(1.0+0.5/(double)mNumCellWidth))||(y>mCryptDepth)||(y<-1e-6))
            {
                mGhostNodeIndices.push_back(i);
            }
        }
    }
    
    void ComputeGhostNodes2()
    {
        assert(mpCylindricalMesh!=NULL);
        
        for (unsigned i=0; i<mpCylindricalMesh->GetNumNodes(); i++)
        {
            double x = mpCylindricalMesh->GetNode(i)->GetPoint().rGetLocation()[0];
            double y = mpCylindricalMesh->GetNode(i)->GetPoint().rGetLocation()[1];
            if ((x<0)||(x>=mCryptWidth)||(y>mCryptDepth)||(y<-1e-6))
            {
               mGhostNodeIndices.push_back(i);
            }
        }
    }
    
    
public:

    ~CryptHoneycombMeshGenerator()
    {
        delete mpMesh;
        delete mpCylindricalMesh;
    }

    /**
     * Crypt Honeycomb Mesh Generator
     * 
     * Creates a honeycomb mesh surrounded by ghost nodes
     * 
     * @param numCellWidth  The number of stem cells you want
     * @param numCellDepth  The number of cells you want along crypt axis
     * @param ghostThickness the thickness of ghost nodes surrounding the mesh (defaults to 2)
     * @param cylindrical whether to generate a cylindrically periodic mesh or not (defaults to false)
     * 
     * Note: this class creates a cancer params instance and sets the crypt width and length
     * accordingly in the parameters class.
     */
    CryptHoneycombMeshGenerator(unsigned numCellWidth, unsigned numCellDepth, unsigned ghostThickness = 3, bool cylindrical=false)
    {
        mCylindrical = cylindrical;
        mNumCellWidth = numCellWidth;
        mCryptWidth = numCellWidth*1.0; //*1 because cells are considered to be size one
        mCryptDepth = sqrt(3)*numCellDepth/2;
        mTop = 0.0; // these set in above function
        mBottom = 0.0;  // these set in above function
        
        mMeshFilename = "2D_temporary_crypt_mesh";
        Make2dPeriodicCryptMesh(numCellWidth,numCellDepth,mCryptWidth, ghostThickness);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<2,2> mesh_reader(testoutput_dir+ mMeshFilename);
        
        if (!mCylindrical)
        {
            mpCylindricalMesh = new Cylindrical2dMesh(0.0,mCryptWidth,mTop,mBottom,mTopNodes,mBottomNodes);// to avoid seg fault when closing
            mpMesh = new ConformingTetrahedralMesh<2,2>;
            mpMesh->ConstructFromMeshReader(mesh_reader);
            ComputeGhostNodes1();
        }
        else
        {   
            mpMesh = new ConformingTetrahedralMesh<2,2>;// to avoid seg fault when closing
            mpCylindricalMesh = new Cylindrical2dMesh(0.0,mCryptWidth,mTop,mBottom,mTopNodes,mBottomNodes);
            mpCylindricalMesh->ConstructFromMeshReader(mesh_reader);
            NodeMap map(mpCylindricalMesh->GetNumNodes());
            mpCylindricalMesh->ReMesh(map); // This makes the mesh cylindrical
            ComputeGhostNodes2();
        }
                
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->SetCryptLength(mCryptDepth);
        p_params->SetCryptWidth(mCryptWidth);
    }
    
    /**
     * Crypt Periodic Honeycomb Mesh Generator
     * 
     * Overwritten constructor for a mesh so mesh can be compressed by changing crypt width
     *
     * @param numNodesAlongWidth  The number of stem cells you want
     * @param numNodesAlongLength  The number of cells you want along crypt axis
     * @param width  The width (circumference) of a crypt in adult cell
     * @param ghosts The thickness of ghost nodes to put around the edge (defaults to 3)
     * @param cylindrical whether the mesh should be cylindrically periodic (defaults to false)
     * 
     * Note: this class creates a cancer params instance and sets the crypt width and length
     * accordingly in the parameters class.
     */
    CryptHoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, double width, unsigned ghosts=3, bool cylindrical = false)
    {
        mCylindrical = cylindrical;
        mNumCellWidth = numNodesAlongWidth;
        mCryptWidth = width; //*1 because cells are considered to be size one
        
        mMeshFilename = "2D_temporary_periodic_crypt_mesh";
        Make2dPeriodicCryptMesh(numNodesAlongWidth,numNodesAlongLength,width,ghosts);
        std::string testoutput_dir;
        OutputFileHandler output_file_handler("");
        testoutput_dir = output_file_handler.GetTestOutputDirectory();
        
        TrianglesMeshReader<2,2> mesh_reader(testoutput_dir+ mMeshFilename);
        
        if (!mCylindrical)
        {
            mpCylindricalMesh = new Cylindrical2dMesh(0.0,mCryptWidth,mTop,mBottom,mTopNodes,mBottomNodes);// to avoid seg fault when closing
            mpMesh = new ConformingTetrahedralMesh<2,2>;
            mpMesh->ConstructFromMeshReader(mesh_reader);
            ComputeGhostNodes1();
        }
        else
        {   
            mpMesh = new ConformingTetrahedralMesh<2,2>;// to avoid seg fault when closing
            mpCylindricalMesh = new Cylindrical2dMesh(0.0,mCryptWidth,mTop,mBottom,mTopNodes,mBottomNodes);
            mpCylindricalMesh->ConstructFromMeshReader(mesh_reader);
            NodeMap map(mpCylindricalMesh->GetNumNodes());
            mpCylindricalMesh->ReMesh(map); // This makes the mesh cylindrical
            ComputeGhostNodes2();
        }
                
        CancerParameters* p_params = CancerParameters::Instance();
        p_params->SetCryptLength(mCryptDepth);
        p_params->SetCryptWidth(mCryptWidth);
        
    }
    
    ConformingTetrahedralMesh<2,2>* GetMesh()
    {
        if (mCylindrical)
        {
            EXCEPTION("A cylindrical mesh was created but a normal mesh is being requested.");
        }
        return mpMesh;
    }
    
    Cylindrical2dMesh* GetCylindricalMesh()
    {
        if (!mCylindrical)
        {
            EXCEPTION("A normal mesh was created but a cylindrical mesh is being requested.");
        }
        return mpCylindricalMesh;
    }
    
    std::vector<unsigned> GetGhostNodeIndices()
    {
        return mGhostNodeIndices;
    }
    
    
    
};
#endif /*CRYPTHONEYCOMBMESHGENERATOR_HPP_*/
