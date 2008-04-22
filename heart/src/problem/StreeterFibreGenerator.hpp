#ifndef STREETERFIBREGENERATOR_HPP_
#define STREETERFIBREGENERATOR_HPP_

#include "ConformingTetrahedralMesh.hpp"
#include "DistanceMapCalculator.hpp"
#include "MemfemMeshReader.hpp"
#include <vector>
#include <set>
#include <fstream>


template<unsigned SPACE_DIM>
class StreeterFibreGenerator
{
private:
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh;
    unsigned mNumNodes, mNumElements;
    
    DistanceMapCalculator<SPACE_DIM>* mpDistanceCalculator;
    
    std::string mEpiFile, mRVFile, mLVFile;
    bool mFilesSet;
    
    std::vector<unsigned> mEpiSurface, mRVSurface, mLVSurface;
    
    std::vector<double> mDistMapEpicardium, mDistMapRightVentricle, mDistMapLeftVentricle;

    enum RegionType_
    {
        LEFT_VENTRICLE_WALL,
        RIGHT_VENTRICLE_WALL,
        LEFT_SEPTUM,
        RIGHT_SEPTUM,
        UNKNOWN
    };
    
    // Area of the septum considered to belong to the each ventricle (relative to 1)
    static const double LEFT_SEPTUM_SIZE = 2.0/3.0;
    static const double RIGHT_SEPTUM_SIZE = 1.0/3.0;

    
    inline RegionType_ GetHeartRegion (unsigned nodeIndex) const
    {      
        
        if (mDistMapRightVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
            mDistMapRightVentricle[nodeIndex] >= mDistMapLeftVentricle[nodeIndex])
        {
            return LEFT_VENTRICLE_WALL;
        }
        
        if (mDistMapLeftVentricle[nodeIndex] >= mDistMapEpicardium[nodeIndex] &&
            mDistMapLeftVentricle[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
        {
            return RIGHT_VENTRICLE_WALL;
        }
        
        if (mDistMapEpicardium[nodeIndex] >= mDistMapLeftVentricle[nodeIndex] &&
            mDistMapEpicardium[nodeIndex] >= mDistMapRightVentricle[nodeIndex])
        {
            if (mDistMapLeftVentricle[nodeIndex] 
                < LEFT_SEPTUM_SIZE*(mDistMapLeftVentricle[nodeIndex] + mDistMapRightVentricle[nodeIndex]))
            {
                return LEFT_SEPTUM;
            }
            else
            {
                return RIGHT_SEPTUM;
            }        
        }
        
        return UNKNOWN;        
    }

    
    inline double GetAveragedThickness(const unsigned nodeIndex, const std::vector<double>& wallThickness) const
    {
        // Initialise the average with the value corresponding to the current node
        double average = wallThickness[nodeIndex];
        unsigned nodes_visited = 1;
        
        // Use a set to store visited nodes
        std::set<unsigned> visited_nodes;
        visited_nodes.insert(nodeIndex);
        
        Node<SPACE_DIM>* p_current_node = mrMesh.GetNode(nodeIndex);
        
        /// \todo The following nested loops appear in DistanceMapCalculator as well, refactor it. Idea: create an iterator over the neighbour nodes in class Node
        
        // Loop over the elements containing the given node
        for(typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
            element_iterator != p_current_node->ContainingElementsEnd();
            ++element_iterator)
        {
            // Get a pointer to the container element
            Element<SPACE_DIM,SPACE_DIM>* p_containing_element = mrMesh.GetElement(*element_iterator);

           // Loop over the nodes of the element
           for(unsigned node_local_index=0; 
               node_local_index<p_containing_element->GetNumNodes();
               node_local_index++)
           {
                Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
                unsigned neighbour_node_index = p_neighbour_node->GetIndex();
            
                // Check if the neighbour node has already been visited
                if (visited_nodes.find(neighbour_node_index) == visited_nodes.end())
                {
                    average += wallThickness[neighbour_node_index];
                    visited_nodes.insert(neighbour_node_index);
                    nodes_visited++;
                }
           }
        }
        
        return average/nodes_visited;
    }

    
    inline void ProcessLine(const std::string& line, std::set<unsigned>& surfaceElements) const
    {
        unsigned num_nodes = 0;
        std::stringstream line_stream(line);                
                    
        while (!line_stream.eof())
        {
            unsigned item;
            line_stream >> item;
            // Shift the nodes, since we are assuming MEMFEM format (numbered from 1 on)
            surfaceElements.insert(item-1);
            
            num_nodes++;
        }        

        if (num_nodes != SPACE_DIM)
        {
            EXCEPTION("Wrong file format");    
        }
        
    }


    void GetNodesAtSurface(const std::string& surfaceFile, std::vector<unsigned>& surfaceVector) const
    {
        // Open the file defining the surface
        std::ifstream file_stream;
        file_stream.open(surfaceFile.c_str());
        if (!file_stream.is_open())
        {
            EXCEPTION("Wrong surface definition file name.");    
        }        
        
        // Temporal storage for the nodes, helps discarting repeated values
        std::set<unsigned> surface_elements;
        
        // Loop over all the triangles and add node indexes to the set
        std::string line;        
        getline(file_stream, line);
        do
        {
            ProcessLine(line, surface_elements);
            
            getline(file_stream, line);
        }
        while(!file_stream.eof());        
        
        // Make vector big enough
        surfaceVector.reserve(surface_elements.size());
        
        // Copy the node indexes from the set to the vector
        for(std::set<unsigned>::iterator element_it=surface_elements.begin();
            element_it != surface_elements.end();
            element_it++)
        {
            surfaceVector.push_back(*element_it);
        }
        
        file_stream.close();
    }


    double GetFibreMaxAngle(const c_vector<RegionType_, SPACE_DIM+1>& nodesRegion) const
    {
        unsigned lv=0, rv=0;
        
        for (unsigned index=0; index<SPACE_DIM+1; index++)
        {
            switch (nodesRegion[index])
            {
                case LEFT_VENTRICLE_WALL:
                case LEFT_SEPTUM:
                    lv++;
                    break;
                    
                case RIGHT_VENTRICLE_WALL:
                case RIGHT_SEPTUM:
                    rv++;
                    break;
                    
                case UNKNOWN:
                    break;                                                                        
            }
        }
        
        // If most of the nodes are in the right ventricle
        if (rv>lv)
        {
            return M_PI/4.0; 
        }

        // Anywhere else
        return M_PI/3.0;        
    }
    


public:   
    StreeterFibreGenerator(ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh):
        mrMesh(rMesh),
        mFilesSet(false)
    {
        mNumNodes = mrMesh.GetNumNodes();
        mNumElements = mrMesh.GetNumElements();         
        mpDistanceCalculator = new DistanceMapCalculator<SPACE_DIM>(mrMesh);                 
    }
    
    ~StreeterFibreGenerator()
    {
        delete mpDistanceCalculator;
    }
    
    /**
     * Sets the files defining the diferent surfaces of the mesh. File format: list of triangles
     * 
     * @param epicardiumFile Epicardium surface
     * @param rightVentricleFile Right Ventricle surface
     * @param rightVentricleFile Left Ventricle surface
     * 
     */
    void SetSurfaceFiles(std::string epicardiumFile,
                         std::string rightVentricleFile,
                         std::string leftVentricleFile)
    {
        mEpiFile = epicardiumFile;        
        mRVFile = rightVentricleFile;
        mLVFile = leftVentricleFile;
        
        mFilesSet = true;
    }
    
    /**
     * Generates an orthotropic fibre orientation model of the ventricular mesh provided. Assumes that the base-apex axis is x. Based on Streeter 1969 and Potse 2006
     * 
     * File format: The first line indicates the number of elements. Each of the following lines contain SPACE_DIM vectors of SPACE_DIM elements for the 
     * direction of the myofibre, the transverse to it in the plane of the myocite laminae and the normal to this laminae.
     * 
     * @param outputDirectory Output directory 
     * @param fibreOrientationFile Output file
     * @param logInfo Tells the method to output extra debug info. To be eliminated once it's fully tested
     * 
     */
    void GenerateOrthotropicFibreOrientation(std::string outputDirectory, std::string fibreOrientationFile, bool logInfo=false)
    {
        if (!mFilesSet)
        {
            EXCEPTION("Files defining the heart surfaces not set");
        }
        
        // Open files
        OutputFileHandler handler(outputDirectory, false);
        out_stream p_fibre_file = handler.OpenOutputFile(fibreOrientationFile);
        out_stream p_regions_file, p_thickness_file, p_ave_thickness_file, p_grad_thickness_file;       
        
        if (logInfo)
        {
            p_regions_file  = handler.OpenOutputFile("node_regions.data");
            p_thickness_file = handler.OpenOutputFile("wall_thickness.data");
            p_ave_thickness_file = handler.OpenOutputFile("averaged_thickness.data");
            p_grad_thickness_file = handler.OpenOutputFile("grad_thickness.data");                       
        }
        
        // First line of the fibre file: number of elements of the mesh
        *p_fibre_file << mNumElements << std::endl;
        
        // Get nodes defining each surface
        GetNodesAtSurface(mEpiFile, mEpiSurface);        
        GetNodesAtSurface(mRVFile, mRVSurface);        
        GetNodesAtSurface(mLVFile, mLVSurface);
        
        // Compute the distance map of each surface        
        mpDistanceCalculator->ComputeDistanceMap(mEpiSurface, mDistMapEpicardium);
        mpDistanceCalculator->ComputeDistanceMap(mRVSurface, mDistMapRightVentricle);
        mpDistanceCalculator->ComputeDistanceMap(mLVSurface, mDistMapLeftVentricle);
              
        /* 
         * Compute wall thickness parameter, e
         */
        std::vector<double> wall_thickness(mNumNodes);
        for (unsigned node_index=0; node_index<mNumNodes; node_index++)
        {
            double dist_epi, dist_endo;
                        
            RegionType_ node_region = GetHeartRegion(node_index);
                        
            switch(node_region)
            {
                case LEFT_VENTRICLE_WALL:
                    dist_epi = mDistMapEpicardium[node_index];
                    dist_endo = mDistMapLeftVentricle[node_index];
                    break;
                    
                case RIGHT_VENTRICLE_WALL:
                    dist_epi = mDistMapEpicardium[node_index];
                    dist_endo = mDistMapRightVentricle[node_index];                
                    break;
                    
                case LEFT_SEPTUM:
                    dist_epi = mDistMapRightVentricle[node_index];
                    dist_endo = mDistMapLeftVentricle[node_index];
                    break;
                    
                case RIGHT_SEPTUM:
                    dist_epi = mDistMapLeftVentricle[node_index];
                    dist_endo = mDistMapRightVentricle[node_index];
                    break;
                
                case UNKNOWN:
                    #define COVERAGE_IGNORE
                    std::cerr << "Wrong distances node: " << node_index << "\t" 
                              << "Epi " << mDistMapEpicardium[node_index] << "\t"
                              << "RV " << mDistMapRightVentricle[node_index] << "\t"
                              << "LV " << mDistMapLeftVentricle[node_index]
                              << std::endl;

                    // Make wall_thickness=0 as in Martin's code                             
                    dist_epi = 1;
                    dist_endo = 0;
                    break;
                    #undef COVERAGE_IGNORE                                                                        
            }
                        
            wall_thickness[node_index] = dist_endo / (dist_endo + dist_epi);           
            
            if (isnan(wall_thickness[node_index]))
            {                
                /*
                 *  A node on both epicardium and lv (or rv) surfaces has wall thickness 0/0.
                 *  By setting its value to 0 we consider it contained only on the lv (or rv) surface.
                 */ 
                wall_thickness[node_index] = 0;                
            }
            
            if (logInfo)
            {
                *p_regions_file << node_region*100 << "\n";
                *p_thickness_file << wall_thickness[node_index] << "\n";                    
            }                                                
        }
        
        /*
         *  For each node, average its value of e with the values of all the neighbours
         */
        std::vector<double> averaged_wall_thickness(mNumNodes);      
        for (unsigned node_index=0; node_index<mNumNodes; node_index++)
        {
            averaged_wall_thickness[node_index] = GetAveragedThickness(node_index, wall_thickness);

            if (logInfo)
            {
                *p_ave_thickness_file << averaged_wall_thickness[node_index] << "\n";                    
            }                                                
            
        }       
        
        /*
         *  Compute the gradient of the averaged wall thickness at the centroid of each tetrahedral
         */
        c_vector<double,SPACE_DIM> grad_ave_wall_thickness;
        
        for (unsigned element_index=0; element_index<mNumElements; element_index++)
        {
            Element<SPACE_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(element_index);

            /*
             *  The gradient of the averaged thickness at the element is:
             *                                                 
             *     grad_ave_wall_thickness[element_index] = ave' * BF * inv(J)
             * 
             *  being : ave, averaged thickness values of the nodes defining the element
             *          J,   the Jacobian of the element as defined in class Element.             
             *                               (-1 -1 -1)
             *          BF,  basis functions ( 1  0  0)
             *                               ( 0  1  0)
             *                               ( 0  0  1)
             * 
             *  Defined as u in Streeter paper
             */
            
            c_vector<double, SPACE_DIM+1> elem_nodes_ave_thickness;
            double element_averaged_thickness = 0.0;
            c_vector<RegionType_, SPACE_DIM+1> elem_nodes_region;
            
            for (unsigned local_node_index=0; local_node_index<SPACE_DIM+1; local_node_index++)
            {
                // Get node's global index
                unsigned global_node_index = p_element->GetNode(local_node_index)->GetIndex();
                
                elem_nodes_ave_thickness[local_node_index] = averaged_wall_thickness[global_node_index];
                elem_nodes_region[local_node_index] = GetHeartRegion(global_node_index);
                
                // Calculate wall thickness averaged value for the element               
                element_averaged_thickness +=  wall_thickness[global_node_index];               
            }
            
            element_averaged_thickness /= SPACE_DIM+1;                
            
            /// \todo basis_functions matrix is 3D specific, work out the generic expression
            assert (SPACE_DIM==3);
            
            c_matrix<double, SPACE_DIM+1, SPACE_DIM> basis_functions( zero_matrix<double>(4u,3u) );
            basis_functions(0,0) = basis_functions(0,1) = basis_functions(0,2) = -1.0;
            basis_functions(1,0) = basis_functions(2,1) = basis_functions(3,2) =  1.0;
                        
            c_matrix<double, SPACE_DIM+1, SPACE_DIM> temp;
            noalias(temp) = prod (basis_functions, *p_element->GetInverseJacobian() );
            noalias(grad_ave_wall_thickness) = prod(elem_nodes_ave_thickness, temp);
            
            if (logInfo)
            {
                *p_grad_thickness_file << grad_ave_wall_thickness[0] << " " << grad_ave_wall_thickness[1] << " " << grad_ave_wall_thickness[2] << std::endl;                    
            }                                                
                          
            /*
             *   Normal to the gradient (v in Streeter paper) which is then the fibre direction (computed
             *  as the cross product with the x-axis)
             */
             /// \todo Assuming the base-apex axis is x
            c_vector<double, SPACE_DIM> fibre_direction = VectorProduct(grad_ave_wall_thickness, Create_c_vector(1.0, 0.0, 0.0));
            
            /*
             *  Longitude direction (w in Streeter paper)  
             */
            c_vector<double, SPACE_DIM> longitude_direction = VectorProduct(grad_ave_wall_thickness, fibre_direction);           
                        
            /*
             *  Compute fibre to v angle: alpha = R*(1-2e)^3
             * 
             *    R is the maximum angle between the fibre and the v axis (heart region dependant)
             *    (1 - 2e)^3 scales it by a value in [-1, 1] defining the rotation of the fibre based 
             *       on the position in the wall   
             */
            double alpha = GetFibreMaxAngle(elem_nodes_region) * pow( (1 - 2*element_averaged_thickness), 3 );
            
            /*
             *  Apply alpha rotation to fibre_direction vector. Solve the system
             * 
             *   ( v(1) v(2) v(3) ) * ( f(1) )   (     ||v||*cos(alpha)    )
             *   ( u(1) u(2) u(3) )   ( f(2) ) = (           0           )
             *   ( w(1) w(2) w(3) )   ( f(3) )   ( ||w||*cos(PI/2 - alpha) )
             * 
             */
            c_matrix<double, SPACE_DIM, SPACE_DIM> fibre_space_matrix;           
            for(unsigned column_index=0; column_index<SPACE_DIM; column_index++)
            {
                fibre_space_matrix(0,column_index) = fibre_direction(column_index);
                fibre_space_matrix(1,column_index) = grad_ave_wall_thickness(column_index);
                fibre_space_matrix(2,column_index) = longitude_direction(column_index);                   
            }
            
            c_matrix<double, SPACE_DIM, SPACE_DIM> inv_fibre_space_matrix = Inverse(fibre_space_matrix);
            
            double norm_v = sqrt(   fibre_direction[0] * fibre_direction[0]
                                  + fibre_direction[1] * fibre_direction[1]
                                  + fibre_direction[2] * fibre_direction[2] );
                                  
            double norm_w = sqrt(   longitude_direction[0] * longitude_direction[0]
                                  + longitude_direction[1] * longitude_direction[1]
                                  + longitude_direction[2] * longitude_direction[2] );                                 
            
            c_vector<double, SPACE_DIM> rotation_rhs(Create_c_vector( norm_v*cos(alpha), 
                                                                      0.0, 
                                                                      norm_w*cos( M_PI/2 - alpha )));

            /// \todo Use LU factorisation to solve this system
            c_vector<double, SPACE_DIM> rotated_fibre_direction = prod( inv_fibre_space_matrix, rotation_rhs);


            /*
             *  Apply alpha rotation to longitude_direction vector. Solve the system
             * 
             *   ( v(1) v(2) v(3) ) * ( l(1) )   ( ||v||*cos(PI/2 + alpha) )
             *   ( u(1) u(2) u(3) )   ( l(2) ) = (         0        )
             *   ( w(1) w(2) w(3) )   ( l(3) )   ( ||w||*cos(alpha) )
             * 
             */           
            rotation_rhs = Create_c_vector( norm_v*cos( M_PI/2 + alpha), 
                                            0.0, 
                                            norm_w*cos( alpha ));

            /// \todo Use LU factorisation to solve this system
            c_vector<double, SPACE_DIM> rotated_longitude_direction = prod( inv_fibre_space_matrix, rotation_rhs);

            assert(rotated_longitude_direction[0] != longitude_direction[0] &&
                   rotated_longitude_direction[1] != longitude_direction[1] &&
                   rotated_longitude_direction[2] != longitude_direction[2] );

            // Output the direction of the myofibre, the transverse to it in the plane of the myocite laminae and the normal to this laminae (in that order)
            *p_fibre_file << rotated_fibre_direction[0]     << " " << rotated_fibre_direction[1]     << " "  << rotated_fibre_direction[2]     << " "
                          << rotated_longitude_direction[0] << " " << rotated_longitude_direction[1] << " "  << rotated_longitude_direction[2] << " "
                          << grad_ave_wall_thickness[0]     << " " << grad_ave_wall_thickness[1]     << " "  << grad_ave_wall_thickness[2]     << " "
                          << std::endl;
            
        }               
        
        p_fibre_file->close();
        
        if (logInfo)
        {
            p_regions_file->close();
            p_thickness_file->close();
            p_ave_thickness_file->close();
            p_grad_thickness_file->close();
        }
    }
    
};

#endif /*STREETERFIBREGENERATOR_HPP_*/
