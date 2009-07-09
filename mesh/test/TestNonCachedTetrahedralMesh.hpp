/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTNONCACHEDTETRAHEDRALMESH_HPP_
#define TESTNONCACHEDTETRAHEDRALMESH_HPP_

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "NonCachedTetrahedralMesh.hpp"
#include "OutputFileHandler.hpp"

#include <ctime>

class TestNonCachedTetrahedralMesh : public CxxTest::TestSuite
{
private:

    unsigned GetMemoryUsage()
    {
#ifdef __linux__        
        OutputFileHandler handler("");                
        std::string file_name = handler.GetOutputDirectoryFullPath() + "memusage.tmp";

        std::stringstream ps_command;
        // Option "-o vsize=" makes ps report process virtual size. "=" after "vsize" prints no column header. 
        ps_command << "ps -o vsize= " << getpid() << " > " << file_name;
        EXPECT0(system, ps_command.str().c_str());

        std::ifstream mem_file;
        mem_file.open(file_name.c_str());
        assert(mem_file.is_open());        
        unsigned vsize;
        mem_file >> vsize;        
        
        std::string rm_command = "rm -rf " + file_name;
        EXPECT0(system, rm_command.c_str());
        
        return vsize;        
#elif
        return 0
#endif
    }
    
    
public:

    void TestConstruct3D()
    {
        // TetrahedralMesh with jacobian caching
        unsigned cached_mem_usage;
        time_t cached_start = std::clock();
        {
            TetrahedralMesh<3,3> mesh;
            mesh.ConstructCuboid(20,20,20);
            cached_mem_usage = GetMemoryUsage();            
        }
        time_t cached_finish = std::clock();
        double cached_construction_time = (cached_finish-cached_start)/(double)CLOCKS_PER_SEC;

        
        // No caching
        unsigned non_cached_mem_usage;        
        time_t non_cached_start = std::clock();
        {
            NonCachedTetrahedralMesh<3,3> mesh;
            mesh.ConstructCuboid(20,20,20);
            non_cached_mem_usage = GetMemoryUsage();                                   
        }
        time_t non_cached_finish = std::clock();
        double non_cached_construction_time = (non_cached_finish-non_cached_start)/(double)CLOCKS_PER_SEC;
        
        // Constructing the non cached object should be quicker
        TS_ASSERT( cached_construction_time > non_cached_construction_time );
               
        // compare mem usage
        TS_ASSERT( cached_mem_usage >= non_cached_mem_usage );
    }
    
    void TestSameJacobianData()
    {
        TetrahedralMesh<3,3> cached_mesh;
        cached_mesh.ConstructCuboid(40,40,40);

        NonCachedTetrahedralMesh<3,3> non_cached_mesh;
        non_cached_mesh.ConstructCuboid(40,40,40);

        TS_ASSERT_EQUALS(cached_mesh.GetNumNodes(), non_cached_mesh.GetNumNodes()); 
        TS_ASSERT_EQUALS(cached_mesh.GetNumElements(), non_cached_mesh.GetNumElements()); 
        //TS_ASSERT_EQUALS(cached_mesh.GetNumFaces(), non_cached_mesh.GetNumFaces()); 

        /*
         *  Check element jacobian data is consistent 
         */
        double cached_access_time = 0.0;
        double non_cached_access_time = 0.0;

        for (unsigned element_index = 0; element_index < cached_mesh.GetNumElements(); element_index++)              
        {
            time_t cached_start = std::clock();        
            c_matrix<double, 3, 3> j_cached;
            c_matrix<double, 3, 3> ij_cached;
            double det_cached;
            cached_mesh.GetInverseJacobianForElement(element_index, j_cached,det_cached,ij_cached);
            cached_access_time += (std::clock()-cached_start)/(double)CLOCKS_PER_SEC;    

            time_t non_cached_start = std::clock();            
            c_matrix<double, 3, 3> j_non_cached;
            c_matrix<double, 3, 3> ij_non_cached;
            double det_non_cached;
            non_cached_mesh.GetInverseJacobianForElement(element_index, j_non_cached,det_non_cached,ij_non_cached);
            non_cached_access_time += (std::clock()-non_cached_start)/(double)CLOCKS_PER_SEC;
            
            TS_ASSERT_EQUALS(det_cached, det_non_cached);
            
            for (unsigned row=0; row<2; row++)
            {
                for (unsigned col=0; col<2; col++)
                {           
                    TS_ASSERT_EQUALS(j_cached(row,col), j_non_cached(row,col));
                    TS_ASSERT_EQUALS(ij_cached(row,col), ij_non_cached(row,col));
                }
            }
        }

        // Retrieving the cached jacobians should be quicker
        TS_ASSERT(cached_access_time < non_cached_access_time);

        /*
         *  Check boundary element jacobian data is consistent 
         */
        for (unsigned boundary_element_index = 0; boundary_element_index < cached_mesh.GetNumBoundaryElements(); boundary_element_index++)              
        {
            c_vector<double, 3> wd_cached;
            double det_cached;
            cached_mesh.GetWeightedDirectionForBoundaryElement(boundary_element_index, wd_cached,det_cached);

            c_vector<double, 3> wd_non_cached;
            double det_non_cached;
            non_cached_mesh.GetWeightedDirectionForBoundaryElement(boundary_element_index, wd_non_cached,det_non_cached);
            
            TS_ASSERT_EQUALS(det_cached, det_non_cached);
            
            for (unsigned row=0; row<2; row++)
            {
                TS_ASSERT_EQUALS(wd_cached(row), wd_non_cached(row));
            }
        }
    }
    
    void TestExceptions()
    {
        NonCachedTetrahedralMesh<3,3> non_cached_mesh;
        non_cached_mesh.ConstructCuboid(1,1,1);

        c_matrix<double, 3, 3> jacobian;
        double det_jacobian;                
        TS_ASSERT_THROWS_ANYTHING(non_cached_mesh.GetJacobianForElement(0u, jacobian, det_jacobian));

        c_vector<double, 3> direction;
        double det_direction;                        
        TS_ASSERT_THROWS_ANYTHING(non_cached_mesh.GetWeightedDirectionForElement(0u, direction, det_direction));                
    }
    
};

#endif /*TESTNONCACHEDTETRAHEDRALMESH_HPP_*/
