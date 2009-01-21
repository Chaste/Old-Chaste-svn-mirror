/*

Copyright (C) University of Oxford, 2008

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
#ifndef PAPILLARYFIBRECALCULATOR_HPP_
#define PAPILLARYFIBRECALCULATOR_HPP_

#include "TetrahedralMesh.hpp"

// Always 3D
class PapillaryFibreCalculator
{
private:

public:
    c_vector<double, 3> GetRadiusVectorForOneElement(TetrahedralMesh<3,3>& rMesh, unsigned elementIndex)
    {
        c_vector<double, 3> centroid = (rMesh.GetElement(elementIndex))->CalculateCentroid();
        // Loops over all papillary face nodes
        c_vector<double,3> coordinates;
        
        double nearest_r_squared=DBL_MAX;
        unsigned nearest_face_node = 0;
          
        TetrahedralMesh<3,3>::BoundaryNodeIterator bound_node_iter = rMesh.GetBoundaryNodeIteratorBegin();
        while (bound_node_iter != rMesh.GetBoundaryNodeIteratorEnd())
        {
            unsigned bound_node_index =  (*bound_node_iter)->GetIndex();           
            coordinates=rMesh.GetNode(bound_node_index)->rGetLocation();

            // Calculates the distance between the papillary face node and the centroid
            double r_squared =  norm_2(centroid-coordinates);
            // Checks to see if it is the smallest so far - if it is, update the current smallest distance
            if (r_squared < nearest_r_squared)
            {
                nearest_r_squared = r_squared;
                nearest_face_node = bound_node_index;
            }
            ++bound_node_iter;
        }
           
        coordinates = rMesh.GetNode(nearest_face_node)->rGetLocation();
        c_vector<double,3> radial_vector = coordinates-centroid;
        return radial_vector;
    }
    
    std::vector< c_vector<double, 3> > GetRadiusVectors(TetrahedralMesh<3,3>& rMesh)
    {
       std::vector < c_vector<double, 3> > radial_vectors;
       
       // Loops over all elements finding radius vector
        TetrahedralMesh<3,3>::ElementIterator iter = rMesh.GetElementIteratorBegin();
        while (iter != rMesh.GetElementIteratorEnd())
        {
            unsigned element_index = (*iter)->GetIndex();
            c_vector<double, 3> vector_for_this_element = GetRadiusVectorForOneElement(rMesh,element_index);

            radial_vectors.push_back(vector_for_this_element);
            
            ++iter;
        }
           
        return radial_vectors;
    }
    
    std::vector< c_matrix<double,3,3> > ConstructStructureTensors(std::vector< c_vector<double, 3> > radiusVectors)
    {
        std::vector< c_matrix<double,3,3> > tensor_i;
        
        for(unsigned i=0;i<radiusVectors.size();i++)
        {
            tensor_i.push_back(outer_prod(radiusVectors[i],radiusVectors[i]));
         
        }
        
        return tensor_i;
     
        
    }
};

#endif /*PAPILLARYFIBRECALCULATOR_HPP_*/
