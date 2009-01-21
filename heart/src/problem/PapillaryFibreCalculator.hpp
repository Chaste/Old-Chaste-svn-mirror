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
// Allow the test class to use the private functions.
friend class TestPapillaryFibreCalculator;
    
private:
//    TetrahedralMesh<3,3>& mrMesh;
//    std::vector< c_vector<double, 3> > mRadiusVectors; 
//    std::vector< c_matrix<double,3,3> > mStructureTensors;    
//    std::vector< c_matrix<double,3,3> > mSmoothedStructureTensors;    

   /**
     * This method calculates the vector from the centroid of an element to all of
     * the boundary nodes. It returns the shortest of the vectors. 
     * 
     * @param rMesh  A reference to the mesh (must be a TetrahedralMesh<3,3>)
     * @param elementIndex  The index of the element we are calculating radial vectors for
     * @return The shortest radial vector
     */
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
    
    /**
     * This method calls GetRadiusVectorForOneElement() for each of the elements and 
     * returns a radial vector for each element of the mesh.
     * 
     * @param rMesh  The mesh
     * @return  A vector (of length # elements) of radial vectors 
     */
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
    
    
    /**
     * This generates structure tensors from the radial vectors by taking 
     * 
     * T = r.r'
     * 
     * @param radiusVectors  A std::vector (of length # elements) that contains the vector for each element that points to the nearest boundary node
     * @return  The structure tensor
     */
    std::vector< c_matrix<double,3,3> > ConstructStructureTensors(std::vector< c_vector<double, 3> > radiusVectors)
    {
        std::vector< c_matrix<double,3,3> > tensor_i;
        
        for(unsigned i=0;i<radiusVectors.size();i++)
        {
            tensor_i.push_back(outer_prod(radiusVectors[i],radiusVectors[i]));
         
        }
        
        return tensor_i;
    }
    
    /**
     * Smoothes the structure tensor components for each papillary element by looping 
     * over all other papillary elements, calculating
     * distance geometric distance between the two elements; 
     * if it is within a certain limit, include this in the Gaussian kernel
     * 
     * @param tensor  The 'rough' tensor for each element
     * @param rMesh  The mesh
     * @return The smoothed tensor for each element
     */
    void SmoothStructureTensors(const std::vector< c_matrix<double,3,3> >& rPreSmoothedTensors, 
                                TetrahedralMesh<3,3>& rMesh, 
                                std::vector<c_matrix<double,3,3> >& rSmoothedTensors)
    {
        assert(rSmoothedTensors.size()==rPreSmoothedTensors.size());

        double g_factor = 0;
        double sigma = 0.05; //cm
        double g_factor_sum = 0;
        double r_max = 0.1; //cm

        for(TetrahedralMesh<3,3>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
            elem_iter != rMesh.GetElementIteratorEnd();
            ++elem_iter)
        {
            rSmoothedTensors[ (*elem_iter)->GetIndex()] = zero_matrix<double>(3,3);
            
            c_vector<double, 3> centroid = (*elem_iter)->CalculateCentroid();  
            g_factor_sum = 0;
            
            for(TetrahedralMesh<3,3>::ElementIterator iter_2 = rMesh.GetElementIteratorBegin();
                iter_2 != rMesh.GetElementIteratorEnd();
                ++iter_2)
            {
                c_vector<double, 3> centroid_2 = (*iter_2)->CalculateCentroid();
                double r = norm_2(centroid-centroid_2);             
                if (r < r_max)
                {
                    g_factor = exp(-r/(2*sigma*sigma));
                
                    g_factor_sum += g_factor;
                
                    rSmoothedTensors[ (*elem_iter)->GetIndex()] += g_factor*rPreSmoothedTensors[ (*iter_2)->GetIndex()];
                }
            }      
        
            rSmoothedTensors[ (*elem_iter)->GetIndex()] /= g_factor_sum;
        }
    }
    
public:

    /// \todo put a constructor in here than can assign the rMesh to a member variable.
 
    /**
     * 
     * 
     */
     std::vector<c_vector<double,3> > CalculateFibreOrientations(TetrahedralMesh<3,3>& rMesh)
     {
        std::vector< c_vector<double, 3> > radial_vectors = GetRadiusVectors(rMesh);
        
        std::vector< c_matrix<double,3,3> > tensors = ConstructStructureTensors(radial_vectors);
        
        std::vector< c_matrix<double,3,3> > smoothed_tensors(tensors.size());
        SmoothStructureTensors(tensors, rMesh, smoothed_tensors);

        // Calculate eigenvalues etc...
        std::vector<c_vector<double,3> > fibre_orientations(tensors.size());
        for(unsigned i=0; i<fibre_orientations.size(); i++)
        {
            fibre_orientations[i] = CalculateEigenvectorForSmallestEigenvalue(smoothed_tensors[i]);
        }

        return fibre_orientations;
     }
};

#endif /*PAPILLARYFIBRECALCULATOR_HPP_*/

