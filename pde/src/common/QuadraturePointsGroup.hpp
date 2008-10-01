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

#ifndef QUADRATUREPOINTSGROUP_HPP_
#define QUADRATUREPOINTSGROUP_HPP_

#include "UblasCustomFunctions.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "GaussianQuadratureRule.hpp"
#include "LinearBasisFunction.hpp"
#include <vector>


/** 
 *  A simple class which takes in a mesh and a quadrature rule, and collects
 *  are the quadrature points (in physical space ie several for each element)
 *  together in one data structure, for access.
 */
template<unsigned DIM>
class QuadraturePointsGroup
{
private :    
    /*< The quadrature points in physical space */
    std::vector<c_vector<double,DIM> > data;
    /*< Number of elements in given mesh */
    unsigned mNumElements;
    /*< Number of quad points per element in given rule */
    unsigned mNumQuadPointsPerElement;

public :
    /** 
     *  Constructor takes in a mesh and a rule and computes and stores all
     *  the quad points in physical space
     */
    QuadraturePointsGroup(ConformingTetrahedralMesh<DIM,DIM>& rMesh,
                          GaussianQuadratureRule<DIM>& rQuadRule)
    {
        mNumElements = rMesh.GetNumElements();
        mNumQuadPointsPerElement = rQuadRule.GetNumQuadPoints();
        data.resize(mNumElements*mNumQuadPointsPerElement, zero_vector<double>(DIM));

        // loop over elements
        for(unsigned elem_index=0; elem_index<rMesh.GetNumElements(); elem_index++)        
        {
            Element<DIM,DIM>& r_elem = *(rMesh.GetElement(elem_index));

            c_vector<double, DIM+1> linear_phi;
            for (unsigned quad_index=0; quad_index<rQuadRule.GetNumQuadPoints(); quad_index++)
            {
                const ChastePoint<DIM>& quadrature_point = rQuadRule.rGetQuadPoint(quad_index);
    
                LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);
    
                // interpolate to calculate quad point
                c_vector<double,DIM> X = zero_vector<double>(DIM);
                for(unsigned node_index=0; node_index<DIM+1; node_index++)
                {
                    X += linear_phi(node_index)*rMesh.GetNode( r_elem.GetNodeGlobalIndex(node_index) )->rGetLocation();
                }
                
                // save the quad point
                assert(elem_index<mNumElements);
                assert(quad_index<mNumQuadPointsPerElement);
                data[ elem_index*mNumQuadPointsPerElement + quad_index ] = X;
            }
        }
    }

    /*< Access the stored quad point by element index and quad index in the element */
    c_vector<double,DIM>& Get(unsigned elementIndex, unsigned quadIndex)
    {
        assert(elementIndex<mNumElements);
        assert(quadIndex<mNumQuadPointsPerElement);
        return data[ elementIndex*mNumQuadPointsPerElement + quadIndex ];
    }

    /*< Get the i-th stored quad point */
    c_vector<double,DIM>& Get(unsigned i)
    {
        assert(i < mNumElements*mNumQuadPointsPerElement);
        return data[i];
    }

    /*< Number of elements in the mesh that was given in the constructor */
    unsigned GetNumElements()
    {
        return mNumElements;
    }
    
    /*< Number of quad points per element in the rule that was given in the constructor */
    unsigned GetNumQuadPointsPerElement()
    {
        return mNumQuadPointsPerElement;
    }
    
    /*< Total size, ie total number of quad points, ie num_elem times num_quad_points_per_elem */
    unsigned Size()
    {
        return mNumElements*mNumQuadPointsPerElement;
    }
};


#endif /*QUADRATUREPOINTSGROUP_HPP_*/
