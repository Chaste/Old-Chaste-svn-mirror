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

#include "VertexCryptBoundaryForce.hpp"

template<unsigned DIM>
VertexCryptBoundaryForce<DIM>::VertexCryptBoundaryForce()
   : AbstractForce<DIM>()
{
}

template<unsigned DIM>
VertexCryptBoundaryForce<DIM>::~VertexCryptBoundaryForce()
{
}

template<unsigned DIM>
void VertexCryptBoundaryForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                         AbstractTissue<DIM>& rTissue)
{
    // Helper variable that is a static cast of the tissue
    VertexBasedTissue<DIM>* p_tissue = static_cast<VertexBasedTissue<DIM>*>(&rTissue);

    /*
     * The boundary force is taken to be a quadratic function (F = (y-a)^2 - b), which passes 
     * through (0,0), but which is fixed at zero for y-values greater than or equal to the 
     * cutoff point. The y-values correspond to the height up the crypt.
     */   

    double cutoff_point = 0.2;    // The y-value at which the boundary force no longer has an effect
    double quadratic_shift = 0.5*cutoff_point;   // Distance to shift the quadratic curve (i.e. 'a')

    // Iterate over vertices in the tissue
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = p_tissue->rGetMesh().GetNodeIteratorBegin();
         node_iter != p_tissue->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // Compute the boundary force on this node
        c_vector<double, DIM> boundary_force = zero_vector<double>(DIM);

        double y_value = node_iter->rGetLocation()[1]; // y-coordinate of node

        if (y_value < cutoff_point)
        {
            boundary_force[1] = 150*(pow((y_value - quadratic_shift),2) - pow(quadratic_shift,2));                   
        }

        rForces[node_iter->GetIndex()] += boundary_force;
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexCryptBoundaryForce<1>;
template class VertexCryptBoundaryForce<2>;
template class VertexCryptBoundaryForce<3>;
