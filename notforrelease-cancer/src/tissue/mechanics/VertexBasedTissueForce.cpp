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
#include "VertexBasedTissueForce.hpp"
#include "VertexBasedTissue.hpp"

template<unsigned DIM>
VertexBasedTissueForce<DIM>::VertexBasedTissueForce(double tissueCellTargetArea)
   : AbstractForce<DIM>(),
     mTissueCellTargetArea(tissueCellTargetArea),
     mDeformationEnergyParameter(1.0),
     mMembraneSurfaceEnergyParameter(1.0),
     mCellCellAdhesionEnergyParameter(1.0)
{
}

template<unsigned DIM>
void VertexBasedTissueForce<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void VertexBasedTissueForce<DIM>::SetMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void VertexBasedTissueForce<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
VertexBasedTissueForce<DIM>::~VertexBasedTissueForce()
{
}

template<unsigned DIM>
void VertexBasedTissueForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                  AbstractTissue<DIM>& rTissue)
{
    // Helper variable that is a static cast of the tissue
    VertexBasedTissue<DIM>* p_tissue = static_cast<VertexBasedTissue<DIM>*>(&rTissue);
    
    // Iterate over vertices in the tissue
    for (unsigned node_index=0; node_index<p_tissue->GetNumNodes(); node_index++)
    {
        // Compute the force on this node
        
        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the Tissue, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all TissueCells in 
         * the tissue. The free energy of each TissueCell is comprised of three
         * parts - a cell deformation energy, a membrant surface-tension energy
         * and a cell-cell adhesion energy.
         * 
         * Note that since the movement of this Node only affects the free energy 
         * of the three TissueCells containing it, we can just consider the 
         * contributions to the free energy gradient from each of these three 
         * TissueCells.
         */

        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
             
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_tissue->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter=containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Find the local index of this node in this element
            unsigned local_index = 0;
            while (p_tissue->GetElement(*iter)->GetNodeGlobalIndex(local_index) != node_index)
            {
                local_index++;
            }
            
            // Compute the area of the element and its gradient at this node
            double element_area = p_tissue->GetElement(*iter)->GetArea();
            c_vector<double, DIM> area_gradient = p_tissue->GetElement(*iter)->GetGradientOfAreaAtNode(local_index);
            
            deformation_contribution += 2*mDeformationEnergyParameter*
                                        (element_area - mTissueCellTargetArea)*area_gradient;
            
            /// \todo Note that in the line above, we'll probably need to add extra code 
            ///       for a variable mTissueCellTargetArea when cell birth is implemented 
            ///       (see #852)

            /// \todo Add contributions from gradients of membrane surface tension energy and 
            ///       cell-cell adhesion energy (see #861)            
        }
         
        c_vector<double, DIM> force_on_node = zero_vector<double>(DIM); /// \todo code up actual forces (see #861)

        rForces[node_index] += force_on_node;
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexBasedTissueForce<1>;
template class VertexBasedTissueForce<2>;
template class VertexBasedTissueForce<3>;
