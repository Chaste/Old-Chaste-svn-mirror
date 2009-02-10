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
#include "VertexBasedTissueForce.hpp"


template<unsigned DIM>
VertexBasedTissueForce<DIM>::VertexBasedTissueForce()
   : AbstractForce<DIM>()
{
}


template<unsigned DIM>
VertexBasedTissueForce<DIM>::~VertexBasedTissueForce()
{
}


template<unsigned DIM>
void VertexBasedTissueForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                       AbstractTissue<DIM>& rTissue)
{
    // Helper instance of CancerParameters
    CancerParameters *p_params = CancerParameters::Instance();

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
         * parts - a cell deformation energy, a membrane surface tension energy
         * and a cell-cell adhesion energy.
         * 
         * Note that since the movement of this Node only affects the free energy 
         * of the three TissueCells containing it, we can just consider the 
         * contributions to the free energy gradient from each of these three 
         * TissueCells.
         */

        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> cell_cell_adhesion_contribution = zero_vector<double>(DIM);
             
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_tissue->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter=containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element and its index
            VertexElement<DIM, DIM>* p_element = p_tissue->GetElement(*iter);
            unsigned element_index = p_element->GetIndex();

            // Find the local index of this node in this element
            unsigned local_index = 0;
            while (p_element->GetNodeGlobalIndex(local_index) != node_index)
            {
                local_index++;
            }


            /******** Start of deformation force calculation ********/

            // Compute the area of this element and its gradient at this node
            double element_area = p_element->GetArea();
            c_vector<double, DIM> element_area_gradient = p_element->GetAreaGradientAtNode(local_index);

            // Get the target area of the cell
            double cell_target_area = p_tissue->GetTargetAreaOfCell(p_tissue->rGetCellUsingLocationIndex(element_index));

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            deformation_contribution += -2*p_params->GetDeformationEnergyParameter()*(element_area - cell_target_area)*element_area_gradient;

            /******** End of membrane force calculation *************/


            /******** Start of membrane force calculation ***********/

            // Compute the perimeter of the element and its gradient at this node
            double element_perimeter = p_element->GetPerimeter();
            c_vector<double, DIM> element_perimeter_gradient = p_element->GetPerimeterGradientAtNode(local_index);

            // Get the target perimeter of the cell
            double cell_target_perimeter = 2*sqrt(M_PI*cell_target_area);

            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            membrane_surface_tension_contribution += -2*p_params->GetMembraneSurfaceEnergyParameter()*(element_perimeter - cell_target_perimeter)*element_perimeter_gradient;
            /******** End of deformation force calculation **********/


            /******** Start of adhesion force calculation ***********/
            /// \todo code this up (#861)
            /******** End of adhesion force calculation *************/
        }
         
        c_vector<double, DIM> force_on_node = deformation_contribution + 
                                              membrane_surface_tension_contribution +
                                              cell_cell_adhesion_contribution;

        rForces[node_index] += force_on_node;
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexBasedTissueForce<1>;
template class VertexBasedTissueForce<2>;
template class VertexBasedTissueForce<3>;
