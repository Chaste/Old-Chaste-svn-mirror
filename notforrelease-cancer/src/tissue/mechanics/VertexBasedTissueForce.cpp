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
c_vector<double, DIM> VertexBasedTissueForce<DIM>::GetDeformationForceContributionAtNode(unsigned localIndex, VertexElement<DIM, DIM>* pElement)
{
    // Compute the area of the element and its gradient at this node
    double element_area = pElement->GetArea();
    c_vector<double, DIM> area_gradient = pElement->GetAreaGradientAtNode(localIndex);
    
    /*
     * Get the target area of the cell. This grows linearly from 0.5*A to A
     * during the G1 phase of the cell cycle, then remains at A for the rest
     * of the cell cycle, where A denotes the cancer parameter mMatureCellTargetArea.
     */

    // Helper variable that is a static cast of the tissue
    //VertexBasedTissue<DIM>* p_tissue = static_cast<VertexBasedTissue<DIM>*>(&rTissue);
    
    double cell_target_area = CancerParameters::Instance()->GetMatureCellTargetArea();

//////// This is to have cells grow once they have divided i think we want to have cells grow before - see #852 //    
//    double cell_age = p_tissue->rGetCellUsingLocationIndex(pElement->GetIndex())->GetAge();
//    double cell_type = p_tissue->rGetCellUsingLocationIndex(pElement->GetIndex())->GetCellType();
//    double g1_duration = p_tissue->rGetCellUsingLocationIndex(pElement->GetIndex())->GetCellCycleModel->GetG1Duration();
//
//    if ( (cell_type!=DIFFERENTIATED) && (cell_age<g1_duration) )
//    {
//        cell_target_area *= 0.5*(1 + cell_age/g1_duration);
//    }
//    
    return 2*CancerParameters::Instance()->GetDeformationEnergyParameter()*(element_area - cell_target_area)*area_gradient;
}


template<unsigned DIM>
c_vector<double, DIM> VertexBasedTissueForce<DIM>::GetMembraneForceContributionAtNode(unsigned localIndex, VertexElement<DIM, DIM>* pElement)
{
    // Compute the perimeter of the element and its gradient at this node
    double element_perimeter = pElement->GetPerimeter();
    
    // \todo this needs to vary like the area in GetDeformationForceContributionAtNode() - see #852
    double target_perimeter = 2*sqrt(M_PI*CancerParameters::Instance()->GetMatureCellTargetArea());
            
    c_vector<double, DIM> perimeter_gradient = pElement->GetPerimeterGradientAtNode(localIndex);
            
    return 2*CancerParameters::Instance()->GetMembraneSurfaceEnergyParameter()*(element_perimeter - target_perimeter)*perimeter_gradient;      
}


template<unsigned DIM>
c_vector<double, DIM> VertexBasedTissueForce<DIM>::GetAdhesionForceContributionAtNode(unsigned localIndex, VertexElement<DIM, DIM>* pElement)
{
    /// \todo code up this method (see #861)
    return zero_vector<double>(DIM);
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
            // Find the local index of this node in this element
            unsigned local_index = 0;
            while (p_tissue->GetElement(*iter)->GetNodeGlobalIndex(local_index) != node_index)
            {
                local_index++;
            }
            
            // Compute force contributions from cell deformation energy, 
            // membrane surface tension energy and a cell-cell adhesion energy
            
            //\todo ticket861 need to uncomment these
            deformation_contribution += GetDeformationForceContributionAtNode(local_index, p_tissue->GetElement(*iter));
            membrane_surface_tension_contribution += GetMembraneForceContributionAtNode(local_index, p_tissue->GetElement(*iter));
            //cell_cell_adhesion_contribution += GetAdhesionForceContributionAtNode(local_index, p_tissue->GetElement(*iter));
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
