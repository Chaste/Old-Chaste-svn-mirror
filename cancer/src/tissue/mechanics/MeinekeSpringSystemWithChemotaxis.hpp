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
#ifndef MEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_
#define MEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Meineke2001SpringSystem.hpp"
#include "CellwiseDataGradient.hpp"

/**
 *  MeinekeSpringSystemWithChemotaxis
 * 
 *  A mechanics system for discrete tissue models based on the Meineke 2001 spring system. 
 *  There is an additional chemotactic force term which couples the tissue to CellwiseData. 
 *  Currently, only LABELLED cells experience this force. 
 * 
 *  Uses Fc = chi(C,|gradC|) gradC/|gradC|, where C is the nutrient concentration
 *  and chi is a specified function. If gradC=0, Fc=0
 */
template<unsigned DIM>
class MeinekeSpringSystemWithChemotaxis  : public Meineke2001SpringSystem<DIM>
{
friend class TestMeinekeSpringSystemWithChemotaxis;
    
private:

    double ChemotacticForceMagnitude(const double nutrientConc, const double nutrientGradientMagnitude);
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<Meineke2001SpringSystem<DIM> >(*this);
    }
    
public:

    MeinekeSpringSystemWithChemotaxis(MeshBasedTissue<DIM>& rTissue);
    
    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM. 
     * The velocities are those that would be returned by the Meineke2001SpringSystem,
     * with a velocity due to the force by chemotaxis added on.
     * 
     * Fc = chi(C,|gradC|) gradC/|gradC|  (if |gradC|>0, else Fc = 0)
     * 
     */
    std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode();
    
};


template<unsigned DIM>
MeinekeSpringSystemWithChemotaxis<DIM>::MeinekeSpringSystemWithChemotaxis(MeshBasedTissue<DIM>& rTissue)
        : Meineke2001SpringSystem<DIM>(rTissue)
{
}

    
template<unsigned DIM>
double MeinekeSpringSystemWithChemotaxis<DIM>::ChemotacticForceMagnitude(const double nutrientConc, const double nutrientGradientMagnitude)
{
    return nutrientConc; // temporary force law - can be changed to something realistic
                         // without tests failing
}


template<unsigned DIM>
std::vector<c_vector<double, DIM> >& MeinekeSpringSystemWithChemotaxis<DIM>::rCalculateVelocitiesOfEachNode()
{
    this->mDrDt = Meineke2001SpringSystem<DIM>::rCalculateVelocitiesOfEachNode();
    
    CellwiseDataGradient<DIM> gradients;        
    gradients.SetupGradients();
    
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mpTissue->Begin();
         cell_iter != this->mpTissue->End();
         ++cell_iter)
    {
        // Only LABELLED cells move chemotactically
        if (cell_iter->GetMutationState()==LABELLED)
        {
            TissueCell& cell = *cell_iter;            
            unsigned node_global_index = cell.GetNodeIndex();

            c_vector<double,DIM>& r_gradient = gradients.rGetGradient(cell.GetNodeIndex());
            double nutrient_concentration = CellwiseData<DIM>::Instance()->GetValue(&cell,0);
            double magnitude_of_gradient = norm_2(r_gradient);
            
            double force_magnitude = ChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);
            
            double damping_constant = this->GetDampingConstant(cell);        
            
            // velocity += viscosity * chi * gradC/|gradC|
            if (magnitude_of_gradient > 0)
            {
                this->mDrDt[node_global_index] += (force_magnitude/(damping_constant*magnitude_of_gradient))*r_gradient;
            }
            // else Fc=0
        }
    }

    return this->mDrDt;
}


#include "TemplatedExport.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MeinekeSpringSystemWithChemotaxis)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MeinekeSpringSystemWithChemotaxis.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MeinekeSpringSystemWithChemotaxis<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<DIM> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MeinekeSpringSystemWithChemotaxis<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<DIM>* p_tissue;

    ar >> p_tissue;
    // Invoke inplace constructor to initialize instance
    ::new(t)MeinekeSpringSystemWithChemotaxis<DIM>(*(static_cast<MeshBasedTissue<DIM>*>(p_tissue)));
}
}
} // namespace ...

#endif /*MEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_*/
