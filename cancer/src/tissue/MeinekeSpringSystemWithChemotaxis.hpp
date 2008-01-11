#ifndef MEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_
#define MEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_


#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Tissue.cpp"
#include "Meineke2001SpringSystem.hpp"
#include "CellwiseData.cpp"
#include "CellwiseDataGradient.hpp"

/**
 *  Meineke2001SystemWithChemotaxis
 * 
 *  A mechanics system for discrete tissue models based on the Meineke 2001 spring system. 
 *  There is an additional chemotactic force term which couples the tissue to CellwiseData.  
 * 
 *  Uses Fc = chi(C,|gradC) gradC/|gradC|, where C is the nutrient concentration
 *  and chi is a specified function. If gradC=0, Fc=0  
 * 
 *  TODO: remove duplicated code calculating damping constant
 */
template<unsigned DIM>
class MeinekeSpringSystemWithChemotaxis  : public Meineke2001SpringSystem<DIM>
{
friend class TestMeinekeSpringSystemWithChemotaxis;
    
private:
    double ChemotacticForceMagnitude(double nutrientConc, double nutrientGradientMagnitude)
    {
        return nutrientConc;
    }
    
public:

    MeinekeSpringSystemWithChemotaxis(Tissue<DIM>& rTissue)
        : Meineke2001SpringSystem<DIM>(rTissue)
    {}
    
    /**
     * Calculates the forces on each node
     *
     * @return the velocity components on each node. Of size NUM_NODES x DIM. 
     * The velocities are those that would be returned by the Meineke2001SpringSystem,
     * with a velocity due to the force by chemotaxis added on.
     * 
     * F_c = chi(C,|gradC|) gradC/|gradC|
     * 
     */
    std::vector<c_vector<double, DIM> >& rCalculateVelocitiesOfEachNode()
    {
        this->mDrDt = Meineke2001SpringSystem<DIM>::rCalculateVelocitiesOfEachNode();
        
        CellwiseDataGradient<DIM> gradients;        
        gradients.SetupGradients();
        
        for (typename Tissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
             cell_iter != this->mrTissue.End();
             ++cell_iter)
        {
            TissueCell& cell = *cell_iter;
            unsigned node_global_index = cell.GetNodeIndex();

            c_vector<double,DIM>& r_gradient = gradients.rGetGradient(cell.GetNodeIndex());            
            double nutrient_concentration = CellwiseData<DIM>::Instance()->GetValue(&cell,0);
            double magnitude_of_gradient = norm_2(r_gradient);

            double force_magnitude = ChemotacticForceMagnitude(nutrient_concentration, magnitude_of_gradient);
        
            double damping_constant = this->GetDampingConstant(cell); 
            
            // velocity += viscosity * chi * gradC/|gradC|
            if(magnitude_of_gradient > 0)
            {
                this->mDrDt[node_global_index] += (force_magnitude/(damping_constant*magnitude_of_gradient))*r_gradient;
            }
            //else Fc=0
        }

        return this->mDrDt;
    }
};

#endif /*MEINEKESPRINGSYSTEMWITHCHEMOTAXIS_HPP_*/
