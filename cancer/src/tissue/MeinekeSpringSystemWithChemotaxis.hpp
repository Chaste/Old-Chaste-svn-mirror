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
        
            // get the damping constant 
            // TODO: this is copied from Meineke2001SpringSystem - needs refactoring!
            double damping_multiplier = 1.0;
            
            if (this->mUseAreaBasedViscosity)
            {
                // use new_damping_const = old_damping_const * (d0+d1*A)
                // where d0,d1 are params and A is the area, and old_damping_const
                // if the damping const if not using mUseAreaBasedViscosity
                
                #define COVERAGE_IGNORE
                assert(DIM==2);
                #undef COVERAGE_IGNORE
                double rest_length = 1.0;
                double d0 = 0.1;
                // this number is such that d0+A*d1=1, where A is the area of a equilibrium
                // cell (=sqrt(3)/4 = a third of the area of a hexagon with edges of size 1)
                double d1 = 2.0*(1.0-d0)/(sqrt(3)*rest_length*rest_length); 
    
                VoronoiTessellation<DIM>& tess = this->mrTissue.rGetVoronoiTessellation();
            
                double area_cell = tess.GetFaceArea(node_global_index);
                
                // the areas should be order 1, this is just to avoid getting infinite areas
                // if an area based viscosity option is chosen without ghost nodes.
                assert(area_cell < 1000);
                
                damping_multiplier = d0 + area_cell*d1;
            }
            
            double damping_constant = CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplier;
            
            if( (cell.GetMutationState()!=HEALTHY) && (cell.GetMutationState()!=APC_ONE_HIT))
            {            
                damping_constant = CancerParameters::Instance()->GetDampingConstantMutant()*damping_multiplier;            
            } 
            
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
