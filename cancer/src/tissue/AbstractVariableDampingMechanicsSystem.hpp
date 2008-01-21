#ifndef ABSTRACTVARIABLEDAMPINGMECHANICSSYSTEM_HPP_
#define ABSTRACTVARIABLEDAMPINGMECHANICSSYSTEM_HPP_

#include "AbstractDiscreteTissueMechanicsSystem.hpp"

/** 
 *  AbstractVariableDampingMechanicsSystem
 *  
 *  An abstract tissue mechanics system that provides a GetDamping() method which
 *  returns a variable damping factor. The two ways that the damping is changed from
 *  the normal value - CancerParameters::GetDampingConstantNormal() - is if
 *  
 *  (a) a cell's is neither HEALTHY nor APC_ONE_HIT
 *  (b) if mUseAreaBasedViscosity is turned on, in which case the damping = old_damping_const*(d0+d1*A),
 *  where A is the cell's voronoi area, and d0,d1 are constants (see code for values), and
 *  old_damping_const is whatever the damping constant would be in (a).
 * 
 *  \todo: make d0,d1 member variables, or allow the user to provide a functional form of
 *  d(A), if this ever becomes needed.
 */
template<unsigned DIM>
class AbstractVariableDampingMechanicsSystem : public AbstractDiscreteTissueMechanicsSystem<DIM>
{
private :
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractDiscreteTissueMechanicsSystem<DIM> >(*this);
        archive & mUseAreaBasedViscosity;
    }

protected :
    /** Whether to use a viscosity that is linear in the cell area, rather than constant */
    bool mUseAreaBasedViscosity;

    /** 
     *  Get the damping constant for this cell - ie d in drdt = F/d
     *  This depends on whether using area-based viscosity has been switched on, and 
     *  on whether the cell is a mutant or not
     */
    double GetDampingConstant(TissueCell& rCell)
    { 
        double damping_multiplier = 1.0;
        
        if (this->mUseAreaBasedViscosity)
        {
            // the subclass had better have said is needs voronoi tessellations..
            assert(this->NeedsVoronoiTessellation());
            
            //  use new_damping_const = old_damping_const * (d0+d1*A)
            //  where d0, d1 are params and A is the area, and old_damping_const
            //  is the damping const if not using mUseAreaBasedViscosity
            #define COVERAGE_IGNORE
            assert(DIM==2);
            #undef COVERAGE_IGNORE

            double rest_length = 1.0;
            double d0 = 0.1;

            // this number is such that d0+A*d1=1, where A is the area of a equilibrium
            // cell (=sqrt(3)/4 = a third of the area of a hexagon with edges of size 1)
            double d1 = 2.0*(1.0-d0)/(sqrt(3)*rest_length*rest_length); 

            VoronoiTessellation<DIM>& tess = this->mrTissue.rGetVoronoiTessellation();
        
            double area_cell = tess.GetFaceArea(rCell.GetNodeIndex());
            
            // the areas should be order 1, this is just to avoid getting infinite areas
            // if an area based viscosity option is chosen without ghost nodes.
            assert(area_cell < 1000);
            
            damping_multiplier = d0 + area_cell*d1;
        }
                
        if( (rCell.GetMutationState()!=HEALTHY) && (rCell.GetMutationState()!=APC_ONE_HIT))
        {            
            return CancerParameters::Instance()->GetDampingConstantMutant()*damping_multiplier;            
        }
        else 
        {
            return CancerParameters::Instance()->GetDampingConstantNormal()*damping_multiplier;
        }
    } 
    

public :
    AbstractVariableDampingMechanicsSystem(Tissue<DIM>& rTissue)
        : AbstractDiscreteTissueMechanicsSystem<DIM>(rTissue)
    {
        mUseAreaBasedViscosity = false;
    }
    

    /**
     * Use an area based viscosity
     */
    void SetAreaBasedViscosity(bool useAreaBasedViscosity)
    {
        assert(DIM == 2);
        mUseAreaBasedViscosity = useAreaBasedViscosity;
    }
};
#endif /*ABSTRACTVARIABLEDAMPINGMECHANICSSYSTEM_HPP_*/
