#ifndef OXYGENBASEDCELLKILLER_HPP_
#define OXYGENBASEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ConformingTetrahedralMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"
#include "CellwiseData.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 *  Kills cells based on the local oxygen concentration.
 *  TODO: need to check whether this takes into account current times or timesteps.
 */
template <unsigned SPACE_DIM>
class OxygenBasedCellKiller : public AbstractCellKiller<SPACE_DIM>
{
private: 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<SPACE_DIM> >(*this);
        
    }
    
public:
    OxygenBasedCellKiller(Tissue<SPACE_DIM>* pTissue)
        : AbstractCellKiller<SPACE_DIM>(pTissue)
    {
    }
    
    void TestAndLabelSingleCellForApoptosis(TissueCell& rCell)
    {        
        if (rCell.GetCellType()!=HEPA_ONE)
        {
            EXCEPTION("OxygenBasedCellKiller is trying to kill a cell that is not of type HEPA_ONE");
        }    
        
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(&rCell);
        
//      TODO: change this line to something like
        // if ( oxygen_concentration < CancerParameters::Instance()->GetHypoxicConcentration() )   
		if ( oxygen_concentration < 0.2 )
        {
        	double hypoxic_duration = rCell.GetHypoxicDuration();
            
            // a little bit of stochasticity here
            double prob_of_death = 1 - 0.5*(oxygen_concentration/0.2); 
            
//      TODO: change this line to something like
		// if (!cell.HasApoptosisBegun() && hypoxic_duration > CancerParameters::Instance()->GetCriticalHypoxicDuration() && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
            if (!rCell.HasApoptosisBegun() && hypoxic_duration > 1.0 && RandomNumberGenerator::Instance()->ranf() < prob_of_death)
            {                     
                rCell.StartApoptosis();
            }
        }           
    }
    
    /**
     * Loops over cells and starts apoptosis if the cell has been hypoxic 
     * for longer than some critical period, and if it is currently hypoxic,
     * with a bit of stochasticity also thrown in for good measure
     */
     virtual void TestAndLabelCellsForApoptosisOrDeath()
     {      
        for( typename Tissue<SPACE_DIM>::Iterator cell_iter = this->mpTissue->Begin();
            cell_iter != this->mpTissue->End();
            ++cell_iter)
        {               
            TestAndLabelSingleCellForApoptosis(*cell_iter);
        }
    }   
};

#include "TemplatedExport.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(OxygenBasedCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TissueSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OxygenBasedCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // save data required to construct instance
    const Tissue<DIM>* const p_tissue = t->GetTissue();
    ar << p_tissue;
}

/**
 * De-serialize constructor parameters and initialise Tissue.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OxygenBasedCellKiller<DIM> * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
    Tissue<DIM>* p_tissue;
    ar >> p_tissue;
    // invoke inplace constructor to initialize instance
    ::new(t)OxygenBasedCellKiller<DIM>(p_tissue);
}
}
} // namespace ...

#endif /*OXYGENBASEDCELLKILLER_HPP_*/
