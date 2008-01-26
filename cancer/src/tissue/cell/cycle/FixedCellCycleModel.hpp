#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleMeinekeCellCycleModel.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellG1Duration + SG2MDuration 
 *  and CancerParameters::TransitCellG1Duration + SG2MDuration)
 */
class FixedCellCycleModel : public AbstractSimpleMeinekeCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleMeinekeCellCycleModel>(*this);
    }
    
    /**
     * Private constructor for identical cells.
     */
    FixedCellCycleModel(double g1Duration, unsigned generation):
    	AbstractSimpleMeinekeCellCycleModel(g1Duration, generation){};
    
protected:    
    
public:

    /**
     * Default constructor - mBirthTime now set in AbstractCellCycleModel()
     * 						 mG1Duration now set in AbstractSimpleCellCycleModel()
     */
    FixedCellCycleModel() {};
    
    AbstractCellCycleModel *CreateDaughterCellCycleModel(); 
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(FixedCellCycleModel)


#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
