#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellG1Duration + SG2MDuration 
 *  and CancerParameters::TransitCellG1Duration + SG2MDuration)
 */
class FixedCellCycleModel : public AbstractSimpleCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
    }
    
    /**
     * Private constructor for identical cells.
     */
    FixedCellCycleModel(double g1Duration):
    	AbstractSimpleCellCycleModel(g1Duration) {};
    
protected:    
    
public:

    /**
     * Default constructor - mBirthTime now set in AbstractCellCycleModel()
     * 						 mG1Duration now set in AbstractSimpleCellCycleModel()
     */
    FixedCellCycleModel() {};
    
    AbstractCellCycleModel *CreateCellCycleModel(); 
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(FixedCellCycleModel)


#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
