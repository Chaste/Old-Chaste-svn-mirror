#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellCycleTime and CancerParameters::TransitCellCycleTime)
 */
class FixedCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    }
public:

    /**
     * Default constructor - mBirthTime now set in AbstractCellCycleModel().
     */
    FixedCellCycleModel() {};

    virtual bool ReadyToDivide();
    
    virtual void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel(); 
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(FixedCellCycleModel)


#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
