#ifndef SIMPLEWNTCELLCYCLEMODEL_HPP_
#define SIMPLEWNTCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CancerParameters.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellG1Duration + SG2MDuration 
 * and CancerParameters::TransitCellG1Duration + SG2MDuration)
 */
class SimpleWntCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mCycleTime;
    }

    /** How long this cell cycle will take, set by constructor or ResetModel() */
    double mCycleTime;
    
public:

    /**
     * Default constructor - mBirthTime now set in AbstractCellCycleModel().
     * 
     * Sets up the mCycleTime randomly
     */
    SimpleWntCellCycleModel()
    : mCycleTime(RandomNumberGenerator::Instance()->
                    NormalRandomDeviate(
                      CancerParameters::Instance()->GetTransitCellG1Duration()
                        + CancerParameters::Instance()->GetSG2MDuration(), 1.0)) {};

    virtual bool ReadyToDivide();
    
    virtual void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel(); 
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleWntCellCycleModel)


#endif /*SIMPLEWNTCELLCYCLEMODEL_HPP_*/
