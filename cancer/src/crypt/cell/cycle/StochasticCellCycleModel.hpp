#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CancerParameters.hpp"

/**
 *  Stochastic cell model
 *
 *  Cell cycle time is deterministic for stem cells and stochastic (normally
 *  distributed with mean CancerParameters::TransitCellCycleTime and variance 1)
 */
class StochasticCellCycleModel : public AbstractCellCycleModel
{
private:
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        // Make sure the singletons we use are archived
        CancerParameters* p_params = CancerParameters::Instance();
        archive & *p_params;
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
    }
    
public:
    StochasticCellCycleModel();
    
    virtual bool ReadyToDivide();
    
    virtual void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel();
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticCellCycleModel)


#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
