#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include <cassert>
#include <iostream>

#include "AbstractSimpleMeinekeCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "Exception.hpp"

/**
 *  Stochastic cell model
 *  
 */
class StochasticCellCycleModel : public AbstractSimpleMeinekeCellCycleModel
{
private:
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleMeinekeCellCycleModel>(*this);
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }
    
    /**
     * Stochastically set the G1 duration.  Called on cell creation at 
     * the start of a simulation, and for both parent and daughter 
     * cells at cell division. 
     */
    void SetG1Duration();
    
    /**
     * Private constructor for identical cells.
     */
    StochasticCellCycleModel(double g1Duration, unsigned generation)
        : AbstractSimpleMeinekeCellCycleModel(g1Duration, generation)
    {}
    
public:
    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    StochasticCellCycleModel() 
    {}
    
    AbstractCellCycleModel *CreateDaughterCellCycleModel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticCellCycleModel)


#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
