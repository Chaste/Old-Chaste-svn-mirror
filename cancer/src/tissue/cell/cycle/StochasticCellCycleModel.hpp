#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 *  Stochastic cell model
 *  
 */
class StochasticCellCycleModel : public AbstractSimpleCellCycleModel
{
private:
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }
    
    /**
     * Private function that should only be called by Reset() and SetCell()
     * this introduces the stochastic element of the model.
     */
    void SetG1Duration();
    
    /**
     * Private constructor for identical cells.
     */
    StochasticCellCycleModel(double g1Duration):
    	AbstractSimpleCellCycleModel(g1Duration) {};
    
public:
    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when SetCell() is called
     */
    StochasticCellCycleModel() {};
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticCellCycleModel)


#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
