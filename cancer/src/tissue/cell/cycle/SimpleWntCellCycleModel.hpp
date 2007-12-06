#ifndef SIMPLEWNTCELLCYCLEMODEL_HPP_
#define SIMPLEWNTCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 *  Simple Wnt-dependent cell cycle model
 *
 */
class SimpleWntCellCycleModel : public AbstractSimpleCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
       
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
    }
    
    /**
     * Private function that should only be called by Reset() and SetCell()
     * this introduces the stochastic element of the model.
     */    
    void SetG1Duration();
    
    /**
     * Private constructor for identical cells.
     */
    SimpleWntCellCycleModel(double g1Duration, unsigned generation):
    	AbstractSimpleCellCycleModel(g1Duration, generation) {};
    
public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when SetCell() is called
     */
    SimpleWntCellCycleModel() {};
    
    bool ReadyToDivide();
        
    AbstractCellCycleModel *CreateCellCycleModel(); 
    
    std::vector<CellType> GetNewCellTypes();
        
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleWntCellCycleModel)


#endif /*SIMPLEWNTCELLCYCLEMODEL_HPP_*/
