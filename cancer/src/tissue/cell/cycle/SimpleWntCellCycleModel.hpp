#ifndef SIMPLEWNTCELLCYCLEMODEL_HPP_
#define SIMPLEWNTCELLCYCLEMODEL_HPP_

#include "FixedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 *  Simple Wnt-dependent cell cycle model
 *
 */
class SimpleWntCellCycleModel : public FixedCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<FixedCellCycleModel>(*this);
       
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
    }
    
    /** Private constructor for creating an identical daughter cell */
    SimpleWntCellCycleModel(double g1Duration)
        : FixedCellCycleModel(g1Duration) {};
        
    /**
     * Private function that should only be called by Reset() and SetCell()
     * this introduces the stochastic element of the model.
     */    
    void SetG1Duration();
    
public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when SetCell() is called
     */
    SimpleWntCellCycleModel()
        : FixedCellCycleModel() {};
    
    virtual bool ReadyToDivide();
        
    AbstractCellCycleModel *CreateCellCycleModel(); 
        
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleWntCellCycleModel)


#endif /*SIMPLEWNTCELLCYCLEMODEL_HPP_*/
