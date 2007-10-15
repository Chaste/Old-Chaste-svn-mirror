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
        archive & mDivisionAge;
    }
    
    /** What time this cell would like to divide (assigned randomly when model is given a cell) */
    double mDivisionAge;
    
    /** Private constructor for creating an identical daughter cell */
    StochasticCellCycleModel(double divisionAge)
        :mDivisionAge(divisionAge) {};
    
    /**
     * Private function that should only be called by Reset() and SetCell()
     * this introduces the stochastic element of the model.
     */
    void SetDivisionAge();
    
public:
    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mDivisionAge is set very high, it is set for the individual cells when SetCell() is called
     */
    StochasticCellCycleModel()
        : mDivisionAge(DBL_MAX) {};
    
    /** 
     * Overridden SetCellMethod - also assigns a DivisionAge based on the cell type.
     * 
     * @param pCell the cell which owns this model.
     */
    void SetCell(TissueCell* pCell);
    
    bool ReadyToDivide();
    
    void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticCellCycleModel)


#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/
