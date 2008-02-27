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
        
        archive & mUseCellTypeDependentG1Duration;
    }

    /**
     * Whether to use different mean G1 durations for different cell types.
     * For use in SetG1Duration().
     */
    bool mUseCellTypeDependentG1Duration;
        
protected:

    /**
     * Stochastically set the G1 duration.  Called on cell creation at 
     * the start of a simulation, and for both parent and daughter 
     * cells at cell division.
     */    
    void SetG1Duration();
    
    /**
     * Private constructor for identical cells.
     */
    SimpleWntCellCycleModel(double g1Duration, unsigned generation, bool useCellTypeDependentG1Duration=false) 
        : AbstractSimpleCellCycleModel(g1Duration, generation),
          mUseCellTypeDependentG1Duration(useCellTypeDependentG1Duration)
    {} 
        
public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    SimpleWntCellCycleModel(bool useCellTypeDependentG1Duration=false)
        : mUseCellTypeDependentG1Duration(useCellTypeDependentG1Duration)
    {}
    
    /** 
     * Overridden UpdateCellCyclePhase() method
     */ 
    void UpdateCellCyclePhase();
        
    AbstractCellCycleModel *CreateDaughterCellCycleModel(); 
        
    void ResetForDivision();
    
    void InitialiseDaughterCell();        
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleWntCellCycleModel)


#endif /*SIMPLEWNTCELLCYCLEMODEL_HPP_*/
