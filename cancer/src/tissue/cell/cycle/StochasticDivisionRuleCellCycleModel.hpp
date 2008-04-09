#ifndef STOCHASTICDIVISIONRULECELLCYCLEMODEL_HPP_
#define STOCHASTICDIVISIONRULECELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 *  Simple cell cycle model for use in crypt projection 
 *  simulations under the non-niche hypothesis. 
 *
 */
class StochasticDivisionRuleCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        
        archive & mDividedSymmetrically;
    }
    
    /** 
     * Whether or not a cell was born out of symmetric division from a stem cell.
     * Note that we do not include the case TRANSIT -> 2 TRANSIT in our definition
     * of 'symmetric division'.
     * 
     * This flag is needed because if a daughter cell is given type STEM 
     * and generation 0 when it is created, we cannot tell whether it was born out 
     * of symmetric division STEM -> 2 STEM or asymmetric division STEM -> STEM + TRANSIT.
     */ 
    bool mDividedSymmetrically;
    
    /**
     *  Stochastically set the G1 duration.  Called on cell creation at 
     *  the start of a simulation, and for both parent and daughter 
     *  cells at cell division.
     * 
     *  The G1 duration is taken from a normal distribution, whose mean is
     *  the G1 duration given in CancerParameters for the cell type, and 
     *  whose standard deviation is 1. 
     */    
    void SetG1Duration();
    
    /**
     * Private constructor for identical cells
     */
    StochasticDivisionRuleCellCycleModel(double g1Duration, unsigned generation, bool dividedSymmetrically) 
        : AbstractSimpleCellCycleModel(g1Duration, generation),
          mDividedSymmetrically(dividedSymmetrically) 
    {}
        
public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    StochasticDivisionRuleCellCycleModel(bool dividedSymmetrically=false)
        : mDividedSymmetrically(dividedSymmetrically) 
    {}
    
    void ResetForDivision();
        
    void InitialiseDaughterCell();
    
    AbstractCellCycleModel* CreateDaughterCellCycleModel();
    
    bool DividedSymmetrically();

};

// Declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticDivisionRuleCellCycleModel)

#endif /*STOCHASTICDIVISIONRULECELLCYCLEMODEL_HPP_*/
