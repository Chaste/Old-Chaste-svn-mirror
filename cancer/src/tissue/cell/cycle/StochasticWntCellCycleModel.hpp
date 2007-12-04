#ifndef STOCHASTICWNTCELLCYCLEMODEL_HPP_
#define STOCHASTICWNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "WntCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

class StochasticWntCellCycleModel : public WntCellCycleModel
{
  private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<WntCellCycleModel>(*this);
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }
    
    /**
     * This is a function which overrides that in WntCellCycleModel and 
     * introduces the stochastic element of this class. 
     * We allow the duration of the G2 phase of the cell cycle to 
     * vary with a mean of its deterministic duration and a standard 
     * deviation of 0.9 hours.
     * 
     * @return the duration of the G2 phases of the cell cycle.
     */
    double GetG2Duration()
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double mean = CancerParameters::Instance()->GetG2Duration();
        double standard_deviation = 0.9;        
        return p_gen->NormalRandomDeviate(mean,standard_deviation);
    }
    
        
  public:
  
    /**
     * The standard constructor called in tests
     */
    StochasticWntCellCycleModel()
      :  WntCellCycleModel()
    {
    }
    
    /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed. Should only be called by the cell itself when it wants to divide.
     */
    StochasticWntCellCycleModel(AbstractOdeSystem* pParentOdeSystem,
                                CellMutationState mutationState,
                                double birthTime, double lastTime,
                                bool inSG2MPhase, bool readyToDivide, double divideTime, unsigned generation)
      : WntCellCycleModel(pParentOdeSystem, mutationState, birthTime, lastTime, 
                          inSG2MPhase, readyToDivide, divideTime, generation)
    {
    }

    /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed. Should only be called by the archiver.
     */
    StochasticWntCellCycleModel(std::vector<double> proteinConcentrations, 
                                CellMutationState mutationState)
      : WntCellCycleModel(proteinConcentrations, mutationState)
    {
    }
    
    /**
     * Returns a new StochasticWntCellCycleModel created with the correct initial conditions.
     *
     * Should be called just after the parent cell cycle model has been .Reset().
     *
     */
    AbstractCellCycleModel* CreateCellCycleModel()
    {
        assert(mpCell!=NULL);
        // calls a cheeky version of the constructor which makes the new cell cycle model
        // the same age as the old one - not a copy at this time.
        return new StochasticWntCellCycleModel(mpOdeSystem, mpCell->GetMutationState(), mBirthTime, mLastTime, mFinishedRunningOdes, mReadyToDivide,mDivideTime, mGeneration);
    }
    
};


// declare identifier for the serializer
BOOST_CLASS_EXPORT(StochasticWntCellCycleModel)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const StochasticWntCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a StochasticWntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, StochasticWntCellCycleModel * t, const unsigned int file_version)
{
    // It doesn't actually matter what values we pass to our standard
    // constructor, provided they are valid parameter values, since the
    // state loaded later from the archive will overwrite their effect in
    // this case.
    // Invoke inplace constructor to initialize instance of my_class

    std::vector<double> state_vars;
    for (unsigned i=0 ; i<9 ; i++)
    {
        state_vars.push_back(0.0);
    }   

    CellMutationState mutation_state = HEALTHY;

    ::new(t)StochasticWntCellCycleModel(state_vars, mutation_state);
}
}
} // namespace ...

#endif /*STOCHASTICWNTCELLCYCLEMODEL_HPP_*/
