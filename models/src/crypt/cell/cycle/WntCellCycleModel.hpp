#ifndef WNTCELLCYCLEMODEL_HPP_
#define WNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"
#include "WntCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
#include "WntGradient.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 *  Wnt-dependent cell cycle model.
 * 
 * Note that this class uses C++'s default copying semantics, and so doesn't implement a copy constructor
 * or operator=.
 */
class WntCellCycleModel : public AbstractCellCycleModel
{
    friend class StochasticWntCellCycleModel;// to allow access to private constructor below.
    friend class boost::serialization::access;   
    
private:
    WntCellCycleOdeSystem* mpOdeSystem;
    static RungeKutta4IvpOdeSolver msSolver;
    double mLastTime;
    double mDivideTime;
    bool mInSG2MPhase;
    bool mReadyToDivide;
    double mInitialWntStimulus;
    WntGradient& mrWntGradient; 
    
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        // reference can be read or written into once mpOdeSystem has been set up
        // mpOdeSystem isn't set up by the first constructor, but is by the second
        // which is now utilised by the load_construct at the bottom of this file.
        archive & mpOdeSystem->rGetStateVariables();   
        archive & mpOdeSystem->rGetMutationState(); 
        archive & mLastTime;
        archive & mDivideTime;
        archive & mInSG2MPhase;
        archive & mReadyToDivide;   
    }
       
protected:    
    virtual double GetWntSG2MDuration();
    
public:

    WntCellCycleModel(double InitialWntStimulus,  WntGradient &rWntGradient);
   
   /**
     * This is needed because a wnt model which is not to be run from the current time is 
     * sometimes needed. Should only be called by the cell itself when it wants to divide.
     * And by archiver to set up a new ODE system.
     */
    WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                      const CryptCellMutationState& rMutationState, 
                      double birthTime, WntGradient &rWntGradient);
     
    virtual ~WntCellCycleModel();
    
    virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>());
    
    virtual void ResetModel();
    
    CryptCellType UpdateCellType();
    
    std::vector< double > GetProteinConcentrations() const;
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    void SetBirthTime(double birthTime);
    
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
    
    void SetCell(MeinekeCryptCell* pCell);
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(WntCellCycleModel)


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
    Archive & ar, const WntCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a WntCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, WntCellCycleModel * t, const unsigned int file_version)
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
    double birth_time = 0.0;
    CryptCellMutationState mutation_state = HEALTHY;
    WntGradient* p_dummy_wnt_gradient = (WntGradient*)NULL; 
    WntGradient& r_wnt_gradient = *p_dummy_wnt_gradient; 
    ::new(t)WntCellCycleModel(state_vars, mutation_state, birth_time, r_wnt_gradient);
}
}
} // namespace ...

#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
