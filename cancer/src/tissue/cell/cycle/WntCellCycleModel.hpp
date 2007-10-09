#ifndef WNTCELLCYCLEMODEL_HPP_
#define WNTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
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
class WntCellCycleModel : public AbstractOdeBasedCellCycleModel
{
private:
    static RungeKutta4IvpOdeSolver msSolver;
    
    friend class boost::serialization::access;   
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL); 
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
        // reference can be read or written into once mpOdeSystem has been set up
        // mpOdeSystem isn't set up by the first constructor, but is by the second
        // which is now utilised by the load_construct at the bottom of this file.
        archive & mpOdeSystem->rGetStateVariables();   
        archive & mpOdeSystem->rGetMutationState(); 
        archive & mInSG2MPhase;
    }
    
    /**
     * Called by ::Initialise() and ::UpdateCellType() only.
     * Updates the mpCell::mCellType to match mpOdeSystem's 
     * beta-catenin levels
     */
    void ChangeCellTypeDueToCurrentBetaCateninLevel();
    
    virtual double GetWntSG2MDuration();
       
protected:   
     
    WntCellCycleOdeSystem* mpOdeSystem;
    /** Is the cell in S->G2->M phase of cell cycle (if so don't run ODEs)*/
    bool mInSG2MPhase;    
    
public:

    WntCellCycleModel();
   
    /**
     * This is needed to create an exact copy of the current cell cycle model
     * (called by CreateCellCycleModel())
     */
    WntCellCycleModel(WntCellCycleOdeSystem* pParentOdeSystem, 
                      const CellMutationState& rMutationState, 
                      double birthTime, double lastTime,
                      bool inSG2MPhase, bool readyToDivide, double divideTime);
   /**
     * This is needed to create an exact copy of the current cell cycle model
     * (called by archiving functions)
     */
    WntCellCycleModel(const std::vector<double>& rParentProteinConcentrations, 
                      const CellMutationState& rMutationState, 
                      double birthTime, double lastTime,
                      bool inSG2MPhase, bool readyToDivide, double divideTime); 
                          
    virtual ~WntCellCycleModel();
    
    virtual bool ReadyToDivide();
    
    virtual void ResetModel();
    
    void UpdateCellType();    
    
    std::vector< double > GetProteinConcentrations() const;
    
    AbstractCellCycleModel *CreateCellCycleModel();
    
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
      
    void Initialise();
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

    CellMutationState mutation_state = HEALTHY;
    ::new(t)WntCellCycleModel(state_vars, mutation_state, 0.0, 0.0, false, false, 0.0);
}
}
} // namespace ...

#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
