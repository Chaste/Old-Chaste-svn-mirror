#ifndef _ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
#define _ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"
#include "AbstractOdeSystem.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to standard cell cycle 
 * ODE models for intracellular protein concentrations. Along the lines
 * of Tyson & Novak etc... 
 */
class AbstractOdeBasedCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        assert(mpOdeSystem!=NULL);
        archive & mpOdeSystem->rGetStateVariables();
        archive & mLastTime;
        archive & mDivideTime;
        archive & mReadyToDivide;
        archive & mFinishedRunningOdes;
    }
    
protected:
    /** The system of ODEs for the cell cycle model */
    AbstractOdeSystem* mpOdeSystem;
    /** The last time the cell cycle ODEs were evaluated.*/
    double mLastTime;
    /** The time at which the cell should divide - Set this to DBL_MAX in constructor.*/
    double mDivideTime;
    /** Whether the cell is ready to divide or not */
    bool mReadyToDivide;
    /** Whether the cell cycle model is currently in a delay (not solving ODEs).*/
    bool mFinishedRunningOdes;
    
public:

    /**
     * Creates an AbstractOdeBasedCellCycleModel, calls SetBirthTime on the 
     * AbstractCellCycleModel to make sure that can be set 'back in time' for
     * cells which did not divide at the current time.
     * 
     * @param lastTime  The birth time of the cell / last time model was evaluated (defaults to the current SimulationTime)
     */
    AbstractOdeBasedCellCycleModel(double lastTime = SimulationTime::Instance()->GetDimensionalisedTime());

    /**
     * This destructor deletes the mpOdeSystem. 
     */
    virtual ~AbstractOdeBasedCellCycleModel();
    
    /**
     * ReadyToDivide() is called by a cell on its model to establish progress through the cell cycle.
     * Can be overridden if something special should happen in the subclass.
     * 
     * @return true if the cell is ready to divide.
     */
    virtual bool ReadyToDivide();
    
    /**
     * This method must be implemented by each subclass - solves the ODEs to a given time and 
     * 
     * @return Whether a stopping event occurred.
     */    
    virtual bool SolveOdeToTime(double currentTime) = 0;
    
    /**
     * This method must be implemented by each subclass
     * FIX COMMENT
     * When the ODEs have reached a stopping event it returns the time at which 
     * the cell should divide, so a delay can be added in for S-G2-M phases if necessary.
     */
    virtual double GetSG2Duration();
    virtual double GetMDuration();
    virtual double GetOdeStopTime() = 0;
    
    /**
     * This overrides the AbstractCellCycleModel::SetBirthTime(double birthTime)
     * because an ODE based cell cycle model has more to reset...
     * 
     * @param birthTime the simulation time when the cell was born
     */
    void SetBirthTime(double birthTime);
    
    /**
     * Returns the protein concentrations at the current time (useful for tests)
     *
     * NB: Will copy the vector - you can't use this to modify the concentrations.
     */
    std::vector<double> GetProteinConcentrations() const;
    
    /**
     * Sets the protein concentrations and time when the model was last evaluated - should only be called by tests
     *
     * @param lastTime the SimulationTime at which the protein concentrations apply
     * @param proteinConcentrations a standard vector of doubles of protein concentrations
     *
     */
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
    
    /**
     * For a naturally cycling model this does not need to be overridden in the 
     * subclasses. But most models should override this function and then 
     * call AbstractOdeBasedCellCycleModel::ResetModel() from inside their version. 
     */
    virtual void ResetModel();
    
};

BOOST_IS_ABSTRACT(AbstractOdeBasedCellCycleModel)

#endif //_ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
