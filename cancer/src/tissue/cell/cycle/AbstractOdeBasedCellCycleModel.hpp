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
    
public:

    /**
     * This destructor deletes mpOdeSystem. Note that some cell cycle models
     * must call Initialise() to set up the ODE system before calling this
     * destructor method, or a seg fault could result.
     */
    virtual ~AbstractOdeBasedCellCycleModel();
    
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
    

};

BOOST_IS_ABSTRACT(AbstractOdeBasedCellCycleModel)

#endif //_ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
