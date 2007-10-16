#ifndef _ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_
#define _ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "CancerParameters.hpp"
#include "SimulationTime.hpp"
#include "AbstractOdeSystem.hpp"

// Needs to be included last
#include <boost/serialization/export.hpp>

/**
 * This class contains all the things common to the Wnt cell cycle ODE based models,
 * the Resetting functions and Updating of cell types etc.
 */
class AbstractWntOdeBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
    }
    
    // no member variables yet - if any are added put them in archive function
    // and add default values in the default constructor.
    
protected:

    static RungeKutta4IvpOdeSolver msSolver;

    AbstractWntOdeBasedCellCycleModel(double lastTime) 
        : AbstractOdeBasedCellCycleModel(lastTime) {};
        
    /**
     * Introduces the delay after ODEs have been solved.
     */
    virtual double GetDivideTime();
    
    /**
     * Introduces the delay after ODEs have been solved,
     * overridden in subclass StochasticWntCellCycleModel
     */
    virtual double GetWntSG2MDuration();
        
public:
    /**
     * Just a default constructor (no member variables)
     */
    AbstractWntOdeBasedCellCycleModel() {};
    
    /**
     * Resets the Wnt Model to the start of the cell cycle (this model does not cycle naturally)
     * Cells are given a new birth time and cell cycle proteins are reset.
     * Note that the wnt pathway proteins maintain their current values.
     *
     * Should only be called by the TissueCell::Divide() method.
     */
    void ResetModel();
    
    /**
     * Updates the current cell type to reflect whether the
     * beta-catenin level has dropped low enough to make it stop dividing.
     * This should only be called when the cell cycle model has been 
     * evaluated to the current time, or it may give misleading results.
     */
    void UpdateCellType();
    
    /**
     * This must be implemented by subclasses to change cell type to reflect
     * current levels of beta-catenin.
     */
    virtual void ChangeCellTypeDueToCurrentBetaCateninLevel() = 0;
    
           
};

BOOST_IS_ABSTRACT(AbstractWntOdeBasedCellCycleModel)

#endif //_ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_


