#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "CellTypes.hpp"
#include "SimulationTime.hpp"
#include "TissueCell.hpp"
#include <vector>

// Needs to be included last
#include <boost/serialization/export.hpp>

class TissueCell; // Circular definition (cells need to know about cycle models and vice-versa).

class AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mBirthTime;
        // Make sure the simulation time gets saved too
        SimulationTime* p_time = SimulationTime::Instance();
        archive & *p_time;
        // DO NOT archive & mpCell; -- The CellCycleModel is only ever archived from the Cell 
        // which knows this and it is handled in the load_construct of TissueCell.
    }
protected:
    TissueCell* mpCell;
    double mBirthTime; // Time to start model from
    
public:

    /**
     * Sets up a new AbstractCellCycleModel, gives it a birth time of the 
     * current simulation time (this is overwritten by some subclasses) 
     */
    AbstractCellCycleModel()
        : mpCell(NULL),
          mBirthTime(SimulationTime::Instance()->GetDimensionalisedTime()) {};

    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractCellCycleModel();
    
    virtual void SetCell(TissueCell* pCell);
    
    virtual void Initialise() {};
    
    TissueCell* GetCell();
    
    /**
     * Refreshes the cell's type using cell cycle information.
     * Doesn't do anything in most classes, overidden in others.
     */
    virtual void UpdateCellType() {};
           
    /**
     * Set the cell's time of birth (usually not required as it should be inside
     * the indivdual cell-cycle-model-constructor, but useful for tests)
     * 
     * @param birthTime the simulation time at this cell's birth.
     * 
     * (This function is overridden in AbstractOdeBasedCellCycleModel).
     */
    virtual void SetBirthTime(double birthTime);
    
    /**
     * Returns the cell's age...
     */
    double GetAge();
    
    /**
     * Returns the cell's birth time...
     */
    double GetBirthTime() const;
    
    /**
     * Determine whether the cell is ready to divide.
     * 
     * @param timeSinceBirth  the elapsed time since the cell was born
     */
    virtual bool ReadyToDivide()=0;
    
    /**
     * Each cell cycle model must be able to be reset after a cell division.
     */
    virtual void ResetModel()=0;
    
    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     * NB It should create an instance which is identical to the host instance.
     * 
     * This method is called in 2 circumstances:
     *  - By the copy constructor and operator= of TissueCell to create a copy of the cell cycle model
     *    when copying a cell. The CreateCellCycleModel just needs to create any instance of the right class,
     *    as operator= on the cell cycle model is then called to ensure the model is copied properly.
     *  - By TissueCell.Divide to create a cell cycle model for the daughter cell. CreateCellCycleModel
     *    must thus produce a cell cycle model in a suitable state for a newly-born cell spawned from the
     *    'current' cell. Note that the parent cell cycle model should be Reset() just before CreateCellCycleModel is
     *    called to copy its state. 
     */
    virtual AbstractCellCycleModel *CreateCellCycleModel()=0;
    
};


BOOST_IS_ABSTRACT(AbstractCellCycleModel)


#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
