#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/base_object.hpp>

#include "MeinekeCryptCellTypes.hpp"
#include "SimulationTime.hpp"
#include <vector>

// Needs to be included last
#include <boost/serialization/export.hpp>

class AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCellType;
        archive & mBirthTime;
        // Make sure the simulation time gets saved too
        SimulationTime* p_time = SimulationTime::Instance();
        archive & *p_time;
    }
protected:
    CryptCellType mCellType;
    double mBirthTime; // Time to start model from
    
public:
    /**
     * Base class with virtual methods needs a virtual destructor.
     */
    virtual ~AbstractCellCycleModel();
    
    /**
     * Set the cell's type.
     * 
     * @param cellType the type of cell defined in MeinekeCryptCellTypes.hpp
     */
    void SetCellType(CryptCellType cellType);
    
    /**
     * Refreshes the cell's type using cell cycle information.
     * 
     * @return CellType
     */
    virtual CryptCellType UpdateCellType();
    
    /**
     * Set the cell's type.
     * 
     * @param cellType the type of cell defined in MeinekeCryptCellTypes.hpp
     */
    CryptCellType GetCellType();
    
    
    /**
     * Set the cell's time of birth (usually not required as it should be inside
     * the indivdual cell-cycle-model-constructor, but useful for tests)
     * 
     * @param birthTime the simulation time at this cell's birth.
     */
    virtual void SetBirthTime(double birthTime)=0;
    
    /**
     * Returns the cell's age...
     */
    double GetAge();
    
    /**
     * Returns the cell's birth time...
     */
    double GetBirthTime();
    
    /**
     * Determine whether the cell is ready to divide.
     * 
     * @param timeSinceBirth  the elapsed time since the cell was born
     */
    virtual bool ReadyToDivide(std::vector<double> cellCycleInfluences = std::vector<double>())=0;
    
    /**
     * Each cell cycle model must be able to be reset after a cell division.
     */
    virtual void ResetModel()=0;
    
    /**
     * Builder method to create new instances of the cell cycle model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     * 
     * This method is called in 2 circumstances:
     *  - By the copy constructor and operator= of MeinekeCryptCell to create a copy of the cell cycle model
     *    when copying a cell. The CreateCellCycleModel just needs to create any instance of the right class,
     *    as operator= on the cell cycle model is then called to ensure the model is copied properly.
     *  - By MeinekeCryptCell.Divide to create a cell cycle model for the daughter cell. CreateCellCycleModel
     *    must thus produce a cell cycle model in a suitable state for a newly-born cell spawned from the
     *    'current' cell. Note that the parent cell cycle model is reset just before CreateCellCycleModel is
     *    called. 
     */
    virtual AbstractCellCycleModel *CreateCellCycleModel()=0;
    
};


// Avoid compiler errors on some systems
// Seems to break things on the Chaste machines, unfortunately :(
BOOST_IS_ABSTRACT(AbstractCellCycleModel)


#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
