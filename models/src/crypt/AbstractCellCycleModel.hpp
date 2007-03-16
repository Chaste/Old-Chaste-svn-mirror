#ifndef ABSTRACTCELLCYCLEMODEL_HPP_
#define ABSTRACTCELLCYCLEMODEL_HPP_

//We want independence from archive type
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/access.hpp>
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
        archive & mpSimulationTime;
    }
protected:
    CryptCellType mCellType;
    double mBirthTime; // Time to start model from
    SimulationTime* mpSimulationTime;
    
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
     */
    virtual AbstractCellCycleModel *CreateCellCycleModel()=0;

};


// Avoid compiler errors on some systems
BOOST_IS_ABSTRACT(AbstractCellCycleModel)


#endif /*ABSTRACTCELLCYCLEMODEL_HPP_*/
