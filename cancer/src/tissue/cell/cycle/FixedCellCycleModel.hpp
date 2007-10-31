#ifndef FIXEDCELLCYCLEMODEL_HPP_
#define FIXEDCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "CancerParameters.hpp"

/**
 *  Fixed cell cycle model
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  CancerParameters::StemCellG1Duration + SG2MDuration 
 *  and CancerParameters::TransitCellG1Duration + SG2MDuration)
 */
class FixedCellCycleModel : public AbstractCellCycleModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    	archive & mG1Duration;
    }
    
    
protected:    

    /** Protected constructor for creating an identical daughter cell */
    FixedCellCycleModel(double g1Duration)
        :mG1Duration(g1Duration) {};
    
    /**
     * The duration of the G1 phase of the cell cycle. This is set once a cell
     * cycle model has been told what cell it belongs to.
     */    
    double mG1Duration;
    
    virtual void SetG1Duration();
    
public:

    /**
     * Default constructor - mBirthTime now set in AbstractCellCycleModel().
     */
    FixedCellCycleModel() :
        mG1Duration(DBL_MAX) {};
    
    /** 
     * Overridden SetCellMethod - also assigns a G1 duration based on the cell type.
     * 
     * @param pCell the cell which owns this model.
     */
    void SetCell(TissueCell* pCell);

    virtual bool ReadyToDivide();
    
    virtual void ResetModel();
    
    AbstractCellCycleModel *CreateCellCycleModel(); 
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(FixedCellCycleModel)


#endif /*FIXEDCELLCYCLEMODEL_HPP_*/
