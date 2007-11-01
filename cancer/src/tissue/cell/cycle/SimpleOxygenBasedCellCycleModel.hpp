#ifndef SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_
#define SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "CellwiseData.cpp"

/**
 *  Simple oxygen-based cell cycle model
 *
 *  \todo: document this class
 */ 
class SimpleOxygenBasedCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        archive & mTimeSpentInG1Phase;
    }
    
    /** Private constructor for creating an identical daughter cell */
    SimpleOxygenBasedCellCycleModel(double g1Duration)
        : AbstractSimpleCellCycleModel(g1Duration) {};
        
    double mTimeSpentInG1Phase;   
        
public:
    SimpleOxygenBasedCellCycleModel() :
          mTimeSpentInG1Phase(0.0) {};
    
    bool ReadyToDivide();
    
    AbstractCellCycleModel* CreateCellCycleModel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
