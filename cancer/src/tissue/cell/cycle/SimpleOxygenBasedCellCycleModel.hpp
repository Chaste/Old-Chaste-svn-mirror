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
        archive & mCurrentHypoxicDuration;
        archive & mCurrentHypoxiaOnsetTime;
    }
    
    double mTimeSpentInG1Phase;   
    double mCurrentHypoxicDuration;
    double mCurrentHypoxiaOnsetTime;
    
    /** Private constructor for creating an identical daughter cell */
    SimpleOxygenBasedCellCycleModel(double g1Duration,
									unsigned generation,
                                    double currentHypoxicDuration,
                                    double currentHypoxiaOnsetTime)
        : AbstractSimpleCellCycleModel(g1Duration,generation),
          mTimeSpentInG1Phase(0.0),
          mCurrentHypoxicDuration(currentHypoxicDuration),
          mCurrentHypoxiaOnsetTime(currentHypoxiaOnsetTime) {};
                
public:
    SimpleOxygenBasedCellCycleModel();
    
    bool ReadyToDivide();
    
    void UpdateHypoxicDuration();
    
    double GetCurrentHypoxicDuration();  
    
    double GetCurrentHypoxiaOnsetTime();  
    
    AbstractCellCycleModel* CreateCellCycleModel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
