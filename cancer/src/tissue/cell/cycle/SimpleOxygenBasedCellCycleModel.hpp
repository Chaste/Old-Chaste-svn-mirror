#ifndef SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_
#define SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "CellwiseData.cpp"

/**
 *  Simple oxygen-based cell cycle model
 *
 *  A simple oxygen-dependent cell cycle model that inherits from 
 *  AbstractSimpleCellCycleModel. The duration of G1 phase depends 
 *  on the local oxygen concentration. A prolonged period of acute
 *  hypoxia leads to the cell being labelled as necrotic. This model
 *  allows for quiescence imposed by transient periods of hypoxia, 
 *  followed by reoxygenation.
 *  
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
    
    /**
     * The time spent in G1 phase so far
     */ 
    double mTimeSpentInG1Phase;
    
    /**
     * How long the current period of hypoxia has lasted
     */    
    double mCurrentHypoxicDuration;
    
    /*
     * The time when the current period of hypoxia began
     */ 
    double mCurrentHypoxiaOnsetTime;
    
    /** 
     * Private constructor for creating an identical daughter cell
     */
    SimpleOxygenBasedCellCycleModel(double g1Duration,
									unsigned generation,
                                    double currentHypoxicDuration,
                                    double currentHypoxiaOnsetTime)
        : AbstractSimpleCellCycleModel(g1Duration,generation),
          mTimeSpentInG1Phase(0.0),
          mCurrentHypoxicDuration(currentHypoxicDuration),
          mCurrentHypoxiaOnsetTime(currentHypoxiaOnsetTime) {};
                
public:

    /**
     * Constructor
     */ 
    SimpleOxygenBasedCellCycleModel();
    
    /** 
     * Overridden UpdateCellCyclePhase() method
     */ 
    void UpdateCellCyclePhase();
    
    /**
     * Method for updating mCurrentHypoxicDuration, 
     * called at the start of ReadyToDivide()
     */ 
    void UpdateHypoxicDuration();
     
    double GetCurrentHypoxicDuration();  
    
    double GetCurrentHypoxiaOnsetTime();  
    
    AbstractCellCycleModel* CreateDaughterCellCycleModel();
    
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
