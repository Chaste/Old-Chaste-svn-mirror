#ifndef SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_
#define SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "CancerParameters.hpp"

/**
 *  Simple oxygen-based cell cycle model
 *
 *  TO DO: document this class
 */
 
class SimpleOxygenBasedCellCycleModel : public FixedCellCycleModel
{
private:

    /** Private constructor for creating an identical daughter cell */
    SimpleOxygenBasedCellCycleModel(double g1Duration)
        : FixedCellCycleModel(g1Duration) {};
        
    double mTimeSpentInG1Phase;   
    double mG1Duration; 
        
public:
    SimpleOxygenBasedCellCycleModel() 
        : FixedCellCycleModel(),
          mTimeSpentInG1Phase(0.0),
          mG1Duration(CancerParameters::Instance()->GetHepaOneCellG1Duration())
    {    
    }
    
    bool ReadyToDivide()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
        bool ready = false;
        
        // get cell's oxygen concentration, which must be positive
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell);
        assert(oxygen_concentration>=0);
                        
        double time_since_birth = GetAge();
               
        if (mpCell->GetCellType()==DIFFERENTIATED)
        {
            mCurrentCellCyclePhase = G_ZERO;   
        }
        else 
        {
            if ( GetAge() < p_params->GetMDuration() )
            {
                mCurrentCellCyclePhase = M;   
            }
            else if ( time_since_birth < p_params->GetMDuration() + mG1Duration )
            {
                mCurrentCellCyclePhase = G_ONE;
                 
                mG1Duration = mG1Duration + (1-oxygen_concentration)*SimulationTime::Instance()->GetTimeStep();
                mTimeSpentInG1Phase = mTimeSpentInG1Phase + SimulationTime::Instance()->GetTimeStep();  
            }
            else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration() )
            {
                mCurrentCellCyclePhase = S; 
            }
            else if ( time_since_birth < p_params->GetMDuration() + mG1Duration + p_params->GetSDuration()  + p_params->GetG2Duration())
            {
                mCurrentCellCyclePhase = G_TWO;   
            }
            else
            {
                ready = true;
                // mCurrentCellCyclePhase = M;
            }
        }
               
        return ready;
        
    }
    
    AbstractCellCycleModel* CreateCellCycleModel()
    {
        return new SimpleOxygenBasedCellCycleModel();
    }
};

// declare identifier for the serializer
BOOST_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
