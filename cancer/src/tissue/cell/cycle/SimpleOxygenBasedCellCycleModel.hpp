#ifndef SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_
#define SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "CancerParameters.hpp"
#include "CellwiseData.cpp"

/**
 *  Simple oxygen-based cell cycle model
 *
 *  TO DO: document this class
 */
 
class SimpleOxygenBasedCellCycleModel : public FixedCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mTimeSpentInG1Phase;
        archive & mG1Duration;   
    }
    
    /** Private constructor for creating an identical daughter cell */
    SimpleOxygenBasedCellCycleModel(double g1Duration)
        : FixedCellCycleModel(g1Duration) {};
        
    double mTimeSpentInG1Phase;   
    double mG1Duration; 
        
public:
    SimpleOxygenBasedCellCycleModel() 
        : FixedCellCycleModel(),
          mTimeSpentInG1Phase(0.0)
    {   
    }
    
    double GetG1Duration()
    {
        return mG1Duration;
    }
    
    bool ReadyToDivide()
    {
        CancerParameters *p_params = CancerParameters::Instance();
        
        // This is rather messy - we want to set the G1 duration at the 
        // earliest opportunity, which will be the first time 
        // ReadyToDivide is called.
        if (mG1Duration==DBL_MAX)
        {   
            SetG1Duration();
        }
        
        bool ready = false;
        
        // get cell's oxygen concentration
        double oxygen_concentration = CellwiseData<2>::Instance()->GetValue(mpCell);
        
        // we want the oxygen concentration to be positive,
        // to within numerical tolerances (hence the -1e-8)
//        assert(oxygen_concentration >= -1e-8);
                        
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
                 
                mG1Duration = mG1Duration + (1-std::max(oxygen_concentration,0.0))*SimulationTime::Instance()->GetTimeStep();
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
