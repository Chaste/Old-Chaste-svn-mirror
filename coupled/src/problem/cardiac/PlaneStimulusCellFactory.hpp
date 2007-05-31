#ifndef PLANESTIMULUSCELLFACTORY_HPP_
#define PLANESTIMULUSCELLFACTORY_HPP_

#include "LuoRudyIModel1991OdeSystem.hpp"
#include "AbstractCardiacCellFactory.hpp"

template<unsigned DIM>
class PlaneStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    PlaneStimulusCellFactory() : AbstractCardiacCellFactory<DIM>(0.01)//Ode timestep
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(-600, 0.5);
    }
    
    PlaneStimulusCellFactory(double timeStep, double stimulusVoltage) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
        // set the new stimulus
        mpStimulus = new InitialStimulus(stimulusVoltage, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (this->mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver,
                                                  this->mTimeStep,
                                                  mpStimulus,
                                                  this->mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(this->mpSolver,
                                                  this->mTimeStep,
                                                  this->mpZeroStimulus,
                                                  this->mpZeroStimulus);
        }
    }
    
    ~PlaneStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};


#endif /*PLANESTIMULUSCELLFACTORY_HPP_*/
