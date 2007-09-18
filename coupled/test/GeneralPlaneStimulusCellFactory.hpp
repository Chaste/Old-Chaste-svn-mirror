#ifndef GENERALPLANESTIMULUSCELLFACTORY_HPP_
#define GENERALPLANESTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"

template <class CELL, unsigned DIM>
class GeneralPlaneStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    // define a new stimulus
    InitialStimulus* mpStimulus;
    
public:
    GeneralPlaneStimulusCellFactory(double timeStep, double numElements) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
        // scale stimulus depending on space_step of elements
        //\todo It looks like the value of the stimulus is specific to 3D

        switch(DIM)
        {
            case 1:
            {
                mpStimulus = new InitialStimulus(-10000000*numElements/64.0, 0.5);
                break;
            }
            case 2:
            {
                mpStimulus = new InitialStimulus(-5000*numElements, 0.5);
                break;
            }
            case 3:
            {
                assert(0);
                break;
            }
            default:
                assert(0);
        }
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        double x = this->mpMesh->GetNode(node)->GetPoint()[0];
        if (x*x<=1e-10)
        {
            return new CELL(this->mpSolver, this->mTimeStep, this->mpStimulus, this->mpZeroStimulus);
        }
        else
        {
            return new CELL(this->mpSolver, this->mTimeStep, this->mpZeroStimulus, this->mpZeroStimulus);
        }
    }
    
    ~GeneralPlaneStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

#endif /*GENERALPLANESTIMULUSCELLFACTORY_HPP_*/
