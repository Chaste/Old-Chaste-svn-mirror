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
    GeneralPlaneStimulusCellFactory(double timeStep, double meshSize, bool useMeshSizeAsMag=false) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
        // scale stimulus depending on space_step of elements
        //\todo It looks like the value of the stimulus is specific to 3D
        if (useMeshSizeAsMag)
        {
            mpStimulus = new InitialStimulus(meshSize, 0.5);
        }
        else
        {
            switch(DIM)
            {
                case 1:
                {
                    mpStimulus = new InitialStimulus(-1e7*meshSize/64.0, 0.5);
                    //wrt mesh 4 which has 64 elements in 1D
                    //Justification: elements go half size with each refinement
                    break;
                }
                case 2:
                {
                    mpStimulus = new InitialStimulus(-1e7*meshSize/64.0, 0.5);
                    //wrt mesh which which has 64 elements across in 2D
                    //Justification: Triangles go quarter size with each refinement, but there are twice as many nodes on 
                    break;
                }
                default:
                    assert(0);
            }
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
