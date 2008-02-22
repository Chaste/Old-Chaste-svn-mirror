#ifndef BIDOMAINFACESTIMULUSCELLFACTORY_HPP_
#define BIDOMAINFACESTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "InitialStimulus.hpp"
#include "RegularStimulus.hpp"

class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
    RegularStimulus *mpRegStimulus;
    
public:
    //Pdetime step is (by default) 0.01
    //Odetime step set below to 0.001 (10:1)
    BidomainFaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.001)
    {
        mpRegStimulus = new RegularStimulus(-900.0*1000, 0.5, 100.0, 0.0);//Same as above, but every 100ms
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpRegStimulus, mpZeroStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus, mpZeroStimulus);
        }
    }
    
    ~BidomainFaceStimulusCellFactory(void)
    {
        delete mpRegStimulus;
    }
};

#endif /*BIDOMAINFACESTIMULUSCELLFACTORY_HPP_*/
