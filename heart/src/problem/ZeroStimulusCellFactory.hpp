#ifndef ZEROSTIMULUSCELLFACTORY_HPP_
#define ZEROSTIMULUSCELLFACTORY_HPP_

#include "AbstractCardiacCellFactory.hpp"

template <class CELL, unsigned DIM>
class ZeroStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{

public:
    ZeroStimulusCellFactory(double timeStep=0.01) : AbstractCardiacCellFactory<DIM>(timeStep)
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new CELL(this->mpSolver, this->mTimeStep, this->mpZeroStimulus, this->mpZeroStimulus);
    }
    
    ~ZeroStimulusCellFactory(void)
    {
    }
};

#endif /*ZEROSTIMULUSCELLFACTORY_HPP_*/
