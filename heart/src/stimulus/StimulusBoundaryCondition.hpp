#ifndef STIMULUSBOUNDARYCONDITION_HPP_
#define STIMULUSBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"
#include "AbstractStimulusFunction.hpp"
#include "PdeSimulationTime.hpp"

/**
 * Boundary condition defined by an AbstractStimlus object.
 */
template<unsigned SPACE_DIM>
class StimulusBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
private:
    AbstractStimulusFunction* mpStimulus;
    
public:
    /**
     * Create a new boundary condition object.
     * 
     * @param pStimulus Stimulus object defining the parameters of the boundary condition
     */
    StimulusBoundaryCondition(AbstractStimulusFunction* pStimulus)
    {
        mpStimulus = pStimulus;
    }
    
    /**
     * @param x The point at which this boundary condition is to be evaluated.
     * @return The constant value given in the constructor.
     */
    double GetValue( const ChastePoint<SPACE_DIM>& ) const
    {
        return mpStimulus->GetStimulus(PdeSimulationTime::GetTime());
    }
};

#endif /*STIMULUSBOUNDARYCONDITION_HPP_*/
