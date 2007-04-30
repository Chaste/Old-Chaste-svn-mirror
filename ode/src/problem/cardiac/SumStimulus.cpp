#include "SumStimulus.hpp"

SumStimulus::SumStimulus(AbstractStimulusFunction* pStimulus1, AbstractStimulusFunction* pStimulus2)
  : mpStimulus1(pStimulus1), mpStimulus2(pStimulus2)
{
}

double SumStimulus::GetStimulus(double time)
{
    return mpStimulus1->GetStimulus(time)+mpStimulus2->GetStimulus(time);
}
