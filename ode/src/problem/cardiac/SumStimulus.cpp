#include "SumStimulus.hpp"

SumStimulus::SumStimulus(AbstractStimulusFunction& Stimulus1, AbstractStimulusFunction& Stimulus2)
  : mpStimulus1(&Stimulus1), mpStimulus2(&Stimulus2)
{
}

double SumStimulus::GetStimulus(double time)
{
    return mpStimulus1->GetStimulus(time)+mpStimulus2->GetStimulus(time);
}
