#include "MultiStimulus.hpp"

void MultiStimulus::AddStimulus(AbstractStimulusFunction* pStimulus)
{
    mStimuli.push_back(pStimulus);
}

double MultiStimulus::GetStimulus(double time)
{
    double total_stimulus = 0.0;
    
    for (unsigned current_stimulus = 0; current_stimulus < mStimuli.size(); ++current_stimulus)
    {
        total_stimulus += mStimuli[current_stimulus]->GetStimulus(time);
    }
    
    return total_stimulus;
}
