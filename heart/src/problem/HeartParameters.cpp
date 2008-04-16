
#include "HeartParameters.hpp"

HeartParameters* HeartParameters::mpInstance = NULL;

HeartParameters* HeartParameters::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new HeartParameters;
    }
    return mpInstance;
}

HeartParameters::HeartParameters()
{
    assert(mpInstance == NULL);
}

void HeartParameters::SetParametersFile(std::string fileName)
{
    mpParameters = new std::auto_ptr<HeartPhysiologicalParametersType>(HeartPhysiologicalParameters(fileName.c_str()));
}

std::auto_ptr<HeartPhysiologicalParametersType> HeartParameters::Parameters()
{
    return *mpParameters;
}

void HeartParameters::Destroy()
{
    delete mpParameters;
}
