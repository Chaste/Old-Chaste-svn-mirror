
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
    // get the parameters using the method 'HeartPhysiologicalParameters(filename)',
    // which returns a std::auto_ptr. We don't want to use a std::auto_ptr because
    // it will delete memory when out of scope, or no longer point when it is copied,
    // so we reallocate memory using a normal pointer and copy the data to there
    std::auto_ptr<HeartPhysiologicalParametersType> x(HeartPhysiologicalParameters(fileName.c_str()));
    mpParameters = new HeartPhysiologicalParametersType(*x);
    assert(mpParameters);
}

HeartPhysiologicalParametersType* HeartParameters::Parameters()
{
    return mpParameters;
}

void HeartParameters::Destroy()
{
    delete mpParameters;
}
