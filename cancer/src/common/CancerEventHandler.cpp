#include "CancerEventHandler.hpp"

double CancerEventHandler::mCpuTime[] = { 0.0 };

const char* CancerEventHandler::EVENT_NAME[] = { "Setup", "Death", "Birth", 
                                           "Remesh", "Tess", "Vel",
                                           "Pos", "Output", "Total" };

void CancerEventHandler::BeginEvent(CancerEventType event)
{
    mCpuTime[event]-= clock()/(CLOCKS_PER_SEC/1000.0); 
}

void CancerEventHandler::EndEvent(CancerEventType event)
{
    mCpuTime[event]+= clock()/(CLOCKS_PER_SEC/1000.0);
}    

void CancerEventHandler::Report()
{
    // Times are in milliseconds
    for (unsigned event=0; event<NUM_EVENTS; event++)
    {
        printf("%2.1e\t", mCpuTime[event]);
        mCpuTime[event]=0.0;
    }
    std::cout << "(milliseconds) \n";
}

void CancerEventHandler::Headings()
{
    for (unsigned event=0; event<NUM_EVENTS; event++)
    {
        printf("%6s\t", EVENT_NAME[event]);
    } 
    std::cout << "\n";
}
