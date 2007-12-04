#include "CancerEventHandler.hpp"

PetscEvent CancerEventHandler::mPetscEvent[] = { 0 };
double CancerEventHandler::mCpuTime[] = { 0.0 };
const char* CancerEventHandler::EVENT_NAME[] = { "Setup", "Death", "Birth", 
                                           "Remesh", "Tess", "Vel",
                                           "Pos", "Output", "Total" };
