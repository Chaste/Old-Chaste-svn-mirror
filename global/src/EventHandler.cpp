#include "EventHandler.hpp"

PetscEvent EventHandler::mPetscEvent[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
double EventHandler::mCpuTime[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
const char* EventHandler::EVENT_NAME[] = { "InMesh", "AssSys", "Ode", 
                                           "Comms", "AssRhs", "NeuBCs", "Ksp",
                                           "Output", "Total" };
