
#include "DebugHelper.hpp"

/* This is a noddy file just to have some code (ie a source file that
 *  isn't templated in the Dealii folder, else the dealii library will 
 *  empty
 */
 
void DebugHelper::Print(std::string message)
{
    std::cout << message << std::endl << std::flush;
}

