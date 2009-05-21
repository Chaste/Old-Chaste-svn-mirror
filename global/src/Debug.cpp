
#include "Debug.hpp"

std::string FormDebugHead()
{
    std::string ret;
    if(PetscTools::NumProcs()==1)
    {
        std::string ret("DEBUG: ");
        return ret;
    }
    else
    {   
        std::stringstream stringstream;
        stringstream << "DEBUG: proc " << PetscTools::GetMyRank() << ": ";
        return stringstream.str();
    }
} 
