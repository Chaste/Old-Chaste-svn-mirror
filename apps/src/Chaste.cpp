/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include <iostream>

// Most of the work is done by this class.  It must be included first.
#include "CardiacSimulation.hpp"

#include "Exception.hpp"
#include "PetscTools.hpp"
#include "Version.hpp"

int main(int argc, char *argv[])
{
    PETSCEXCEPT(PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL) );
    //Only show one copy of copyright/header
    if (PetscTools::AmMaster())
    {
        std::cout << "Copyright (C) University of Oxford, 2005-2009 \n\n\
\
Chaste is free software: you can redistribute it and/or modify \n\
it under the terms of the Lesser GNU General Public License as published by \n\
the Free Software Foundation, either version 2.1 of the License, or \n\
(at your option) any later version. \n\n\
\
Chaste is distributed in the hope that it will be useful, \n\
but WITHOUT ANY WARRANTY; without even the implied warranty of \n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n\
Lesser GNU General Public License for more details. \n\n\
\
You should have received a copy of the Lesser GNU General Public License \n\
along with Chaste.  If not, see <http://www.gnu.org/licenses/>.\n\n";

        //Compilation information
        std::cout << "This version of Chaste was compiled on:\n";
        std::cout << ChasteBuildInfo::GetBuildTime() << " by " << ChasteBuildInfo::GetBuilderUnameInfo() << " (uname)\n";
        std::cout << "from revision number " << ChasteBuildInfo::GetRevisionNumber() << " with build type " << ChasteBuildInfo::GetBuildInformation() << ".\n\n";
    }

    if (!PetscTools::IsSequential())
    {
        ///Information to show that Chaste is being run in parallel
        for (unsigned i=0; i<PetscTools::GetNumProcs(); i++)
        {
            if (i==PetscTools::GetMyRank())
            {
                std::cout << "Chaste launched on process " << PetscTools::GetMyRank()
                    << " of " << PetscTools::GetNumProcs() << "." << std::endl << std::flush;
            }
            PetscTools::Barrier();
        }
    }    
    try
    {
        if (argc<2)
        {
            if (PetscTools::AmMaster())
            {
                std::cout << "Usage: Chaste parameters_file\n";
            }
            return -1;
        }
        std::string xml_file_name(argv[1]);
        // Creates & runs the simulation
        CardiacSimulation simulation(xml_file_name);

        return 0;
    }
    catch (Exception& e)
    {
        std::cout << e.GetMessage() << std::endl;
        return 1;
    }

    PetscFinalize();
}
