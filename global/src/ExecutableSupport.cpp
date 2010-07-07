/*

Copyright (C) University of Oxford, 2005-2010

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

#include "ExecutableSupport.hpp"

#include <iostream>

#include "CommandLineArguments.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"
#include "Version.hpp"
#include "OutputFileHandler.hpp"
#include <sys/utsname.h>

void ExecutableSupport::InitializePetsc(int* pArgc, char*** pArgv)
{
    // Store the arguments in case other code needs them
    CommandLineArguments::Instance()->p_argc = pArgc;
    CommandLineArguments::Instance()->p_argv = pArgv;
    // Initialise PETSc
    PETSCEXCEPT(PetscInitialize(pArgc, pArgv, PETSC_NULL, PETSC_NULL));
}

void ExecutableSupport::ShowCopyright()
{
    //Only show one copy of copyright/header
    if (PetscTools::AmMaster())
    {
        std::cout << "Copyright (C) University of Oxford, 2005-2010 \n\n\
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
}

void ExecutableSupport::ShowParallelLaunching()
{
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
}

void ExecutableSupport::WriteMachineInfoFile(std::string fileBaseName)
{
    OutputFileHandler out_file_handler("",false);
    std::stringstream file_name;
    file_name << fileBaseName << "." <<  PetscTools::GetMyRank();
    out_stream out_file = out_file_handler.OpenOutputFile(file_name.str());
    *out_file << "Process " << PetscTools::GetMyRank() << " of " 
        << PetscTools::GetNumProcs() << "." << std::endl << std::flush;
    
    
        
    struct utsname uts_info;
    uname(&uts_info);
    
    *out_file << "uname sysname  = " << uts_info.sysname << std::endl << std::flush;
    *out_file << "uname nodename = " << uts_info.nodename << std::endl << std::flush;
    *out_file << "uname release  = " << uts_info.release << std::endl << std::flush;
    *out_file << "uname version  = " << uts_info.version << std::endl << std::flush;
    *out_file << "uname machine  = " << uts_info.machine << std::endl << std::flush;
    char buffer[100];
    FILE * system_info;
    
    *out_file << "\nInformation on number and type of processors:\n";
    system_info = popen("grep ^model.name /proc/cpuinfo" ,"r");    
    while ( fgets(buffer, 100, system_info) != NULL )
    {
        *out_file << buffer;
    }    
    fclose(system_info);

    *out_file << "\nInformation on processor caches, in the same order as above:\n";
    system_info = popen("grep ^cache.size /proc/cpuinfo" ,"r");    
    while ( fgets(buffer, 100, system_info) != NULL )
    {
        *out_file << buffer;
    }    
    fclose(system_info);
    
    *out_file << "\nInformation on system memory:\n";
    system_info = popen("grep ^MemTotal /proc/meminfo" ,"r");    
    while ( fgets(buffer, 100, system_info) != NULL )
    {
        *out_file << buffer;
    }    
    fclose(system_info);
    
    out_file->close();
}

void ExecutableSupport::StandardStartup(int* pArgc, char*** pArgv)
{
    InitializePetsc(pArgc, pArgv);
    ShowCopyright();
    ShowParallelLaunching();
}

void ExecutableSupport::PrintError(const std::string& rMessage, bool masterOnly)
{
    if (!masterOnly || PetscTools::AmMaster())
    {
        std::cerr << rMessage << std::endl;
    }
}

void ExecutableSupport::FinalizePetsc()
{
    PetscFinalize();
}
