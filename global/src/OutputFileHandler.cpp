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


#include "OutputFileHandler.hpp"

#include <cstdlib>
#include <sys/stat.h>

#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ArchiveLocationInfo.hpp"

OutputFileHandler::OutputFileHandler(const std::string &rDirectory,
                                     bool cleanOutputDirectory)
{
    //Is it a valid request for a directory?
    if (rDirectory.find("..") != std::string::npos)
    {
        EXCEPTION("Will not create directory: " + rDirectory +
                " due to it potentially being above, and cleaning, CHASTE_TEST_OUTPUT.");
    }

    mDirectory = MakeFoldersAndReturnFullPath(rDirectory);

    // Clean the directory (default)
    if (rDirectory != "" && cleanOutputDirectory) // Don't clean CHASTE_TEST_OUTPUT
    {
        std::string command = "test -e " + mDirectory + ".chaste_deletable_folder";
        int return_value = system(command.c_str());
        if (return_value!=0)
        {
            EXCEPTION("Cannot delete " + mDirectory + " because signature file \".chaste_deletable_folder\" is not present.");
        }

        // Are we the master process?  Only the master should delete files
        if (PetscTools::AmMaster())
        {
            //Remove whatever was there before
            //Note that the /* part prevents removal of hidden files (.filename), which is useful in NFS systems
            EXPECT0(system,"rm -rf " + mDirectory + "/*");
        }
        // Wait for master to finish before going on to use the directory.
        PetscTools::Barrier();
    }
}

std::string OutputFileHandler::GetChasteTestOutputDirectory()
{
    char *chaste_test_output = getenv("CHASTE_TEST_OUTPUT");
    std::string directory_root;
    if (chaste_test_output == NULL || *chaste_test_output == 0)
    {
        // Default to 'testoutput' folder within the current directory
        directory_root = "./testoutput/";
    }
    else
    {
        directory_root = std::string(chaste_test_output);
    }
    AddTrailingSlash(directory_root);

    return directory_root;
}


std::string OutputFileHandler::MakeFoldersAndReturnFullPath(const std::string& rDirectory)
{
    std::string directory_root = GetChasteTestOutputDirectory();
    std::string directory = directory_root + rDirectory;
    AddTrailingSlash(directory);

    // Are we the master process?  Only the master should make any new directories
    if (PetscTools::AmMaster())
    {
        // Re-create the output directory structure:
        std::string command = "test -d " + directory;
        int return_value = system(command.c_str());
        if (return_value!=0)
        {
            // We make as many folders as necessary here.
            EXPECT0(system,"mkdir -p " + directory);

            // Put the Chaste signature file in all folders we have created
            EXPECT0(system,"touch " + directory + ".chaste_deletable_folder");
        }
    }
    // Wait for master to finish before going on to use the directory.
    PetscTools::Barrier();

    return directory;
}


std::string OutputFileHandler::GetOutputDirectoryFullPath()
{
    return mDirectory;
}

out_stream OutputFileHandler::OpenOutputFile(const std::string& rFileName,
                                             std::ios_base::openmode mode)
{
    out_stream p_output_file(new std::ofstream((mDirectory+rFileName).c_str(), mode));
    if (!p_output_file->is_open())
    {
        EXCEPTION("Could not open file " + rFileName + " in " + mDirectory);
    }
    return p_output_file;
}


out_stream OutputFileHandler::OpenOutputFile(const std::string& rFileName,
                                             unsigned number,
                                             const std::string& rFileFormat,
                                             std::ios_base::openmode mode)
{
    std::stringstream string_stream;
    string_stream << rFileName << number << rFileFormat;
    return OpenOutputFile(string_stream.str(), mode);
}

void OutputFileHandler::SetArchiveDirectory()
{
    ArchiveLocationInfo::SetArchiveDirectory(GetOutputDirectoryFullPath());
}

void OutputFileHandler::AddTrailingSlash(std::string& rDirectory)
{
    // Add a trailing slash if not already there
    if (! ( *(rDirectory.end()-1) == '/'))
    {
        rDirectory = rDirectory + "/";
    }
}
