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


#define CHECK_SYSTEM(cmd) EXPECT0(system, cmd)


OutputFileHandler::OutputFileHandler(const std::string &rDirectory,
                                     bool cleanOutputDirectory)
{
    // Are we the master process?  Only the master should make any new directories
    mAmMaster = PetscTools::AmMaster();
    mDirectory = GetOutputDirectoryFullPath(rDirectory);
    //Is it a valid request for a directory?
    if (rDirectory != "" && rDirectory.find("..") == std::string::npos)
    {
        // Are we the master process?  Only the master should make any new directories
        if (cleanOutputDirectory && mAmMaster)
        {
            std::string directory_to_move_to = GetOutputDirectoryFullPath("last_cleaned_directory");
            IGNORE_RET(system, "rm -rf " + directory_to_move_to);
            // Re-create the special directory
            mkdir(directory_to_move_to.c_str(), 0775);
            CHECK_SYSTEM("mv " + mDirectory + " " + directory_to_move_to);
            //system(("rm -rf " + mDirectory).c_str());
            // Re-create the output directory
            mkdir(mDirectory.c_str(), 0775);
        }
//This not always collective so we can't have a barrier       PetscTools::Barrier();
        struct stat st;
        while (stat(mDirectory.c_str(), &st) != 0)
        {
            //Wait until directory becomes available
        }
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
        // Add a trailing slash if not already there
        if (! ( *(directory_root.end()-1) == '/'))
        {
            directory_root = directory_root + "/";
        }
    }

    return directory_root;
}


std::string OutputFileHandler::GetOutputDirectoryFullPath(const std::string& rDirectory)
{
    std::string directory_root = GetChasteTestOutputDirectory();
    std::string directory = directory_root + rDirectory;
    // Make sure it exists (ish)
    if (mAmMaster)
    {
        CHECK_SYSTEM("mkdir -p " + directory);
    }

    // Add a trailing slash if not already there
    if (! ( *(directory.end()-1) == '/'))
    {
        directory = directory + "/";
    }
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

bool OutputFileHandler::IsMaster()
{
    return mAmMaster;
}

void OutputFileHandler::SetArchiveDirectory()
{
    ArchiveLocationInfo::SetArchiveDirectory(GetOutputDirectoryFullPath());
}
