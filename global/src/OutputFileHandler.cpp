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

#define CHECK_SYSTEM(cmd) EXPECT0(system, cmd)


OutputFileHandler::OutputFileHandler(const std::string &rDirectory,
                                     bool rCleanOutputDirectory)
{
    // Are we the master process?  Only the master should do any writing to disk
    mAmMaster = PetscTools::AmMaster();
    mDirectory = GetOutputDirectoryFullPath(rDirectory);

    // Clean the output dir?
    if (rCleanOutputDirectory && mAmMaster &&
        rDirectory != "" && rDirectory.find("..") == std::string::npos)
    {
        std::string directory_to_move_to = GetOutputDirectoryFullPath("last_cleaned_directory");
        CHECK_SYSTEM("rm -rf " + directory_to_move_to);
        // Re-create the special directory
        mkdir(directory_to_move_to.c_str(), 0775);
        CHECK_SYSTEM("mv " + mDirectory + " " + directory_to_move_to);
        //system(("rm -rf " + mDirectory).c_str());
        // Re-create the output directory
        mkdir(mDirectory.c_str(), 0775);
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


std::string OutputFileHandler::GetOutputDirectoryFullPath(std::string directory)
{
    std::string directory_root = GetChasteTestOutputDirectory();
    directory = directory_root + directory;
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


out_stream OutputFileHandler::OpenOutputFile(std::string filename,
                                             std::ios_base::openmode mode)
{
    out_stream p_output_file(new std::ofstream((mDirectory+filename).c_str(), mode));
    if (!p_output_file->is_open())
    {
        EXCEPTION("Could not open file " + filename + " in " + mDirectory);
    }
    return p_output_file;
}


out_stream OutputFileHandler::OpenOutputFile(std::string fileName,
                                             unsigned number,
                                             std::string fileFormat,
                                             std::ios_base::openmode mode)
{
    std::stringstream string_stream;
    string_stream << fileName << number << fileFormat;
    return OpenOutputFile(string_stream.str(), mode);
}

bool OutputFileHandler::IsMaster()
{
    return mAmMaster;
}
