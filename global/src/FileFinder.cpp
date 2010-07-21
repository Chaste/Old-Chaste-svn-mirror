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

#include "FileFinder.hpp"

#include "ChasteBuildRoot.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "PetscTools.hpp"
#include <fstream>
#include <cassert>
#include <sys/stat.h>


FileFinder::FileFinder()
    : mAbsPath("UNSET!")
{
}

FileFinder::FileFinder(const std::string& rRelativePath, RelativeTo::Value relativeTo)
{
    SetAbsolutePath(rRelativePath, relativeTo);
}


void FileFinder::SetAbsolutePath(const std::string& rRelativePath, RelativeTo::Value relativeTo)
{
    switch (relativeTo)
    {
        case RelativeTo::ChasteSourceRoot:
            mAbsPath = ChasteBuildRootDir() + rRelativePath;
            break;

        case RelativeTo::ChasteTestOutput:
            mAbsPath = OutputFileHandler::GetChasteTestOutputDirectory() + rRelativePath;
            break;

        case RelativeTo::CWD:
            mAbsPath = GetCurrentWorkingDirectory() + "/" + rRelativePath;
            break;

        case RelativeTo::Absolute:
            mAbsPath = rRelativePath;
            break;

        default:
            // Getting here is impossible
            NEVER_REACHED;
            break;
    }
    
    if (IsDir())
    {
        if (*(mAbsPath.end()-1) != '/')
        {
            mAbsPath = mAbsPath + "/";
        }
    }
}

bool FileFinder::Exists() const
{
    struct stat our_stats;
    int retcode = stat(GetAbsolutePath().c_str(), &our_stats);
    return (retcode == 0);
}

bool FileFinder::IsFile() const
{
    bool result;
    struct stat our_stats;
    int retcode = stat(GetAbsolutePath().c_str(), &our_stats);
    if (retcode == 0)
    {
        result = S_ISREG(our_stats.st_mode);
    }
    else
    {
        // If it doesn't exist, it isn't a file
        result = false;
    }
    return result;
}

bool FileFinder::IsDir() const
{
    bool result;
    struct stat our_stats;
    int retcode = stat(GetAbsolutePath().c_str(), &our_stats);
    if (retcode == 0)
    {
        result = S_ISDIR(our_stats.st_mode);
    }
    else
    {
        // If it doesn't exist, it isn't a directory
        result = false;
    }
    return result;
}

std::string FileFinder::GetAbsolutePath() const
{
    return mAbsPath;
}

bool FileFinder::IsNewerThan(const FileFinder& rOtherEntity) const
{
    assert(Exists());
    assert(rOtherEntity.Exists());
    struct stat our_stats, other_stats;
    stat(GetAbsolutePath().c_str(), &our_stats);
    stat(rOtherEntity.GetAbsolutePath().c_str(), &other_stats);
    return our_stats.st_mtime > other_stats.st_mtime;
}

bool FileFinder::IsAbsolutePath(const std::string& rPath)
{
    return rPath[0]=='/';
}
