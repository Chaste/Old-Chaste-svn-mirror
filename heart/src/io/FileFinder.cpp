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
#include <fstream>
#include <sys/stat.h>

FileFinder::FileFinder(const cp::path_type& rPath)
{
    SetAbsolutePath(rPath);
}

FileFinder::FileFinder(const std::string& rPath, cp::relative_to_type relativeTo)
{
    cp::path_type path(rPath);
    path.relative_to(relativeTo);
    SetAbsolutePath(path);
}


void FileFinder::SetAbsolutePath(const cp::path_type& rPath)
{
    std::string leaf_path(rPath);

    switch (rPath.relative_to())
    {
        case cp::relative_to_type::chaste_source_root:
            mAbsPath = ChasteBuildRootDir() + leaf_path;
            break;

        case cp::relative_to_type::chaste_test_output:
            mAbsPath = OutputFileHandler::GetChasteTestOutputDirectory() + leaf_path;
            break;

        case cp::relative_to_type::cwd:
            mAbsPath = GetCurrentWorkingDirectory() + "/" + leaf_path;
            break;

        case cp::relative_to_type::absolute:
            mAbsPath = leaf_path;
            break;

        default:
            // Getting here is impossible due to the schema
            NEVER_REACHED;
            break;
    }
}

bool FileFinder::Exists() const
{
    std::ifstream file(mAbsPath.c_str());
    bool exists = file.is_open();
    file.close();

    return exists;
}

std::string FileFinder::GetAbsolutePath() const
{
    return mAbsPath;
}

bool FileFinder::IsNewerThan(const FileFinder& rOtherFile) const
{
    assert(Exists());
    assert(rOtherFile.Exists());
    struct stat our_stats, other_stats;
    stat(GetAbsolutePath().c_str(), &our_stats);
    stat(rOtherFile.GetAbsolutePath().c_str(), &other_stats);
    return our_stats.st_mtime > other_stats.st_mtime;
}
