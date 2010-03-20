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

// Must be included before any other serialisation headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <sstream>
#include <iostream>

#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ProcessSpecificArchive.hpp"
#include "Exception.hpp"

/**
 * Specialization for input archives.
 * @param rDirectory
 * @param rFileName
 * @param relativeToChasteTestOutput
 * @param procId
 */
template<>
ArchiveOpener<boost::archive::text_iarchive, std::ifstream>::ArchiveOpener(
        const std::string& rDirectory,
        const std::string& rFileName,
        bool relativeToChasteTestOutput,
        unsigned procId)
    : mpCommonStream(NULL),
      mpPrivateStream(NULL),
      mpCommonArchive(NULL),
      mpPrivateArchive(NULL)
{
    // Figure out where things live
    if (relativeToChasteTestOutput)
    {
        OutputFileHandler handler(rDirectory, false);
        handler.SetArchiveDirectory();
    }
    else
    {
        ArchiveLocationInfo::SetArchiveDirectory(rDirectory, relativeToChasteTestOutput);
    }
    std::string private_path = ArchiveLocationInfo::GetProcessUniqueFilePath(rFileName, procId);
    std::stringstream common_path;
    common_path << ArchiveLocationInfo::GetArchiveDirectory() << rFileName;

    // Try to open the main archive for replicated data
    mpCommonStream = new std::ifstream(common_path.str().c_str(), std::ios::binary);
    if (!mpCommonStream->is_open())
    {
        delete mpCommonStream;
        EXCEPTION("Cannot load main archive file: " + common_path.str());
    }

    try
    {
        mpCommonArchive = new boost::archive::text_iarchive(*mpCommonStream);
    }
    catch (boost::archive::archive_exception& boost_exception)
    {
        if (boost_exception.code == boost::archive::archive_exception::unsupported_version)
        {
            // This is forward compatibility issue.  We can't open the archive because it's been written by a more recent Boost.
            delete mpCommonArchive;
            delete mpCommonStream;
            EXCEPTION("Could not open Boost archive '" + common_path.str() + "' because it was written by a more recent Boost.  Check process-specific archives too");
        }
        else
        {
        // We don't understand the exception, so we shouldn't continue
#define COVERAGE_IGNORE
            throw boost_exception;
#undef COVERAGE_IGNORE
        }
    }
    // Try to open the secondary archive for distributed data
    mpPrivateStream = new std::ifstream(private_path.c_str(), std::ios::binary);
    if (!mpPrivateStream->is_open())
    {
        delete mpPrivateStream;
        delete mpCommonArchive;
        delete mpCommonStream;
        EXCEPTION("Cannot load secondary archive file: " + private_path);
    }
    mpPrivateArchive = new boost::archive::text_iarchive(*mpPrivateStream);
    ProcessSpecificArchive<boost::archive::text_iarchive>::Set(mpPrivateArchive);
}

template<>
ArchiveOpener<boost::archive::text_iarchive, std::ifstream>::~ArchiveOpener()
{
    ProcessSpecificArchive<boost::archive::text_iarchive>::Set(NULL);
    delete mpPrivateArchive;
    delete mpPrivateStream;
    delete mpCommonArchive;
    delete mpCommonStream;
}

/**
 * Specialization for output archives.
 * @param rDirectory
 * @param rFileName
 * @param relativeToChasteTestOutput
 * @param procId
 */
template<>
ArchiveOpener<boost::archive::text_oarchive, std::ofstream>::ArchiveOpener(
        const std::string& rDirectory,
        const std::string& rFileName,
        bool relativeToChasteTestOutput,
        unsigned procId)
    : mpCommonStream(NULL),
      mpPrivateStream(NULL),
      mpCommonArchive(NULL),
      mpPrivateArchive(NULL)
{
    // Check for user error
    if (procId != PetscTools::GetMyRank())
    {
        EXCEPTION("Specifying the secondary archive file ID doesn't make sense when writing.");
    }
    // Figure out where things live
    if (relativeToChasteTestOutput)
    {
        OutputFileHandler handler(rDirectory, false);
        handler.SetArchiveDirectory();
    }
    else
    {
        ArchiveLocationInfo::SetArchiveDirectory(rDirectory, relativeToChasteTestOutput);
    }
    std::string private_path = ArchiveLocationInfo::GetProcessUniqueFilePath(rFileName);
    std::stringstream common_path;
    common_path << ArchiveLocationInfo::GetArchiveDirectory() << rFileName;

    // Create master archive for replicated data
    if (PetscTools::AmMaster())
    {
        mpCommonStream = new std::ofstream(common_path.str().c_str());
        if (!mpCommonStream->is_open())
        {
            delete mpCommonStream;
            EXCEPTION("Failed to open main archive file for writing: " + common_path.str());
        }
    }
    else
    {
        // Non-master processes need to go through the serialization methods,
        // but not write any data.
        mpCommonStream = new std::ofstream("/dev/null");
        #define COVERAGE_IGNORE
        if (!mpCommonStream->is_open())
        {
            delete mpCommonStream;
            EXCEPTION("Failed to open dummy archive file '/dev/null' for writing");
        }
        #undef COVERAGE_IGNORE
    }
    mpCommonArchive = new boost::archive::text_oarchive(*mpCommonStream);

    // Create secondary archive for distributed data
    mpPrivateStream = new std::ofstream(private_path.c_str());
    if (!mpPrivateStream->is_open())
    {
        delete mpPrivateStream;
        delete mpCommonArchive;
        delete mpCommonStream;
        EXCEPTION("Failed to open secondary archive file for writing: " + private_path);
    }
    mpPrivateArchive = new boost::archive::text_oarchive(*mpPrivateStream);
    ProcessSpecificArchive<boost::archive::text_oarchive>::Set(mpPrivateArchive);
}

template<>
ArchiveOpener<boost::archive::text_oarchive, std::ofstream>::~ArchiveOpener()
{
    ProcessSpecificArchive<boost::archive::text_oarchive>::Set(NULL);
    delete mpPrivateArchive;
    delete mpPrivateStream;
    delete mpCommonArchive;
    delete mpCommonStream;

    /* In a parallel setting, make sure all processes have finished writing before
     * continuing, to avoid nasty race conditions.
     * For example, many tests will write an archive then immediately read it back
     * in, which could easily break without this.
     */
    PetscTools::Barrier("~ArchiveOpener");
}



