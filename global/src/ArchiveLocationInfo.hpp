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
#ifndef ARCHIVELOCATIONINFO_HPP_
#define ARCHIVELOCATIONINFO_HPP_

#include <string>
#include <sstream>
#include <cassert>
#include <iostream>

#include "Exception.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"

/**
 * Mini-class to help with 'archiving' various classes that don't write their
 * data directly to the archive file.  They thus need to know information
 * about where the archive is being written to, in order to write their own
 * files into the same folder.  The main methods are GetArchiveDirectory and
 * SetArchiveDirectory.
 *
 * This functionality is used by the meshes, LinearSystem and HeartConfig.
 *
 * For the benefit of the meshes (and the cell_based code), there are also
 * shortcut methods SetMeshPathname and GetMeshFilename, allowing you to
 * specify the base file name for the mesh.  This is needed because the cell_based
 * code adds timestamp information to the file name.
 */
class ArchiveLocationInfo
{
private:

    /** Directory that archives are being written to. */
    static std::string mDirPath;

    /** Whether mDirPath is relative to CHASTE_TEST_OUTPUT
     *  true if the directory provided lives in CHASTE_TEST_OUTPUT, false if it's relative to CWD 
     */
    static bool mDirIsRelativeToChasteTestOutput;

    /** Mesh filename (relative to #mDirPath). */
    static std::string mMeshFilename;
    

public:
    /**
     * Set the location to write mesh files.
     * @param rDirectory  the directory to write to.
     * @param rFilename  the base name (minus extension) for the mesh files.
     */
    static void SetMeshPathname(const std::string& rDirectory, const std::string& rFilename)
    {
        SetArchiveDirectory(rDirectory);
        mMeshFilename = rFilename;
    }

    /**
     * Set the filename for mesh files.
     * @param rFilename  the base name (minus extension) for the mesh files, used to put on a timestamp.
     */
    static void SetMeshFilename(const std::string& rFilename)
    {
        mMeshFilename = rFilename;
    }

    /**
     * Get the filename that the mesh should be written to
     * @return mesh filename
     */
    static std::string GetMeshFilename()
    {
        if (mMeshFilename == "")
        {
            EXCEPTION("ArchiveLocationInfo::mMeshFilename has not been set");
        }
        return mMeshFilename;
    }

    /**
     * Get the directory that archives are being written to.
     * Will always end in a '/'.
     * @return full path to directory
     */
    static std::string GetArchiveDirectory()
    {
        if (mDirPath == "")
        {
            EXCEPTION("ArchiveLocationInfo::mDirPath has not been set");
        }
        return mDirPath;
    }

    /**
     * Get the full path to an output file which has a name unique to the current
     * process.  Useful for ensuring that each process writes to / reads from a
     * separate file when running in parallel.
     *
     * The path will have the form "path_to_output_dir/rFileName.process_rank"
     *
     * @param rFileName  the base file name
     * @param procId  the process id number (defaults to current process)
     * @return  a full path to the file for this process
     */
    static std::string GetProcessUniqueFilePath(const std::string& rFileName,
                                                unsigned procId=PetscTools::GetMyRank())
    {
        std::stringstream filepath;
        filepath << GetArchiveDirectory() << rFileName << "." << procId;
        return filepath.str();
    }

    /**
     * Set the directory that archives are being written to.
     *
     * @param rDirectory  the directory in question.
     * @param relativeToChasteTestOutput  Whether to convert directory to an absolute path using the
     *                      OutputFileHandler (defaults to true)
     */
    static void SetArchiveDirectory(const std::string& rDirectory, bool relativeToChasteTestOutput=true)
    {
        mDirPath = rDirectory;
        if (! ( *(mDirPath.end()-1) == '/'))
        {
            mDirPath = mDirPath + "/";
        }
        
        mDirIsRelativeToChasteTestOutput = relativeToChasteTestOutput;        
    }

    /**
     * Get the directory to which the archives are being written.
     * Remove CHASTE_TEST_OUTPUT prefix (assuming it exists)
     * Will always end in a '/'.
     * @return relative path to directory
     */
    static std::string GetArchiveRelativePath()
    {
        if (mDirIsRelativeToChasteTestOutput)
        {
            std::string chaste_output=OutputFileHandler::GetChasteTestOutputDirectory();
            //Find occurrence of CHASTE_TEST_OUTPUT in string
            std::string::size_type pos=mDirPath.find(chaste_output, 0);
            if (pos == std::string::npos)
            {
                EXCEPTION("Full path doesn't give a directory relative to CHASTE_TEST_OUTPUT");
            }
            assert(pos == 0); //Expect this as a prefix, not in the middle of the string
    
            return mDirPath.substr(chaste_output.length());
        }
        else
        {
            return mDirPath;
        }
    }
    
    /**
     *  Get wheter the directory provided is relative to CHASTE_TEST_OUTPUT
     *  @return true if the directory provided lives in CHASTE_TEST_OUTPUT, false if it's relative to CWD 
     */
    static bool GetIsDirRelativeToChasteTestOutput()
    {
        return mDirIsRelativeToChasteTestOutput;
    }
};

#endif /*ARCHIVELOCATIONINFO_HPP_*/
