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

#ifndef FILEFINDER_HPP_
#define FILEFINDER_HPP_

#include <string>

/**
 * Structure encapsulating the enumeration of path 'types', i.e. what a path
 * can be relative to.  This allows us to write things like RelativeTo::ChasteTestOutput
 * for readability.
 */
struct RelativeTo
{
    /**
     * What things a path can be relative to.
     */
    enum Value
    {
        CWD, /**< The current working directory */
        ChasteTestOutput, /**< The CHASTE_TEST_OUTPUT folder */
        ChasteSourceRoot, /**< The root of the Chaste source tree */
        Absolute /**< This is an absolute path */
    };
};

/**
 * A helper class for finding files or directories, given paths which can be relative
 * to various locations (e.g. the Chaste source tree root, the current directory, the
 * Chaste test output directory, or an absolute path).
 */
class FileFinder
{
private:
    /** The absolute path to our file */
    std::string mAbsPath;

protected:
    /**
     * Determine the absolute path to this file/dir.
     * Used by constructor and HeartFileFinder.
     *
     * @param rRelativePath  the relative path to the file to find
     * @param relativeTo  what it's relative to
     */
    void SetAbsolutePath(const std::string& rRelativePath,
                         RelativeTo::Value relativeTo);

    /**
     * Default constructor for subclasses to use.  They @b must call
     * SetAbsolutePath() in their constructor.
     */
    FileFinder();
public:

    /**
     * Main constructor.
     * @param rPath  the path to the file/dir to find
     * @param relativeTo  how to interpret this path
     */
    FileFinder(const std::string& rPath, RelativeTo::Value relativeTo);

    /**
     * Test whether we exist.
     */
    bool Exists() const;
    
    /**
     * Are we pointing at a file?
     */
    bool IsFile() const;
    
    /**
     * Are we pointing at a directory?
     */
    bool IsDir() const;

    /**
     * Get the absolute path to this file/dir.
     */
    std::string GetAbsolutePath() const;

    /**
     * Test whether this file/dir is newer than another file/dir.
     * Compares modification times.
     * @param rOtherEntity  the entity to test against.
     */
    bool IsNewerThan(const FileFinder& rOtherEntity) const;
};

#endif /*FILEFINDER_HPP_*/
