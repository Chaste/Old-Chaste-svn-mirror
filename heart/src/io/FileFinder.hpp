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

#include "HeartConfig.hpp"

/**
 * A helper class for finding files, given paths which can be relative to various
 * locations (e.g. the Chaste source tree root, the current directory, the Chaste
 * test output directory, or an absolute path).
 */
class FileFinder
{
private:
    /** The absolute path to our file */
    std::string mAbsPath;
    
    /**
     * Determine the absolute path to this file.
     * Used by constructors.
     * 
     * @param rPath  the path to the file to find
     */
    void SetAbsolutePath(const cp::path_type& rPath);
public:
    /**
     * Create a file finder for the given path.
     * This type includes both a path name, and an attribute specifying how this
     * should be interpreted.  See the XML schema for details.
     * 
     * @param rPath  the path to the file to find
     */
    FileFinder(const cp::path_type& rPath);
    
    /**
     * Alternative constructor, taking the components of a cp::path_type.
     * @param rPath  the path to the file to find
     * @param relativeTo  how to interpret this path
     */
    FileFinder(const std::string& rPath, cp::relative_to_type relativeTo);
    
    /**
     * Test whether our file exists.
     */
    bool Exists() const;
    
    /**
     * Get the absolute path to this file.
     */
    std::string GetAbsolutePath() const;
};

#endif /*FILEFINDER_HPP_*/
