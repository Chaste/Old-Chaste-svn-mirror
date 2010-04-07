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


#ifndef INVENTORVORONOIWRITER_HPP_
#define INVENTORVORONOIWRITER_HPP_

#include <string>

#include "VoronoiTessellation.hpp"
#include "OutputFileHandler.hpp"

/**
 * A concrete InventorVoronoiWriter class. Used to write VoronoiTessellations to file.
 */
class InventorVoronoiWriter
{
protected:

    OutputFileHandler* mpOutputFileHandler; /**< Output file handler */
    std::string mBaseName; /**< Base name for the input files */

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    InventorVoronoiWriter(const std::string& rDirectory,
                          const std::string& rBaseName,
                          const bool clearOutputDir=true);

    /**
     * Destructor.
     */
    ~InventorVoronoiWriter();

    /**
     *  Write the Voronoi tessellation in Inventor format.
     *
     * @param rTessellation the Voronoi tessellation
     */
    void Write(const VoronoiTessellation<3>& rTessellation);

    /**
     *  Scale the vertex of each cell toward the centre of that cell by the given scaleFactor
     *  and write.
     *
     * @param rTessellation the Voronoi tessellation
     * @param scaleFactor the scale factor
     */
    void ScaleAndWrite(VoronoiTessellation<3>& rTessellation, double scaleFactor);

};


#endif /*INVENTORVORONOIWRITER_HPP_*/
