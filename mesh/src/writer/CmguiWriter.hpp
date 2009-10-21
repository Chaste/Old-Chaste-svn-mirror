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

#ifndef CMGUIWRITER_HPP_
#define CMGUIWRITER_HPP_

#include "AbstractTetrahedralMeshWriter.hpp"

/**
 * Header for node base file (.exnode)
 */
static const char CmguiNodeFileHeader[] = " #Fields=1\n\
 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  Value index= 1, #Derivatives= 0\n\
   y.  Value index= 2, #Derivatives= 0\n\
   z.  Value index= 3, #Derivatives= 0\n";

/**
 * Header for element base file (.exelem)
 */
static const char CmguiElementFileHeader[] = "Shape.  Dimension=3, simplex(2;3)*simplex*simplex\n\
 #Scale factor sets= 0\n\
 #Nodes= 4\n";
 
/**
 * Header for element base file (.exelem), this coems after the definition of the number of fields
 */
static const char CmguiCoordinatesFileHeader[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n";
       
/**
 * Header for additional fields in the element base file (.exelem), 
 * Here we assume all additional fields will be interpolated by cmgui in the same way
 */
static const char CmguiAdditonalFieldHeader[] = " field, rectangular cartesian, #Components=1\n\
   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n";
/**
 *  CmguiWriter
 *
 *  Writes a mesh in Cmgui (the visualisation frontend of CMISS) format. Creates an exnode
 *  file and a exelem file. Note that the lines and faces are not written in the exelem
 *  file, so to load the data in Cmgui, you must use 'generate_faces_and_lines', i.e.
 *
 *  gfx read node base_file
 *  gfx read elem base_file generate_faces_and_lines
 *  gfx cr win
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CmguiWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM,SPACE_DIM>
{
private:
    
    /**
     * For storage of names of additional fields.
     */
    std::vector<std::string> mAdditionalFieldNames;
    
public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param rCleanDirectory  whether to clean the directory (defaults to true) \todo make this cleanDirectory to be consistent with other writers? (#991)
     */
    CmguiWriter(const std::string& rDirectory,
                const std::string& rBaseName,
                const bool& rCleanDirectory=true);

    /**
     * Write mesh data to files.
     */
    void WriteFiles();
    
    /**
     * Set any additional field that we want cmgui to visualize (interpolated over) elements and surfaces
     * 
     * @param rFieldNames is a reference to a vector of string containing the names of each additional field name
     */
    void SetAdditionalFieldNames(std::vector<std::string>& rFieldNames);

    /**
     * Destructor.
     */
    virtual ~CmguiWriter()
    {}
};

#endif /*CMGUIWRITER_HPP_*/
