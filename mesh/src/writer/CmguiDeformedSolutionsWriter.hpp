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


#ifndef CMGUIDEFORMEDSOLUTIONSWRITER_HPP_
#define CMGUIDEFORMEDSOLUTIONSWRITER_HPP_

#include "CmguiMeshWriter.hpp"
#include "QuadraticMesh.hpp"


/** 
 *  CmguiDeformedSolutionsWriter
 * 
 *  A class for writing a mesh, and solutions from a solid mechanics (ie deformed meshes) problem, 
 *  in Cmgui output format.
 * 
 *  Inherits from CmguiMeshWriter
 */
template<unsigned DIM>
class CmguiDeformedSolutionsWriter : public CmguiMeshWriter<DIM, DIM>
{
private:
    /** 
     *  The quadratic mesh used in the mechanics simulation. solution_0.exnode and solution_0.exelem
     *  (the only exelem file written) will be written using this mesh
     */
    QuadraticMesh<DIM>* mpQuadraticMesh;
    
    /** 
     *  A counter is given whenever WriteDeformationPositions() is called, this variable
     *  stores the last one used
     */
    unsigned mFinalCounter;

public:
    /**
     *  Constructor 
     *  @param outputDirectory The output directory for the Cmgui files
     *  @param baseName The base name for the Cmgui output files - the files written will be
     *   [basename_0.exnode, [basename]_0.exelem; [basename]_1.exnode, [basename]_2.exnode, ..
     *  @param rQuadraticMesh The quadratic mesh used in the mechanics simulation
     */
    CmguiDeformedSolutionsWriter(std::string outputDirectory,
                                 std::string baseName,
                                 QuadraticMesh<DIM>& rQuadraticMesh);
    
    /**
     *  Write [basename]_0.exnode, [basename]_0.exelem using the quadratic mesh
     */
    void WriteInitialMesh();
    
    /**
     *  Write [basename]_i.exnode using the given deformed positions
     *  @param rDeformedPositions std::vector of deformed positions to be used, must have size equal to number
     *  of nodes in the mesh
     *  @param counter the value "i" in "[basename]_i.exnode" to be used.
     */
    void WriteDeformationPositions(std::vector<c_vector<double,DIM> >& rDeformedPositions,
                                   unsigned counter);
    
    /**
     *  Writes a small cmgui script called LoadSolutions.com, for loading the output that has been written. 
     *  Assumes the output was solution_0.exnode .. solution_N.exnode, where N is the counter that was 
     *  given in the last call to WriteDeformationPositions()
     */
    void WriteCmguiScript();
    
    /**
     *  For a simulation that has already been run, convert the chaste output to cmgui format.
     *  @inputDirectory The directory the chaste output is in
     *  @inputFileBaseName The base name for the chaste output
     *  @finalCounter   The final counter, ie the value N for the final file <basename>_N.nodes.
     *  The first file is assumed to be <basename>_0.nodes. The files <basename>_0.nodes up to
     *  <basename>_N.nodes will be converted to cmgui format and put in the output directory
     *  given in the constructor.
     */   
    void ConvertOutput(std::string inputDirectory, 
                       std::string inputFileBaseName, 
                       unsigned finalCounter);
};



#endif /*CMGUIDEFORMEDSOLUTIONSWRITER_HPP_*/

