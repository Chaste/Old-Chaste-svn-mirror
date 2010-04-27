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

#ifndef _GENERICMESHREADER_HPP_
#define _GENERICMESHREADER_HPP_

#include "AbstractMeshReader.hpp"
//Possible subclasses to delegate
#include "TrianglesMeshReader.hpp"
#include "MemfemMeshReader.hpp"

/**
 * A generic mesh reader
 * Uses a delegated member variable of type AbstractMeshReader to probe for files
 * which may be read via
 *  - TrianglesMeshReader
 *  - MemfemMeshReader
 * 
 * Probing is done during construction.
 * 
 * Thereafter all public methods use the public methods of the delegated class and 
 * exceptions are passed back to the caller without being caught locally
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class GenericMeshReader : public AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Delegated mesh reader.
     *  Used to probe various types of mesh reader
     */
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>* mpMeshReader;
public:
    /**
     * Constructor.
     *
     * @param pathBaseName  the base name of the files from which to read the mesh data
     */
    GenericMeshReader(std::string pathBaseName)
    {
        try
        {
            mpMeshReader=new TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>(pathBaseName);
        }
        catch (Exception triangles_exception)
        {
            try
            {
                mpMeshReader=new MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>(pathBaseName);
            }
            catch (Exception memfem_exception)
            {
                std::string combined_message = "Could not open appropriate mesh files for "+pathBaseName+"\n";
                combined_message += "Triangle format: "+triangles_exception.GetShortMessage()+"\n"; 
                combined_message += "Memfem format: "+memfem_exception.GetShortMessage()+"\n";
                EXCEPTION(combined_message);
            }
        }
        
    }
    
    
    /**
     * Destructor
     */
    ~GenericMeshReader()
    {
        delete mpMeshReader;
    }


    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline unsigned GetNumNodes() const
    {
        return mpMeshReader->GetNumNodes();
    }

    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline unsigned GetNumElements() const
    {
        return mpMeshReader->GetNumElements();
    }

    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline unsigned GetNumFaces() const
    {
        return mpMeshReader->GetNumFaces();
    }
    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline unsigned GetNumElementAttributes() const
    {
        return mpMeshReader->GetNumElementAttributes();
    }

/** \todo #1323 COVER
 */

//    /** 
//     * Method uses the public method of the delegated mesh reader
//     */
//    inline unsigned GetNumFaceAttributes() const
//    {
//        return mpMeshReader->GetNumFaceAttributes();
//    }

    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline void Reset()
    {
        mpMeshReader->Reset();
    }
    
    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline std::vector<double> GetNextNode()
    {
        return mpMeshReader->GetNextNode();
    }

    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline ElementData GetNextElementData()
    {
        return mpMeshReader->GetNextElementData();
    }
    /** 
     * Method uses the public method of the delegated mesh reader
     */
    inline ElementData GetNextFaceData()
    {
        return mpMeshReader->GetNextFaceData();
    }

  
};

#endif //_GENERICMESHREADER_HPP_
