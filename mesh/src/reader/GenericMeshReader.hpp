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
#include "VtkMeshReader.hpp"

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
     *    (either absolute, or relative to the current directory)
     */
    GenericMeshReader(std::string pathBaseName)
    {
        try
        {
            mpMeshReader = new TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>(pathBaseName);
        }
        catch (const Exception& r_triangles_exception)
        {
            try
            {
                mpMeshReader = new MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>(pathBaseName);
            }
            catch (const Exception& r_memfem_exception)
            {
#ifdef CHASTE_VTK
                try
                {
                    mpMeshReader = new VtkMeshReader<ELEMENT_DIM, SPACE_DIM>(pathBaseName);
                }
                catch (const Exception& r_vtk_exception)
                {
#endif // CHASTE_VTK
                    std::string eol("\n");
                    std::string combined_message = "Could not open appropriate mesh files for " + pathBaseName + eol;
                    combined_message += "Triangle format: " + r_triangles_exception.GetShortMessage() + eol;
                    combined_message += "Memfem format: " + r_memfem_exception.GetShortMessage() + eol;
#ifdef CHASTE_VTK
                    combined_message += "Vtk format: " + r_vtk_exception.GetShortMessage() + eol;
#endif // CHASTE_VTK
                    EXCEPTION(combined_message);
#ifdef CHASTE_VTK
                }
#endif // CHASTE_VTK
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
    unsigned GetNumNodes() const
    {
        return mpMeshReader->GetNumNodes();
    }

    /**
     * Method uses the public method of the delegated mesh reader
     */
    unsigned GetNumElements() const
    {
        return mpMeshReader->GetNumElements();
    }

    /**
     * Method uses the public method of the delegated mesh reader
     */
    unsigned GetNumFaces() const
    {
        return mpMeshReader->GetNumFaces();
    }
    /**
     * Method uses the public method of the delegated mesh reader
     */
    unsigned GetNumElementAttributes() const
    {
        return mpMeshReader->GetNumElementAttributes();
    }

    /**
     * Method uses the public method of the delegated mesh reader
     */
    unsigned GetNumFaceAttributes() const
    {
        return mpMeshReader->GetNumFaceAttributes();
    }

    /**
     * Method uses the public method of the delegated mesh reader
     */
    void Reset()
    {
        mpMeshReader->Reset();
    }

    /**
     * Method uses the public method of the delegated mesh reader
     */
    std::vector<double> GetNextNode()
    {
        return mpMeshReader->GetNextNode();
    }

    /**
     * Method uses the public method of the delegated mesh reader
     */
    ElementData GetNextElementData()
    {
        return mpMeshReader->GetNextElementData();
    }
    /**
     * Method uses the public method of the delegated mesh reader
     */
    ElementData GetNextFaceData()
    {
        return mpMeshReader->GetNextFaceData();
    }

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global node index
     * @return a vector of the coordinates of the node
     */
    std::vector<double> GetNode(unsigned index)
    {
        return mpMeshReader->GetNode(index);
    }


    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global element index
     * @return a vector of the node indices of the element (and any attribute infomation, if there is any)
     */
    ElementData GetElementData(unsigned index)
    {
        return mpMeshReader->GetElementData(index);
    }


    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global face index
     * @return a vector of the node indices of the face (and any attribute/containment infomation, if there is any)
     */
    ElementData GetFaceData(unsigned index)
    {
        return mpMeshReader->GetFaceData(index);
    }


    /**
     *  Normally throws an exception.  When a NCL file is available, returns a list of the elements
     *  that contain the node (only available for binary files).
     *
     * @param index  The global node index
     * @return a vector of the node indices of the face (and any attribute/containment infomation, if there is any)
     */
    std::vector<unsigned> GetContainingElementIndices(unsigned index)
    {
        return mpMeshReader->GetContainingElementIndices(index);
    }

    /**
     * Get the base name (less any extension) for mesh files.  Only implemented for some mesh types.
     * Method uses the public method of the delegated mesh reader
     */
    std::string GetMeshFileBaseName()
    {
        return mpMeshReader->GetMeshFileBaseName();
    }

    /**
     * Returns true if reading binary files, false if reading ascii files.
     */
    bool IsFileFormatBinary()
    {
        return mpMeshReader->IsFileFormatBinary();
    }

    /**
     * Returns true if there is a node connectivity list (NCL) file available.
     */
    bool HasNclFile()
    {
        return mpMeshReader->HasNclFile();
    }


};

#endif //_GENERICMESHREADER_HPP_
