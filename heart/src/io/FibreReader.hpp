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

#ifndef FIBREREADER_HPP_
#define FIBREREADER_HPP_

#include <string>
#include <fstream>
#include <vector>

#include "UblasIncludes.hpp"
#include "FileFinder.hpp"

/**
 * TODO:DOX
 * 
 */
template<unsigned DIM>
class FibreReader
{

private:
    /** File stream to use for GetTokensAtNextLine */
    std::ifstream mDataFile;

    /** Absolute path of the file being read */
    std::string mFilePath;

    /** Number of lines of data in the file, read from the first line of the file */
    unsigned mNumLinesOfData;

    /** How many items we expect to find per line: DIM for axisymmetric, DIM*DIM for orthotropic */
    unsigned mNumItemsPerLine;

    /** Vector which entries read from a line in a file is put into. */
    std::vector<double> mTokens;

    /**
     *  Read a line of numbers from #mDataFile.
     *  Sets up the member variable #mTokens with the data in the next line.
     *  @return the number of data entries put into #mTokens
     */
    unsigned GetTokensAtNextLine();

    /**
     *  Read number of elements from #mDataFile.
     *  Note: Must be called before GetTokensAtNextLine (it assumes that
     *  it's reading the first line).
     */
    void ReadNumLinesOfDataFromFile();

public:
    /**
     * Create a new FibreReader.
     * 
     * @param rFileFinder  the path to the fibre direction file
     * @param axiOrOrtho  1 for axisymmetric, 2 for orthotropic
     */
    FibreReader(FileFinder& rFileFinder, unsigned axiOrOrtho);
    
    /**
     *  Destructor closes file.
     */
    ~FibreReader();

    /**
     * Read the next fibre direction matrix from the file.  Must only be used when
     * reading an orthotropic file.  These have lines of the form
     * \code
     *  fibre0 fibre1 fibre2 sheet0 sheet1 sheet2 normal0 normal1 normal2
     * \endcode
     * which are converted to the matrix
     * \code
     *     [ fibre0   sheet0   normal0  ]
     *     [ fibre1   sheet1   normal1  ]
     *     [ fibre2   sheet2   normal2  ]
     * \endcode
     * 
     * @param rFibreMatrix  matrix to be filled in
     */
    void GetNextFibreSheetAndNormalMatrix(c_matrix<double,DIM,DIM>& rFibreMatrix);

    /**
     * Read the next fibre direction vector from the file.  Must only be used when
     * reading an axisymmetric file.  These have lines of the form
     * \code
     *  fibre0 fibre1 fibre2
     * \endcode
     *
     * @param rFibreVector  vector to be filled in
     */
    void GetNextFibreVector(c_vector<double,DIM>& rFibreVector);

    /**
     *  Get the number of lines of data in the file - this is the value read from
     *  the first line.
     */
    unsigned GetNumLinesOfData()
    {
        return mNumLinesOfData;
    }
};

#endif /*FIBREREADER_HPP_*/
