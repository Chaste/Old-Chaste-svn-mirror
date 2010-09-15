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

#include "HeartConfig.hpp"
#include "FileFinder.hpp"
#include "UblasIncludes.hpp"
#include "HeartFileFinder.hpp"

/**
 * TODO:DOX
 * 
 */
template<unsigned DIM>
class FibreReader
{

private:
    /** File stream to use for GetTokensAtNextLine*/
    std::ifstream mDataFile;

    /** Absolute path of the file being read */
    std::string mFilePath;

    /** Number of lines of data in the file, read from the first line of the file */
    unsigned mNumLinesOfData;

    /** Vector which entries read from a line in a file is put into. */
    std::vector<double> mTokens;

    /**
     *  Read a line of numbers from mDataFile.
     *  Sets up the member variable mTokens with the data in the next line.
     *  @return the number of data entries put into mTokens
     */
    unsigned GetTokensAtNextLine();

    /**
     *  Read number of elements from mDataFile
     *  Note: Must be called before GetTokensAtNextLine (it assumes that
     *  it's reading the first line).
     */
    void ReadNumLinesOfDataFromFile();

public:
    /**
     * TODO:DOX
     * 
     * @param rFileFinder the path to the fibre direction file
     */
    FibreReader(HeartFileFinder& rFileFinder);
    
    /**
     *  Destructor closes file
     */
    ~FibreReader();

    /**
     * TODO:DOX
     * 
     * 
     * @param rFibreMatrix Matrix to be filled in.
     */
    void GetNextFibreSheetAndNormalMatrix(c_matrix<double,DIM,DIM>& rFibreMatrix);

    /**
     *  Get the number of lines of data in the matrix - this is the value read from
     *  the first line.
     */
    unsigned GetNumLinesOfData()
    {
        return mNumLinesOfData;
    }
};

#endif /*FIBREREADER_HPP_*/
