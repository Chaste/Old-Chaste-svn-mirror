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

#include "FibreReader.hpp"

#include <sstream>
#include "Exception.hpp"

template<unsigned DIM>
FibreReader<DIM>::FibreReader(FileFinder& rFileFinder, unsigned axiOrOrtho)
{
    if ((axiOrOrtho != 1) && (axiOrOrtho != 2))
    {
        EXCEPTION("Fibres must be axisymmetric or orthotropic");
    }
    if (axiOrOrtho == 1)
    {
        mNumItemsPerLine = DIM;
    }
    else
    {
        mNumItemsPerLine = DIM*DIM;
    }
    mTokens.resize(mNumItemsPerLine);

    mFilePath = rFileFinder.GetAbsolutePath();
    mDataFile.open(mFilePath.c_str());
    if (!mDataFile.is_open())
    {
        EXCEPTION("Failed to open fibre file " + rFileFinder.GetAbsolutePath());
    }

    // Note: this method will close the file on error
    ReadNumLinesOfDataFromFile();
}

template<unsigned DIM>
FibreReader<DIM>::~FibreReader()
{
    mDataFile.close();
}

template<unsigned DIM>
void FibreReader<DIM>::GetNextFibreSheetAndNormalMatrix(c_matrix<double,DIM,DIM>& rFibreMatrix)
{
    if (mNumItemsPerLine != DIM*DIM)
    {
        EXCEPTION("Use GetNextFibreVector when reading axisymmetric fibres.");
    }

    unsigned num_entries = GetTokensAtNextLine();
    if(num_entries < mNumItemsPerLine)
    {
        std::stringstream string_stream;
        string_stream << "A line is incomplete in " << mFilePath
                      << " - each line should contain " << DIM*DIM << " entries";
        EXCEPTION(string_stream.str());
    }

    for(unsigned i=0; i<DIM; i++)
    {
        for(unsigned j=0; j<DIM; j++)
        {
            rFibreMatrix(j,i) = mTokens[DIM*i + j];
        }
    }
}

template<unsigned DIM>
void FibreReader<DIM>::GetNextFibreVector(c_vector<double,DIM>& rFibreVector)
{
    if (mNumItemsPerLine != DIM)
    {
        EXCEPTION("Use GetNextFibreSheetAndNormalMatrix when reading orthotropic fibres.");
    }

    unsigned num_entries = GetTokensAtNextLine();
    if(num_entries < mNumItemsPerLine)
    {
        std::stringstream string_stream;
        string_stream << "A line is incomplete in " << mFilePath
                      << " - each line should contain " << DIM << " entries";
        EXCEPTION(string_stream.str());
    }

    for(unsigned i=0; i<DIM; i++)
    {
        rFibreVector(i) = mTokens[i];
    }
}


template<unsigned DIM>
unsigned FibreReader<DIM>::GetTokensAtNextLine()
{
    std::string line;

    bool comment_line;
    bool blank_line;

    assert(mTokens.size() == mNumItemsPerLine);

    do
    {
        getline(mDataFile, line);

        if (line.empty() && mDataFile.eof())
        {
            mDataFile.close();
            EXCEPTION("Fibre orientation file contains less fibre definitions than the number of lines of data defined at the beginning of the file");
        }

        comment_line = (line.find('#',0) != std::string::npos);
        blank_line = (line.find_first_not_of(" \t",0) == std::string::npos);
    }
    while(comment_line || blank_line);

    std::string::iterator  iter=line.end();
    iter--;
    unsigned nchars2delete = 0;
    while(*iter == ' ')
    {
        nchars2delete++;
        iter--;
    }
    line.erase(line.length()-nchars2delete);

    std::stringstream line_stream(line);

    unsigned index = 0;
    while (!line_stream.eof())
    {
        double item;
        line_stream >> item;
        if(index >= mNumItemsPerLine)
        {
            EXCEPTION("Too many entries in a line in " + mFilePath);
        }
        mTokens[index++] = item;
    }

    return index; // the number of entries put into mTokens
}


template<unsigned DIM>
void FibreReader<DIM>::ReadNumLinesOfDataFromFile()
{
    if (GetTokensAtNextLine() != 1)
    {
        mDataFile.close();
        EXCEPTION("First (non comment) line of the fibre orientation file should contain the number of lines of data in the file (and nothing else)");
    }

    mNumLinesOfData = (unsigned) mTokens[0];
}



template class FibreReader<1>;
template class FibreReader<2>;
template class FibreReader<3>;





