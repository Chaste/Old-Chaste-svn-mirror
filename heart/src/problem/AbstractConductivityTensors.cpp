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

#include "AbstractConductivityTensors.hpp"
#include "Exception.hpp"
#include <sstream>

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::OpenFibreOrientationFile()
{
    mDataFile.open(mFibreOrientationFilename.c_str());
    if (!mDataFile.is_open())
    {
        EXCEPTION("Wrong fibre orientation file name "+mFibreOrientationFilename);
    }
}

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::CloseFibreOrientationFile()
{
    mDataFile.close();
}

template<unsigned SPACE_DIM>
unsigned AbstractConductivityTensors<SPACE_DIM>::GetTokensAtNextLine(std::vector<double>& rTokens)
{
    std::string line;

    bool comment_line;
    bool blank_line;
    //We've assuming this is a fresh vector.  Why?
    assert(rTokens.size() == 0);
    do
    {
        getline(mDataFile, line);

        if (mDataFile.eof())
        {
            CloseFibreOrientationFile();
            EXCEPTION("Fibre orientation file contains less fibre definitions than the number of elements in the mesh");
        }

        comment_line = (line.find('#',0) != std::string::npos);
        blank_line = (line.find_first_not_of(" \t",0) == std::string::npos);
    }
    while(comment_line || blank_line);


    std::stringstream line_stream(line);

    // Read all the numbers from the line
    while (!line_stream.eof())
    {
        double item;
        line_stream >> item;
        rTokens.push_back(item);
    }

    return rTokens.size();
}

template<unsigned SPACE_DIM>
unsigned AbstractConductivityTensors<SPACE_DIM>::GetNumElementsFromFile()
{
    std::vector<double> tokens;

    if (GetTokensAtNextLine(tokens) != 1)
    {
        CloseFibreOrientationFile();
        EXCEPTION("First (non comment) line of the fibre orientation file should contain the number of elements of the mesh (and nothing else)");
    }

    return (unsigned) tokens[0];
}


template<unsigned SPACE_DIM>
AbstractConductivityTensors<SPACE_DIM>::AbstractConductivityTensors()
    : mNumElements(1),
      mUseNonConstantConductivities(false),
      mUseFibreOrientation(false),
      mInitialised(false)
{
    double init_data[]={DBL_MAX, DBL_MAX, DBL_MAX};

    for (unsigned dim=0; dim<SPACE_DIM; dim++)
    {
         mConstantConductivities[dim] = init_data[dim];
    }
}

template<unsigned SPACE_DIM>
AbstractConductivityTensors<SPACE_DIM>::~AbstractConductivityTensors()
{
}

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::SetFibreOrientationFile(const std::string &rFibreOrientationFilename)
{
    mUseFibreOrientation = true;
    mFibreOrientationFilename = rFibreOrientationFilename;
}

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::SetConstantConductivities(c_vector<double, 1> constantConductivities)
{
    if (SPACE_DIM != 1)
    {
        EXCEPTION("Wrong number of conductivities provided");
    }

    mUseNonConstantConductivities = false;
    mConstantConductivities = constantConductivities;
}

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::SetConstantConductivities(c_vector<double, 2> constantConductivities)
{
    if (SPACE_DIM != 2)
    {
        EXCEPTION("Wrong number of conductivities provided");
    }

    mUseNonConstantConductivities = false;
    mConstantConductivities = constantConductivities;
}

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::SetConstantConductivities(c_vector<double, 3> constantConductivities)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("Wrong number of conductivities provided");
    }

    mUseNonConstantConductivities = false;
    mConstantConductivities = constantConductivities;
}

template<unsigned SPACE_DIM>
void AbstractConductivityTensors<SPACE_DIM>::SetNonConstantConductivities(std::vector<c_vector<double, SPACE_DIM> >* pNonConstantConductivities)
{
    mUseNonConstantConductivities = true;
    mpNonConstantConductivities = pNonConstantConductivities;
}

template<unsigned SPACE_DIM>
c_matrix<double,SPACE_DIM,SPACE_DIM>& AbstractConductivityTensors<SPACE_DIM>::operator[](const unsigned index)
{
    assert(mInitialised);

    if (!mUseNonConstantConductivities && !mUseFibreOrientation)
    {
        return mTensors[0];
    }
    else
    {
        assert(index < mNumElements);
        return mTensors[index];
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractConductivityTensors<1>;
template class AbstractConductivityTensors<2>;
template class AbstractConductivityTensors<3>;
