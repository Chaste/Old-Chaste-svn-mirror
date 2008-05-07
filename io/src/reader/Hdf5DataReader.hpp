/*

Copyright (C) University of Oxford, 2008

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
#ifndef HDF5DATAREADER_HPP_
#define HDF5DATAREADER_HPP_

#include <hdf5.h>
#include <cassert>
#include <vector>
#include <map>
#include <petscvec.h>
#include <iostream>
#include "Exception.hpp"
#include "OutputFileHandler.hpp"

class Hdf5DataReader
{
private:
    // defined in writer too, todo: define it once
    const static unsigned MAX_DATASET_RANK=3;

    std::string mDirectory; /**< Directory output files will be stored in. */
    std::string mBaseName; /**< The base name for the output data files. */
    
    hid_t mFileId;
    
    hid_t mVariablesDatasetId;    
    unsigned mVariablesDatasetRank;
    hsize_t mVariablesDatasetSizes[MAX_DATASET_RANK];

    bool mIsUnlimitedDimensionSet; /**< Is the unlimited dimension set */
    hid_t mTimeDatasetId;
    hsize_t mNumberTimesteps;
    
    std::map<std::string, unsigned>    mVariableToColumnIndex;  
    std::vector<std::string>    mVariableNames;  
    
    std::map<std::string, std::string> mVariableToUnit;
    
    bool mIsDataComplete;
    std::vector<unsigned> mIncompleteNodeIndices; 
    
public:

    Hdf5DataReader(std::string directory, std::string baseName, bool make_absolute=true);
    
    std::vector<double> GetVariableOverTime(std::string variableName, unsigned nodeIndex);
    
    void GetVariableOverNodes(Vec data, std::string variableName, unsigned timestep=0);
    
    std::vector<double> GetUnlimitedDimensionValues();
    
    unsigned GetNumberOfRows()
    {
        return mVariablesDatasetSizes[1];
    }
    
    std::vector<std::string> GetVariableNames()
    {
        return mVariableNames;
    }
    std::string GetUnit(std::string variableName)
    {
        return mVariableToUnit[variableName];
    }
    
    bool IsDataComplete()
    {
        return mIsDataComplete;
    }
    
    std::vector<unsigned> GetIncompleteNodeMap()
    {
        return mIncompleteNodeIndices;
    }
    
    void Close();

};

#endif /*HDF5DATAREADER_HPP_*/
