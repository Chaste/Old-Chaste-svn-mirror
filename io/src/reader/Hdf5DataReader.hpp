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
    
public:

    Hdf5DataReader(std::string directory, std::string baseName, bool make_absolute=true);
    
    std::vector<double> GetVariableOverTime(std::string variableName, unsigned nodeIndex);
    
    void GetVariableOverNodes(Vec data, std::string variableName, unsigned timestep=0);
    
    std::vector<double> GetUnlimitedDimensionValues();
    
    std::vector<std::string> GetVariableNames()
    {
        return mVariableNames;
    }
    std::string GetUnit(std::string variableName)
    {
        return mVariableToUnit[variableName];
    }
    
    void Close();

};

#endif /*HDF5DATAREADER_HPP_*/
