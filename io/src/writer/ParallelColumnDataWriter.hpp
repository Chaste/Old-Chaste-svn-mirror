#ifndef PARALLELCOLUMNDATAWRITER_HPP_
#define PARALLELCOLUMNDATAWRITER_HPP_


#include "ColumnDataWriter.hpp"
#include <petscvec.h>
class ParallelColumnDataWriter  : public ColumnDataWriter
{
private:
    bool mIsParallel;        
    bool mAmMaster;
    Vec mConcentrated;
    VecScatter mToMaster;
    int mNumProcs; //\todo remove
    int mMyRank;
public:
	ParallelColumnDataWriter(std::string directory, std::string baseName);
    virtual ~ParallelColumnDataWriter();
    void PutVector(int variableID, Vec PetscVector);
    void PutVariable(int variableID, double variableValue,long dimensionPosition = -1);
    void EndDefineMode();
    void AdvanceAlongUnlimitedDimension(); 
    void Close();
};

#endif /*PARALLELCOLUMNDATAWRITER_HPP_*/
