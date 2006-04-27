#ifndef PARALLELCOLUMNDATAWRITER_HPP_
#define PARALLELCOLUMNDATAWRITER_HPP_


#include "ColumnDataWriter.hpp"
#include <petscvec.h>
class ParallelColumnDataWriter  : public ColumnDataWriter
{
private:
    bool mIsParallel;        
    bool mAmMaster;
    int mNumProcs;
    int mMyRank;
public:
	ParallelColumnDataWriter(std::string directory, std::string baseName);
    virtual ~ParallelColumnDataWriter();
    void PutVector(int variableID, Vec PetscVector);
    
    void EndDefineMode();
};

#endif /*PARALLELCOLUMNDATAWRITER_HPP_*/
