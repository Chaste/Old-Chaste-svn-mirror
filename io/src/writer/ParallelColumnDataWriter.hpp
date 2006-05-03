#ifndef PARALLELCOLUMNDATAWRITER_HPP_
#define PARALLELCOLUMNDATAWRITER_HPP_


#include "ColumnDataWriter.hpp"
#include <petscvec.h>
class ParallelColumnDataWriter  : public ColumnDataWriter
{
private:
    bool mIsParallel;        /**< set to true in constructor if running in parallel*/
    bool mAmMaster;          /**< set to true in constructor for process is the rank 0 process*/
    Vec mConcentrated;       /**< Vector to hold concentrated copy of distributed vector on the master process*/
    VecScatter mToMaster;    /**< variable holding information for concentrating a vector*/
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
