#ifndef PARALLELCOLUMNDATAWRITER_HPP_
#define PARALLELCOLUMNDATAWRITER_HPP_


#include "ColumnDataWriter.hpp"
#include <petscvec.h>
class ParallelColumnDataWriter  : public ColumnDataWriter
{
public:
	ParallelColumnDataWriter(std::string directory, std::string baseName);
    virtual ~ParallelColumnDataWriter();
    void PutVector(int variableID, Vec PetscVector);
};

#endif /*PARALLELCOLUMNDATAWRITER_HPP_*/
