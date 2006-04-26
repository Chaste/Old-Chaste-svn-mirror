#include "ParallelColumnDataWriter.hpp"

ParallelColumnDataWriter::ParallelColumnDataWriter(std::string directory, std::string baseName)
: ColumnDataWriter::ColumnDataWriter(directory, baseName)
{    
}

void
ParallelColumnDataWriter::PutVector(int variableID, Vec petscVector)
{
    //\todo Doesn't work in parallel 
    int lo, hi;
    double *petsc_vector_array;
    VecGetArray(petscVector, &petsc_vector_array); 
    VecGetOwnershipRange(petscVector,&lo,&hi);
        
    for (int global_index=lo ; global_index<hi ; global_index++)
    {
        PutVariable(variableID, petsc_vector_array[global_index - lo], global_index);
    }
    
    VecRestoreArray(petscVector, &petsc_vector_array);      
}

ParallelColumnDataWriter::~ParallelColumnDataWriter()
{
}
