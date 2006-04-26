#include "ParallelColumnDataWriter.hpp"

ParallelColumnDataWriter::ParallelColumnDataWriter(std::string directory, std::string baseName)
: ColumnDataWriter::ColumnDataWriter(directory, baseName)
{    
    int num_procs;
    MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
    if (num_procs==1){
        is_parallel=false;
    } else {
        is_parallel=true;
    } 
    
    int my_rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
    if (my_rank==0){
        am_master=true;
    } else {
        am_master=false;
    } 
    
    if (is_parallel)
    {
        std::stringstream suffix;
        suffix << std::setfill('0') << std::setw(3) << my_rank;
                   mBaseName = mBaseName+suffix.str(); 
    }    
}

void
ParallelColumnDataWriter::PutVector(int variableID, Vec petscVector)
{
    int lo, hi;
    double *petsc_vector_array;
    VecGetArray(petscVector, &petsc_vector_array); 
    VecGetOwnershipRange(petscVector,&lo,&hi);
        
    for (int global_index=hi-1 ; global_index>=lo ; global_index--)
    {
        PutVariable(variableID, petsc_vector_array[global_index - lo], global_index);
    }
    
    VecRestoreArray(petscVector, &petsc_vector_array);      
}

ParallelColumnDataWriter::~ParallelColumnDataWriter()
{
    //\todo Doesn't work in parallel 

    //Concatenate files...
}
