#include "ParallelColumnDataWriter.hpp"
#include "global/src/Exception.hpp"
#include <iostream>

ParallelColumnDataWriter::ParallelColumnDataWriter(std::string directory, std::string baseName)
: ColumnDataWriter::ColumnDataWriter(directory, baseName)
{    
    
    MPI_Comm_size(PETSC_COMM_WORLD, &mNumProcs);
    if (mNumProcs==1){
        mIsParallel=false;
    } else {
        mIsParallel=true;
    } 
    
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &mMyRank);
    if (mMyRank==0){
        mAmMaster=true;
    } else {
        mAmMaster=false;
    } 
    
//    if (mIsParallel)
//    {
//        std::stringstream suffix;
//        suffix << std::setfill('0') << std::setw(3) << my_rank;
//        mBaseName = mBaseName+suffix.str(); 
//    }    
}

void
ParallelColumnDataWriter::PutVector(int variableID, Vec petscVector)
{
    int lo, hi, size;
    double *petsc_vector_array;
    VecGetArray(petscVector, &petsc_vector_array); 
    VecGetOwnershipRange(petscVector,&lo,&hi);
    VecGetSize(petscVector,&size);
    
    if(size != mFixedDimensionSize)
    {
        //std::cout << "fixed_dim: " << mFixedDimensionSize << ", vec_size " << size << "\n";
        throw Exception("Size of vector does not match FixedDimensionSize.");
    }
      
    for (int i=0; i< mNumProcs; i++)
    { 
        if (mMyRank == i)
        {   
            for (int global_index=lo ; global_index<hi ; global_index++)
            {
                 PutVariable(variableID, petsc_vector_array[global_index - lo], global_index);
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
    
     VecRestoreArray(petscVector, &petsc_vector_array);      
}

ParallelColumnDataWriter::~ParallelColumnDataWriter()
{
    //\todo Doesn't work in parallel 

    //Concatenate files...
}

void ParallelColumnDataWriter::EndDefineMode()
{
//    if(mAmMaster)
//    {
        ColumnDataWriter::EndDefineMode();
//    }
//    else
//    {
//        mIsInDefineMode = false;
//    }
}
