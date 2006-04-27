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
      
           
    if (mAmMaster)
    {
        for (int global_index=lo ; global_index<hi ; global_index++)
        {
            PutVariable(variableID, petsc_vector_array[global_index - lo], global_index);
        }
    }
    
    VecRestoreArray(petscVector, &petsc_vector_array);      
}



void ParallelColumnDataWriter::EndDefineMode()
{
   //std::cout<<"In EndDefineMode mMyRank="<< mMyRank<<
   //   " mpCurrentOutputFile="<<mpCurrentOutputFile<<"\n";
     if(mAmMaster)
    {
        ColumnDataWriter::EndDefineMode();
    }
    else
    {
        mIsInDefineMode = false;
    }
}

void ParallelColumnDataWriter::PutVariable(int variableID, double variableValue,long dimensionPosition)
{
    //std::cout<<"In PutVariable mMyRank="<< mMyRank<<
    //  " mpCurrentOutputFile="<<mpCurrentOutputFile<<"\n";
    if (mMyRank == 0) 
    {
       //Master process is allowed to write
       ColumnDataWriter::PutVariable(variableID,  variableValue, dimensionPosition);
    } else {
       if(variableID == UNLIMITED_DIMENSION_VAR_ID)
       {
            //The unlimited dimension is written to the ancillary file by the master
            return;
       }
       
        assert(0);
        // We're not going to let individual bits of data be written by
        // slave processes    
    }     
}


ParallelColumnDataWriter::~ParallelColumnDataWriter()
{
}



void ParallelColumnDataWriter::AdvanceAlongUnlimitedDimension()
{
    //Make sure that everyone has queued their messages
    MPI_Barrier(PETSC_COMM_WORLD);
    
    //std::cout<<"In AdvanceAlongUnlimitedDimension mMyRank="<< mMyRank<<
    //  " mpCurrentOutputFile="<<mpCurrentOutputFile<<"\n";
    
    if (mAmMaster){
        //\todo
        //This is where the master is going to take messages from the 
        //slaves and write them
        ColumnDataWriter::DoAdvanceAlongUnlimitedDimension(); 
    }
}

void ParallelColumnDataWriter::Close()
{
    //std::cout<<"In Close mMyRank="<< mMyRank<<"\n";
    
    //\todo.. we may still have queued messages at this point.
    if (mAmMaster){
        ColumnDataWriter::Close();
    }
}


