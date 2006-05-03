#include "ParallelColumnDataWriter.hpp"
#include "global/src/Exception.hpp"
#include <iostream>

ParallelColumnDataWriter::ParallelColumnDataWriter(std::string directory, std::string baseName)
: ColumnDataWriter::ColumnDataWriter(directory, baseName)
{    
    mConcentrated=NULL;
    
    int num_procs, my_rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
    if (num_procs==1){
        mIsParallel=false;
    } else {
        mIsParallel=true;
    } 
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
    if (my_rank==0){
        mAmMaster=true;
    } else {
        mAmMaster=false;
    } 
    

}

void
ParallelColumnDataWriter::PutVector(int variableID, Vec petscVector)
{
    int size;
    VecGetSize(petscVector,&size);
    
    if(size != mFixedDimensionSize)
    {
        //std::cout << "fixed_dim: " << mFixedDimensionSize << ", vec_size " << size << "\n";
        throw Exception("Size of vector does not match FixedDimensionSize.");
    }
    
    //Construct the appropriate "scatter" object to concentrate the vector on the master
    if (mConcentrated==NULL){
        VecScatterCreateToZero(petscVector, &mToMaster, &mConcentrated);
    }  
    
    VecScatterBegin(petscVector, mConcentrated, INSERT_VALUES, SCATTER_FORWARD, mToMaster);
    VecScatterEnd(petscVector, mConcentrated, INSERT_VALUES, SCATTER_FORWARD, mToMaster);
           
    if (mAmMaster)
    {
        double *concentrated_vector;
        VecGetArray(mConcentrated, &concentrated_vector); 
        for (int i=0 ; i<size; i++)
        {
            ColumnDataWriter::PutVariable(variableID, concentrated_vector[i], i);
        }
        VecRestoreArray(mConcentrated, &concentrated_vector);      
    }
    
}



void ParallelColumnDataWriter::EndDefineMode()
{
    if(mAmMaster)
    {
        ColumnDataWriter::EndDefineMode();
    }
    else
    {
        mIsInDefineMode = false;
    }
    
/*    int num_procs, my_rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
    for (int i=0; i<num_procs; i++)
    {
        if (i==my_rank)
        {
            std::cerr<<"On "<<my_rank<<":\n";
            system("ls -l testoutput/ParallelColumnWriter.info");
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
*/    
}

/**
 * There are two ways of calling PutVariable
 * 1) All processes call it as a collective operation from the user's code.
 *    This only makes sense if they are writing the unlimited dimension (time) variable.
 *    It is an error if any non-master process writes anything other than
 *      unlimited dimension 
 * 2) The master calls the equivalent method in the parent class after concentrating 
 *      the data into a single Vec (ie. from the method PutVector() above).
 */
void ParallelColumnDataWriter::PutVariable(int variableID, double variableValue,long dimensionPosition)
{
    if(variableID != UNLIMITED_DIMENSION_VAR_ID)
    {
       throw Exception("Non-master processes cannot write to disk.");
    }
    
    if (mAmMaster) 
    {
       //Master process is allowed to write
       ColumnDataWriter::PutVariable(variableID,  variableValue, dimensionPosition);
    }
}


ParallelColumnDataWriter::~ParallelColumnDataWriter()
{
    if (mConcentrated != NULL) {
        VecScatterDestroy(mToMaster);
        VecDestroy(mConcentrated);
    }
        
}



void ParallelColumnDataWriter::AdvanceAlongUnlimitedDimension()
{
    //Make sure that everyone has queued their messages
    MPI_Barrier(PETSC_COMM_WORLD);
    
//    std::cout<<"In AdvanceAlongUnlimitedDimension mAmMaster="<< mAmMaster<<
//     " mpCurrentOutputFile="<<mpCurrentOutputFile<<"\n"<<std::flush;
    
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


