#ifndef _TESTCOLUMNDATAREADERWRITER_HPP_
#define _TESTCOLUMNDATAREADERWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "ColumnDataWriter.hpp"
#include "ParallelColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include <petsc.h>
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"

class TestParallelColumnDataReaderWriter : public CxxTest::TestSuite
{

private:
    ParallelColumnDataWriter *mpParallelWriter;
    ColumnDataReader *mpReader;
    
    
    const static int num_nodes=10;
    
    
public:

    void TestParallelColumnWriter(void)
    {
    
    
        int time_var_id, var1_id, var2_id;
        
        //Make a parallel data writer
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter = new ParallelColumnDataWriter("TestParallelColumnDataWriter","ParallelColumnWriter"));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpParallelWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->DefineFixedDimension("Node","dimensionless", num_nodes));
        
        TS_ASSERT_THROWS_NOTHING(var1_id = mpParallelWriter->DefineVariable("Var1","LightYears"));
        TS_ASSERT_THROWS_NOTHING(var2_id = mpParallelWriter->DefineVariable("Var2","Angstroms"));
        
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->EndDefineMode());
        
        MPI_Barrier(PETSC_COMM_WORLD);
        
        std::string output_dir = mpParallelWriter->GetOutputDirectory();
        
        // Test that the ouput directory and the .info file was created
        TS_ASSERT_EQUALS(system(("test -f " + output_dir + "ParallelColumnWriter.info").c_str()), 0);
        
        
        //Set up some data in PETSc vectors
        Vec var1, var2, var3;
        VecCreate(PETSC_COMM_WORLD, &var1);
        VecCreate(PETSC_COMM_WORLD, &var2);
        VecCreate(PETSC_COMM_WORLD, &var3);
        VecSetSizes(var1, PETSC_DECIDE, num_nodes);
        VecSetSizes(var2, PETSC_DECIDE, num_nodes);
        VecSetSizes(var3, PETSC_DECIDE, num_nodes+1);
        VecSetFromOptions(var1);
        VecSetFromOptions(var2);
        VecSetFromOptions(var3);
        
        double* var1_array,* var2_array;
        VecGetArray(var1, &var1_array);
        VecGetArray(var2, &var2_array);
        int lo, hi;
        VecGetOwnershipRange(var1,&lo,&hi);
        
        for (int global_index=lo ; global_index<hi ; global_index++)
        {
            var1_array[global_index-lo]=global_index;
            var2_array[global_index-lo]=global_index+100;
        }
        
        VecRestoreArray(var1, &var1_array);
        VecAssemblyBegin(var1);
        VecAssemblyEnd(var1);
        VecRestoreArray(var2, &var2_array);
        VecAssemblyBegin(var2);
        VecAssemblyEnd(var2);
        
        //Write out the data (Conventional)
        /*VecGetArray(var1, &var1_array);
        VecGetArray(var2, &var2_array); 
        mpWriter->PutVariable(time_var_id, 0); 
        for (int global_index=lo ; global_index<hi ; global_index++)
        {        
            mpWriter->PutVariable(var1_id, var1_array[global_index - lo], global_index);
            mpWriter->PutVariable(var2_id, var2_array[global_index - lo], global_index);
        }
        mpWriter->AdvanceAlongUnlimitedDimension();
        VecRestoreArray(var1, &var1_array); 
        VecRestoreArray(var2, &var2_array);
        */
        
        //Write out the data (Parallel)
        mpParallelWriter->PutVariable(time_var_id, 0.1);
        TS_ASSERT_THROWS_NOTHING( mpParallelWriter->PutVector(var1_id, var1) );
        TS_ASSERT_THROWS_NOTHING( mpParallelWriter->PutVector(var2_id, var2) );
        // Throws since var3 is the wrong size
        TS_ASSERT_THROWS_ANYTHING(mpParallelWriter->PutVector(var1_id, var3) );
        
        // No-op if not master, writes anyway if we are
        TS_ASSERT_THROWS_NOTHING(mpParallelWriter->PutVariable(var1_id, 0.0, 0));
        
        TS_ASSERT_THROWS_NOTHING( mpParallelWriter->AdvanceAlongUnlimitedDimension() );
        
        //Change the data
        VecSqrt(var1);
        VecSqrt(var2);
        
        //Write out the data again (Conventional)
        /*
        VecGetArray(var1, &var1_array); 
        VecGetArray(var2, &var2_array); 
        mpWriter->PutVariable(time_var_id, 1); 
        for (int global_index=lo ; global_index<hi ; global_index++)
        {        
            mpWriter->PutVariable(var1_id, var1_array[global_index - lo], global_index);
            mpWriter->PutVariable(var2_id, var2_array[global_index - lo], global_index);
        }
        */
                
        
        //Write out the data (Parallel)
        mpParallelWriter->PutVariable(time_var_id, 0.2);
        mpParallelWriter->PutVector(var1_id, var1);
        mpParallelWriter->PutVector(var2_id, var2);
        //mpParallelWriter->AdvanceAlongUnlimitedDimension();
        
        //delete mpParallelWriter;
        delete mpParallelWriter;
        
        MPI_Barrier(PETSC_COMM_WORLD);
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"ParallelColumnWriter.info io/test/data/ColumnWriter.info").c_str()),
                         0);
                         
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"ParallelColumnWriter_000000.dat io/test/data/ColumnWriter_000000.dat").c_str()),
                         0);
                         
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"ParallelColumnWriter_000001.dat io/test/data/ColumnWriter_000001.dat").c_str()),
                         0);
                         
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"ParallelColumnWriter_unlimited.dat io/test/data/ColumnWriter_unlimited.dat").c_str()),
                         0);
                         
        VecDestroy(var1);
        VecDestroy(var2);
        VecDestroy(var3);
    }
    
    
    void TestPutSlice(void)
    {
        // create a vector slice
        
        const unsigned problem_size=10;
        DistributedVector::SetProblemSize(problem_size);
        Vec striped=DistributedVector::CreateVec(2);
        
        DistributedVector distributed_vector(striped);
        DistributedVector::Stripe zeros(distributed_vector,0);
        DistributedVector::Stripe ones(distributed_vector,1);
        // write some values
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            zeros[index] =  0;
            ones[index] =  1;
        }
        
        // write to file with parallel data writer
        
        ParallelColumnDataWriter* p_parallel_writer = new ParallelColumnDataWriter("TestParallelColumnDataWriterStripe","Stripe");
        unsigned time_var_id = p_parallel_writer->DefineUnlimitedDimension("Time","msecs");
        unsigned var1_id = p_parallel_writer->DefineVariable("Var1","LightYears");
        p_parallel_writer->DefineFixedDimension("Node","dimensionless", problem_size);
        p_parallel_writer->EndDefineMode();
        MPI_Barrier(PETSC_COMM_WORLD);
        std::string output_dir = p_parallel_writer->GetOutputDirectory();
        
        p_parallel_writer->PutVariable(time_var_id, 0.1);
        p_parallel_writer->PutVectorStripe(var1_id, ones);
        p_parallel_writer->AdvanceAlongUnlimitedDimension();
        
        // check file
        
        MPI_Barrier(PETSC_COMM_WORLD);
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"Stripe.info io/test/data/Stripe.info").c_str()),
                         0);
                         
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"Stripe_000000.dat io/test/data/Stripe_000000.dat").c_str()),
                         0);
                         
        TS_ASSERT_EQUALS(system(
                             ("diff "+output_dir+"Stripe_unlimited.dat io/test/data/Stripe_unlimited.dat").c_str()),
                         0);
        
        VecDestroy(striped);
        delete p_parallel_writer;
    }
    
    // Read back the data written in the test above
    void TestColumnReader(void)
    {
        //There is no *Parallel* ColumnDataReader.  Since everyone might
        //need to know everything there's no point in only one processor opening the
        //file
        
        //Make a parallel data writer
        TS_ASSERT_THROWS_NOTHING(mpReader = new ColumnDataReader("TestParallelColumnDataWriter","ParallelColumnWriter"));
        
        //Check that there's the correct number of files
        std::vector<double> time_stamps;
        time_stamps=mpReader->GetUnlimitedDimensionValues();
        TS_ASSERT_EQUALS(time_stamps[0],0.1);
        TS_ASSERT_EQUALS(time_stamps[1],0.2);
        TS_ASSERT_EQUALS(time_stamps.size(),2u);
        
        //Check that some of the data is correct
        std::vector<double> var1_node4;
        var1_node4 = mpReader->GetValues("Var1",4);
        TS_ASSERT_EQUALS(var1_node4[0],4.0); //First time step
        TS_ASSERT_EQUALS(var1_node4[1],2.0); //Second time step
        std::vector<double> var2_node4;
        var2_node4 = mpReader->GetValues("Var2",4);
        TS_ASSERT_EQUALS(var2_node4[0],104.0); //First time step
        TS_ASSERT_DELTA(var2_node4[1],sqrt(104.0),1e-4); //Second time step
        
        TS_ASSERT_THROWS_ANYTHING(mpReader->GetValues("LifeSigns",4));
        TS_ASSERT_THROWS_ANYTHING(mpReader->GetValues("Var1",10));
        
        //Delete the reader: makes sure that files are closed
        delete mpReader;
    }
};

#endif //_TESTCOLUMNDATAREADERWRITER_HPP_
