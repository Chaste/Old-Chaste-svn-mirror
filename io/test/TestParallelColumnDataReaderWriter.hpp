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
#include "global/src/Exception.hpp"
#include <petsc.h>
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
         
class TestParallelColumnDataReaderWriter : public CxxTest::TestSuite 
{

private: 
	ParallelColumnDataWriter *mpParallelWriter;
	//ColumnDataWriter *mpWriter;
    //ColumnDataReader *mpReader;
	
    const static int num_nodes=10;
	
	 
public:

  
    
    void testParallelColumnWriter(void)
    {
        
        int my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
//        if (my_rank==0)
//        {
//            system("rm -f /tmp/testoutput/ParallelColumnWriter*");
//        }
        
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

        TS_ASSERT_EQUALS(system(("test -f " + output_dir + "ParallelColumnWriter.info").c_str()), 0);
 
 
        //Make a conventional data writer
        /*TS_ASSERT_THROWS_NOTHING(mpWriter = new ColumnDataWriter("/tmp/testoutput","ColumnWriter"));
        TS_ASSERT_THROWS_NOTHING(time_var_id = mpWriter->DefineUnlimitedDimension("Time","msecs"));
        TS_ASSERT_THROWS_NOTHING(mpWriter->DefineFixedDimension("Node","dimensionless", num_nodes));
 
        
        TS_ASSERT_THROWS_NOTHING(var1_id = mpWriter->DefineVariable("Var1","LightYears"));
        TS_ASSERT_THROWS_NOTHING(var2_id = mpWriter->DefineVariable("Var2","Angstroms"));
        TS_ASSERT_THROWS_NOTHING(mpWriter->EndDefineMode());
        */        

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
        TS_ASSERT_THROWS_ANYTHING(mpParallelWriter->PutVector(var1_id, var3) );
        TS_ASSERT_THROWS_ANYTHING(mpParallelWriter->PutVariable(var1_id, 0.0, 0));

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
        
        //delete mpWriter;      
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
    
};

#endif //_TESTCOLUMNDATAREADERWRITER_HPP_
