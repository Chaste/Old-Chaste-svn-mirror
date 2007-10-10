#ifndef TESTPARALLELWRITERPERFORMANCE_HPP_
#define TESTPARALLELWRITERPERFORMANCE_HPP_

#include <cxxtest/TestSuite.h>
#include "ParallelColumnDataWriter.hpp"
#include "DistributedVector.hpp"
#include <petsc.h>
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"

class TestParallelWriterPerformance : public CxxTest::TestSuite
{
public:
    void Test1()
    {
        const unsigned SIZE=1000;
        const unsigned REPETITIONS=10;
        // create a distibuted vector
        DistributedVector::SetProblemSize(SIZE);
        Vec petsc_vec=DistributedVector::CreateVec();
        DistributedVector distributed_vector(petsc_vec);
        for (DistributedVector::Iterator index = DistributedVector::Begin();
             index!= DistributedVector::End();
             ++index)
        {
            distributed_vector[index] =  -(double)(index.Local*index.Global);
        }
        distributed_vector.Restore();
        
        // set up a parallel writer
        ParallelColumnDataWriter parallel_writer("TestParallelWriterPerformance","ParallelColumnWriter");
        unsigned time_var_id = parallel_writer.DefineUnlimitedDimension("Time","msecs");
        parallel_writer.DefineFixedDimension("Node","dimensionless", SIZE);
        unsigned var1_id = parallel_writer.DefineVariable("Var1","LightYears");
        parallel_writer.EndDefineMode();
        
        // write multiple times
        for (unsigned i=0; i<REPETITIONS; i++)
        {
            double time=(double)i;
            parallel_writer.PutVariable(time_var_id, time);
            parallel_writer.PutVector(var1_id, petsc_vec);
            parallel_writer.AdvanceAlongUnlimitedDimension();
        }
        
        VecDestroy(petsc_vec);
    }
};
#endif /*TESTPARALLELWRITERPERFORMANCE_HPP_*/
