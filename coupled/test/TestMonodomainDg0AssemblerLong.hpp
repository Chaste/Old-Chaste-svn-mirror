#ifndef _TESTMONODOMAINDG0ASSEMBLERLONG_HPP_
#define _TESTMONODOMAINDG0ASSEMBLERLONG_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"

#include <cxxtest/TestSuite.h>
#include "petscvec.h"
#include <vector>
//#include <iostream>

#include "ConformingTetrahedralMesh.cpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"


#include <time.h>

class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    InitialStimulus *mpStimulus;
    int mNodeNum;
public:
    PointStimulus2dCellFactory(int nodeNum) : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new InitialStimulus(-6000.0, 0.5);
        mNodeNum = nodeNum;
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (node == mNodeNum)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~PointStimulus2dCellFactory(void)
    {
        delete mpStimulus;
    }
};


class TestMonodomainDg0AssemblerLong : public CxxTest::TestSuite 
{
public:

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    // We run for 500 ms and then check that all the voltages at the final time
    // have returned to the resting potential of -84.5
    // test should take about 30mins (or less)
    void TestMonodomainDg02DWithPointStimulusInTheVeryCentreOfTheMesh( void )
    {   
        // To time the solve
        time_t start,end;
        double dif;
        time (&start);
        
        PointStimulus2dCellFactory cell_factory(60); // Central node

        MonodomainProblem<2> monodomain_problem(&cell_factory);

        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(500);   // 500 ms
        monodomain_problem.SetOutputDirectory("testoutput/MonoDg02dWithPointStimulusLong");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_2dWithPointStimulusLong");
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
        
        // To time the solve
        time (&end);
        dif = difftime (end,start);
 //       printf ("\nSolve took %.2lf seconds. \n", dif );
        
        double* p_voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&p_voltage_array, lo, hi); 
    
        // test whether voltages and gating variables are in correct ranges
        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;
            double Ek    = -77.0;
            
            TS_ASSERT_LESS_THAN_EQUALS( p_voltage_array[local_index] , Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-p_voltage_array[local_index] + (Ek-30), 0);
                
            std::vector<double> odeVars = monodomain_problem.GetMonodomainPde()->
                                           GetCardiacCell(global_index)->rGetStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1 
                if ((j!=4) && (j!=3))
                {
                    TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);        
                    TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);        
                }
            }
        }
        
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            /*
             * Test that corners are 'equal', and centres of sides.
             * Irregularities in which way the triangles are oriented make
             * this rather difficult, especially since the edges are sampled
             * during the upstroke.
             */
             
            // corners
            TS_ASSERT_DELTA(p_voltage_array[0], p_voltage_array[10],  0.1);
            TS_ASSERT_DELTA(p_voltage_array[0], p_voltage_array[110], 0.1);
            TS_ASSERT_DELTA(p_voltage_array[0], p_voltage_array[120], 0.1);
            
            // centres of edges
            TS_ASSERT_DELTA(p_voltage_array[5], p_voltage_array[55],  0.1);
            TS_ASSERT_DELTA(p_voltage_array[5], p_voltage_array[65],  0.1);
            TS_ASSERT_DELTA(p_voltage_array[5], p_voltage_array[115], 0.1);
            
            int num_nodes = monodomain_problem.rGetMesh().GetNumNodes();
            // test final voltages have returned to the resting potential
            for (int i=0; i<num_nodes; i++)
            {
                TS_ASSERT_DELTA(p_voltage_array[i], -84.5, 1);
            }
                        
        }
        monodomain_problem.RestoreVoltageArray(&p_voltage_array);
    }   
};
#endif //_TESTMONODOMAINDG0ASSEMBLERLONG_HPP_
