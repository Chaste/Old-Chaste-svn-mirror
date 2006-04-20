#ifndef _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_

// Element.hpp includes the Boost ublas objects - these need to
// be included early...  We think.  We're not that sure.
#include "Element.hpp"


#include <cxxtest/TestSuite.h>
#include <petsc.h>
#include "SimpleLinearSolver.hpp"
#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <iostream>
#include "Node.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleDg0ParabolicAssembler.hpp"  
#include "TrianglesMeshReader.hpp"

//#include "MonodomainPdeFitzHughNagumo.hpp"
#include "MonodomainPde.hpp"
#include "MonodomainProblem.hpp"

#include "MonodomainDg0Assembler.hpp"
#include "ColumnDataWriter.hpp"
#include <cmath>
#include "EulerIvpOdeSolver.hpp"
#include "FitzHughNagumo1961OdeSystem.hpp"
#include "PetscSetupAndFinalize.hpp"

class FhnEdgeStimulusCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    InitialStimulus *mpStimulus;
public:
    FhnEdgeStimulusCellFactory() : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new InitialStimulus(-600.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(int node)
    {
        if (mpMesh->GetNodeAt(node)->GetPoint()[0] == 0.0)
        {
            return new FitzHughNagumo1961OdeSystem(mpSolver, mpStimulus, mTimeStep);
        }
        else
        {
            return new FitzHughNagumo1961OdeSystem(mpSolver, mpZeroStimulus, mTimeStep);
        }
    }
    
    ~FhnEdgeStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};

class TestMonodomainFitzHughNagumoWithDg0Assembler : public CxxTest::TestSuite 
{   
public:
    
    /** \todo Not yet fully working...
     * May not be necessary any more, since the Luo-Rudy model seems to work.
     */


   // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    void TestMonodomainFitzHughNagumoWithEdgeStimulus( void )
    {   
        FhnEdgeStimulusCellFactory cell_factory;
        
        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(2);   // 2 ms
        monodomain_problem.SetOutputDirectory("testoutput/FhnWithEdgeStimulus");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainFhn_2dWithEdgeStimulus");

        monodomain_problem.Initialise();
        
        monodomain_problem.Solve();
        
        double* voltage_array;
        int lo, hi;
        monodomain_problem.GetVoltageArray(&voltage_array, lo, hi); 
    
 
        int num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);

        if (num_procs == 1)
        {
            /*
             * Test the top right node against the right one in the 1D case, 
             * comparing voltage, and then test all the nodes on the right hand 
             * side of the square against the top right one, comparing voltage.
             */
            bool need_initialisation = true;
            double voltage;

            need_initialisation = true;

            // Test the RHS of the mesh
            for (int i = 0; i < monodomain_problem.rGetMesh().GetNumNodes(); i++)
            {
                if (monodomain_problem.rGetMesh().GetNodeAt(i)->GetPoint()[0] == 0.1)
                {
                    // x = 0 is where the stimulus has been applied
                    // x = 0.1cm is the other end of the mesh and where we want to 
                    //       to test the value of the nodes
                    
                    if (need_initialisation)
                    {
                        voltage = voltage_array[i];
                        need_initialisation = false;
                    }
                    else
                    {
                        // Tests the final voltages for all the RHS edge nodes
                        // are close to each other.
                        // This works as we are using the 'criss-cross' mesh,
                        // the voltages would vary more with a mesh with all the
                        // triangles aligned in the same direction.
      //
        //\todo  - Make this test something sensible
                       //   TS_ASSERT_DELTA(voltage_array[i], voltage, 0.01);

                       // std::cout << "y=" << monodomain_problem.mMesh.GetNodeAt(i)->GetPoint()[1] << std::endl;
                    }
                    
                    
                  
        //
        //\todo  - Make this test something sensible
                //    TS_ASSERT_DELTA(voltage_array[i], -65.0087, 0.01);
                }
            }
        }
        monodomain_problem.RestoreVoltageArray(&voltage_array);
        
     
    }   



}; 
#endif //_TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_
