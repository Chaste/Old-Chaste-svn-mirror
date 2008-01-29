#ifndef _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_
#define _TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>


#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"

#include "FitzHughNagumo1961OdeSystem.hpp"
#include "ReplicatableVector.hpp"


class FhnEdgeStimulusCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    InitialStimulus *mpStimulus;
public:
    FhnEdgeStimulusCellFactory() : AbstractCardiacCellFactory<2>(0.01)
    {
        mpStimulus = new InitialStimulus(-10.0, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new FitzHughNagumo1961OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new FitzHughNagumo1961OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
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

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    void TestMonodomainFitzHughNagumoWithEdgeStimulus( void )
    {
        FhnEdgeStimulusCellFactory cell_factory;
        
        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/2D_0_to_1mm_400_elements");
        monodomain_problem.SetEndTime(1.2);   // 1.2 ms
        monodomain_problem.SetOutputDirectory("FhnWithEdgeStimulus");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainFhn_2dWithEdgeStimulus");
        
        monodomain_problem.Initialise();
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);
        monodomain_problem.GetMonodomainPde()->SetIntracellularConductivityTensor(0.01*identity_matrix<double>(2));
        
        monodomain_problem.Solve();
        
        
        /*
        * Test the top right node against the right one in the 1D case, 
        * comparing voltage, and then test all the nodes on the right hand 
        * side of the square against the top right one, comparing voltage.
        */
        bool need_initialisation = true;
        double probe_voltage=0.0;
        Vec voltage=monodomain_problem.GetVoltage();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);
        need_initialisation = true;
        
        // Test the RHS of the mesh
        for (unsigned i = 0; i < monodomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (monodomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes
                
                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[i];
                    need_initialisation = false;
                }
                else
                {
                    // Tests the final voltages for all the RHS edge nodes
                    // are close to each other.
                    // This works as we are using the 'criss-cross' mesh,
                    // the voltages would vary more with a mesh with all the
                    // triangles aligned in the same direction.
                    
                    TS_ASSERT_DELTA(voltage_replicated[i], probe_voltage, 1e-10);
                }
                
                TS_ASSERT_DELTA(voltage_replicated[i], 0.1394, 1e-4);
            }
        }
        
        
        
    }
    
    
    
};
#endif //_TESTMONODOMAINFITZHUGHNAGUMOWITHDG0ASSEMBLER_HPP_
