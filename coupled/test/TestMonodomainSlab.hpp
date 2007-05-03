#ifndef _TESTMONODOMAINSLAB_HPP_
#define _TESTMONODOMAINSLAB_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "ReplicatableVector.hpp"

class FaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    InitialStimulus *mpStimulus;
public:
    FaceStimulusCellFactory() : AbstractCardiacCellFactory<3>(0.01)
    {
        mpStimulus = new InitialStimulus(-600.0*1000, 0.5);
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        if (mpMesh->GetNode(node)->GetPoint()[0] == 0.0)
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpStimulus);
        }
        else
        {
            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        }
    }
    
    ~FaceStimulusCellFactory(void)
    {
        delete mpStimulus;
    }
};



class TestMonodomainSlab : public CxxTest::TestSuite
{

public:

    // Solve on a 3D 1mm by 1mm by 1mm mesh (space step = 0.1mm), stimulating
    // the left face.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp
    void TestMonodomainDg03DWithFaceStimulus( void )
    {
        FaceStimulusCellFactory cell_factory;
        
        MonodomainProblem<3> monodomain_problem(&cell_factory);
        
        monodomain_problem.SetMeshFilename("mesh/test/data/3D_0_to_1mm_6000_elements");
        monodomain_problem.SetEndTime(4);   // 4 ms
        monodomain_problem.SetOutputDirectory("MonoDg03dWithFaceStimulus");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_3dWithFaceStimulus");
        monodomain_problem.Initialise();
        
        monodomain_problem.Solve();
        
        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);
        
        /*
         * Test the top right node against the right one in the 1D case, 
         * comparing voltage, and then test all the nodes on the right hand 
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage;
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        need_initialisation = true;
        
        // Test the RHF of the mesh
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
                    TS_ASSERT_DELTA(voltage_replicated[i], probe_voltage, 1.1);
                    // std::cout << "y=" << monodomain_problem.mMesh.GetNode(i)->GetPoint()[1] << std::endl;
                }
                
                // Check against 1d case - the TestMonodomainDg01D test is run
                // for 4ms 
                
                TS_ASSERT_DELTA(voltage_replicated[i], 21.3, 1);
                
            }
        }
    }
};


#endif //_TESTMONODOMAINSLAB_HPP_
