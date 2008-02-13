#ifndef _TESTNEUMANNHEARTSTIMULUS_HPP_
#define _TESTNEUMANNHEARTSTIMULUS_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"
#include "ReplicatableVector.hpp"


class NoStimulus2dCellFactory : public AbstractCardiacCellFactory<1>
{

public:
    NoStimulus2dCellFactory() : AbstractCardiacCellFactory<1>(0.01)
    {
    }
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {

            return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);

    }
    
    ~NoStimulus2dCellFactory(void)
    {
    }
};



class TestNeumannHeartStimulus : public CxxTest::TestSuite
{
public:

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    void TestMonodomainConstantStimulus() throw(Exception)
    {
        NoStimulus2dCellFactory cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1mm_10_elements");
        monodomain_problem.SetEndTime(2);   // ms
        monodomain_problem.SetOutputDirectory("MonoNeumanConst");
        monodomain_problem.SetOutputFilenamePrefix("NewMonodomainLR91_1d");        
        monodomain_problem.SetIntracellularConductivities(0.0005);        
        monodomain_problem.Initialise();
        
        // create boundary conditions container
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1> *p_bc_stim = new ConstBoundaryCondition<1>(2);
        ConstBoundaryCondition<1> *p_bc_no_flux = new ConstBoundaryCondition<1>(0);
                
        // get mesh
        ConformingTetrahedralMesh<1,1> &mesh = monodomain_problem.rGetMesh();
        // loop over boundary elements
        ConformingTetrahedralMesh<1, 1>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            if (((*iter)->GetNodeLocation(0))[0]==0.0)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_stim);
            }
            else
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bc_no_flux);
            }
            
            //AddDirichletBoundaryCondition(*iter, p_boundary_condition, indexOfUnknown);
            iter++;
        }
            // if the element is on the left of the mesh, add a stimulus to the bcc
        
        // pass the bcc to the monodomain problem
        monodomain_problem.SetBoundaryConditionsContainer(&bcc);
        
        monodomain_problem.GetMonodomainPde()->SetSurfaceAreaToVolumeRatio(1.0);
        monodomain_problem.GetMonodomainPde()->SetCapacitance(1.0);
        
        monodomain_problem.Solve();
        
        // check some voltages    
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double atol=5e-3;
        
        TS_ASSERT_DELTA(voltage_replicated[1], 94.6426, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 49.7867, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 30.5954, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 21.6782, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -33.9983, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -52.2396, atol);
              
    }
};

#endif //_TESTNEUMANNHEARTSTIMULUS_HPP_
