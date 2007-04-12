#ifndef _TESTMONODOMAINNOSTIMULUS_HPP_
#define _TESTMONODOMAINNOSTIMULUS_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "petscvec.h"
#include <vector>

#include "PropagationPropertiesCalculator.hpp"
#include "ColumnDataReader.hpp"

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class ZeroStimulusCellFactory : public AbstractCardiacCellFactory<1>
{
public:
    ZeroStimulusCellFactory(double timeStep) : AbstractCardiacCellFactory<1>(timeStep)
    {}
    
    AbstractCardiacCell* CreateCardiacCellForNode(unsigned node)
    {
        return new LuoRudyIModel1991OdeSystem(mpSolver, mTimeStep, mpZeroStimulus);
        
    }
};


/* TestMonodomainNoStimulus - based on TestMonodomainConductionVelocity
 *
 * No initial stimulus applied.
 * Check that the voltage of all cells is constant thoroughout the mesh
 * at any point in time and never lower than the resting potential
 * of the LR cell = -85.0 mV
 *
 * Best run with optimisation on.
 */
class TestMonodomainNoStimulus : public CxxTest::TestSuite
{
public:

    void TestZeroStimulus()
    {
        ZeroStimulusCellFactory cell_factory(0.01); // ODE time step (ms)
        MonodomainProblem<1> monodomain_problem(&cell_factory);
        
        monodomain_problem.SetMeshFilename("mesh/test/data/1D_0_to_1_20_elements");
        monodomain_problem.SetEndTime(30);   // 30 ms
        monodomain_problem.SetOutputDirectory("MonoNoStim");
        monodomain_problem.SetOutputFilenamePrefix("MonodomainNoStimLR91_1d");
        monodomain_problem.Initialise();

        monodomain_problem.Solve();
        
        ReplicatableVector voltage_replicated(monodomain_problem.GetVoltage());
        double constant_voltage = voltage_replicated[0];
        TS_ASSERT_LESS_THAN(-85.0, constant_voltage);
        
        for (unsigned index=0; index<voltage_replicated.size() ; index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index] , constant_voltage, 1E-5);
        }
    }
};
#endif //_TESTMONODOMAINNOSTIMULUS_HPP_
