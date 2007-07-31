#ifndef TESTCARDIACELECTROMECHANICS_HPP_
#define TESTCARDIACELECTROMECHANICS_HPP_


#include <cxxtest/TestSuite.h>
#include "CardiacMechanicsAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "FiniteElasticityTools.hpp"
#include "DealiiMonodomainAssembler.hpp" // to be replaced by [Mono/Bi]Dg0Assmebler...
#include "InitialStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "AbstractCardiacCell.hpp"
#include "LuoRudyIModel1991OdeSystem.hpp"

class TestCardiacElectroMechanics : public CxxTest::TestSuite
{
public:
    // basic test of the dealii monodomain assembler (which is only going to be used
    // until dealii can be connected to chaste)
    void TestMonodomainDealiiAssembler() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 0.1);
        mesh.refine_global(4);
        
        AbstractStimulusFunction* p_stimulus = new InitialStimulus(-1e6, 0.5);;
        AbstractStimulusFunction* p_zero_stim = new ZeroStimulus;
        EulerIvpOdeSolver euler_solver;
        
        double time_step = 0.01;
        
        std::vector< AbstractCardiacCell* > cells(mesh.n_vertices());
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.ReachedEnd())
        {
            if(vertex_iter.GetVertex()[0]==0)
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new LuoRudyIModel1991OdeSystem(&euler_solver, time_step, p_stimulus);
            }
            else
            {
                cells[vertex_iter.GetVertexGlobalIndex()] = new LuoRudyIModel1991OdeSystem(&euler_solver, time_step, p_zero_stim);
            }
            
            vertex_iter.Next();
        }
        
        DealiiMonodomainAssembler<2> monodomain_assembler(&mesh, cells);
        
        double end_time = 2;
        unsigned num_timesteps = (unsigned)floor(end_time/time_step+0.5);
        monodomain_assembler.Solve(0, end_time, num_timesteps);
        
        for(unsigned i=0; i<cells.size(); i++)
        {
            delete cells[i];
        }
        
        TS_ASSERT_DELTA(monodomain_assembler.GetVoltage()[1], 41.8787, 1e-3);
    }
};
#endif /*TESTCARDIACELECTROMECHANICS_HPP_*/
