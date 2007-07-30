#ifndef TESTCARDIACMECHANICSASSEMBLER_HPP_
#define TESTCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "CardiacMechanicsAssembler.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "FiniteElasticityTools.hpp"


class TestCardiacMechanicsAssembler : public CxxTest::TestSuite
{
public :
    void TestExceptions() throw(Exception)
    {
    }
    
    
    void TestWithZeroActiveTension() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/ZeroActiveTension");
        
        std::vector<double> active_tensions(mesh.n_vertices(), 0.0);
        cardiac_mech_assembler.SetActiveTensions(active_tensions);

        cardiac_mech_assembler.Solve();
        
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumNewtonIterations(), 0u);
    }
};
#endif /*TESTCARDIACMECHANICSASSEMBLER_HPP_*/
