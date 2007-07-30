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
    

    void TestCompareJacobians() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/ZeroActiveTension");
        
        std::vector<double> active_tension(mesh.n_vertices(), 0.0);
        for(TriangulationVertexIterator<2> vertex_iter(&mesh); !vertex_iter.ReachedEnd(); vertex_iter.Next())
        {
            double x = vertex_iter.GetVertex()[0];
            active_tension[vertex_iter.GetVertexGlobalIndex()] = 0.1*x;
        }
        cardiac_mech_assembler.SetActiveTension(active_tension);

        TS_ASSERT_THROWS_NOTHING( cardiac_mech_assembler.CompareJacobians(1.5e-7) );
    }

    
    void TestWithZeroActiveTension() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/ZeroActiveTension");
        
        std::vector<double> active_tension(mesh.n_vertices(), 0.0);
        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();
        
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumNewtonIterations(), 0u);
    }
    

    void TestSpecifiedActiveTension() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);
        
        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw<2> material_law(0.02);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, 
                                                            "CardiacMech/SpecifiedActiveTension",
                                                            &material_law);
        
        std::vector<double> active_tension(mesh.n_vertices(), 0.0);
        for(TriangulationVertexIterator<2> vertex_iter(&mesh); !vertex_iter.ReachedEnd(); vertex_iter.Next())
        {
            double x = vertex_iter.GetVertex()[0];
            active_tension[vertex_iter.GetVertexGlobalIndex()] = 0.1*x;
        }
        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();
        
        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node 1 is the bottom-right corner node.
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1), 1.3793, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), 0.2051, 1e-3);
    }

    
};
#endif /*TESTCARDIACMECHANICSASSEMBLER_HPP_*/
