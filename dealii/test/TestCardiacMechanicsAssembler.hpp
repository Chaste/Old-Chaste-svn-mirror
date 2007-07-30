#ifndef TESTCARDIACMECHANICSASSEMBLER_HPP_
#define TESTCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "CardiacMechAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "FiniteElasticityTools.hpp"


class TestCardiacMechAssembler : public CxxTest::TestSuite
{

private:
    // little helper method
    // set up a active tension that is constant along any fibre (indep of x), but grows linearly with y
    template<unsigned DIM>
    void SetUpLinearActiveTension(Triangulation<DIM>& rMesh, double value, std::vector<double>& rActiveTension)
    {
        unsigned current = 0;   
        for(typename Triangulation<DIM>::cell_iterator element_iter = rMesh.begin_active(); 
            element_iter!=rMesh.end();
            element_iter++)
        {
            double y = element_iter->vertex(0)[1];
            for(unsigned q=0; q<pow(3,DIM); q++) // assumes there's 3 quad points in each direction
            {
                rActiveTension[current++] = value*y;
            }
        }
    }
    
public :    
    void TestCompareJacobians() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw<2> material_law(0.02);

        CardiacMechAssembler<2> cardiac_mech_assembler(&mesh, 
                                                       "CardiacMech/SpecifiedActiveTensionStretching",
                                                       &material_law);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        
        for(unsigned i=0; i<active_tension.size(); i++)
        {
            active_tension[i] = 0.1;
        }

        cardiac_mech_assembler.SetActiveTension(active_tension);

        TS_ASSERT_THROWS_NOTHING( cardiac_mech_assembler.CompareJacobians(2e-6) );
    }

    
    void TestWithZeroActiveTension() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/ZeroActiveTension");
        
        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();
        
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumNewtonIterations(), 0u);
    }
    

    void TestSpecifiedActiveTensionStretching() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw<2> material_law(0.02);

        CardiacMechAssembler<2> cardiac_mech_assembler(&mesh, 
                                                       "CardiacMech/SpecifiedActiveTensionStretching",
                                                       &material_law);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        SetUpLinearActiveTension<2>(mesh, 0.1, active_tension); 
        

        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();
        
        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 1 is the bottom-right corner node, 
        // and the deformation is quite large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1), 0.9782, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), 0.2777, 1e-3);
        
        // FIXME: try running solve again here - should be instantaneous but isn't
    }

    void TestSpecifiedActiveTensionSquashing() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw<2> material_law(0.02);

        CardiacMechAssembler<2> cardiac_mech_assembler(&mesh, 
                                                       "CardiacMech/SpecifiedActiveTensionSquashing",
                                                       &material_law);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        SetUpLinearActiveTension<2>(mesh, -0.05, active_tension);
        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 1 is the bottom-right corner node, 
        // and the deformation is quite large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1),  0.9872, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), -0.1000, 1e-3);
    }
};
#endif /*TESTCARDIACMECHANICSASSEMBLER_HPP_*/
