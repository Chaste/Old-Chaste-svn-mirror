#ifndef TESTCARDIACMECHANICSASSEMBLER_HPP_
#define TESTCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "CardiacMechanicsAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "FiniteElasticityTools.hpp"


class TestCardiacMechanicsAssemblerLong : public CxxTest::TestSuite
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
    void TestSpecifiedActiveTensionCompression3d() throw(Exception)
    {
        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<3> zero;
        FiniteElasticityTools<3>::FixFacesContainingPoint(mesh, zero);
        
        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw<3> material_law(0.02, 0.0);

        CardiacMechanicsAssembler<3> cardiac_mech_assembler(&mesh, 
                                                            "CardiacMech/SpecifiedActiveTensionCompression3d",
                                                            &material_law);

        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumQuadPointsPerElement(), 27u);
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetTotalNumQuadPoints(), 27u*mesh.n_active_cells());

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        SetUpLinearActiveTension<3>(mesh, 0.01, active_tension); 
        

        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        cardiac_mech_assembler.StaticSolve();
        
        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 1 is a corner node, 
        // The deformation is not that large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1), 1.0012, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), 0.0020, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[2](1), 0.0009, 1e-3);
        
        std::vector<double>& lambda = cardiac_mech_assembler.rGetLambda();
        std::vector<std::vector<double> > quad_points 
            = FiniteElasticityTools<3>::GetQuadPointPositions(mesh,3);
        
        
        // the lambdas should be less than 1 (positive T_a => compression), and also
        // should be near the same for any particular value of y, ie the same along any 
        // fibre. Lambda should decrease approx linearly with y. Uncomment trace and 
        // view in matlab (plot y against lambda) to observe this. The parameters 
        // 0.07, 0.01 etc were obtained by looking at the plot. 
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[i][1];
            double mid = 1 - 0.07*y;
            double range;
            if(y>0.3)
            {
                range = 0.01; // lambdas very close together away from y=0
            }
            else
            {
                range = 0.04;  // more spread in lambda nearer y=0
            }
            
            TS_ASSERT_LESS_THAN(lambda[i], mid + range);
            TS_ASSERT_LESS_THAN(mid - range, lambda[i]);
            
            // don't delete:
            //std::cout << quad_points[i][0] << " " << quad_points[i][1] << " " << quad_points[i][2] << " " << lambda[i] << "\n";
        }
        
    }
};
#endif /*TESTCARDIACMECHANICSASSEMBLER_HPP_*/
