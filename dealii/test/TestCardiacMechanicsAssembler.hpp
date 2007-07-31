#ifndef TESTCARDIACMECHANICSASSEMBLER_HPP_
#define TESTCARDIACMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "CardiacMechanicsAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "FiniteElasticityTools.hpp"


class TestCardiacMechanicsAssembler : public CxxTest::TestSuite
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
    
    template<unsigned DIM>
    std::vector<std::vector<double> >  GetQuadPointPositions(Triangulation<DIM>& rMesh)
    {
        QGauss<DIM> quadrature_formula(3);
        unsigned n_q_points = quadrature_formula.n_quadrature_points;

        std::vector<std::vector<double> > ret(DIM);
        for(unsigned i=0; i<DIM; i++)
        {
            ret[i].resize(n_q_points * rMesh.n_active_cells());
        }
        
        
        FE_Q<DIM> fe(2);
        FEValues<DIM> fe_values(fe, quadrature_formula,
                                UpdateFlags(update_q_points));
        
        unsigned current=0;                        
        for(typename Triangulation<DIM>::cell_iterator element_iter = rMesh.begin_active(); 
            element_iter!=rMesh.end();
            element_iter++)
        {
            fe_values.reinit(element_iter);
            
            std::vector<Point<DIM> > quad_points = fe_values.get_quadrature_points();

            for(unsigned q=0; q<quad_points.size(); q++)
            {
                for(unsigned i=0; i<DIM; i++)
                {
                    ret[i][current] = quad_points[q][i];
                }
                current++;
            }
        }
        
        return ret;
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

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, 
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

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/ZeroActiveTension");
        
        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();
        
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumNewtonIterations(), 0u);
    }
    

    void TestSpecifiedActiveTensionCompression() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        
        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);
        
        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw<2> material_law(0.02);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, 
                                                            "CardiacMech/SpecifiedActiveTensionCompression",
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
        
        std::vector<double>& lambda = cardiac_mech_assembler.GetLambda();
        std::vector<std::vector<double> > quad_points = GetQuadPointPositions<2>(mesh);
        
        // the lambdas should be less than 1 (positive T_a => compression), and also
        // should be near the same for any particular value of y, ie the same along any 
        // fibre. Lambda should decrease approx linearly with y. Uncomment trace and 
        // view in matlab (plot y against lambda) to observe this. The parameters 
        // 0.35,0.05 etc were obtained by looking at the plot. 
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[1][i];
            double mid = 1 - 0.35*y;
            double range = 0.05 + 0.07*y;
            
            TS_ASSERT_LESS_THAN(lambda[i], mid + range);
            TS_ASSERT_LESS_THAN(mid - range, lambda[i]);
            
            // don't delete:
            //std::cout << quad_points[0][i] << " " << quad_points[1][i] << " " << lambda[i] << "\n";
        }
        
        // FIXME: try running solve again here - should be instantaneous but isn't
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

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, 
                                                           "CardiacMech/SpecifiedActiveTensionStretching",
                                                            &material_law);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        SetUpLinearActiveTension<2>(mesh, -0.05, active_tension); // doesn't converge if -0.1
        cardiac_mech_assembler.SetActiveTension(active_tension);

        cardiac_mech_assembler.Solve();

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 1 is the bottom-right corner node, 
        // and the deformation is quite large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1),  0.9872, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), -0.1000, 1e-3);
        
        std::vector<double>& lambda = cardiac_mech_assembler.GetLambda();
        std::vector<std::vector<double> > quad_points = GetQuadPointPositions<2>(mesh);
        
        // the lambdas should be greater than 1 (negative T_a => stretch), and also
        // should be near the same for any particular value of y, ie the same along any 
        // fibre. Lambda should decrease approx linearly with y. Uncomment trace and 
        // view in matlab (plot y against lambda) to observe this. The parameters 
        // 0.29, 0.04 were obtained by looking at the plot. 
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[1][i];
            double mid = 1 + 0.29*y;
            double range = 0.04;
            
            TS_ASSERT_LESS_THAN(lambda[i], mid + range);
            TS_ASSERT_LESS_THAN(mid - range, lambda[i]);
            
            // don't delete:
            //std::cout << quad_points[0][i] << " " << quad_points[1][i] << " " << lambda[i] << "\n";
        }
    }
};
#endif /*TESTCARDIACMECHANICSASSEMBLER_HPP_*/
