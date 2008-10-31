/*

Copyright (C) University of Oxford, 2008

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


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

public :
    void TestCompareJacobians() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/SpecifiedActiveTensionStretching");

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);

        for(unsigned i=0; i<active_tension.size(); i++)
        {
            active_tension[i] = 0.1;
        }

        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        // tol for the default pole-zero law
        TS_ASSERT_THROWS_NOTHING( cardiac_mech_assembler.CompareJacobians(2e-4) );
        // tol for mooney-riv(0.02)
        //TS_ASSERT_THROWS_NOTHING( cardiac_mech_assembler.CompareJacobians(2e-6) );
    }


    void TestWithZeroActiveTension() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/ZeroActiveTension");

        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumQuadPointsPerElement(), 9u);
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetTotalNumQuadPoints(), mesh.n_active_cells()*9u);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        cardiac_mech_assembler.StaticSolve(); // the times are unused as explicit

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
        MooneyRivlinMaterialLaw2<2> material_law(2);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh,
                                                            "CardiacMech/SpecifiedActiveTensionCompression",
                                                            &material_law);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        SetUpLinearActiveTension<2>(mesh, 5, active_tension);


        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        cardiac_mech_assembler.StaticSolve();

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 1 is the bottom-right corner node,
        // and the deformation is reasonably large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1), 0.9994, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), 0.1154, 1e-3);

        std::vector<double>& lambda = cardiac_mech_assembler.rGetLambda();
        std::vector<std::vector<double> > quad_points
           = FiniteElasticityTools<2>::GetQuadPointPositions(mesh,3);

        // the lambdas should be less than 1 (positive T_a => compression), and also
        // should be near the same for any particular value of y, ie the same along any
        // fibre. Lambda should decrease approx linearly with y. Uncomment trace and
        // view in matlab (plot y against lambda) to observe this. The parameters
        // 0.35,0.05 etc were obtained by looking at the plot.
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[i][1];
            double mid = 1 - 0.25*y;
            double range = 0.03;

            TS_ASSERT_LESS_THAN(lambda[i], mid + range);
            TS_ASSERT_LESS_THAN(mid - range, lambda[i]);

            // don't delete:
            //std::cout << quad_points[i][0] << " " << quad_points[i][1] << " " << lambda[i] << "\n";
        }

        // hardcoded test
        TS_ASSERT_DELTA(lambda[34], 0.9665, 1e-4);
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
        MooneyRivlinMaterialLaw2<2> material_law(2);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh,
                                                           "CardiacMech/SpecifiedActiveTensionStretching",
                                                            &material_law);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.0);
        SetUpLinearActiveTension<2>(mesh, -2.5, active_tension); // doesn't converge if -0.1
        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        cardiac_mech_assembler.StaticSolve();

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 1 is the bottom-right corner node,
        // and the deformation is reasonably large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](1),  0.9895, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](1), -0.0694, 1e-3);

        std::vector<double>& lambda = cardiac_mech_assembler.rGetLambda();
        std::vector<std::vector<double> > quad_points
           = FiniteElasticityTools<2>::GetQuadPointPositions(mesh,3);


        // the lambdas should be greater than 1 (negative T_a => stretch), and also
        // should be near the same for any particular value of y, ie the same along any
        // fibre. Lambda should increase approx linearly with y. Uncomment trace and
        // view in matlab (plot y against lambda) to observe this. The parameters
        // 0.29, 0.04 were obtained by looking at the plot.
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[i][1];
            double mid = 1 + 0.15*y;
            double range = 0.02;

            TS_ASSERT_LESS_THAN(lambda[i], mid + range);
            TS_ASSERT_LESS_THAN(mid - range, lambda[i]);

            // don't delete:
            //std::cout << quad_points[i][0] << " " << quad_points[i][1] << " " << lambda[i] << "\n";
        }

        // hardcoded test
        TS_ASSERT_DELTA(lambda[34], 1.0179, 1e-4);
    }


    // constant active tension and compression of fibres set to be in the y-direction
    void TestConstActiveTensionCompressionFibreAt90Deg() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        // fix the x=0 surface
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 1, 0.0);

        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw2<2> material_law(2);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh,
                                                           "CardiacMech/SpecifiedActiveTensionComp90Deg",
                                                            &material_law);

        Tensor<2,2> zero_mat;
        // not orthogonal
        TS_ASSERT_THROWS_ANYTHING(cardiac_mech_assembler.SetFibreSheetMatrix(zero_mat));

        Tensor<2,2> P;
        P[0][0] =  0;
        P[0][1] =  1;
        P[1][0] = -1;
        P[1][1] =  0;

        cardiac_mech_assembler.SetFibreSheetMatrix(P);


        // constant active tension this time
        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 5);
        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        cardiac_mech_assembler.StaticSolve();

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 2 is the top-right corner node,
        // and the deformation is reasonably large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](2), 1.1952, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](2), 0.7428, 1e-3);

        std::vector<double>& lambda = cardiac_mech_assembler.rGetLambda();
        std::vector<std::vector<double> > quad_points
           = FiniteElasticityTools<2>::GetQuadPointPositions(mesh,3);

        // the lambdas should be less than 1 and approx constant away from the fixed boundary
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[i][1];

            if (y>0.5) // ie away from the fixed boundary y=0
            {
                TS_ASSERT_LESS_THAN(lambda[i], 0.8);
                TS_ASSERT_LESS_THAN(0.69, lambda[i]);
            }

            // don't delete:
            //std::cout << quad_points[i][0] << " " << quad_points[i][1] << " " << lambda[i] << "\n";
        }

        // hardcoded test
        TS_ASSERT_DELTA(lambda[34], 0.7925, 1e-4);
    }



    // constant active tension and compression of fibres set to be at 45 degrees
    void TestConstActiveTensionCompressionFibreAt45Deg() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        // fix the x=0 surface
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 1, 0.0);

        // specify this material law so the test continues to pass when the default
        // material law is changed.
        MooneyRivlinMaterialLaw2<2> material_law(2);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh,
                                                           "CardiacMech/SpecifiedActiveTensionComp45Deg",
                                                            &material_law);

        Tensor<2,2> P;
        double theta = 0.785398163; //pi/4;
        P[0][0] =  cos(theta);
        P[0][1] =  sin(theta);
        P[1][0] = -sin(theta);
        P[1][1] =  cos(theta);

        cardiac_mech_assembler.SetFibreSheetMatrix(P);


        // constant active tension this time
        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 5);
        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        cardiac_mech_assembler.StaticSolve();

        // have visually checked the answer and seen that it looks ok, so have
        // a hardcoded test here. Node that 2 is the top-right corner node,
        // and the deformation is reasonably large
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[0](2), 1.6276, 1e-3);
        TS_ASSERT_DELTA( cardiac_mech_assembler.rGetDeformedPosition()[1](2), 0.8998, 1e-3);

        std::vector<double>& lambda = cardiac_mech_assembler.rGetLambda();
        std::vector<std::vector<double> > quad_points
           = FiniteElasticityTools<2>::GetQuadPointPositions(mesh,3);

        // the lambdas should be less than 1 and approx constant away from the fixed boundary
        for(unsigned i=0; i<lambda.size(); i++)
        {
            double y = quad_points[1][i];

            if (y>0.5) // ie away from the fixed boundary y=0
            {
                TS_ASSERT_LESS_THAN(lambda[i], 0.751);  // varies much less than previous test for some reason..?
                TS_ASSERT_LESS_THAN(0.73, lambda[i]);
            }

            // don't delete:
            //std::cout << quad_points[0][i] << " " << quad_points[1][i] << " " << lambda[i] << "\n";
        }

        // hardcoded test
        TS_ASSERT_DELTA(lambda[34], 0.7773, 1e-4);
    }

    void TestWithPoleZeroLaw() throw(Exception)
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechanicsAssembler<2> cardiac_mech_assembler(&mesh, "CardiacMech/PoleZeroConstActiveTension");

        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumQuadPointsPerElement(), 9u);
        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetTotalNumQuadPoints(), mesh.n_active_cells()*9u);

        std::vector<double> active_tension(cardiac_mech_assembler.GetTotalNumQuadPoints(), 0.1);
        cardiac_mech_assembler.SetForcingQuantity(active_tension);

        // just test it run ok
        cardiac_mech_assembler.StaticSolve();

        TS_ASSERT_EQUALS(cardiac_mech_assembler.GetNumNewtonIterations(), 3u);
    }

    void TestWithScalingOfNashHunterPoleZeroLaw2()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        Point<2> zero;
        FiniteElasticityTools<2>::FixFacesContainingPoint(mesh, zero);

        CardiacMechanicsAssembler<2> assembler(&mesh,"");

        std::vector<double> active_tension(assembler.GetTotalNumQuadPoints(), 0.1);
        assembler.SetForcingQuantity(active_tension);
        assembler.Solve(0,0.01,0.01);

        CardiacMechanicsAssembler<2> assembler_with_scaling(&mesh,"");

        assembler_with_scaling.SetForcingQuantity(active_tension);
        assembler_with_scaling.SetScaling(0.5);
        assembler_with_scaling.Solve(0,0.01,0.01);

        std::vector<Vector<double> >& position1
            = assembler.rGetDeformedPosition();

        std::vector<Vector<double> >& position2
            = assembler_with_scaling.rGetDeformedPosition();

        for(unsigned i=0; i<position1[0].size(); i++)
        {
            TS_ASSERT_DELTA(position1[0](i), position2[0](i), 1.1e-5);
            TS_ASSERT_DELTA(position1[1](i), position2[1](i), 1.1e-5);
        }
    }
};
#endif /*TESTCARDIACMECHANICSASSEMBLER_HPP_*/
