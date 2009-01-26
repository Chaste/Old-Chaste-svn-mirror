/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef TESTFINITEELASTICITYASSEMBLERLONG_HPP_
#define TESTFINITEELASTICITYASSEMBLERLONG_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"
#include "MooneyRivlinMaterialLaw2.hpp"
#include "PolynomialMaterialLaw3d2.hpp"
#include "ExponentialMaterialLaw2.hpp"
#include "FiniteElasticityTools.hpp"


double MATERIAL_PARAM = 0.05;
double ALPHA = 0.2;

double SHEARS_ALPHA = 0.1;
double SHEARS_C1 = 0.01;

Vector<double> TractionForMixedShears(const Point<2>& X)
{
    assert(X[0]==1 || X[1]==0 || X[1]==1);
    
    Vector<double> traction(2);
    if(X[1]==1)
    {
        traction(0) = -2*SHEARS_ALPHA * SHEARS_C1;
        traction(1) =  0;
    }
    else if (X[1]==0)
    {
        traction(0) = 2*SHEARS_ALPHA*SHEARS_C1;
        traction(1) = 0;
    }
    else if (X[0]==1)
    {
        traction(0) = 0;
        traction(1) = -2*SHEARS_ALPHA*SHEARS_C1; // also equal to -2*c2*beta;
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}

Vector<double> MyBodyForce(const Point<2>& X)
{
    assert(X[0]>=0 && X[0]<=1 && X[1]>=0 && X[1]<=1);

    Vector<double> body_force(2);
    double lam = 1+ALPHA*X[0];
    body_force(0) = -2*MATERIAL_PARAM * ALPHA;
    body_force(1) = -2*MATERIAL_PARAM * 2*ALPHA*ALPHA*X[1]/(lam*lam*lam);
    return body_force;
}

Vector<double> MyTraction(const Point<2>& X)
{
    Vector<double> traction(2);
    
    double lam = 1+ALPHA*X[0];
    if(X[0]==1)
    {
        traction(0) =  2*MATERIAL_PARAM * (lam - 1.0/lam);
        traction(1) = -2*MATERIAL_PARAM * X[1]*ALPHA/(lam*lam);
    }
    else if(X[1]==0)
    {
        traction(0) =  2*MATERIAL_PARAM * X[1]*ALPHA/(lam*lam);
        traction(1) = -2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else if(X[1]==1)
    {
        traction(0) = -2*MATERIAL_PARAM * X[1]*ALPHA/(lam*lam);
        traction(1) =  2*MATERIAL_PARAM * (-lam + 1.0/lam);
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}


// longer tests, including 3d tests, for nightly runs
class TestFiniteElasticityAssemblerLong : public CxxTest::TestSuite
{
public :

    // Run same simulation on two meshes (one more refined than the other)
    // and test they agree on shared gridpoints
    void Test2dProblemOnSquareForConvergence() throw(Exception)
    {
        ////////////////////////////////////////////////
        // run 1: on a mesh which is refined 3 times..
        ////////////////////////////////////////////////

        // note small value of body force and mooney-rivlin const (as if rescaled)
        // - speeds up GMRES
        Vector<double> body_force(2);
        body_force(0) = 0.06;
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(0.02);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);


        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/simple2d");


        finite_elasticity.StaticSolve();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = finite_elasticity.rGetDeformedPosition();

        //////////////////////////////////////////////////////////////
        // run 2: same problem, on a mesh which is refined 4 times..
        //////////////////////////////////////////////////////////////
        Triangulation<2> mesh_refined;
        GridGenerator::hyper_cube(mesh_refined, 0.0, 1.0);
        mesh_refined.refine_global(4);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh_refined, 0, 0.0);

        FiniteElasticityAssembler<2> finite_elasticity_ref(&mesh_refined,
                                                           &mooney_rivlin_law,
                                                           body_force,
                                                           1.0,
                                                           "finite_elas/simple2d");
        finite_elasticity_ref.StaticSolve();

        // get deformed position
        std::vector<Vector<double> >& deformed_position_ref
            = finite_elasticity_ref.rGetDeformedPosition();


        //////////////////////////////////////////////////////////////
        // compare the solution with that on the previous mesh
        //////////////////////////////////////////////////////////////
        for (unsigned i=0; i<deformed_position[0].size(); i++)
        {
            TS_ASSERT_DELTA(deformed_position[0](i), deformed_position_ref[0](i), 1e-2);
            TS_ASSERT_DELTA(deformed_position[1](i), deformed_position_ref[1](i), 1e-2);
        }

        // check nothing has changed
        TS_ASSERT_DELTA( deformed_position[0](6), 1.2158, 1e-3);
        TS_ASSERT_DELTA( deformed_position[1](6), 0.5,    1e-3);
    }


    // 3d simulation with result compared to alternative simulation that
    // uses linear bases.
    // TODO: Ideally, we would want to compare to an identical simulation
    // which uses quadratics, so that the answers should be identical
    void Test3dProblemOnCube() throw(Exception)
    {
        Vector<double> body_force(3);
        body_force(1) = 20;

        MooneyRivlinMaterialLaw2<3> mooney_rivlin_law(1,2);

        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 0.1);
        mesh.refine_global(3);
        FiniteElasticityTools<3>::SetFixedBoundary(mesh,0,0.0);

        FiniteElasticityAssembler<3> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/simple3d");
        finite_elasticity.StaticSolve();

        // get undeformed position
        std::vector<Vector<double> >& undeformed_position
            = finite_elasticity.rGetUndeformedPosition();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = finite_elasticity.rGetDeformedPosition();

        TS_ASSERT_EQUALS(deformed_position.size(),3u);


        // compare against a simulation using my (pras) phd finite
        // elasticity code (linear/piecewise constant) basis functions for
        // displacement/pressure, 8*8*8 nodes in the cube, same material law
        // gravity. That code had been compared with cmiss and gave exactly the
        // same answers when cmiss was run with the same mesh,bases etc.

        // note: would seem sensible to compare this directly with cmiss but,
        // for technical reasons cmiss is a complete utter *^%&*$ when it comes
        // to large meshes. my code was compared to cmiss on a small mesh, they
        // agreed exactly. we can't use a small problem here as we are not
        // solving identical problems (different bases) and require some amount
        // of numerical convergence

        // compare the deformation at the unfixed corner nodes to the linear
        // simulation

        // assume the difference between the two simulations is due to the
        // different bases used/numerical convergence not being reached.

        // been visually verified that the two simulations agree quite close
        // qualitatively (note that it is quite a large deformation)

        // verify node 1 is the point corresponding undeformed position (0.1,0,0)
        TS_ASSERT_DELTA(undeformed_position[0](1), 0.1, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[1](1),   0, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[2](1),   0, 1e-12);
        // the linear simulation moves this node to (0.1077, 0.0324, 0.0000)
        TS_ASSERT_DELTA(deformed_position[0](1), 0.1077, 1e-2);
        TS_ASSERT_DELTA(deformed_position[1](1), 0.0324, 1e-2);
        TS_ASSERT_DELTA(deformed_position[2](1), 0.0000, 1e-2);

        // verify node 2 is the point corresponding undeformed position (0.1,0,0.1)
        TS_ASSERT_DELTA(undeformed_position[0](2), 0.1, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[1](2),   0, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[2](2), 0.1, 1e-12);
        // the linear simulation moves this node to (0.1077, 0.0324, 0.1000)
        TS_ASSERT_DELTA(deformed_position[0](2), 0.1077, 1e-2);
        TS_ASSERT_DELTA(deformed_position[1](2), 0.0324, 1e-2);
        TS_ASSERT_DELTA(deformed_position[2](2), 0.1000, 1e-2);

        // verify node 5 is the point corresponding undeformed position (0.1,0.1,0)
        TS_ASSERT_DELTA(undeformed_position[0](5), 0.1, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[1](5), 0.1, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[2](5),   0, 1e-12);
        // the linear simulation moves this node to (0.0887, 0.1310, 0.0002)
        TS_ASSERT_DELTA(deformed_position[0](5), 0.0887, 1e-2);
        TS_ASSERT_DELTA(deformed_position[1](5), 0.1310, 1e-2);
        TS_ASSERT_DELTA(deformed_position[2](5), 0.0002, 1e-2);

        // verify node 6 is the point corresponding undeformed position (0.1,0.1,0.1)
        TS_ASSERT_DELTA(undeformed_position[0](6), 0.1, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[1](6), 0.1, 1e-12);
        TS_ASSERT_DELTA(undeformed_position[2](6), 0.1, 1e-12);
        // the linear simulation moves this node to (0.0887, 0.1310, 0.0998)
        TS_ASSERT_DELTA(deformed_position[0](6), 0.0887, 1e-2);
        TS_ASSERT_DELTA(deformed_position[1](6), 0.1310, 1e-2);
        TS_ASSERT_DELTA(deformed_position[2](6), 0.0998, 1e-2);
    }



    void Test3dProblemOnCubeFixedDisplacement() throw(Exception)
    {
        Vector<double> body_force(3); // zero vector
        body_force(2)=0.05;

        MooneyRivlinMaterialLaw2<3> mooney_rivlin_law(0.02,0.02);

        Triangulation<3> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(1);

        Triangulation<3>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<3>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double z = element_iter->face(face_index)->center()(2);

                    if (fabs(z)<1e-4)
                    {
                        // z=0, label as fixed boundary
                        element_iter->face(face_index)->set_boundary_indicator(FIXED_BOUNDARY);
                    }
                    else if (fabs(z-1)<1e-4)
                    {
                        // z=1, label as dirichlet boundary
                        element_iter->face(face_index)->set_boundary_indicator(DIRICHLET_BOUNDARY);
                    }
                    else
                    {
                        // z!=0 or 1, label as neumann boundary
                        element_iter->face(face_index)->set_boundary_indicator(NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }

        FiniteElasticityAssembler<3> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/fixeddisp3d");

        DoFHandler<3>& dof_handler = finite_elasticity.rGetDofHandler();


        std::map<unsigned,double> boundary_values;

        std::vector<bool> component_mask(3+1); // dim+1
        component_mask[0] = true;
        component_mask[1] = true;
        component_mask[2] = true;
        component_mask[3] = false;

        VectorTools::interpolate_boundary_values(dof_handler,
                                                 FIXED_BOUNDARY,
                                                 ZeroFunction<3>(3+1),  // note the "+1" here! - number of components
                                                 boundary_values,
                                                 component_mask);

        assert(!boundary_values.empty());


        VectorTools::interpolate_boundary_values(dof_handler,
                                                 DIRICHLET_BOUNDARY,
                                                 ComponentSelectFunction<3>(2,-0.1,4),
                                                 boundary_values,
                                                 component_mask);

        finite_elasticity.SetBoundaryValues(boundary_values);


        finite_elasticity.StaticSolve();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = finite_elasticity.rGetDeformedPosition();

        // UPDATE THE NODE POSITIONS:
        // GetVertex returns a reference to a Point<DIM>, so this changes the mesh
        // directly. Do this so the new volume can be computed
        TriangulationVertexIterator<3> vertex_iter(&mesh);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            vertex_iter.GetVertex()[0] = deformed_position[0](vertex_index);
            vertex_iter.GetVertex()[1] = deformed_position[1](vertex_index);
            vertex_iter.GetVertex()[2] = deformed_position[2](vertex_index);

            vertex_iter.Next();
        }

        // compute the deformed volume
        // NOTE: this aren't very accurate volumes, since we have lost the
        // positions of the extra nodes (those used with quadratic basis functions)
        // and the measure() function below must use linear interpolation. Hence
        // the high tolerances
        double deformed_volume = 0.0;
        element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double element_volume = element_iter->measure();
            TS_ASSERT_DELTA(element_volume, 1.0/mesh.n_active_cells(), 1e-2);

            deformed_volume += element_volume;
            element_iter++;
        }

        TS_ASSERT_DELTA(deformed_volume, 1.0, 0.05);
    }


/* NOTES:
   OLD TEST: see test with MIXED SHEARS BELOW

   - heterogeneity with gravity doesn't work
   - without gravity, the solution should be zero disp, discontinuous pw const
      pressure. BUT, bases don't allow this (unlike if we use linear bases for
      displacement and pw constant for pressure)
   - instead solution is nearly zero displacement, with crinkling at interface
   - this a modelling problem? trying to specify pw const material laws
   - a functional analysis thing?
   - converges with zero gravity, so should converge with small non-zero gravity,
      the fact that it doesn't may indicate a bug
*/
    void xTestOnHeterogeneousProblem()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            if (element_iter->center()[0] < 0.5)
            {
                element_iter->set_material_id(5);
            }
            else
            {
                element_iter->set_material_id(6);
            }
            element_iter++;
        }

        std::vector<unsigned> material_ids;
        material_ids.push_back(5);
        material_ids.push_back(6);

        Vector<double> body_force(2);
        body_force(1) = 0.01;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law_stiff(1);
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law_weak(0.2);

        std::vector<AbstractIncompressibleMaterialLaw2<2>*> material_laws;
        material_laws.push_back(&mooney_rivlin_law_stiff);
        material_laws.push_back(&mooney_rivlin_law_weak);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       NULL,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/heterogeneous2d");

        finite_elasticity.SetMaterialLawsForHeterogeneousProblem(material_laws, material_ids);
        finite_elasticity.StaticSolve();

        // get undeformed position
        std::vector<Vector<double> >& undeformed_position
            = finite_elasticity.rGetUndeformedPosition();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = finite_elasticity.rGetDeformedPosition();

        for (unsigned vertex_index=0; vertex_index<deformed_position[0].size(); vertex_index++)
        {
            // todo: TEST THESE!!
            double X = undeformed_position[0](vertex_index);
            double Y = undeformed_position[1](vertex_index);
            double x = deformed_position[0](vertex_index);
            double y = deformed_position[1](vertex_index);

            std::cout << vertex_index << " " << X << " " << Y
                                      << " " << x << " " << y << "\n";
        }
    }


    void TestWithMixedShears() throw(Exception)
    {
        double c1 = SHEARS_C1;
        double c2 = 0.02;
        double alpha = SHEARS_ALPHA;
        double beta = (c1/c2)*alpha;

        Vector<double> body_force(2);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(4);

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);


        MooneyRivlinMaterialLaw2<2> law1(c1);
        MooneyRivlinMaterialLaw2<2> law2(c2);

        std::vector<AbstractIncompressibleMaterialLaw2<2>*> laws;
        laws.push_back(&law1);
        laws.push_back(&law2);

        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            if (element_iter->center()[0] <= 0.5)
            {
                element_iter->set_material_id(5);
            }
            else
            {
                element_iter->set_material_id(6);
            }
            element_iter++;
        }

        std::vector<unsigned> material_ids;
        material_ids.push_back(5);
        material_ids.push_back(6);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       NULL,
                                                       body_force,
                                                       1.0,
                                                       "dealii_finite_elas/mixed_shears");

        finite_elasticity.SetMaterialLawsForHeterogeneousProblem(laws, material_ids);
        finite_elasticity.SetFunctionalTractionBoundaryCondition(TractionForMixedShears);

        // solve
        finite_elasticity.StaticSolve();
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        for(unsigned i=0; i < r_deformed_position[0].size(); i++)
        {
            double X = r_undeformed_position[0](i);
            double Y = r_undeformed_position[1](i);
    
            double exact_y = X<0.5 ? Y - alpha*X : Y - beta*X + (beta-alpha)/2;
             
            double tol = 1e-2;
            if(fabs(X)<1e-6)
            {
                tol = 1e-9;
            }
            
            TS_ASSERT_DELTA( r_deformed_position[0](i), X, tol );
            TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, tol );
        }
        
        // don't check the final pressure
        //  if quad-hex is like quad-tet they won't be anywhere near...
    }

    //
    // Test using the body force and surface tractions corresponding
    // to x = X+0.5*alpha*X^2, y=Y/(1+alpha*X)
    //
    void TestWithFunctionalData() throw(Exception)
    {
        Vector<double> body_force(2);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(4);

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        MooneyRivlinMaterialLaw2<2> law(MATERIAL_PARAM);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &law,
                                                       body_force,
                                                       1.0,
                                                       "dealii_finite_elas/functional");

        finite_elasticity.SetFunctionalTractionBoundaryCondition(MyTraction);
        finite_elasticity.SetFunctionalBodyForce(MyBodyForce);

        // solve
        finite_elasticity.StaticSolve();


        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        for(unsigned i=0; i < r_deformed_position[0].size(); i++)
        {
            double X = r_undeformed_position[0](i);
            double Y = r_undeformed_position[1](i);
    
            double exact_x = X + 0.5*ALPHA*X*X;
            double exact_y = Y/(1+ALPHA*X);
             
            TS_ASSERT_DELTA( r_deformed_position[0](i), exact_x, 1e-3 );
            TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, 1e-3 );
        }

        // check the final pressure
        Vector<double>& full_solution = finite_elasticity.rGetCurrentSolution();
        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        while (!vertex_iter.ReachedEnd())
        {
            double pressure = full_solution(vertex_iter.GetDof(2));
            TS_ASSERT_DELTA(pressure/(2*MATERIAL_PARAM), 1.0, 1.3e-3);
            vertex_iter.Next();
        }
    }

};
#endif /*TESTFINITEELASTICITYASSEMBLERLONG_HPP_*/
