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


#ifndef TESTFINITEELASTICITYASSEMBLER_HPP_
#define TESTFINITEELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteElasticityAssembler.cpp"

#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

#include "MooneyRivlinMaterialLaw2.hpp"
#include "PolynomialMaterialLaw3d2.hpp"
#include "ExponentialMaterialLaw2.hpp"

#include "FiniteElasticityTools.hpp"

// class for non-zero dirichlet boundary conditions
class MyFunction : public Function<2>
{
private:
    double mLambda;
public:
    MyFunction(double lambda)
     : Function<2>(3),
       mLambda(lambda)
    {
    }
    
    void vector_value(const Point<2>& p, Vector<double> &values) const
    {
        assert(values.size()==3);
        values(0) = 0.0;
        values(1) = (mLambda-1)*p(1);
        values(2) = 0.0;
    }
};
    


class TestFiniteElasticityAssembler : public CxxTest::TestSuite
{
public :
    void TestExceptions() throw(Exception)
    {
        Vector<double> body_force(2);
        Vector<double> bad_body_force(3);
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

        // should throw because the mesh has no surface elements set as the fixed boundary
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssembler<2> bad_fe1(&mesh,&mooney_rivlin_law,body_force,1.0,""));

        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        // should throw because of the negative density
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssembler<2> bad_fe2(&mesh,&mooney_rivlin_law,body_force,-1.0,""));
        // should throw because the body force is the wrong size
        TS_ASSERT_THROWS_ANYTHING(FiniteElasticityAssembler<2> bad_fe3(&mesh,&mooney_rivlin_law,bad_body_force,1.0,""));

        // should be ok now
        TS_ASSERT_THROWS_NOTHING(FiniteElasticityAssembler<2> fe3(&mesh,&mooney_rivlin_law,body_force,1.0,""));

        std::vector<unsigned> material_ids;
        material_ids.push_back(0);

        std::vector<AbstractIncompressibleMaterialLaw2<2>*> material_laws;
        material_laws.push_back(&mooney_rivlin_law);
        material_laws.push_back(&mooney_rivlin_law);

        FiniteElasticityAssembler<2> fe(&mesh,&mooney_rivlin_law,body_force,1.0,"");

        // should thrown because material_laws and material_ids are not the same size
        TS_ASSERT_THROWS_ANYTHING(fe.SetMaterialLawsForHeterogeneousProblem(material_laws,material_ids));

        // check for exception is mesh contains elements whose material id is not
        // equal to either of those passed in in material_ids
        material_ids.clear();
        material_ids.push_back(5);
        material_ids.push_back(6);
        TS_ASSERT_THROWS_ANYTHING(fe.SetMaterialLawsForHeterogeneousProblem(material_laws,material_ids));
    }


    void TestCompareJacobians() throw(Exception)
    {
        Vector<double> body_force(2);
        body_force(0) = 6.0;

        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);


        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");
        TS_ASSERT_THROWS_NOTHING( finite_elasticity.CompareJacobians() );
    }


    // just tests whether the method rGetUndeformedPosition() returns a data structure
    // that consistent with the mesh..
    void TestGetUndeformedPosition() throw(Exception)
    {
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(2.0);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");


        // get undeformed position
        std::vector<Vector<double> >& undeformed_position = finite_elasticity.rGetUndeformedPosition();
        TS_ASSERT_EQUALS(undeformed_position.size(), 2u);
        TS_ASSERT_EQUALS(undeformed_position[0].size(), mesh.n_vertices());
        TS_ASSERT_EQUALS(undeformed_position[1].size(), mesh.n_vertices());

        TriangulationVertexIterator<2> vertex_iter(&mesh);

        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();

            TS_ASSERT_DELTA(undeformed_position[0](vertex_index), posn(0), 1e-12);
            TS_ASSERT_DELTA(undeformed_position[1](vertex_index), posn(1), 1e-12);

            vertex_iter.Next();
        }
    }

    // A test where the solution should be zero displacement
    // It mainly tests that the initial guess was set up correctly to
    // the final correct solution, ie u=0, p=zero_strain_pressure (!=0)
    void TestWithZeroDisplacement() throw(Exception)
    {
        Vector<double> body_force(2); //zero
        double c1 = 3.0;
        MooneyRivlinMaterialLaw2<2> mooney_rivlin_law(c1);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);
        FiniteElasticityTools<2>::SetFixedBoundary(mesh, 0, 0.0);

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &mooney_rivlin_law,
                                                       body_force,
                                                       1.0,
                                                       "");

        finite_elasticity.StaticSolve();

        TS_ASSERT_EQUALS(finite_elasticity.GetNumNewtonIterations(), 0u);

        // get undeformed position
        std::vector<Vector<double> >& undeformed_position
            = finite_elasticity.rGetUndeformedPosition();

        // get deformed position
        std::vector<Vector<double> >& deformed_position
            = finite_elasticity.rGetDeformedPosition();

        for (unsigned i=0; i<deformed_position[0].size(); i++)
        {
            TS_ASSERT_DELTA(undeformed_position[0](i), deformed_position[0](i), 1e-8);
            TS_ASSERT_DELTA(undeformed_position[1](i), deformed_position[1](i), 1e-8);
        }

        // check the final pressure
        Vector<double>& full_solution = finite_elasticity.rGetCurrentSolution();
        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        while (!vertex_iter.ReachedEnd())
        {
            // get the pressure at this node
            double pressure = full_solution(vertex_iter.GetDof(2));
            TS_ASSERT_DELTA(pressure, 2*c1, 1e-6);
            vertex_iter.Next();
        }
    }


    void Test2dProblemOnSquare() throw(Exception)
    {
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


        // Solve again and check that no newton iterations were needed.
        finite_elasticity.StaticSolve();
        TS_ASSERT_EQUALS(finite_elasticity.GetNumNewtonIterations(), 0u);


        // get deformed position
        std::vector<Vector<double> >& deformed_position = finite_elasticity.rGetDeformedPosition();
        TS_ASSERT_EQUALS(deformed_position.size(), 2u);
        TS_ASSERT_EQUALS(deformed_position[0].size(), mesh.n_vertices());
        TS_ASSERT_EQUALS(deformed_position[1].size(), mesh.n_vertices());

        // some hardcoded tests
        TS_ASSERT_DELTA(deformed_position[0](15),0.32141,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](15),0.86898,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](32),0.00000,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](32),0.87500,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](38),0.47755,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](38),0.87779,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](62),0.31559,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](62),0.77431,1e-4);

        TS_ASSERT_DELTA(deformed_position[0](80),0.14713,1e-4);
        TS_ASSERT_DELTA(deformed_position[1](80),0.79705,1e-4);

        // todo: TEST THESE!!

        // also get the solution vector directly to check the deformed position
        // object was set up correctly...
        Vector<double>& solution = finite_elasticity.rGetCurrentSolution();
        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            Point<2> old_posn = vertex_iter.GetVertex();

            Point<2> new_posn;
            new_posn(0) = old_posn(0)+solution(vertex_iter.GetDof(0));
            new_posn(1) = old_posn(1)+solution(vertex_iter.GetDof(1));

            TS_ASSERT_DELTA(deformed_position[0](vertex_index), new_posn(0), 1e-12);
            TS_ASSERT_DELTA(deformed_position[1](vertex_index), new_posn(1), 1e-12);

            //// UPDATE THE NODE POSITIONS
            // GetVertex returns a reference to a Point<DIM>, so this changes the mesh
            // directly. Do this so the new volume can be computed
            vertex_iter.GetVertex()[0] = new_posn(0);
            vertex_iter.GetVertex()[1] = new_posn(1);

            vertex_iter.Next();
        }

        // compute the deformed volume
        // NOTE: this aren't very accurate volumes, since we have lost the
        // positions of the extra nodes (those used with quadratic basis functions)
        // and the measure() function below must use linear interpolation. Hence
        // the high tolerances
        double deformed_volume = 0.0;
        Triangulation<2>::active_cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            double element_volume = element_iter->measure();
            TS_ASSERT_DELTA(element_volume, 1.0/mesh.n_active_cells(), 1e-2);

            deformed_volume += element_volume;
            element_iter++;
        }

        TS_ASSERT_DELTA(deformed_volume, 1.0, 1e-2);
    }
    
    
    /**
     *  Solve a problem with non-zero dirichlet boundary conditions 
     *  and non-zero tractions. THIS TEST COMPARES AGAINST AN EXACT SOLUTION.
     * 
     *  Choosing the deformation x=X/lambda, y=lambda*Y, with a 
     *  Mooney-Rivlin material, then
     *   F = [1/lam 0; 0 lam], T = [2*c1-p*lam^2, 0; 0, 2*c1-p/lam^2], 
     *   sigma = [2*c1/lam^2-p, 0; 0, 2*c1*lam^2-p].
     *  Choosing p=2*c1*lam^2, then sigma = [2*c1/lam^2-p 0; 0 0].
     *  The surface tractions are then 
     *   TOP and BOTTOM SURFACE: 0
     *   RHS: s = SN = J*invF*sigma*N = [lam 0; 0 1/lam]*sigma*[1,0] 
     *          = [2*c1(1/lam-lam^3), 0]
     * 
     *  So, we have to specify displacement boundary conditions (y=lam*Y) on 
     *  the LHS (X=0), and traction bcs (s=the above) on the RHS (X=1), and can
     *  compare the computed displacement and pressure against the true solution.  
     * 
     */
    void TestSolveWithNonZeroBoundaryConditions() throw(Exception)
    {
        double lambda = 0.85;
        double c1 = 0.02;
        Vector<double> body_force(2);
        MooneyRivlinMaterialLaw2<2> law(c1);

        Triangulation<2> mesh;
        GridGenerator::hyper_cube(mesh, 0.0, 1.0);
        mesh.refine_global(3);

  
        ////////////////////////////////////////////////////////
        // define dirichlet (X=0) and Neumann (X=1) boundaries
        ////////////////////////////////////////////////////////
        unsigned component = 0;
        double value = 1.0;
        Triangulation<2>::cell_iterator element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())
        {
            for (unsigned face_index=0; face_index<GeometryInfo<2>::faces_per_cell; face_index++)
            {
                if (element_iter->face(face_index)->at_boundary())
                {
                    double component_val = element_iter->face(face_index)->center()(component);
                    if (fabs(component_val)<1e-4)
                    {
                        // X=0, label as dirichlet boundary
                        element_iter->face(face_index)->set_boundary_indicator(DIRICHLET_BOUNDARY);
                    }
                    else if (fabs(component_val - value)<1e-4)
                    {
                        // X=1, label as neumann boundary
                        element_iter->face(face_index)->set_boundary_indicator(NEUMANN_BOUNDARY);
                    }
                }
            }
            element_iter++;
        }

        FiniteElasticityAssembler<2> finite_elasticity(&mesh,
                                                       &law,
                                                       body_force,
                                                       1.0,
                                                       "finite_elas/dealii_non_zero_bcs");

        // apply traction
        Vector<double> traction(2);
        traction(0) = 2*c1*(pow(lambda,-1) - lambda*lambda*lambda);
        traction(1) = 0;
        finite_elasticity.SetConstantSurfaceTraction(traction);

        //////////////////////////////////////////////////
        // define non-zero dirichlet boundary conditions
        //////////////////////////////////////////////////
        std::map<unsigned,double> boundary_values;

        std::vector<bool> component_mask(2+1); // dim+1
        component_mask[0] = true;
        component_mask[1] = true;
        component_mask[2] = false;

        DoFHandler<2>& dof_handler = finite_elasticity.rGetDofHandler();
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 DIRICHLET_BOUNDARY,
                                                 MyFunction(lambda),
                                                 boundary_values,
                                                 component_mask);

        assert(!boundary_values.empty());
        finite_elasticity.SetBoundaryValues(boundary_values);


        // solve
        finite_elasticity.StaticSolve();
                
                
        // compare                            
        std::vector<Vector<double> >& r_deformed_position = finite_elasticity.rGetDeformedPosition();
        std::vector<Vector<double> >& r_undeformed_position = finite_elasticity.rGetUndeformedPosition();
        
        for(unsigned i=0; i < r_deformed_position[0].size(); i++)
        {
            double X = r_undeformed_position[0](i);
            double Y = r_undeformed_position[1](i);
            double exact_x = (1.0/lambda)*X;
            double exact_y = lambda*Y;
            
            double tol = 1e-6;
            if(fabs(X)<1e-6)
            {
                tol = 1e-9;
            }
            
            TS_ASSERT_DELTA( r_deformed_position[0](i), exact_x, tol );
            TS_ASSERT_DELTA( r_deformed_position[1](i), exact_y, tol );
        }
        
        // check the final pressure
        Vector<double>& full_solution = finite_elasticity.rGetCurrentSolution();
        DofVertexIterator<2> vertex_iter(&mesh, &dof_handler);

        while (!vertex_iter.ReachedEnd())
        {
            // get the pressure at this node
            double pressure = full_solution(vertex_iter.GetDof(2));
            TS_ASSERT_DELTA(pressure, 2*c1*lambda*lambda, 1e-6);
            vertex_iter.Next();
        }
    }
};
#endif /*TESTFINITEELASTICITYASSEMBLER_HPP_*/
