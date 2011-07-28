/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTSTOKESFLOWSOLVER_HPP_
#define TESTSTOKESFLOWSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "StokesFlowSolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"


class TestStokesFlow : public CxxTest::TestSuite
{
public:

	/*
	 * Test that the mesh is calculated correctly on the cannonical triangle.
	 * Tests against the analytical solution calculated by hand.
	 */
	void noTestAssembleOnElementStokes()	throw(Exception)
	{
        EXIT_IF_PARALLEL; // defined in PetscTools

		QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("pde/test/data/canonical_triangle_quadratic",2,2,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<unsigned> dirichlet_nodes;
        dirichlet_nodes.push_back(0);

        double mu = 1.0;
        c_vector<double,2> body_force = zero_vector<double>(2);

        StokesFlowSolver<2> solver(mu,
                                   &mesh,
                                   body_force,
                                   "",
                                   dirichlet_nodes);

        c_matrix<double, 15, 15 > A_elem;
        c_matrix<double, 15, 15 > A_elem_precond;
        c_vector<double, 15 > b_elem;

        solver.AssembleOnElement(*(mesh.GetElement(0)), A_elem, A_elem_precond, b_elem);


		double A[6][6]={
                            {      1.0,  1.0/6.0,  1.0/6.0,      0.0, -2.0/3.0, -2.0/3.0},
                            {  1.0/6.0,  1.0/2.0,      0.0,      0.0,      0.0, -2.0/3.0},
                            {  1.0/6.0,      0.0,  1.0/2.0,      0.0, -2.0/3.0,      0.0},
                            {      0.0,      0.0,      0.0,  8.0/3.0, -4.0/3.0, -4.0/3.0},
                            { -2.0/3.0,      0.0, -2.0/3.0, -4.0/3.0,  8.0/3.0,      0.0},
                            { -2.0/3.0, -2.0/3.0,      0.0, -4.0/3.0,      0.0,  8.0/3.0}
                       };

		double Bx[6][3]={
                            { -1.0/6.0,      0.0,      0.0},
                            {      0.0,  1.0/6.0,      0.0},
                            {      0.0,      0.0,      0.0},
                            {  1.0/6.0,  1.0/6.0,  1.0/3.0},
                            { -1.0/6.0, -1.0/6.0, -1.0/3.0},
                            {  1.0/6.0, -1.0/6.0,      0.0},
                        };

		double By[6][3]={
                            { -1.0/6.0,      0.0,      0.0},
                            {      0.0,      0.0,      0.0},
                            {      0.0,      0.0,  1.0/6.0},
                            {  1.0/6.0,  1.0/3.0,  1.0/6.0},
                            {  1.0/6.0,      0.0, -1.0/6.0},
                            { -1.0/6.0, -1.0/3.0, -1.0/6.0},
                        };


        c_matrix<double,15,15> exact_ael = zero_matrix<double>(15);

        // The diagonal 6x6 blocks
        for (int i=0; i<6; i++)
        {
            for (int j=0; j<6; j++)
            {
                exact_ael(2*i,  2*j)   = A[i][j];
                exact_ael(2*i+1,2*j+1) = A[i][j];
            }
        }

         // The 6x3 Blocks
        for (int i=0; i<6; i++)
        {
            for (int j=0; j<3; j++)
            {
                exact_ael(2*i,12+j)   = -Bx[i][j];
                exact_ael(2*i+1,12+j) = -By[i][j];
                //- as -Div(U)=0
                exact_ael(12+j,2*i)   = -Bx[i][j];
                exact_ael(12+j,2*i+1) = -By[i][j];
            }
        }

		for (int i=0; i<15; i++)
		{
			for (int j=0; j<15; j++)
			{
			    if (fabs(A_elem(i,j))<1e-9)
                {
                    A_elem(i,j) = 0.0;
                }

				TS_ASSERT_DELTA(A_elem(i,j),exact_ael(i,j),1e-9);
			}
			TS_ASSERT_DELTA(b_elem(i),0.0,1e-9);
		}

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Preallocation failure: requested number of nonzeros per row greater than number of columns");
        Warnings::QuietDestroy();
	}

	/*
	 * solution is u = [x, -y], p = const (= 1 as applying zero-Neumann on RHS)
	 * Dirichlet BC applied on three sides, zero-Neumann on the other (so pressure is fully defined)
	 * Just 2 elements
	 */
	void TestStokesWithDirichletVerySimple() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // set up a mesh on [0 1]x[0 1]
        unsigned num_elem = 1;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        // material params
        double mu = 1.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x=mesh.GetNode(i)->rGetLocation()[0];
            double y=mesh.GetNode(i)->rGetLocation()[1];
            // Only apply on Top Bottom and LHS
            if ( x ==0.0  || y ==0.0 || y == 1.0 )
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = x;
                flow(1) = -y;
                dirichlet_flow.push_back(flow);
            }
        }
        assert(dirichlet_flow.size()==7);


        StokesFlowSolver<2> solver(mu,
                                   &mesh,
                                   zero_vector<double>(2),
                                   "SimpleStokesFlow",
                                   dirichlet_nodes,
                                   &dirichlet_flow);

        solver.Solve();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x=mesh.GetNode(i)->rGetLocation()[0];
            double y=mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), x,  1e-5);
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), -y, 1e-5);
        }

        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA(solver.rGetPressures()[i], 1.0, 1e-5 );
        }

        // Test Warnings
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Preallocation failure: requested number of nonzeros per row greater than number of columns");
        Warnings::QuietDestroy();
    }


	/*
	 * Solution is u = [y(1-y), 0], p = 2(1-x)
	 * Dirichlet BC applied on three sides, zero-Neumann on the other (so pressure is fully defined)
	 */
	void TestStokesWithImposedPipeCondition() throw(Exception)
    {
		EXIT_IF_PARALLEL; // defined in PetscTools

		// set up a mesh on [0 1]x[0 1]
		unsigned num_elem = 5;
		QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

		// material params
		double mu = 1.0;

		// Boundary flow
		std::vector<unsigned> dirichlet_nodes;
		std::vector<c_vector<double,2> > dirichlet_flow;
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			double x=mesh.GetNode(i)->rGetLocation()[0];
			double y=mesh.GetNode(i)->rGetLocation()[1];
			if ( x ==0.0 || y == 0.0 || y == 1.0)
			{
				dirichlet_nodes.push_back(i);
				c_vector<double,2> flow = zero_vector<double>(2);

				flow(0) = y*(1-y);
				flow(1) = 0.0;
				dirichlet_flow.push_back(flow);
			}
		}

		assert(dirichlet_flow.size()== 6*num_elem +1);

		c_vector<double,2> body_force = zero_vector<double>(2);

		StokesFlowSolver<2> solver(mu,
								   &mesh,
								   body_force,
								   "PipeStokesFlow",
								   dirichlet_nodes,
								   &dirichlet_flow);

		//Uncomment to make errors smaller
		//solver.SetKspAbsoluteTolerance(1e-12);

		solver.Solve();

		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			double y=mesh.GetNode(i)->rGetLocation()[1];

			double exact_flow_x = y*(1-y);
			double exact_flow_y = 0.0;

			TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), exact_flow_x, 1e-3);
			TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), exact_flow_y, 1e-3);
		}

		for (unsigned i=0; i<mesh.GetNumVertices(); i++)
		{
			double x=mesh.GetNode(i)->rGetLocation()[0];

			double exact_pressure = 2*(1-x);

			TS_ASSERT_DELTA( solver.rGetPressures()[i], exact_pressure, 1e-3);
		}
    }

	/*
	 * Solution is u = [20xy^3, 5x^4-5y^4], p = 60x^2y-20y^3+const
	 * Dirichlet BC applied on all 4 sides so pressure is not fully defined.
	 * This might not be a great test problem discuss with Dave Kay.
	 * Just 2 elements
	 */
	void TestStokesWithAnalyticSolution() throw(Exception)
    {
        EXIT_IF_PARALLEL; // defined in PetscTools

        // set up a mesh on [-1 1]x[-1 1]
        unsigned num_elem = 20;
        QuadraticMesh<2> mesh(2.0/num_elem, 2.0, 2.0);
        mesh.Translate(-1.0,-1.0);

        // material params
        double mu = 1.0;

        // Boundary flow
        std::vector<unsigned> dirichlet_nodes;
        std::vector<c_vector<double,2> > dirichlet_flow;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
        	double x=mesh.GetNode(i)->rGetLocation()[0];
        	double y=mesh.GetNode(i)->rGetLocation()[1];
            if ( x ==-1.0  || x ==1.0 || y == -1.0 || y == 1.0)
            {
                dirichlet_nodes.push_back(i);
                c_vector<double,2> flow = zero_vector<double>(2);

                flow(0) = 20.0*x*y*y*y;
                flow(1) = 5.0*x*x*x*x - 5.0*y*y*y*y;
                dirichlet_flow.push_back(flow);
            }
        }

        assert(dirichlet_flow.size()== 8*num_elem);

        c_vector<double,2> body_force = zero_vector<double>(2);

        StokesFlowSolver<2> solver(mu,
                                   &mesh,
                                   body_force,
                                   "AnalyticalStokesFlow",
                                   dirichlet_nodes,
                                   &dirichlet_flow);

        // Change tolerance for coverage
        solver.SetKspAbsoluteTolerance(1e-8);

        solver.Solve();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x=mesh.GetNode(i)->rGetLocation()[0];
            double y=mesh.GetNode(i)->rGetLocation()[1];

            double exact_flow_x = 20.0*x*y*y*y;
            double exact_flow_y = 5.0*x*x*x*x - 5.0*y*y*y*y;

            TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), exact_flow_x, 1e-2);
            TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), exact_flow_y, 1e-2);
        }

    	double x0=mesh.GetNode(0)->rGetLocation()[0];
		double y0=mesh.GetNode(0)->rGetLocation()[1];
		double exact_pressure_0 = 60.0*x0*x0*y0 -20.0*y0*y0*y0;
		double pressure_diff = solver.rGetPressures()[0] - exact_pressure_0;

        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
        	double x=mesh.GetNode(i)->rGetLocation()[0];
			double y=mesh.GetNode(i)->rGetLocation()[1];

			double exact_pressure = 60.0*x*x*y -20.0*y*y*y;

            TS_ASSERT_DELTA( solver.rGetPressures()[i], exact_pressure + pressure_diff, 1e-0 );//TODO This error is masive!
        }
    }

	/*
	 * simulation with regularised lid driven cavity u=1-x^4 on the top
	 */
	void TestStokesWithLidCavity() throw(Exception)
	{
        EXIT_IF_PARALLEL; // defined in PetscTools

		// set up a mesh on [-1 1]x[-1 1]
		unsigned num_elem = 5;
		QuadraticMesh<2> mesh(2.0/num_elem, 2.0, 2.0);
	    mesh.Translate(-1.0,-1.0);

		// material params
		double mu = 1.0;

		// Boundary flow
		std::vector<unsigned> dirichlet_nodes;
		std::vector<c_vector<double,2> > dirichlet_flow;
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			double x=mesh.GetNode(i)->rGetLocation()[0];
			double y=mesh.GetNode(i)->rGetLocation()[1];
			if ( x == -1.0 || x == 1.0 || y == -1.0)
			{
				dirichlet_nodes.push_back(i);
				c_vector<double,2> flow = zero_vector<double>(2);
				dirichlet_flow.push_back(flow);
			}
			else if (y == 1.0) // this doesnt include corners
			{
				dirichlet_nodes.push_back(i);
				c_vector<double,2> flow = zero_vector<double>(2);

				flow(0) = 1-x*x*x*x;
				flow(1) = 0.0;
				dirichlet_flow.push_back(flow);
			}

		}

		assert(dirichlet_flow.size()== 8*num_elem);

		c_vector<double,2> body_force = zero_vector<double>(2);

		StokesFlowSolver<2> solver(mu,
								   &mesh,
								   body_force,
								   "LidCavityStokesFlow",
								   dirichlet_nodes,
								   &dirichlet_flow);

		//Uncomment to make errors smaller
		//solver.SetKspAbsoluteTolerance(1e-12);

		solver.Solve();

	    //TODO Test something
	}

	/*
	 * Solution is u = [y(1-y), 0], p = 2(2-x)
	 * Dirichlet BC applied on top and bottom, Neumann on the ends p(x=0)=3, p(x=1)=1,
	 */
	void TestPoiseuilleFlow() throw(Exception)
	{
		EXIT_IF_PARALLEL; // defined in PetscTools

		// set up a mesh on [0 1]x[0 1]
		unsigned num_elem = 5;
		QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

		// material params
		double mu = 1.0;

		// Boundary flow
		std::vector<unsigned> dirichlet_nodes;
		std::vector<c_vector<double,2> > dirichlet_flow;
		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			double y=mesh.GetNode(i)->rGetLocation()[1];
			// Fix top and bottom
			if ( y == 0.0 || y == 1.0)
			{
				dirichlet_nodes.push_back(i);
				c_vector<double,2> flow = zero_vector<double>(2);
				dirichlet_flow.push_back(flow);
			}
		}
		assert(dirichlet_flow.size()== 4*num_elem+2);

		// Normal stress boundary

		std::vector<BoundaryElement<1,2>*> boundary_elems;
		std::vector<c_vector<double,2> > normal_stresses;

		for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0]) < 1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);

                c_vector<double,2> normal_stress = zero_vector<double>(2);
                normal_stress[0] = 3;

                normal_stresses.push_back(normal_stress);
            }
            else if(fabs((*iter)->CalculateCentroid()[0] - 1.0) < 1e-4)
            {
            	BoundaryElement<1,2>* p_element = *iter;
				boundary_elems.push_back(p_element);

				c_vector<double,2> normal_stress = zero_vector<double>(2);
				//This is negative because the outward pointing normal at this edge is opposite to the other edge
				normal_stress[0] = -1;

				normal_stresses.push_back(normal_stress);
            }
        }
        assert(boundary_elems.size()==2.0*num_elem);



		c_vector<double,2> body_force = zero_vector<double>(2);

		StokesFlowSolver<2> solver(mu,
								   &mesh,
								   body_force,
								   "PoiseuilleFlow",
								   dirichlet_nodes,
								   &dirichlet_flow);

		solver.SetSurfaceNormalStressBoundaryConditions(boundary_elems, normal_stresses);

		//Uncomment to make errors smaller
		//solver.SetKspAbsoluteTolerance(1e-12);

		solver.Solve();

		for (unsigned i=0; i<mesh.GetNumNodes(); i++)
		{
			double y=mesh.GetNode(i)->rGetLocation()[1];

			double exact_flow_x = y*(1-y);
			double exact_flow_y = 0.0;

			TS_ASSERT_DELTA(solver.rGetVelocities()[i](0), exact_flow_x, 1e-3);
			TS_ASSERT_DELTA(solver.rGetVelocities()[i](1), exact_flow_y, 1e-3);
		}

		for (unsigned i=0; i<mesh.GetNumVertices(); i++)
		{
			double x=mesh.GetNode(i)->rGetLocation()[0];

			double exact_pressure = 2*(1-x) + 1;

			TS_ASSERT_DELTA( solver.rGetPressures()[i], exact_pressure, 1e-3);
		}
	}



};

#endif // TESTSTOKESFLOWSOLVER_HPP_
