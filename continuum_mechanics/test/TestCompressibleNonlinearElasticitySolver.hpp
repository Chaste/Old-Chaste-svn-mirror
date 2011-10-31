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

#ifndef TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_
#define TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ToyCompressibleMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "NonlinearElasticityTools.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleExponentialLaw.hpp"

/*
 * All these are for the MyBodyForce and MySurfaceTraction functions below.
 * See TestAgainstExactNonlinearSolution().
 */
static const double C_PARAM = 1.0;
static const double D_PARAM = 0.1;
static const double A_PARAM = 0.1;
static const double Q_PARAM = 0.9;
static const double m = -0.5; // -1/DIM
static const double w1 = C_PARAM*pow(Q_PARAM,2*m);

double ComputeLambda(double X)
{
    return 1+A_PARAM*X;
}

double ComputeI1(double X, double Y)
{
    double lam = ComputeLambda(X);
    return Q_PARAM*Q_PARAM*lam*lam + A_PARAM*A_PARAM*Y*Y/(lam*lam*lam*lam) + 1.0/(lam*lam);
}

double ComputeW3(double X, double Y)
{
    return m*C_PARAM*ComputeI1(X,Y)*pow(Q_PARAM,2*m-2) + D_PARAM*(1-1.0/Q_PARAM);
}

double ComputeDW3dX(double X, double Y)
{
    double lam = ComputeLambda(X);
    return m*C_PARAM*pow(Q_PARAM,2*m-2)*(2*Q_PARAM*Q_PARAM*lam - 4*A_PARAM*A_PARAM*Y*Y*pow(lam,-5) - 2/(lam*lam*lam))*A_PARAM;
}

double ComputeDW3dY(double X, double Y)
{
    double lam = ComputeLambda(X);
    return m*C_PARAM*pow(Q_PARAM,2*m-2)*A_PARAM*A_PARAM*2*Y/(lam*lam*lam*lam);
}

c_matrix<double,2,2> Compute1stPkStress(double X, double Y)
{
    c_matrix<double,2,2> S;

    double w3 = ComputeW3(X,Y);
    double lam = ComputeLambda(X);

    S(0,0) = 2*( w1*Q_PARAM*lam + w3*Q_PARAM/(lam) );
    S(0,1) = -2* w1*Y*A_PARAM/(lam*lam);
    S(1,0) = 2* w3*Y*A_PARAM*Q_PARAM/(lam*lam);
    S(1,1) = 2*( w1/lam + w3*lam*Q_PARAM*Q_PARAM );

    return S;
}

c_vector<double,2> MyBodyForce(c_vector<double,2>& rX, double t)
{
    assert(rX(0)>=0 && rX(0)<=1 && rX(1)>=0 && rX(1)<=1);

    double lam = ComputeLambda(rX(0));
    double a = A_PARAM;
    double q = Q_PARAM;
    double w3 = ComputeW3(rX(0),rX(1));
    double dw3dX = ComputeDW3dX(rX(0),rX(1));
    double dw3dY = ComputeDW3dY(rX(0),rX(1));

    double dS00dX = 2*(w1*q*a - w3*q*a/(lam*lam) + dw3dX*q/(lam));
    double dS01dX = 2*(2*w1*a*a*rX(1)/(lam*lam*lam));
    double dS10dY = 2*(w3*a*q/(lam*lam) + a*rX(1)*q*dw3dY/(lam*lam));
    double dS11dY = 2*lam*dw3dY*q*q;

    c_vector<double,2> body_force;
    body_force(0) = -dS00dX-dS10dY;
    body_force(1) = -dS01dX-dS11dY;
    return body_force;
}

c_vector<double,2> MyTraction(c_vector<double,2>& rX, double t)
{
    c_matrix<double,2,2> S = Compute1stPkStress(rX(0), rX(1));

    c_vector<double,2> traction = zero_vector<double>(2);

    if (fabs(rX(0)-1.0) <= 1e-12) //Right edge
    {
        traction(0) = S(0,0);
        traction(1) = S(0,1);
    }
    else if (fabs(rX(1))  <= 1e-12) //Bottom edge
    {
        traction(0) = -S(1,0);
        traction(1) = -S(1,1);
    }
    else if (fabs(rX(1) - 1.0) <= 1e-12)//Top edge
    {
        traction(0) = S(1,0);
        traction(1) = S(1,1);
    }
    else
    {
        NEVER_REACHED;
    }
    return traction;
}

class TestCompressibleNonlinearElasticitySolver : public CxxTest::TestSuite
{
public:

    /*
     * This is purely for coverage of assembling a 3D system (and also uses
     * alternative, heterogeneous constructor, also for coverage).
     */
    void TestAssembleSystem3D() throw (Exception)
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);

        mesh.ConstructFromMeshReader(mesh_reader1);

        ToyCompressibleMaterialLaw<3> law(1.0, 0.0, -1.0);
        std::vector<AbstractMaterialLaw<3>*> laws;
        laws.push_back(&law);

        std::vector<unsigned> fixed_nodes;
        fixed_nodes.push_back(0);

        SolidMechanicsProblemDefinition<3> problem_defn(mesh);

        CompressibleNonlinearElasticitySolver<3> solver(mesh,
                                                        problem_defn,
                                                        laws,
                                                        "");


        solver.AssembleSystem(true, true);

        TS_ASSERT(PetscMatTools::CheckSymmetry(solver.mJacobianMatrix));
    }

	// compare computed Jacobian against a numerically computed
	// Jacobian
    void TestAssembleSystem() throw (Exception)
    {
        QuadraticMesh<2> mesh(1.0/2, 1.0, 1.0);
        CompressibleExponentialLaw<2> law;


        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "");
        solver.AssembleSystem(true, true);

        ///////////////////////////////////////////////////////////////////
        // test whether residual vector is currently zero (as
        // current solution should have been initialised to u=0
        ///////////////////////////////////////////////////////////////////
        ReplicatableVector rhs_vec(solver.mResidualVector);
        TS_ASSERT_EQUALS( rhs_vec.GetSize(), 2U*25U );
        for (unsigned i=0; i<rhs_vec.GetSize(); i++)
        {
            TS_ASSERT_DELTA(rhs_vec[i], 0.0, 1e-12);
        }

        ///////////////////////////////////////////////////////////////////
        // compute numerical Jacobian and compare with analytic jacobian
        // (about u=0)
        ///////////////////////////////////////////////////////////////////
        unsigned num_dofs = rhs_vec.GetSize();
        double h = 1e-6;

        int lo, hi;
        MatGetOwnershipRange(solver.mJacobianMatrix, &lo, &hi);

        for (unsigned j=0; j<num_dofs; j++)
        {
            solver.rGetCurrentSolution().clear();
            solver.rGetCurrentSolution().resize(num_dofs, 0.0);
            solver.rGetCurrentSolution()[j] += h;

            solver.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( solver.mResidualVector );

            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = PetscMatTools::GetElement(solver.mJacobianMatrix,i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec[i])/h;
                    if ((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-3);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                    }

//                    double diff =  analytic_matrix_val - numerical_matrix_val;
//                    if(fabs(diff)<1e-6)
//                    {
//                        diff = 0.0;
//                    }
//                    std::cout << diff << " ";
                }
            }
//            std::cout << "\n";
        }
        PetscTools::Barrier();


        //////////////////////////////////////////////////////////
        // compare numerical and analytic jacobians again, this
        // time using a non-zero displacement, u=lambda x, v = mu y
        // (lambda not equal to 1/nu)
        //////////////////////////////////////////////////////////
        double lambda = 1.2;
        double mu = 1.0/1.3;

        solver.rGetCurrentSolution().clear();
        solver.rGetCurrentSolution().resize(num_dofs, 0.0);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            solver.rGetCurrentSolution()[2*i]   = (lambda-1)*mesh.GetNode(i)->rGetLocation()[0];
            solver.rGetCurrentSolution()[2*i+1] = (mu-1)*mesh.GetNode(i)->rGetLocation()[1];
        }

        solver.AssembleSystem(true, true);
        ReplicatableVector rhs_vec2(solver.mResidualVector);

        h=1e-8; // needs to be smaller for this one [COMMENT COPIED FROM INCOMPRESSIBLE VERSION OF THIS TEST]

        for (unsigned j=0; j<num_dofs; j++)
        {
            solver.rGetCurrentSolution()[j] += h;
            solver.AssembleSystem(true, false);

            ReplicatableVector perturbed_rhs( solver.mResidualVector );

            for (unsigned i=0; i<num_dofs; i++)
            {
                if ((lo<=(int)i) && ((int)i<hi))
                {
                    double analytic_matrix_val = PetscMatTools::GetElement(solver.mJacobianMatrix,i,j);
                    double numerical_matrix_val = (perturbed_rhs[i] - rhs_vec2[i])/h;
                    if ((fabs(analytic_matrix_val)>1e-6) && (fabs(numerical_matrix_val)>1e-6))
                    {
                        // relative error
                        TS_ASSERT_DELTA( (analytic_matrix_val-numerical_matrix_val)/analytic_matrix_val, 0.0, 1e-3);
                    }
                    else
                    {
                        // absolute error
                        TS_ASSERT_DELTA(analytic_matrix_val, numerical_matrix_val, 1e-4);
                    }

//                    double diff =  analytic_matrix_val - numerical_matrix_val;
//
//                    if(fabs(diff)<1e-5)
//                    {
//                        diff = 0.0;
//                    }
//                    std::cout << diff << " ";
                }
            }
//            std::cout << "\n";

            solver.rGetCurrentSolution()[j] -= h;
        }


    }


    // It just tests that nothing happens if zero force and tractions are given
    void TestWithZeroDisplacement() throw(Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        ToyCompressibleMaterialLaw<2> law(1.0, 0.0, -1.0);

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "");

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);

        // Get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = solver.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }

        // Coverage of exceptions
        MooneyRivlinMaterialLaw<2> incompressible_law(1.0,1.0);
        TS_ASSERT_THROWS_CONTAINS(CompressibleNonlinearElasticitySolver<2> bad_solver(mesh,problem_defn,&incompressible_law,""),  "ompressibleNonlinearElasticitySolver must take in a compressible material law");

        std::vector<AbstractMaterialLaw<2>*> incompressible_laws;
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            incompressible_laws.push_back(&incompressible_law);
        }
        TS_ASSERT_THROWS_CONTAINS(CompressibleNonlinearElasticitySolver<2> bad_solver(mesh,problem_defn,incompressible_laws,""),  "CompressibleNonlinearElasticitySolver must take in a compressible material law");
    }

    /**
     * Test against an exact solution.
     *
     * Suppose the deformation is given by  x = (alpha X, beta Y), with a nonlinear Mooney-Rivlin material law
     * W(I1,I2,I3) = c(I1*I3^{-1/2} -3) - d(I3^{1/2} - 1)^2
     *
     * On the unit square we specify displacement boundaries on the X=0 which match this deformation, we assume
     * zero body force, zero traction boundary conditions on the top/bottom surfaces, and fixed traction value, s, on
     * the X=1 surface. Using the above deformation and material law we can compute S by hand, in terms of alpha and
     * beta, and then S11 defines s, and S22=0 gives a relationship between alpha and beta.
     *
     * In this case (writing a for alpha, etc), 0.5*S22 = c/a + a^2 b(-0.5c*(a^2+b^2)/(ab)^3 + d(1 - 1/(ab))
     *
     * which gives a cubic equation to determine beta given alpha, c and d. Let D=d/c, then:
     *  (2Da^3)b^3 + (1-2Da^2)b^2  - a^2 = 0
     *
     * For a given alpha we can use matlab to get the solution (choosing the positive real root):
     *
     *  >> a=0.9; D=0.5;
     *  >> roots([2*D*a*a*a, 1-2*D*a*a, 0.0, -a*a])
     *  ans =
     *    -0.608190204001744 + 0.890314286611269i
     *    -0.608190204001744 - 0.890314286611269i
     *     0.955749406631746
     */
    void TestSolveForSimpleDeformationWithCompMooneyRivlin() throw(Exception)
    {
        double c = 2.2;
        double d = 1.1;
        double alpha = 0.9;
        double beta = 0.955749406631746;

        double w1 = c/(alpha*beta); // dW_dI1
        double w3 = -0.5*c*(alpha*alpha+beta*beta)*pow(alpha*beta,-3) + d*(1.0 - 1.0/(alpha*beta)); // dW_dI3

        double traction_value = 2*w1*alpha + 2*w3*alpha*beta*beta;

        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleMooneyRivlinMaterialLaw<2> law(c, d);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if ( fabs(mesh.GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
                c_vector<double,2> new_position;
                new_position(0) = 0;
                new_position(1) = beta*mesh.GetNode(i)->rGetLocation()[1];
                locations.push_back(new_position);
            }
        }

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = traction_value;
        traction(1) = 0;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[0] - 1.0)<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==num_elem);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetFixedNodes(fixed_nodes, locations);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);


        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "comp_nonlin_compMR_simple");

        // Coverage
        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<fixed_nodes.size(); i++)
        {
            unsigned index = fixed_nodes[i];
            TS_ASSERT_DELTA(r_solution[index](0), locations[i](0), 1e-8);
            TS_ASSERT_DELTA(r_solution[index](1), locations[i](1), 1e-8);
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double exact_x = alpha*mesh.GetNode(i)->rGetLocation()[0];
            double exact_y = beta*mesh.GetNode(i)->rGetLocation()[1];

            TS_ASSERT_DELTA( r_solution[i](0), exact_x, 1e-5 );
            TS_ASSERT_DELTA( r_solution[i](1), exact_y, 1e-5 );
        }

        MechanicsEventHandler::Headings();
        MechanicsEventHandler::Report();
    }

    /**
     * Test using a nonlinear material law and for a nonlinear deformation
     *
     * We take the unit square as before, fix one side, use the compressible Mooney-Rivlin
     * law,
     *  W(I1,I2,I3) = c(I1*I3^{-1/2} -3) - d(I3^{1/2} - 1)^2,
     * and want to prescribe tractions on the remaining sides and a body force such that
     *
     *  x = [ q(X + aX^2/2), Y/(1+aX) ]
     *
     * Note that when q=1, this is an incompressible nonlinear deformation (see similar test of
     * incompressible solver). q adds some compressibility
     *
     * Then after a page of algebra, we can derive what the 1st PK stress is, which allows us to
     * determine the required traction and body force.
     *
     * The calculation is written out fully in the FiniteElementImplementations document.
     */
    void TestAgainstExactNonlinearSolution() throw(Exception)
    {
        unsigned num_elem = 10;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        CompressibleMooneyRivlinMaterialLaw<2> law(C_PARAM,D_PARAM);

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        std::vector<BoundaryElement<1,2>*> boundary_elems;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            // get all boundary elems except those on X=0
            if (fabs((*iter)->CalculateCentroid()[0])>1e-6)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
            }
        }
        assert(boundary_elems.size()==3*num_elem);


        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetBodyForce(MyBodyForce);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, MyTraction);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "comp_nonlin_elas_exact_soln");


        solver.Solve();

        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double X = mesh.GetNode(i)->rGetLocation()[0];
            double Y = mesh.GetNode(i)->rGetLocation()[1];

            double exact_x = Q_PARAM*(X + 0.5*A_PARAM*X*X);
            double exact_y = Y/(1+A_PARAM*X);

            TS_ASSERT_DELTA(r_solution[i](0), exact_x, 1e-4);
            TS_ASSERT_DELTA(r_solution[i](1), exact_y, 1e-4);
        }

        // Check that the last matrix was symmetric
        TS_ASSERT(PetscMatTools::CheckSymmetry(solver.mJacobianMatrix));

        ////////////////////////////////////////////////////////////////////
        // Completely separately, we now test the AssembleOnBoundaryElement
        // method for the situation where a normal pressure on the deformed surface
        // is chosen. This is tested in the incompressible case by solving
        // a full problem and testing against an exact solution, here we just
        // check that AssembleOnBoundaryElement in the compressible solver matches
        // AssembleOnBoundaryElement in the incompressible solve
        ////////////////////////////////////////////////////////////////////
        double pressure = 12.32423;
        problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure);

        c_matrix<double,6,6> a_elem;
        c_vector<double,6> b_elem;
        solver.AssembleOnBoundaryElement(*(boundary_elems[0]), a_elem, b_elem, true, false, 0);

        MooneyRivlinMaterialLaw<2> mooney_rivlin_incompressible(1.0);

        IncompressibleNonlinearElasticitySolver<2> incompressible_solver(mesh,
                                                                         problem_defn,
                                                                         &mooney_rivlin_incompressible,
                                                                         "");

        c_matrix<double,8,8> a_elem_incompressible;
        c_vector<double,8> b_elem_incompressible;

        incompressible_solver.mCurrentSolution = solver.mCurrentSolution;

        incompressible_solver.AssembleOnBoundaryElement(*(boundary_elems[0]), a_elem_incompressible, b_elem_incompressible, true, false, 0);

        // incompressible_b = [compressible_b  0 0 ]
        for(unsigned i=0; i<6; i++)
        {
            TS_ASSERT_DELTA( b_elem_incompressible(i), b_elem(i), 1e-12 );
        }
    }


    void TestCheckPositiveDefinitenessOfJacobianMatrix() throw(Exception)
    {
        unsigned num_elem = 10;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleExponentialLaw<2> law;

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "");

        solver.AssembleSystem(false,true);
        unsigned N = solver.mNumDofs;

        Vec test_vec = PetscTools::CreateAndSetVec(N, 0.0);
        Vec product_vec = PetscTools::CreateAndSetVec(N, 0.0);

        for(unsigned i=0; i<N; i++)
        {
            PetscVecTools::SetElement(test_vec, i, 1.0);

            MatMult(solver.mJacobianMatrix,test_vec,product_vec);
            double vT_J_v;
            VecDot(product_vec, test_vec, &vT_J_v);
            //std::cout << vT_J_v << " ";
            TS_ASSERT_LESS_THAN(0.0, vT_J_v);

            PetscVecTools::SetElement(test_vec, i, 0.0);

        }

        VecDestroy(test_vec);
        VecDestroy(product_vec);
    }

    // Solve using an exponential material law. Doesn't test against an exact solution, just that check that the
    // solver converges. Doesn't seem very robust.
    void TestSolveForSimpleDeformationWithExponentialLaw() throw(Exception)
    {
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        CompressibleExponentialLaw<2> law;

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);

        // works ok with g1 = -1,1,2,3. Doesn't newton converge for g1=-2, and gets worse as num_elem increases..
        c_vector<double,2> gravity;
        gravity(0) = 2.0;
        gravity(1) = 0.0;
        problem_defn.SetBodyForce(gravity);

        CompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                        problem_defn,
                                                        &law,
                                                        "CompressibleExponentialLawSolve");


        solver.Solve();

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 5u);

        // check node 5 is (1,0)
        assert( fabs(mesh.GetNode(5)->rGetLocation()[0] - 1.0) < 1e-8 );
        assert( fabs(mesh.GetNode(5)->rGetLocation()[1] - 0.0) < 1e-8 );


        std::vector<c_vector<double,2> >& r_solution = solver.rGetDeformedPosition();

        TS_ASSERT_DELTA(r_solution[5](0), 1.06156, 1e-4);
        TS_ASSERT_DELTA(r_solution[5](1), 0.00510, 1e-4);
    }
};

#endif /* TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_ */
