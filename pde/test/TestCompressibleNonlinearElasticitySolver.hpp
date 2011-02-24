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
#include "PetscSetupAndFinalize.hpp"
#include "ToyCompressibleMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "NonlinearElasticityTools.hpp"



/// todos: #1699:
///   stop using linear systems, two matrices in compressible case, matrix memory allocation, avoid function pointers, write tutorial
///   possible material law refactor




// all these are for the MyBodyForce and MySurfaceTraction functions below. See TestAgainstExactNonlinearSolution
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

    S(0,0) = 2*( w1*Q_PARAM*lam + w3/(lam*Q_PARAM) );
    S(0,1) = -2* w1*Y*A_PARAM/(lam*lam);
    S(1,0) = 2* w3*Y*A_PARAM/(Q_PARAM*lam*lam);
    S(1,1) = 2*( w1/lam + w3*lam );

    return S;
}


c_vector<double,2> MyBodyForce(c_vector<double,2>& X)
{
    assert(X(0)>=0 && X(0)<=1 && X(1)>=0 && X(1)<=1);

    double lam = ComputeLambda(X(0));
    double a = A_PARAM;
    double q = Q_PARAM;
    double w3 = ComputeW3(X(0),X(1));
    double dw3dX = ComputeDW3dX(X(0),X(1));
    double dw3dY = ComputeDW3dY(X(0),X(1));

    double dS00dX = 2*(w1*q*a - w3*a/(q*lam*lam) + dw3dX/(q*lam));
    double dS01dX = 2*(2*w1*a*a*X(1)/(lam*lam*lam));
    double dS10dY = 2*(w3*a/(q*lam*lam) + a*X(1)*dw3dY/(q*lam*lam));
    double dS11dY = 2*lam*dw3dY;

    c_vector<double,2> body_force;
    body_force(0) = -dS00dX-dS10dY;
    body_force(1) = -dS01dX-dS11dY;
    return body_force;
}

c_vector<double,2> MyTraction(c_vector<double,2>& X)
{
    c_matrix<double,2,2> S = Compute1stPkStress(X(0), X(1));

    c_vector<double,2> traction = zero_vector<double>(2);

    if (fabs(X(0)-1.0) <= 1e-12) //Right edge
    {
        traction(0) = S(0,0);
        traction(1) = S(0,1);
    }
    else if (fabs(X(1))  <= 1e-12) //Bottom edge
    {
        traction(0) = -S(1,0);
        traction(1) = -S(1,1);
    }
    else if (fabs(X(1) - 1.0) <= 1e-12)//Top edge
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
    // This is purely for coverage of assembling a 3D system (and also uses alternative, heterogeneous
    // constructor, also for coverage)
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

        CompressibleNonlinearElasticitySolver<3> solver(&mesh,
                                                        laws,
                                                        zero_vector<double>(3),
                                                        1.0,
                                                        "",
                                                        fixed_nodes);


        solver.AssembleSystem(true, true);
    }

    // It just tests that nothing happens if zero force and tractions are given
    void TestWithZeroDisplacement() throw(Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        ToyCompressibleMaterialLaw<2> law(1.0, 0.0, -1.0);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

        CompressibleNonlinearElasticitySolver<2> solver(&mesh,
                                                        &law,
                                                        zero_vector<double>(2),
                                                        1.0,
                                                        "",
                                                        fixed_nodes);

        solver.Solve();
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 0u);

        // get deformed position
        std::vector<c_vector<double,2> >& r_deformed_position
            = solver.rGetDeformedPosition();

        for (unsigned i=0; i<r_deformed_position.size(); i++)
        {
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[0], r_deformed_position[i](0), 1e-8);
            TS_ASSERT_DELTA(mesh.GetNode(i)->rGetLocation()[1], r_deformed_position[i](1), 1e-8);
        }
    }

    /**
     *  Suppose the deformation is given to be x = (alpha X, beta Y), and the material law is the toy
     *  law W(I1,I2,I3) = c1(I1-3) - c3(I3-1)     (where c3=-c1)
     *  On the unit square we specify displacement boundaries on the X=0 which match this deformation, we assume
     *  zero body force, traction boundary conditions on the top/bottom surfaces, and fixed traction value, s, on
     *  the X=1 surface.
     *
     *  The given deformation corresponds to constant F=diag(alpha,beta) => constant 1st PK stress
     *  S=diag(2*c1*alpha + 2*c3/alpha, 2*c1*beta + 2*c3/beta), which satisfies dS/dX_M = 0, so the equilibrium
     *  equation is trivially satisfied with zero body force. To match the traction boundary conditions we need
     *  S11 = s   -- defines the traction
     *  S22 = 0   -- this gives a relationship between alpha and beta, although in this case beta = sqrt(-c3/c1)
     *  is indep of alpha
     */
    void TestSolveForSimpleDeformation() throw(Exception)
    {
        double c = 2.0;
        double alpha = 0.9;
        double beta = 1; // = sqrt(-c3/c1);
        double traction_value = 2*c*alpha - 2*c/alpha;

        c_vector<double,2> body_force = zero_vector<double>(2);
        unsigned num_elem = 5;

        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);
        ToyCompressibleMaterialLaw<2> law(c, 0.0, -c);

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

        CompressibleNonlinearElasticitySolver<2> solver(&mesh,
                                                        &law,
                                                        body_force,
                                                        1.0,
                                                        "comp_nonlin_elas_non_zero_bcs",
                                                        fixed_nodes,
                                                        &locations);

        solver.SetSurfaceTractionBoundaryConditions(boundary_elems, tractions);

        // coverage
        solver.SetKspAbsoluteTolerance(1e-10);

        solver.Solve();

        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 3u); // 'hardcoded' answer, protects against jacobian getting messed up

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
     *  Same deformation as TestSolveForSimpleDeformation (see description for this). The deformation is given by
     *  x = (alpha X, beta Y), but with a nonlinear material law
     *  W(I1,I2,I3) = c(I1*I3^{-1/2} -3) - d(I3^{1/2} - 1)^2
     *  On the unit square we specify displacement boundaries on the X=0 which match this deformation, we assume
     *  zero body force, traction boundary conditions on the top/bottom surfaces, and fixed traction value, s, on
     *  the X=1 surface. As in the above test, we have to compute S by hand, in terms of alpha and beta, then S11
     *  defines s, and S22=0 gives a relationship between alpha and beta.
     *
     *  In this case (writing a for alpha, etc), 0.5*S22 = c/a + (1/b)(-0.5c*(a^2+b^2)/(ab)^3 + d(1 - 1/(ab))
     *
     *  which gives a quartic equation to determine beta:
     *   (a^2)b^4 + (Da^3)b^3 - (Da^2+1/2)b^2 - a^2/2
     *  where D=d/c.
     *  For a given alpha we can use matlab to get the solution (choosing the positive real root):
     *
     *  >> a=0.9; D=0.5;
     *  >> roots([a*a, D*a*a*a, -D*a*a-0.5,0,-0.5*a*a])
     *  ans =
     *   -1.415560175004011
     *    1.048769433904755
     *   -0.041604629450372 + 0.578844505760885i
     *   -0.041604629450372 - 0.578844505760885i
     */
    void TestSolveForSimpleDeformationWithCompMooneyRivlin() throw(Exception)
    {
        double c = 2.2;
        double d = 1.1;
        double alpha = 0.9;
        double beta = 1.048769433904755;

        //// another choice
        //double d = 0.0;
        //double alpha = 0.9;
        //double beta = 1.039313636950019;

        double w1 = c/(alpha*beta);
        double w3 = -0.5*c*(alpha*alpha+beta*beta)*pow(alpha*beta,-3) + d*(1.0 - 1.0/(alpha*beta));

        double traction_value = 2*w1*alpha + 2*w3/alpha;

        c_vector<double,2> body_force = zero_vector<double>(2);
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

        CompressibleNonlinearElasticitySolver<2> solver(&mesh,
                                                        &law,
                                                        body_force,
                                                        1.0,
                                                        "comp_nonlin_compMR_simple",
                                                        fixed_nodes,
                                                        &locations);

        solver.SetSurfaceTractionBoundaryConditions(boundary_elems, tractions);

        // coverage
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
     *  Test using a nonlinear material law and for a nonlinear deformation
     *
     *  We take the unit square as before, fix one side, use the compressible Mooney-Rivlin
     *  law,
     *  W(I1,I2,I3) = c(I1*I3^{-1/2} -3) - d(I3^{1/2} - 1)^2,
     *  and want to prescribe tractions on the remaining sides and a body force such that
     *
     *  x = [ q(X + aX^2/2), Y/(1+aX) ]
     *
     *  Note that when q=1, this is an incompressible nonlinear deformation (see similar test of
     *  incompressible solver). q adds so compressibility
     *
     *  The after a page of algebra, we can derive what the 1st PK stress is, which allows us to
     *  determine the required traction and body force.
     *
     */
    void TestAgainstExactNonlinearSolution() throw(Exception)
    {
        unsigned num_elem = 10;
        QuadraticMesh<2> mesh(1.0/num_elem, 1.0, 1.0);

        CompressibleMooneyRivlinMaterialLaw<2> law(C_PARAM,D_PARAM);

        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh,0,0);

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

        CompressibleNonlinearElasticitySolver<2> solver(&mesh,
                                                        &law,
                                                        zero_vector<double>(2),
                                                        1.0,
                                                        "comp_nonlin_elas_exact_soln",
                                                        fixed_nodes);

        solver.SetFunctionalBodyForce(MyBodyForce);
        solver.SetFunctionalTractionBoundaryCondition(boundary_elems, MyTraction);

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
    }
};

#endif /* TESTCOMPRESSIBLENONLINEARELASTICITYSOLVER_HPP_ */
