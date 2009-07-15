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
#ifndef TESTSOLVINGCOUPLEDPDES_HPP_
#define TESTSOLVINGCOUPLEDPDES_HPP_

#include <cxxtest/TestSuite.h>

#include <vector>
#include <cmath>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "AbstractLinearAssembler.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "ReplicatableVector.hpp"

#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>


//////////////////////////////////////////////////////////////////////////////
// a simple pde : u_xx + u_yy + x = 0
//////////////////////////////////////////////////////////////////////////////
class MySimplePde : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& x)
    {
        return x[0];
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>* )
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& x)
    {
        return identity_matrix<double>(2);
    }
};


//////////////////////////////////////////////////////////////////////////////
// an assembler to solve the 'coupled' 2-unknown problem
//    u_xx + u_yy +         x  = 0
//    v_xx + v_yy + \lambda x  = 0
//
//   \lambda is taken in in the constructor
//////////////////////////////////////////////////////////////////////////////
template<class CONCRETE = boost::mpl::void_>
class MySimpleCoupledAssembler : public AbstractLinearAssembler<2,2,2, true, MySimpleCoupledAssembler<CONCRETE> >
{
public:
    static const unsigned E_DIM = 2u;
    static const unsigned S_DIM = 2u;
    static const unsigned P_DIM = 2u;

    typedef MySimpleCoupledAssembler<CONCRETE> SelfType;
    typedef AbstractLinearAssembler<2,2,2, true, SelfType> BaseClassType;
    friend class AbstractStaticAssembler<2, 2, 2, true, SelfType>;

private:
    double mLambda;

    virtual c_matrix<double,2*(2+1),2*(2+1)> ComputeMatrixTerm(c_vector<double, 2+1>& rPhi,
                                                               c_matrix<double, 2, 2+1>& rGradPhi,
                                                               ChastePoint<2>& rX,
                                                               c_vector<double,2>& rU,
                                                               c_matrix<double,2,2>& rGradU,
                                                               Element<2,2>* pElement)
    {
        c_matrix<double,2*(2+1),2*(2+1)> ret = zero_matrix<double>(2*(2+1), 2*(2+1));

        // the following can be done more efficiently using matrix slices and prods
        // and so on (see BidomainDg0Assembler) - efficiency not needed for this
        // test though
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    ret(2*i,  2*j)   += rGradPhi(k,i)*rGradPhi(k,j);
                    ret(2*i+1,2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
                }
            }
        }
        return ret;
    }


    virtual c_vector<double,2*(2+1)> ComputeVectorTerm(c_vector<double, 2+1>& rPhi,
                                                       c_matrix<double, 2, 2+1>& rGradPhi,
                                                       ChastePoint<2>& rX,
                                                       c_vector<double,2>& rU,
                                                       c_matrix<double,2,2>& rGradU,
                                                       Element<2,2>* pElement)
    {
        c_vector<double,2*(2+1)> ret;

        for (unsigned i=0; i<3; i++)
        {
            ret(2*i)   =         rX[0]*rPhi(i);
            ret(2*i+1) = mLambda*rX[0]*rPhi(i);
        }
        return ret;
    }


    virtual c_vector<double, 2*2> ComputeVectorSurfaceTerm(const BoundaryElement<2-1,2>& rSurfaceElement,
                                                           c_vector<double,2>& rPhi,
                                                           ChastePoint<2>& rX )
    {
        // D_times_grad_u_dot_n  = (D gradu) \dot n
        double D_times_grad_u_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 0);
        double D_times_grad_v_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, rX, 1);

        c_vector<double, 2*2> ret;
        for (int i=0; i<2; i++)
        {
            ret(2*i)   = rPhi(i)*D_times_grad_u_dot_n;
            ret(2*i+1) = rPhi(i)*D_times_grad_v_dot_n;
        }

        return ret;
    }


public:
    MySimpleCoupledAssembler(TetrahedralMesh<2,2>* pMesh,
                             BoundaryConditionsContainer<2,2,2>* pBoundaryConditions,
                             double lambda) :
            AbstractAssembler<2,2,2>(),
            BaseClassType()
    {
        this->mpMesh = pMesh;
        this->mpBoundaryConditions = pBoundaryConditions;
        mLambda = lambda;
    }
};

template<class CONCRETE>
struct AssemblerTraits<MySimpleCoupledAssembler<CONCRETE> >
{
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     MySimpleCoupledAssembler<CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CVT_CLS>::type
            CVT_CLS;
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     MySimpleCoupledAssembler<CONCRETE>,
                                     typename AssemblerTraits<CONCRETE>::CMT_CLS>::type
            CMT_CLS;
    typedef typename boost::mpl::if_<boost::mpl::is_void_<CONCRETE>,
                                     AbstractStaticAssembler<2u,2u,2u,true,MySimpleCoupledAssembler<CONCRETE> >,
                                     typename AssemblerTraits<CONCRETE>::CMT_CLS>::type
            INTERPOLATE_CLS;
};


//////////////////////////////////////////////////////////////////////////////
// an assembler to solve the coupled 2-unknown problem
//    u_xx + u_yy + v = f(x,y)
//    v_xx + v_yy + u = g(x,y)
//
//   where f and g are chosen so that (with zero-dirichlet boundary conditions)
//   the solution is
//       u = sin(pi*x)sin(pi*x),   v = sin(2*pi*x)sin(2*pi*x)
//
// ComputeMatrixTerm() and ComputeSurfaceRhsTerm() are identical to MySimpleCoupledAssembler
//  (above), so this class inherits from MySimpleCoupledAssembler
//////////////////////////////////////////////////////////////////////////////
class AnotherCoupledAssembler : public MySimpleCoupledAssembler<AnotherCoupledAssembler>
{
public:
    static const unsigned E_DIM = 2u;
    static const unsigned S_DIM = 2u;
    static const unsigned P_DIM = 2u;

    friend class AbstractStaticAssembler<2, 2, 2, true, MySimpleCoupledAssembler<AnotherCoupledAssembler> >;

private:
    double f(double x,double y)
    {
        return -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y) + sin(2*M_PI*x)*sin(2*M_PI*y);
    }

    double g(double x,double y)
    {
        return -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y) + sin(M_PI*x)*sin(M_PI*y);
    }

    virtual c_vector<double,2*(2+1)> ComputeVectorTerm(c_vector<double, 2+1>& rPhi,
                                                       c_matrix<double, 2, 2+1>& rGradPhi,
                                                       ChastePoint<2>& rX,
                                                       c_vector<double,2>& rU,
                                                       c_matrix<double,2,2>& rGradU,
                                                       Element<2,2>* pElement)
    {
        c_vector<double,2*(2+1)> ret;

        for (unsigned i=0; i<3; i++)
        {
            ret(2*i)   = ( rU(1) - f(rX[0],rX[1]) )*rPhi(i);   // = (v-f(x,y))*phi_i
            ret(2*i+1) = ( rU(0) - g(rX[0],rX[1]) )*rPhi(i);   // = (u-g(x,y))*phi_i
        }
        return ret;
    }


public :
    AnotherCoupledAssembler(TetrahedralMesh<2,2>* pMesh,
                            BoundaryConditionsContainer<2,2,2>* pBoundaryConditions) :
            AbstractAssembler<2,2,2>(),
            MySimpleCoupledAssembler<AnotherCoupledAssembler>(pMesh, pBoundaryConditions,0.0)
    {}
};

template<>
struct AssemblerTraits<AnotherCoupledAssembler>
{
    typedef AnotherCoupledAssembler CVT_CLS;
    typedef MySimpleCoupledAssembler<AnotherCoupledAssembler> CMT_CLS;
    typedef AbstractStaticAssembler<2u,2u,2u,true,AnotherCoupledAssembler> INTERPOLATE_CLS;
};


//////////////////////////////////////////////////////////////////////////////
// test class
//////////////////////////////////////////////////////////////////////////////
class TestSolvingCoupledPdes : public CxxTest::TestSuite
{
public:

    /*  Solve:
     *     u_xx + u_yy + x  = 0
     *     v_xx + v_yy + 2x = 0
     *  with zero dirichlet on boundary
     *
     *  This is obviously really just two virtually identical uncoupled
     *  problems
     */
    void TestSimpleCoupledPde() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ////////////////////////////////////////////////////////////////
        // Solve coupled system using assembler defined above
        ////////////////////////////////////////////////////////////////

        // boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns;
        bcc_2unknowns.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc_2unknowns.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // lambda in MySimpleCoupledAssembler = 2
        MySimpleCoupledAssembler<> assembler_2unknowns(&mesh,&bcc_2unknowns,2.0);
        Vec result_2unknowns = assembler_2unknowns.Solve();
        ReplicatableVector result_2unknowns_repl(result_2unknowns);


        ///////////////////////////////////////////////////////////////////
        // Now solve u_xx + u_yy + x = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////

        // Instantiate PDE object
        MySimplePde pde;  //defined above

        // boundary conditions for 1-unknown problem
        BoundaryConditionsContainer<2,2,1> bcc_1unknown;
        bcc_1unknown.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Assembler
        SimpleLinearEllipticAssembler<2,2> assembler_1unknown(&mesh,&pde,&bcc_1unknown);

        Vec result_1unknown = assembler_1unknown.Solve();
        ReplicatableVector result_1unknown_repl(result_1unknown);

        // check the u solutions (result_2unknowns_repl[2*i]) is equal to the
        // solution of the 1-unknown problem and the v solutions
        // (result_2unknowns_repl[2*i+1]) are equal to two times the 1-unknown
        // solution
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i]  ,   result_1unknown_repl[i], 1e-10);
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i+1], 2*result_1unknown_repl[i], 1e-10);
            //std::cout << result_1unknown_repl[i] << " ";
        }

        VecDestroy(result_2unknowns);
        VecDestroy(result_1unknown);
    }


    /*  Solve:
     *     u_xx + u_yy + x = 0
     *     v_xx + v_yy + x = 0
     *  with neumann boundary conditions (the same on both u and v)
     *  on part of the boundary
     *
     *  This is obviously two identical uncoupled problems
     */
    void TestSimpleCoupledPdeWithNeumannBoundaryConditions() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ////////////////////////////////////////////////////////////////
        // Solve coupled system using assembler defined above
        ////////////////////////////////////////////////////////////////

        // boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns;

        // du/dn = -0.5 on r=1
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2> *p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        ConstBoundaryCondition<2> *p_boundary_condition1 = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition,0);
            bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition1,1);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        p_boundary_condition1 = new ConstBoundaryCondition<2>(2.0);
        bcc_2unknowns.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition,0);
        bcc_2unknowns.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition1,1);

        // use assembler to solve (with lambda in MySimpleCoupledAssembler = 1)

        MySimpleCoupledAssembler<> assembler_2unknowns(&mesh,&bcc_2unknowns,1.0);

        Vec result_2unknowns = assembler_2unknowns.Solve();
        ReplicatableVector result_2unknowns_repl(result_2unknowns);


        ///////////////////////////////////////////////////////////////////
        // Now solve u_xx + u_yy + x = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////

        // Instantiate PDE object
        MySimplePde pde;  //defined above

        // boundary conditions for 1-unknown problem
        BoundaryConditionsContainer<2,2,1> bcc_1unknown;

        iter = mesh.GetBoundaryElementIteratorBegin();
        p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc_1unknown.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        bcc_1unknown.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition);

        // Assembler
        SimpleLinearEllipticAssembler<2,2> assembler_1unknown(&mesh,&pde,&bcc_1unknown);

        Vec result_1unknown = assembler_1unknown.Solve();
        ReplicatableVector result_1unknown_repl(result_1unknown);

        // check the u solutions (result_2unknowns_repl[2*i]) is equal to the
        // solution of the 1-unknown problem and the v solutions
        // (result_2unknowns_repl[2*i+1]) are equal to the 1-unknown
        // solution
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i]  , result_1unknown_repl[i], 1e-6);
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i+1], result_1unknown_repl[i], 1e-6);
            //std::cout << result_1unknown_repl[i] << " ";
        }

        VecDestroy(result_2unknowns);
        VecDestroy(result_1unknown);
    }


    /*
     *  Solve a real coupled problem:
     *     u_xx + u_yy  + v = f(x,y)
     *     v_xx + v_yy  + u = g(x,y)
     *
     *  where f and g are chosen so that (with zero-dirichlet boundary conditions)
     *  the solution is
     *
     *     u = sin(pi*x)sin(pi*x),   v = sin(2*pi*x)sin(2*pi*x)
     */
    void TestRealCoupledPde() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // boundary conditions for 2-unknown problem
        BoundaryConditionsContainer<2,2,2> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // purpose-made assembler for this problem:
        AnotherCoupledAssembler assembler(&mesh,&bcc);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];

            double u = sin(M_PI*x)*sin(M_PI*y);
            double v = sin(2*M_PI*x)*sin(2*M_PI*y);

            // need lower tolerance for v because v is higher frequency
            // and not captured very well on this mesh
            TS_ASSERT_DELTA( result_repl[2*i]  , u, 0.02);
            TS_ASSERT_DELTA( result_repl[2*i+1], v, 0.1);
        }

        VecDestroy(result);
    }
};
#endif /*TESTSOLVINGCOUPLEDPDES_HPP_*/
