#ifndef TESTSOLVINGCOUPLEDPDES_HPP_
#define TESTSOLVINGCOUPLEDPDES_HPP_


#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "PetscSetupAndFinalize.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "TrianglesMeshReader.cpp"
#include "AbstractLinearStaticProblemAssembler.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "ReplicatableVector.hpp"


#define PI 3.14159265


//////////////////////////////////////////////////////////////////////////////
// a simple pde : u_xx + u_yy + x = 0
//////////////////////////////////////////////////////////////////////////////
class MySimplePde : public AbstractLinearEllipticPde<2>
{
public:
    double ComputeLinearSourceTerm(Point<2> x) 
    {
        return x[0];
    }
    
    double ComputeNonlinearSourceTerm(Point<2>, double )
    {
        return 0.0;
    }
    
    c_matrix<double,2,2> ComputeDiffusionTerm(Point<2> x)
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
class MySimpleCoupledAssembler : public AbstractLinearStaticProblemAssembler<2,2,2>
{
    double mLambda;
    
    virtual c_matrix<double,2*(2+1),2*(2+1)> ComputeLhsTerm(const c_vector<double, 2+1> &rPhi,
                                                            const c_matrix<double, 2, 2+1> &rGradPhi,
                                                            const Point<2> &rX,
                                                            const c_vector<double,2> &u)
    {
        c_matrix<double,2*(2+1),2*(2+1)> ret;
        
        // the following can be done more efficiently using matrix slices and prods
        // and so on (see BidomainDg0Assembler) - efficiency not needed for this 
        // test though
        for(unsigned i=0; i<3; i++)
        {
            for(unsigned j=0; j<3; j++)
            {
                for(unsigned k=0; k<2; k++)
                {
                    ret(2*i,  2*j)   += rGradPhi(k,i)*rGradPhi(k,j);
                    ret(2*i+1,2*j+1) += rGradPhi(k,i)*rGradPhi(k,j);
                }
            }
        }
        return ret;
    }    

    
    virtual c_vector<double,2*(2+1)> ComputeRhsTerm(const c_vector<double, 2+1> &rPhi,
                                                    const Point<2> &rX,
                                                    const c_vector<double,2> &u)
    {
        c_vector<double,2*(2+1)> ret;
        
        for(unsigned i=0; i<3; i++)
        {
            ret(2*i)   =         rX[0]*rPhi(i);
            ret(2*i+1) = mLambda*rX[0]*rPhi(i);
        }
        return ret;
    }
  

    virtual c_vector<double, 2*2> ComputeSurfaceRhsTerm(const BoundaryElement<2-1,2> &rSurfaceElement,
                                                        const c_vector<double,2> &phi,
                                                        const Point<2> &x )
    {
        // D_times_grad_u_dot_n  = (D gradu) \dot n
        double D_times_grad_u_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, x, 0);
        double D_times_grad_v_dot_n = this->mpBoundaryConditions->GetNeumannBCValue(&rSurfaceElement, x, 1);

        c_vector<double, 2*2> ret;
        for(int i=0; i<2; i++)
        {
            ret(2*i)   = phi(i)*D_times_grad_u_dot_n;
            ret(2*i+1) = phi(i)*D_times_grad_v_dot_n;
        }
    
        return ret;
    }
  
  
public:
    MySimpleCoupledAssembler(ConformingTetrahedralMesh<2,2>* pMesh,
                             BoundaryConditionsContainer<2,2,2>* pBoundaryConditions,
                             double lambda) :
            AbstractLinearStaticProblemAssembler<2,2,2>()            
    {
        mpMesh = pMesh;
        mpBoundaryConditions = pBoundaryConditions;
        mLambda = lambda;
    }    
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
// ComputeLhsTerm() and ComputeSurfaceRhsTerm() are identical to MySimpleCoupledAssembler
//  (above), so this class inherits from MySimpleCoupledAssembler
//////////////////////////////////////////////////////////////////////////////
class AnotherCoupledAssembler : public MySimpleCoupledAssembler
{
    double f(double x,double y)
    {
        return -2*PI*PI*sin(PI*x)*sin(PI*y) + sin(2*PI*x)*sin(2*PI*y);
    }
        
    double g(double x,double y)
    {
        return -8*PI*PI*sin(2*PI*x)*sin(2*PI*y) + sin(PI*x)*sin(PI*y);
    }
        
    virtual c_vector<double,2*(2+1)> ComputeRhsTerm(const c_vector<double, 2+1> &rPhi,
                                                    const Point<2> &rX,
                                                    const c_vector<double,2> &u)
    {
        c_vector<double,2*(2+1)> ret;
     
        for(unsigned i=0; i<3; i++)
        {
            ret(2*i)   = ( u(1) - f(rX[0],rX[1]) )*rPhi(i);   // = (v-f(x,y))*phi_i
            ret(2*i+1) = ( u(0) - g(rX[0],rX[1]) )*rPhi(i);   // = (u-g(x,y))*phi_i
        }
        return ret;
    }


public :
    AnotherCoupledAssembler(ConformingTetrahedralMesh<2,2>* pMesh,
                            BoundaryConditionsContainer<2,2,2>* pBoundaryConditions) :
            MySimpleCoupledAssembler(pMesh, pBoundaryConditions,0.0)            
    {
    }    
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
    void TestSimpleCoupledPde( void ) throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ////////////////////////////////////////////////////////////////
        // Solve coupled system using assembler defined above
        ////////////////////////////////////////////////////////////////        

        // boundary conditions for 2-unknown problem       
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns(mesh.GetNumNodes());
        bcc_2unknowns.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc_2unknowns.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // lambda in MySimpleCoupledAssembler = 2 
        MySimpleCoupledAssembler assembler_2unknowns(&mesh,&bcc_2unknowns,2.0);
        Vec result_2unknowns = assembler_2unknowns.Solve();
        ReplicatableVector result_2unknowns_repl(result_2unknowns);
        
        
        ///////////////////////////////////////////////////////////////////        
        // Now solve u_xx + u_yy + x = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////        
        
        // Instantiate PDE object
        MySimplePde pde;  //defined above

        // boundary conditions for 1-unknown problem       
        BoundaryConditionsContainer<2,2,1> bcc_1unknown(mesh.GetNumNodes());
        bcc_1unknown.DefineZeroDirichletOnMeshBoundary(&mesh);
                
        // Assembler
        SimpleLinearEllipticAssembler<2,2> assembler_1unknown(&mesh,&pde,&bcc_1unknown);
                
        Vec result_1unknown = assembler_1unknown.Solve();
        ReplicatableVector result_1unknown_repl(result_1unknown);
        
        // check the u solutions (result_2unknowns_repl[2*i]) is equal to the 
        // solution of the 1-unknown problem and the v solutions
        // (result_2unknowns_repl[2*i+1]) are equal to two times the 1-unknown 
        // solution
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i]  ,   result_1unknown_repl[i], 1e-10);
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i+1], 2*result_1unknown_repl[i], 1e-10);
            //std::cout << result_1unknown_repl[i] << " ";
        }
    }


    /*  Solve:
     *     u_xx + u_yy + x = 0
     *     v_xx + v_yy + x = 0
     *  with neumann boundary conditions (the same on both u and v)
     *  on part of the boundary
     * 
     *  This is obviously two identical uncoupled problems
     */ 
    void TestSimpleCoupledPdeWithNeumannBoundaryConditions( void ) throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        ////////////////////////////////////////////////////////////////
        // Solve coupled system using assembler defined above
        ////////////////////////////////////////////////////////////////        

        // boundary conditions for 2-unknown problem       
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns(mesh.GetNumNodes());

        // du/dn = -0.5 on r=1
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        ConstBoundaryCondition<2>* p_boundary_condition1 = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition,0);
            bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition1,1);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        p_boundary_condition1 = new ConstBoundaryCondition<2>(2.0);
        bcc_2unknowns.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), p_boundary_condition,0);
        bcc_2unknowns.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), p_boundary_condition1,1);

        // use assembler to solve (with lambda in MySimpleCoupledAssembler = 1) 
        
        MySimpleCoupledAssembler assembler_2unknowns(&mesh,&bcc_2unknowns,1.0);

        Vec result_2unknowns = assembler_2unknowns.Solve();
        ReplicatableVector result_2unknowns_repl(result_2unknowns);
        
        
        ///////////////////////////////////////////////////////////////////        
        // Now solve u_xx + u_yy + x = 0 as an uncoupled 1-unknown problem
        ///////////////////////////////////////////////////////////////////        
        
        // Instantiate PDE object
        MySimplePde pde;  //defined above

        // boundary conditions for 1-unknown problem       
        BoundaryConditionsContainer<2,2,1> bcc_1unknown(mesh.GetNumNodes());
        
        iter = mesh.GetBoundaryElementIteratorBegin();
        p_boundary_condition = new ConstBoundaryCondition<2>(-0.5);
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            bcc_1unknown.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
            iter++;
        }
        // u = 2 at some point on the boundary, say node 1
        p_boundary_condition = new ConstBoundaryCondition<2>(2.0);
        bcc_1unknown.AddDirichletBoundaryCondition(mesh.GetNodeAt(1), p_boundary_condition);
                
        // Assembler
        SimpleLinearEllipticAssembler<2,2> assembler_1unknown(&mesh,&pde,&bcc_1unknown);
                
        Vec result_1unknown = assembler_1unknown.Solve();
        ReplicatableVector result_1unknown_repl(result_1unknown);
        
        // check the u solutions (result_2unknowns_repl[2*i]) is equal to the 
        // solution of the 1-unknown problem and the v solutions
        // (result_2unknowns_repl[2*i+1]) are equal to the 1-unknown 
        // solution
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i]  , result_1unknown_repl[i], 1e-6);
            TS_ASSERT_DELTA(result_2unknowns_repl[2*i+1], result_1unknown_repl[i], 1e-6);
            //std::cout << result_1unknown_repl[i] << " ";
        }
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
    void TestRealCoupledPde( void ) throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4096_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
      
        // boundary conditions for 2-unknown problem       
        BoundaryConditionsContainer<2,2,2> bcc(mesh.GetNumNodes());
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,0); // zero dirichlet for u
        bcc.DefineZeroDirichletOnMeshBoundary(&mesh,1); // zero dirichlet for v

        // purpose-made assembler for this problem:
        AnotherCoupledAssembler assembler(&mesh,&bcc);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNodeAt(i)->GetPoint()[0];
            double y = mesh.GetNodeAt(i)->GetPoint()[1];
            
            double u = sin(PI*x)*sin(PI*y);
            double v = sin(2*PI*x)*sin(2*PI*y);

            // need lower tolerance for v because v is higher frequency 
            // and not captured very well on this mesh
            TS_ASSERT_DELTA( result_repl[2*i]  , u, 0.02);
            TS_ASSERT_DELTA( result_repl[2*i+1], v, 0.1); 
        }
    }
};
#endif /*TESTSOLVINGCOUPLEDPDES_HPP_*/
