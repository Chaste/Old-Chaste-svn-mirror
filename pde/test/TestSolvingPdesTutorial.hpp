/*
 * 
 *  Chaste Pde Tutorial - this page get automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 * 
 * 
 */  
#ifndef TESTSOLVINGPDESTUTORIAL_HPP_
#define TESTSOLVINGPDESTUTORIAL_HPP_

/* 
 * In this test we will solve the PDE: div(D grad u) + u = x^2^+y^2^, in 2D, where
 * D is the diffusion tensor (1 1; 0 1) (ie D11=D12=D22=1, D21=0), on a square
 * domain, with boundary conditions u=0 on x=0 or y=0, (D gradu).n = 0 on x=1 and y=1,
 * where n is the surface normal.
 * 
 * EMPTYLINE
 * 
 * First we include the header needed to define this class as a test suite */
#include <cxxtest/TestSuite.h>
/* This is the class that is needed to solve a linear elliptic pde */
#include "SimpleLinearEllipticAssembler.hpp"
/* This is needed to read mesh datafiles of the 'Triangles' format */
#include "TrianglesMeshReader.cpp"
/* !PetscSetupAndFinalize.hpp must be included in every test that uses Petsc. It
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"



/* We need to create a class representing the PDE we want to solve, which will be 
 * passed into the solver. The PDE we are solving is of the type 
 * {{{AbstractLinearEllipticPde}}}, which is an abstract class with 3 pure methods
 * which have to implemented. The template variable in the dimension of the space.
 */
class MyPde : public AbstractLinearEllipticPde<2>
{
private:
    /* For efficiency, we will save the diffusion that will be returned by one of the
     * class' methods as a member variable */
    c_matrix<double,2,2> mDiffusionTensor;

public:
    /* The constructor just sets up the diffusion tensor. */
    MyPde()
    {
        mDiffusionTensor(0,0) = 1.0;
        mDiffusionTensor(0,1) = 1.0;
        mDiffusionTensor(1,0) = 0.0;
        mDiffusionTensor(1,1) = 1.0;
    }

    /* The first method which has to be implemented returns the constant 
     * (not dependent on u) part of the source term, which for our PDE is 
     * x^2^ + y^2^ */
    double ComputeConstantInUSourceTerm(const ChastePoint<2>& rX)
    {
        return rX[0]*rX[0] + rX[1]*rX[1];
    }

    /* The second  which has to be implemented returns the coefficient in the linear-in-u 
     * part of the source term, which for our PDE is just 1.0 */
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>& )
    {
        return 1.0;
    }    

    /* The third method returns the diffusion tensor D */
    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return mDiffusionTensor;
    }
};


/* Next, we define the test suite (a class). It is sensible to name it the same
 * as the filename. The class should inherit from {{{CxxTest::TestSuite}}} */
class TestSolvingPdesTutorial : public CxxTest::TestSuite
{
/* All individual test defined in this test suite '''must''' be declared as public */ 
public:
    /* Define a particular test */
    void TestExampleWithLinearEllipticAssembler()
    {
        /* First we declare a mesh reader which reads mesh data files of the 'Triangles'
         * format. The path given is the relative to the main Chaste directory. The reader
         * will look for three datafiles, [name].nodes, [name].ele and (in 2d or 3d) 
         * [name].edge. Note that the first template argument here is the dimension of the 
         * elements in the mesh ({{{ELEM_DIM}}}), the second the dimension of the nodes, 
         * i.e. the dimension of the space the mesh lives in ({{{SPACE_DIM}}}). Usually 
         * {{{ELEM_DIM}}} and {{{SPACE_DIM}}} will be equal. */ 
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        /* Now declare a tetrahedral mesh with the same dimensions */
        ConformingTetrahedralMesh<2,2> mesh;
        /* Construct the mesh using the mesh reader */
        mesh.ConstructFromMeshReader(mesh_reader);
        
        /* Next we instantiate an instance of our PDE we which to solve */ 
        MyPde pde;

        /* A set of boundary conditions are stored in a {{{BoundaryConditionsContainer}}}. The 
         * three template arguments are ELEM_DIM, SPACE_DIM and PROBLEM_DIM, the latter being
         * the number of unknowns we are solving for. In this case it is 1. */
        BoundaryConditionsContainer<2,2,1> bcc;
        
        /* Defining the boundary conditions is the only particularly fiddly part of solving PDEs,
         * unless they are very simple, such as u=0 on the boundary, which could be done
         * as follows */
        //bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        /* We want to specify u=0 on x=0 and y=0. To do this, get a boundary node iterator
         * from the mesh */
        ConformingTetrahedralMesh<2,2>::BoundaryNodeIterator iter 
           = mesh.GetBoundaryNodeIteratorBegin();
        /* Then loop over the boundary nodes, getting the x and y value */                
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            /* if x=0 or y=0.. */
            if((x==0) || (y==0))
            {
                /* ..create a new {{{ConstBoundaryConditions}}} object. This is a subclass of
                 * {{{AbstractBoundaryCondition}}}, and tells the caller what value to return
                 * given a particular point in space. We say that value should be 0.0, then
                 * associate it as a boundary condition to be associated with this node (ie. {{{*iter}}}
                 */
                ConstBoundaryCondition<2>* p_dirichlet_boundary_condition = new ConstBoundaryCondition<2>(0.0);
                bcc.AddDirichletBoundaryCondition(*iter, p_dirichlet_boundary_condition);
            }
            iter++;
        }
        
        /* Now we create Neumann boundary conditions for the ''surface elements'' on x=1 and y=1. Note that
         * Dirichlet boundary conditions are defined on nodes, whereas Neumann boundary conditions are 
         * defined on surface elements. Note also that the natural boundary condition statement for this
         * PDE is (D grad u).n = g(x) (where n is the surface normal), and g a prescribed function, 
         * ''not'' something like du/dn=g(x). Hence the boundary condition we are specifying is
         * (D grad u).n = 0.
         * 
         * EMPTYLINE
         * 
         * '''Important note for 1D:''' This means that if we were solving 2u_xx_=f(x) in 1D, and
         * wanted to specify du/dx=1 on the RHS boundary, the boundary value we have to specify is
         * -2, as (D gradu).n = -2 when du/dx=1
         * 
         * EMPTYLINE
         * 
         * To define Neumann bcs, we have to loop over surface elements, using the iterator
         * provided by the mesh class. 
         */ 
         
         ConformingTetrahedralMesh<2,2>::BoundaryElementIterator surf_iter 
           = mesh.GetBoundaryElementIteratorBegin();
         ConstBoundaryCondition<2>* p_neumann_boundary_condition = new ConstBoundaryCondition<2>(0.0);
         while (surf_iter < mesh.GetBoundaryElementIteratorEnd())
         {  
            /* Get the x and y values of any node (here, the 0th) in the element */
            unsigned node_index = (*surf_iter)->GetNodeGlobalIndex(0);
            double x = mesh.GetNode(node_index)->GetPoint()[0];
            double y = mesh.GetNode(node_index)->GetPoint()[1];
            
            /* if x=1 or y=1.. */
            if( (fabs(x-1.0) < 1e-6) || (fabs(y-1.0) < 1e-6) )
            {
                /* associate the boundary condition with the surface element */
                bcc.AddNeumannBoundaryCondition(*surf_iter, p_neumann_boundary_condition);
            }
            
            /* and increment the iterator */
            surf_iter++;
        }

        
        
        /* Next we define the assembler - the solver of the PDE. (Assembler is a bit of a 
         * misnomer - assemblers both assemble the finite element equations, and solve them.
         * To solve {{{AbstractLinearEllipticPde}}} (which is the type of pde {{{MyPde}}} is),
         * we use a {{{SimpleLinearEllipticAssembler}}}. The assembler, again templated over
         * {{{ELEM_DIM}}} and {{{SPACE_DIM}}} needs to be given (pointers to) the mesh, 
         * pde and boundary conditions.
         */ 
        SimpleLinearEllipticAssembler<2,2> assembler(&mesh,&pde,&bcc);
        
        /* To solve, just call {{{Solve()}}}. A Petsc vector is returned. */
        Vec result = assembler.Solve();
        
        /* It is a pain to access the individual components of a Petsc vector, even in
         * sequential. A helper class called {{{ReplicatableVector}}} has been created. Create
         * an instance of one of these, using the petsc {{{Vec}}} as the data. The ith
         * component of {{{result}}} can now be obtained by simply doing {{{result_repl[i]}}}.
         */
        ReplicatableVector result_repl(result);
        
        /* Let us write out the solution to a file. To do this, create an 
         * {{{OutputFileHandler}}}, passing in the directory we want files written to.
         * This is relative to the directory defined by the CHASTE_TEST_OUTPUT environment
         * variable - usually /tmp/chaste/testout. Note by default the output directory
         * passed in is cleaned. To avoid this, {{{false}}} can be passed in as a second 
         * parameter
         */
        OutputFileHandler output_file_handler("TestSolvingPdeTutorial");
        
        /* Create an {{{out_stream}}}, which is a stream to a particular file. An {{{out_stream}}}
         * is a pointer to a ofstream */
        out_stream p_file = output_file_handler.OpenOutputFile("linear_solution.txt");

        /* Loop over the entries of the solution */
        for (unsigned i=0; i<result_repl.size(); i++)
        {
            /* Get the x and y-values of the node corresponding to this entry. The method 
             * {{{GetNode}}} on the mesh class returns a pointer to a {{{Node}}} */
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];

            /* Get the computed solution at this node from the {{{ReplicatableVector}}} */
            double u = result_repl[i];

            /* Finally, write x, y and u to the output file. The solution can be 
             * visualised in (eg) matlab, using the commands: 
             * {{{sol=Load('linear_solution.txt'); plot3(sol(:,1),sol(:,2),sol(:,3);}}}*/
            (*p_file) << x << " " << y << " " << u << "\n"; 
        }

        /* All Petsc {{{Vec}}}s should be destroyed when they are no longer needed */
        VecDestroy(result);
    }
    
/* ''Remember the semi-colon at the end of every class you define!'' */
};
    


#endif /*TESTSOLVINGPDESTUTORIAL_HPP_*/
