#ifndef TESTABSTRACTDEALIIASSEMBLER_HPP_
#define TESTABSTRACTDEALIIASSEMBLER_HPP_


#include <cxxtest/TestSuite.h>
#include "AbstractDealiiAssembler.hpp"
#include "grid/tria_boundary_lib.h"

// simple concrete assembler for laplace's equation 
template<unsigned DIM>
class LaplacesAssembler : public AbstractDealiiAssembler<DIM>
{
private:    
    FE_Q<DIM> mFe;    
    
    void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                           Vector<double>&        elementRhs,
                           FullMatrix<double>&    elementMatrix,
                           bool                   assembleVector,
                           bool                   assembleMatrix)
    {
        static QGauss<DIM>   quadrature_formula(3);
        const unsigned n_q_points = quadrature_formula.n_quadrature_points;
        

        FEValues<DIM> fe_values(mFe, quadrature_formula,
                                UpdateFlags(update_values    |
                                            update_gradients |
                                            update_JxW_values));
                                            
                                            
        const unsigned dofs_per_element = mFe.dofs_per_cell;
        
        std::vector<unsigned> local_dof_indices(dofs_per_element);
        
        elementMatrix = 0;
        elementRhs = 0;
        
        elementIter->get_dof_indices(local_dof_indices);
        
        fe_values.reinit(elementIter); // compute fe values for this element
        
        for (unsigned q_point=0; q_point<n_q_points; q_point++)
        {
            for (unsigned i=0; i<dofs_per_element; i++)
            {
                if(assembleMatrix)
                {
                    for (unsigned j=0; j<dofs_per_element; j++)
                    {
                        elementMatrix(i,j) +=   fe_values.shape_grad(i,q_point)
                                              * fe_values.shape_grad(j,q_point)
                                              * fe_values.JxW(q_point);
                    }
                }
                if(assembleVector)
                {
                    elementRhs(i) +=  - fe_values.shape_value(i,q_point)
                                      * fe_values.JxW (q_point);
                }
            }
        }
    }

    void ApplyDirichletBoundaryConditions()
    {
        std::map<unsigned int,double> boundary_values;
        VectorTools::interpolate_boundary_values(this->mDofHandler,
                                                 0,
                                                 ZeroFunction<DIM>(),
                                                 boundary_values);
                                                 
        MatrixTools::apply_boundary_values(boundary_values,
                                           this->mSystemMatrix,
                                           this->mCurrentSolution,
                                           this->mRhsVector);
    }

    void DistributeDofs()
    {
        this->mDofHandler.distribute_dofs(mFe);
    }


public :  
    LaplacesAssembler(Triangulation<DIM>* pMesh) 
        : AbstractDealiiAssembler<DIM>(pMesh),
          mFe(1)
    {  
        // distribute dofs
        this->mDofHandler.distribute_dofs(mFe);
        this->InitialiseMatricesVectorsAndConstraints();
        
        this->mDofsPerElement = mFe.dofs_per_cell;
        
        assert(pMesh!=NULL);
    }
    
    
    void Solve()
    {
        this->AssembleSystem(true, true);
        
        SolverControl solver_control(1000, 1e-12, false, false);
        PrimitiveVectorMemory<> vector_memory;
        SolverCG<>              cg(solver_control, vector_memory);
        
        cg.solve(this->mSystemMatrix, this->mCurrentSolution, this->mRhsVector, PreconditionIdentity());

        // remember this!
        this->mHangingNodeConstraints.distribute(this->mCurrentSolution);
    }
};



class TestAbstractDealiiAssembler : public CxxTest::TestSuite
{
public:
    // solve laplaces equation on a circle and check the result
    void TestWithLaplacesEquation2d()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(3);

        LaplacesAssembler<2> laplaces(&mesh);
        laplaces.Solve();
        
        Vector<double> solution;
        laplaces.GetSolutionAtVertices(solution);
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);

        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            // to derive this, write laplacian is cylindrical coords 
            double expected_val = 0.25*(r*r-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 1e-2);
            vertex_iter.Next();
        }
    }


    // solve laplaces equation on a circle, using a mesh with hanging nodes,
    // and check the result, to verify the constraints code is working properly
    void TestWithLaplacesEquation2dWithHangingNodes()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(2);

        unsigned counter=0;
        Triangulation<2>::active_cell_iterator  element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            if(counter%2==0)
            {
                element_iter->set_refine_flag();
            }
            
            element_iter++;
            counter++;
        }
        mesh.execute_coarsening_and_refinement();

        LaplacesAssembler<2> laplaces(&mesh);
        laplaces.Solve();
        
        Vector<double> solution;
        laplaces.GetSolutionAtVertices(solution);
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);

        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            // to derive this, write laplacian is cylindrical coords 
            double expected_val = 0.25*(r*r-1);
            // std::cout << posn[0] << " " << posn[1] << " " << solution(index) << "\n";
            TS_ASSERT_DELTA(expected_val, solution(index), 1e-2);
            vertex_iter.Next();
        }
    }


    // solve laplaces equation on a sphere and check the result
    void TestWithLaplacesEquation3d()
    {
        Triangulation<3> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<3> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(3);

        LaplacesAssembler<3> laplaces(&mesh);
        laplaces.Solve();
        
        Vector<double> solution;
        laplaces.GetSolutionAtVertices(solution);
        
        TriangulationVertexIterator<3> vertex_iter(&mesh);

        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<3> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            // to derive this, write laplacian is spherical coords 
            double expected_val = (1.0/6.0)*(r*r-1); 
            TS_ASSERT_DELTA(expected_val, solution(index), 1e-2);
            vertex_iter.Next();
        }
    }


    // solve laplaces equation on a circle, then refine every other element, check 
    // the solution vector has been interpolated correctly, check can solve on the
    // refined mesh correctly, then repeat.
    void TestWithLaplacesEquation2dWithRefinement()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(2);

        LaplacesAssembler<2> laplaces(&mesh);
        laplaces.Solve();
        
        Vector<double> solution;
        laplaces.GetSolutionAtVertices(solution);

        // refine every other element
        unsigned counter=0;
        Triangulation<2>::active_cell_iterator  element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            if(counter%2==0)
            {
                element_iter->set_refine_flag();
            }
            
            element_iter++;
            counter++;
        }

        laplaces.RefineCoarsen();

        // get the interpolated current solution
        laplaces.GetSolutionAtVertices(solution);
        // the size of the new interpolated solution should be the size of the mesh        
        TS_ASSERT_EQUALS(mesh.n_vertices(),solution.size());

        TriangulationVertexIterator<2> vertex_iter(&mesh);
        
        // the solution will have been interpolated onto the new nodes, should
        // be approx correct
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            double expected_val = 0.25*(r*r-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 5e-2); //reduced tol
            vertex_iter.Next();
        } 
        

        // solve again
        laplaces.Solve();
        laplaces.GetSolutionAtVertices(solution);

        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            // to derive this, write laplacian is cylindrical coords 
            double expected_val = 0.25*(r*r-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 1e-2);
            vertex_iter.Next();
        } 


        ////////////////////////
        // and again...
        ////////////////////////
        counter=0;
        element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            if(counter%2==0)
            {
                element_iter->set_refine_flag();
            }
            
            element_iter++;
            counter++;
        }

        laplaces.RefineCoarsen();
        
        // get the interpolated current solution
        laplaces.GetSolutionAtVertices(solution);
        // the size of the new interpolated solution should be the size of the mesh        
        TS_ASSERT_EQUALS(mesh.n_vertices(),solution.size());

        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            double expected_val = 0.25*(r*r-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 5e-2); //reduced tol
            vertex_iter.Next();
        } 

        // solve again
        laplaces.Solve();
        laplaces.GetSolutionAtVertices(solution);

        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            // to derive this, write laplacian is cylindrical coords 
            double expected_val = 0.25*(r*r-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 1e-2);
            vertex_iter.Next();
        } 
    }

    // test the RefineCoarsen method on the AbstractDealiiAssembler interpolates
    // extra vectors correctlu
    void TestInterpolationWithLaplacesEquation2d()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(2);

        LaplacesAssembler<2> laplaces(&mesh);

        Vector<double> some_vector(mesh.n_vertices());
        Vector<double> another_vector(mesh.n_vertices());

        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            // linear data
            some_vector(index) = posn[0] + 2*posn[1];
            another_vector(index) = 5*posn[0] - posn[1];
            vertex_iter.Next();
        } 

        // refine every other element
        unsigned counter=0;
        Triangulation<2>::active_cell_iterator  element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            if(counter%2==0)
            {
                element_iter->set_refine_flag();
            }
            
            element_iter++;
            counter++;
        }

        unsigned num_vertices_before = mesh.n_vertices();

        laplaces.AddVectorForInterpolation(&some_vector);
        laplaces.AddVectorForInterpolation(&another_vector);
        laplaces.RefineCoarsen();
    
        // check something was refined
        TS_ASSERT_LESS_THAN(num_vertices_before,mesh.n_vertices());
        
        // check the vector has grown
        TS_ASSERT_EQUALS(some_vector.size(), mesh.n_vertices());
        TS_ASSERT_EQUALS(another_vector.size(), mesh.n_vertices());

        // check its been interpolated properly
        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            // the vectors should have new data, which should be 
            // interpolated. note the high tolerance - the interpolation
            // works but isn't that accurate, for some reason
            TS_ASSERT_DELTA(some_vector(index), posn[0]+2*posn[1], 1e-1);
            TS_ASSERT_DELTA(another_vector(index), 5*posn[0]-posn[1], 1e-1);
            vertex_iter.Next();
        } 
    }

    // coarsening every element and test that extra vectors have correct values 
    // on the new mesh
    void TestPureCoarsening()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(2);

        LaplacesAssembler<2> laplaces(&mesh);

        Vector<double> some_vector(mesh.n_vertices());

        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            some_vector(index) = posn[0] + 2*posn[1];
            vertex_iter.Next();
        } 

        // coarsen every element
        Triangulation<2>::active_cell_iterator  element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            element_iter->set_coarsen_flag();
            element_iter++;
        }

        unsigned num_vertices_before = mesh.n_vertices();

        laplaces.AddVectorForInterpolation(&some_vector);
        laplaces.RefineCoarsen();
    
        // check mesh was coarsened
        // !! note the call of n_USED_vertices(), not n_vertices() !!
        TS_ASSERT_LESS_THAN(mesh.n_used_vertices(), num_vertices_before);
        
        // this vector should still have the size n_vertices. however,
        // only the active vertices will be filled in.
        TS_ASSERT_EQUALS(some_vector.size(), mesh.n_vertices());

        // check its been interpolated properly
        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            TS_ASSERT_DELTA(some_vector(index), posn[0]+2*posn[1], 1e-12);
            vertex_iter.Next();
        } 
    }


    // solve laplace's eqn and then coarsen mesh. check the current solution
    // has the correct values and size following this. 
    void TestPureCoarsenInterpolateSoln()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(3);

        LaplacesAssembler<2> laplaces(&mesh);
        laplaces.Solve();

        // coarsen elements
        Triangulation<2>::active_cell_iterator  element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            element_iter->set_coarsen_flag();
            element_iter++;
        }

        laplaces.RefineCoarsen();

        // get the interpolated current solution
        Vector<double> solution;
        laplaces.GetSolutionAtVertices(solution);

        TS_ASSERT_EQUALS(mesh.n_vertices(),solution.size());

        TriangulationVertexIterator<2> vertex_iter(&mesh);

        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            double expected_val = 0.25*(r*r-1);

            TS_ASSERT_DELTA(expected_val, solution(index), 1e-1); 
            vertex_iter.Next();
        } 
    }

    // solve, refine part of the mesh and coarsen another, check solution
    // and extra vec interpolated correctly, then check can solve correctly
    // on the new mesh.
    void TestRefineCoarsenAndThenSolveWithLaplacesEqn2d()
    {
        Triangulation<2> mesh;
        GridGenerator::hyper_ball(mesh);
        HyperBallBoundary<2> boundary;
        mesh.set_boundary(0, boundary);
        mesh.refine_global(3);

        LaplacesAssembler<2> laplaces(&mesh);

        // solve to set up the current_solution internally
        laplaces.Solve();

        Vector<double> some_vector(mesh.n_vertices());
        
        TriangulationVertexIterator<2> vertex_iter(&mesh);
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
  
            Point<2> posn = vertex_iter.GetVertex();
            some_vector(index) = posn[0] + 2*posn[1];
            vertex_iter.Next();
        }

        // refine and coarsen elements
        unsigned counter=0;
        Triangulation<2>::active_cell_iterator  element_iter = mesh.begin_active();
        while (element_iter!=mesh.end())  
        {
            if(counter<10)
            {
                element_iter->set_refine_flag();
            }
            else
            {
                element_iter->set_coarsen_flag();
            }
            element_iter++;
            counter++;
        }

        // add vector for interpolation and refine-coarsen the mesh
        laplaces.AddVectorForInterpolation(&some_vector);

        laplaces.RefineCoarsen();

        // get the interpolated current solution
        Vector<double> solution;
        laplaces.GetSolutionAtVertices(solution);

        // the current solution should have been interpolated in the
        // RefineCoarsen, as should have some_vector. The size should
        // be n_vvertices for both, even if it  
        TS_ASSERT_EQUALS(solution.size(),    mesh.n_vertices());
        TS_ASSERT_EQUALS(some_vector.size(), mesh.n_vertices());

        // the solution will have been interpolated onto the new nodes, should
        // be approx correct
        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();

            Point<2> posn = vertex_iter.GetVertex();
            double r = std::sqrt(posn.square());
            double expected_val = 0.25*(r*r-1);

            TS_ASSERT_DELTA(some_vector(index), posn[0]+2*posn[1], 1e-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 1e-1); 
            vertex_iter.Next();
        } 

        // solve on the RefineCoarsen-ed mesh
        laplaces.Solve();
        
        // get the new current solution
        laplaces.GetSolutionAtVertices(solution);

        // check its size and values
        TS_ASSERT_EQUALS(solution.size(), mesh.n_vertices());
        vertex_iter.Begin();
        while(!vertex_iter.End())
        {
            unsigned index = vertex_iter.GetVertexGlobalIndex();
            Point<2> posn = vertex_iter.GetVertex();

            double r = std::sqrt(posn.square());
            double expected_val = 0.25*(r*r-1);
            TS_ASSERT_DELTA(expected_val, solution(index), 5e-2); 

            vertex_iter.Next();
        } 
    }
};


#endif /*TESTABSTRACTDEALIIASSEMBLER_HPP_*/
