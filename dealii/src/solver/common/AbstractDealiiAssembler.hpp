#ifndef ABSTRACTDEALIIASSEMBLER_HPP_
#define ABSTRACTDEALIIASSEMBLER_HPP_

#include <grid/grid_out.h>
#include <grid/grid_in.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>

#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>

#include <fe/fe_q.h>
#include <fe/fe_system.h>

#include <dofs/dof_tools.h>

#include <fe/fe_values.h>
#include <base/quadrature_lib.h>

#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_gmres.h>

#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>

// for dealing with hanging nodes..
#include <dofs/dof_constraints.h>


#include <numerics/solution_transfer.h>

#include <math.h>
#include <iostream>

#include "Exception.hpp"
#include "TriangulationVertexIterator.hpp"
#include "DofVertexIterator.hpp"

/**
 *  Abstract assembler with common functionality for most assemblers.
 * 
 *  The concrete class needs to implement 
 *  1. a constructor, which sets up an FE_Q or FESystem object, calls DistributeDofs
 *  and sets up mDofsPerElement
 *  2. have a method called DistributeDofs() which just calls
 *  this->mDofHandler.distribute_dofs(fe), where fe is the FE_Q or FESystem object
 *  3. Implements AssembleOnElement
 *  4. Implements ApplyDirichletBoundaryConditions
 *  5. A Solve method which actually calls AssembleSystem and then does the linear
 *  system solving/Newton's method, time-looping etc..
 */
template<unsigned DIM>
class AbstractDealiiAssembler
{
protected:
    /*< The mesh to be solve on */
    Triangulation<DIM>*   mpMesh;
    
    /*< Degrees of freedom handler */
    DoFHandler<DIM>       mDofHandler;
    
    /*< A structure needed for handling hanging nodes */
    ConstraintMatrix      mHangingNodeConstraints;
    
    /*< A structure needed for setting up a sparse matrix */
    SparsityPattern       mSparsityPattern;
    
    /**
     *  The main system matrix. Eg, the global stiffness matrix in a static linear problem
     *  or the Jacobian in a nonlinear problem
     */
    SparseMatrix<double> mSystemMatrix;
    
    /**
     *  The main rhs vector. Eg the global load vector in a static linear problem, or
     *  the residual in a nonlinear problem
     */
    Vector<double>       mRhsVector;
    
    /**
     *  The current solution in a time-dependent problem, or current guess in a nonlinear
     *  static problem (or both in a time-dependent nonlinear problem), or final solution 
     *  in a static linear problem
     */
    Vector<double>       mCurrentSolution;
    
    /**
     *  The size of element stiffness vector. Depends on the FeSystem object (which
     *  depends on the concrete class). Must be set up to the correct value in the 
     *  constructor of the concrete class
     */
    unsigned             mDofsPerElement;
    
    /** 
     *  The method RefineCoarsen in this class refines/coarsens the mesh according to whether 
     *  elements have been labelled for refinement/coarsening. It also interpolates the
     *  current solution vector onto the new mesh. The user may want other data 
     *  interpolated, in which case they can pass in those vectors using
     *  AddVectorForInterpolation(). Such vectors are stored here. They should be 
     *  vertex-wise vectors, ie if x is the vector to be interpolated onto the new mesh
     *  x(i)=value for x at vertex i. It will be linearly interpolated
     */
    std::vector<Vector<double>*> mVectorsToInterpolate;
 
    /**
     *  The main function to be implemented in the concrete class
     * 
     *  Assemble of the element matrix and/or element vector for the given element. Called
     *  by AssembleSystem
     * 
     *  @elementIter Iterator pointing at current element
     *  @elementRhs Small vector to be filled in. Should be of size AbstractDealiiAssembler::mDofsPerElement
     *  @elementMatrix Small matrix to be filled in. Should be of square, of size AbstractDealiiAssembler::mDofsPerElement
     *  @assembleVector Whether to assemble the small vector
     *  @assembleMatrix Whether to assemble the small matrix
     */
    virtual void AssembleOnElement(typename DoFHandler<DIM>::active_cell_iterator  elementIter,
                                   Vector<double>&        elementRhs,
                                   FullMatrix<double>&    elementMatrix,
                                   bool                   assembleVector,
                                   bool                   assembleMatrix)=0;
                                   
    /**
     *  A pure method which needs to implemented in the concrete class which 
     *  applies the dirichlet boundary conditions to the system matrix and 
     *  system vector. Called at the end of AssembleSystem()
     */
    virtual void ApplyDirichletBoundaryConditions()=0;
    
    
    /**
     *  A pure method which needs to be implemented in the concrete class which
     *  distributes the dofs using which fe object is being used in the concrete 
     *  class. For example if the concrete class uses a simple FE_Q<DIM> object 
     *  called mFe, this method should just be 
     *  this->mDofHandler.distribute_dofs(mFe);
     */
    virtual void DistributeDofs()=0;


    /**
     *  Initialise the system matrix, system rhs vector, current solution vector, and 
     *  hanging nodes constraints objects
     *  
     *  This should only be called when the DofHandler knows how big it should be, ie what
     *  the FE system is, ie after DistributeDofs() has been called
     */
    void InitialiseMatricesVectorsAndConstraints()
    {
        // HANGING NODES, SEE DEALII TUTORIAL 2
        // clear the constrait matrix
        mHangingNodeConstraints.clear();
        // form constraints
        DoFTools::make_hanging_node_constraints(mDofHandler, mHangingNodeConstraints);
        // some postprocessing
        mHangingNodeConstraints.close();
        
        // form sparsity pattern
        mSparsityPattern.reinit(mDofHandler.n_dofs(),
                                mDofHandler.n_dofs(),
                                mDofHandler.max_couplings_between_dofs());
                                
        DoFTools::make_sparsity_pattern(mDofHandler, mSparsityPattern);
        
        // see dealii tutorial 2
        mHangingNodeConstraints.condense(mSparsityPattern);
        
        
        mSparsityPattern.compress();
        
        // initialise vectors and matrices
        mSystemMatrix.reinit(mSparsityPattern);
        mCurrentSolution.reinit(mDofHandler.n_dofs());
        mHangingNodeConstraints.distribute(mCurrentSolution);

        mRhsVector.reinit(mDofHandler.n_dofs());
    }
    
    
    /**
     *  AssembleSystem
     * 
     *  Loops over the elements and assembles the system matrix (global stiffness matrix, 
     *  or jacobian etc, depending on problem) and system vector. It calls the pure methods
     *  AssembleOnElement() and ApplyDirichletBoundaryConditions() which need to be 
     *  implemented in the concrete class.
     * 
     *  @param assembleVector A boolean saying whether to assemble the system vector
     *  @param assembleMatrix A boolean saying whether to assemble the system matrix
     */
    void AssembleSystem(bool assembleVector, bool assembleMatrix)
    {
        // if this fails, mDofsPerElement hasn't been set. It should
        // have been set in the constructor of the concrete class, using
        // something like
        // mDofsPerElement = mFeSystem.dofs_per_cell
        assert(mDofsPerElement>0);
        
        // if this fails, InitialiseMatricesVectorsAndConstraints() probably
        // hasn't been called..
        assert(mRhsVector.size()!=0);

        FullMatrix<double>   element_matrix(mDofsPerElement, mDofsPerElement);
        Vector<double>       element_rhs(mDofsPerElement);
        // the dofs associated with the nodes of an element
        std::vector<unsigned> local_dof_indices(mDofsPerElement);
        
        typename DoFHandler<DIM>::active_cell_iterator  element_iter = mDofHandler.begin_active();
        
        if (assembleVector)
        {
            mRhsVector = 0;
        }

        if (assembleMatrix)
        {
            mSystemMatrix = 0;
        }
                
        //unsigned elem_counter = 0;
        
        // loop over elements
        while (element_iter!=mDofHandler.end()) 
        {
            // zero the small matrix and vector
            element_matrix = 0;
            element_rhs = 0;
            
            element_iter->get_dof_indices(local_dof_indices);
            
            // assemble the small matrix and vector
            AssembleOnElement(element_iter,
                              element_rhs,
                              element_matrix,
                              assembleVector,
                              assembleMatrix);
            
            // add to the full matrix and vector
            for (unsigned i=0; i<mDofsPerElement; i++)
            {
                if (assembleMatrix)
                {
                    for (unsigned j=0; j<mDofsPerElement; j++)
                    {
                        mSystemMatrix.add(local_dof_indices[i],
                                          local_dof_indices[j],
                                          element_matrix(i,j));
                    }
                }

                if (assembleVector)
                {
                    mRhsVector(local_dof_indices[i]) += element_rhs(i);                    
                }
            }
            
            element_iter++;
        }
        
        // note this has to be done before applying dirichlet bcs
        if (assembleMatrix)
        {
            mHangingNodeConstraints.condense(mSystemMatrix);
        }
        if (assembleVector)
        {
            mHangingNodeConstraints.condense(mRhsVector);
        }

        ApplyDirichletBoundaryConditions();
        
        // stupid thing won't quit if variables become NaN (and says norm_rhs_vec=0 too!), 
        // so have to check here
        for(unsigned i=0; i<mRhsVector.size(); i++)
        {
            if( isnan(mRhsVector(i)))
            {
                EXCEPTION("Component of the system rhs vector became NaN - check for division by zero."); 
            }
        }
    }
    
    /**
     *  Compute the L2 norm of the current residual vector divided by it's length.
     * 
     *  This method obviously only makes sense if the assembler is for a nonlinear
     *  PDE, and assumes mRhsVector is the residual vector. 
     *  
     */
    double CalculateResidualNorm()
    {
        return mRhsVector.norm_sqr()/mDofHandler.n_dofs();
    }
    
    /**
     *  Take one Newton step.
     * 
     *  This method obviously only makes sense if the assembler is for a nonlinear
     *  PDE, and assumes mSystemMatrix is the jacobian matrix and mRhsVector is the 
     *  residual vector. It solves for the update vector, and determines best damping 
     *  value.
     * 
     *  NOTE: gmres, identity preconditioning, num iterations etc are all hardcoded 
     *  in here at the moment.
     */
    void TakeNewtonStep()
    {
        // compute Jacobian
        AssembleSystem(true, true);
        
        // solve the linear system
        SolverControl  solver_control(200000, 1e-6, false, true);
        PrimitiveVectorMemory<> vector_memory;
        
        Vector<double> update;
        update.reinit(mDofHandler.n_dofs());
        
        SolverGMRES<>::AdditionalData gmres_additional_data(200);
        SolverGMRES<>  gmres(solver_control, vector_memory, gmres_additional_data);
        
        gmres.solve(mSystemMatrix, update, mRhsVector, PreconditionIdentity());
    
        // deal with hanging nodes - form a continuous solutions
        mHangingNodeConstraints.distribute(update);
        
        // save the old current solution
        Vector<double> old_solution = mCurrentSolution;
        
        double best_norm_resid = 1e10;
        double best_damping_value = 0.0;
        
        std::vector<double> damping_values;
        damping_values.reserve(12);
        damping_values.push_back(0.0);
        damping_values.push_back(0.05);
        for (unsigned i=1; i<=10; i++)
        {
            damping_values.push_back((double)i/10.0);
        }
        
        for (unsigned i=0; i<damping_values.size(); i++)
        {
            mCurrentSolution.equ(1.0, old_solution, -damping_values[i], update);
            
            // compute residual
            AssembleSystem(true, false);
            double norm_resid = CalculateResidualNorm();
            
            std::cout << "\tTesting s = " << damping_values[i] << ", |f| = " << norm_resid << "\n" << std::flush;
            if (norm_resid < best_norm_resid)
            {
                best_norm_resid = norm_resid;
                best_damping_value = damping_values[i];
            }
        }
        
        
        if (best_damping_value == 0.0)
        {
            #define COVERAGE_IGNORE
            std::cout << "\nResidual does not decrease in newton direction, quitting\n" << std::flush;
            assert(0); // got a weird error once with an exception here, hence the assert(0);
            #undef COVERAGE_IGNORE
        }
        else
        {
            std::cout << "\tBest s = " << best_damping_value << "\n"  << std::flush;
        }
        // implement best update and recalculate residual
        mCurrentSolution.equ(1.0, old_solution, -best_damping_value, update);
    }    
    
    
public :
    AbstractDealiiAssembler(Triangulation<DIM>* pMesh) :
            mDofHandler(*pMesh)  // associate the mesh with the dof handler
    {
        // probably will fail in the mDofHandler(*pMesh) line above before here if
        // pMesh==NULL
        assert(pMesh!=NULL);
        mpMesh = pMesh;
               
        // initially set mDofsPerElement to be zero so can check it
        // has been set in AssembleSystem(). It should be set in the
        // constructor of a concrete class using something like
        // mDofsPerElement = mFeSystem.dofs_per_cell;
        mDofsPerElement = 0;
    }
    
    
    /**
     *  Refine and coarsen a mesh, depending on whether the refine_flag and coarsen_flag
     *  has been set on the elements of the mesh. This method calls 
     *  execute_coarsening_and_refinement() on the mesh in order to do the refinement
     *  and coarsening, but also handles interpolation as well. The current solution
     *  vector is automatically interpolated onto the new mesh. The user can also add 
     *  extra vector to be interpolated as well, by calling AddVectorForInterpolation()
     *  before this method. These vectors should have vertex-wise data (ie x(i) is the 
     *  value of x at vertex i), and will be linearly interpolated onto the new mesh.
     *  The vector will be resized and it's values changed in this method
     * 
     *  //\todo: possibly incredibly inefficient - make efficient
     */
    void RefineCoarsen()
    {
        // a linear fe object for linear interpolation of any extra vectors
        FE_Q<DIM> linear_fe(1);
        
        // shouldn't use more than one SolutionTransfer (interpolation doesn't work
        // if two or more soluation transfer objects are used on the same mesh), 
        // but we have to (we want to use one solution transfer for the current solution,
        // with interpolation based on whatever basis functions are used in the concrete
        // version of this class, and another with linear interpolation on the extra 
        // vectors. To get round this, we make a copy of the mesh(!)
        // TODO: make this efficient
        Triangulation<DIM> copy_of_mesh;
        copy_of_mesh.copy_triangulation(*mpMesh);
        
        DoFHandler<DIM> linear_dof_handler(copy_of_mesh);
        linear_dof_handler.distribute_dofs(linear_fe);

        // prepare both meshes
        mpMesh->prepare_coarsening_and_refinement();
        copy_of_mesh.prepare_coarsening_and_refinement();

        // solution transfer for the current solution
        SolutionTransfer<DIM,double> transfer_for_cur_soln(mDofHandler);
        transfer_for_cur_soln.prepare_for_coarsening_and_refinement(mCurrentSolution);
        
        // if there are any other vecs to interpolate, prepare a 
        // solution transfer for them
        unsigned num_vecs_to_interpolate = mVectorsToInterpolate.size();
        SolutionTransfer<DIM,double> transfer_for_other_vecs(linear_dof_handler);
        std::vector<Vector<double> > vecs_by_dofs(num_vecs_to_interpolate);

        if(num_vecs_to_interpolate>0)
        {
            for(unsigned i=0; i<num_vecs_to_interpolate; i++)
            {
                // convert the vertex-wise data to dof-wise data
                vecs_by_dofs[i].reinit(mpMesh->n_used_vertices());
                vecs_by_dofs[i]=0;

                DofVertexIterator<DIM> dof_vertex_iter(mpMesh, &linear_dof_handler);
                while(!dof_vertex_iter.ReachedEnd())
                {
                    unsigned index = dof_vertex_iter.GetVertexGlobalIndex();
                    unsigned dof = dof_vertex_iter.GetDof(0);
                    vecs_by_dofs[i](dof) = (*(mVectorsToInterpolate[i]))(index);
                    
                    dof_vertex_iter.Next();
                }
            }
    
            transfer_for_other_vecs.prepare_for_coarsening_and_refinement(vecs_by_dofs);
        }

        // refine coarsen
        mpMesh->execute_coarsening_and_refinement();
        copy_of_mesh.execute_coarsening_and_refinement();
    
        // redistribute dofs
        DistributeDofs();
        linear_dof_handler.distribute_dofs(linear_fe);

        // interpolate to get the new current solution
        Vector<double> new_current_soln(this->mDofHandler.n_dofs());
        transfer_for_cur_soln.interpolate(mCurrentSolution, new_current_soln); 

        // resize matrices, vectors etc, re-setup hanging node constraints... 
        InitialiseMatricesVectorsAndConstraints();

        mCurrentSolution = new_current_soln;

        // if there were other vecs to interpolate
        if(num_vecs_to_interpolate>0)
        {
            // set up new interpolated vectors
            std::vector<Vector<double> > new_vecs_by_dofs(num_vecs_to_interpolate);

            for(unsigned i=0; i<num_vecs_to_interpolate; i++)
            {
                new_vecs_by_dofs[i].reinit(linear_dof_handler.n_dofs());
            }

            // interpolate and save in this new vectors
            transfer_for_other_vecs.interpolate(vecs_by_dofs, new_vecs_by_dofs);
                        
            // convert from dof-wise data back to vertex-wise data and put back
            // in the original objects
            for(unsigned i=0; i<num_vecs_to_interpolate; i++)
            {
                (*(mVectorsToInterpolate[i])).reinit(mpMesh->n_vertices());
                DofVertexIterator<DIM> dof_vertex_iter(mpMesh, &linear_dof_handler);
                while(!dof_vertex_iter.ReachedEnd())
                {
                    unsigned index = dof_vertex_iter.GetVertexGlobalIndex();
                    unsigned dof = dof_vertex_iter.GetDof(0);

                    (*(mVectorsToInterpolate[i]))(index) = new_vecs_by_dofs[i](dof);
                    dof_vertex_iter.Next();
                } 
            }
        }

        ApplyDirichletBoundaryConditions();
        
        mVectorsToInterpolate.clear();
    }    
    
    /**
     *  Get a reference to the current solution vector. This is a dof-wise vector, is
     *  x(i) is the solution for degree of freedom i. GetSolutionAtVertices() will
     *  probably be more useful.
     */
    Vector<double>& rGetCurrentSolution()
    {
        return mCurrentSolution;
    }    
    
    /**
     *  Get a reference to dof handler. Won't generally be needed. 
     */
    DoFHandler<DIM>& rGetDofHandler()
    {
        return mDofHandler;
    }

    /**
     *  Get the mesh
     */
    Triangulation<DIM>* GetMesh()
    {
        return mpMesh;
    }
        
    /** 
     *  Get the solution (for a particular unknown) as a vertex-wise vector (ie 
     *  call GetSolutionAtVertices(x, j) and then x(i) will be the solution at vertex 
     *  i, for unknown j
     */
    void GetSolutionAtVertices(Vector<double>& rSolutionAtVertices, unsigned unknown=0)
    {
        rSolutionAtVertices.reinit(mpMesh->n_vertices());
        DofVertexIterator<DIM> vertex_iter(mpMesh, &mDofHandler);
        while (!vertex_iter.ReachedEnd())
        {
            unsigned vertex_index = vertex_iter.GetVertexGlobalIndex();
            rSolutionAtVertices(vertex_index) = mCurrentSolution(vertex_iter.GetDof(unknown));
            vertex_iter.Next();
        }
    }
    
    /**
     *  Add a vector for interpolation when RefineCoarsen is called. See
     *  RefineCoarsen
     */
    void AddVectorForInterpolation(Vector<double>* pVector)
    {
        mVectorsToInterpolate.push_back(pVector);
    }
    
    virtual ~AbstractDealiiAssembler()
    {}
};

#endif /*ABSTRACTDEALIIASSEMBLER_HPP_*/


