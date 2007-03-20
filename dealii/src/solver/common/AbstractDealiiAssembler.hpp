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


/**
 *  Abstract assembler with common functionality for most assemblers.
 */
template<unsigned DIM>
class AbstractDealiiAssembler
{
protected:
    /*< The mesh to be solve on */
    Triangulation<DIM>*   mpMesh;
    
    // an FE_Q object seems to be equivalent to our basis functions
    // It is templated over dimension, with the order of the bases taken
    // in the contructor
    // note that this must be defined before mDofHandler!
    /*FE_Q<DIM>            mFe;*/
    
    /*< Degrees of freedom handler */
    DoFHandler<DIM>      mDofHandler;
    
    /*< A structure needed for handling hanging nodes */
    ConstraintMatrix     mHangingNodeConstraints;
    
    /*< A structure needed for setting up a sparse matrix
     */
    SparsityPattern      mSparsityPattern;
    
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
     *  Initialise the system matrix, system rhs vector, current solution vector, and 
     *  hanging nodes constraints objects
     *  
     *  This should only be called when the DofHandler knows how big it should be, ie what
     *  the FE system is, ie after something like "mDofHandler.distribute_dofs(mFeSystem);"
     *  has been called
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
        
        unsigned elem_counter = 0;
        
        while (element_iter!=mDofHandler.end())  // huh? mDof.end() returns an element iterator?
        {
            // zero the small matrix and vector
            element_matrix = 0;
            element_rhs = 0;
            
            element_iter->get_dof_indices(local_dof_indices);
            
            AssembleOnElement(element_iter,
                              element_rhs,
                              element_matrix,
                              assembleVector,
                              assembleMatrix);
                              
            /*if(assembleMatrix)
            {
                std::cout << elem_counter++ << " of " << mpMesh->n_active_cells() << "\n" << std::flush;
            }*/
            
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
        //if(assembleMatrix) { std::cout << "\n"; }
        
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
};

#endif /*ABSTRACTDEALIIASSEMBLER_HPP_*/
