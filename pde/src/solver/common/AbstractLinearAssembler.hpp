#ifndef _ABSTRACTLINEARASSEMBLER_HPP_
#define _ABSTRACTLINEARASSEMBLER_HPP_


#include <vector>
#include <iostream>
#include <petscvec.h>

#include "AbstractStaticAssembler.hpp"


/**
 *  AbstractLinearAssembler. See AbstractAssembler for usage.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
class AbstractLinearAssembler : public AbstractStaticAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>
{
private:
    bool mMatrixIsConstant;
    double mLinearSolverRelativeTolerance;
    double mLinearSolverAbsoluteTolerance;
    bool mUseLinearSolverAbsoluteTolerance;

protected:
    
    /** Hack for dynamic mixin */
    void SetMatrixIsConst(bool matrixIsConstant = true)
    {
         mMatrixIsConstant = matrixIsConstant;
     }
    
    /**
     * Apply Dirichlet boundary conditions to the linear system.
     */
    void ApplyDirichletConditions(Vec /* unused */, bool applyToMatrix)
    {
        this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), applyToMatrix);
    }
    
    /**
     * Create the linear system object if it hasn't been already.
     * 
     * Can use an initial solution as PETSc template, or base it on the mesh size.
     */
    void InitialiseForSolve(Vec initialSolution)
    {
        if (this->mpLinearSystem == NULL)
        {
            if (initialSolution == NULL)
            {
                // Static problem, create linear system using the size
                unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
                this->mpLinearSystem = new LinearSystem(size);
            }
            else
            {
                // Use the currrent solution (ie the initial solution)
                // as the template in the alternative constructor of
                // LinearSystem. This is to avoid problems with VecScatter.
                this->mpLinearSystem = new LinearSystem(initialSolution);
            }

            this->mpLinearSystem->SetMatrixIsConstant(mMatrixIsConstant);
            if(mUseLinearSolverAbsoluteTolerance)
            {
                this->mpLinearSystem->SetAbsoluteTolerance(mLinearSolverAbsoluteTolerance);
            }
            else
            {
                this->mpLinearSystem->SetRelativeTolerance(mLinearSolverRelativeTolerance);
            }
        }
    }
        
    bool ProblemIsNonlinear()
    {
        return false;
    }
    
    /**
     * Solve a static pde, or a dynamic pde for 1 timestep.
     * 
     * The mesh, pde and boundary conditions container must be set first.
     */
    virtual Vec StaticSolve(Vec currentSolutionOrGuess=NULL,
                            double currentTime=0.0,
                            bool assembleMatrix=true)
    {
        this->AssembleSystem(true, assembleMatrix, currentSolutionOrGuess, currentTime);
        PetscInt vec_size;
        if (currentSolutionOrGuess)
        {
            VecGetSize(currentSolutionOrGuess, &vec_size);
        }
        if (currentSolutionOrGuess && (unsigned)vec_size == this->mpLinearSystem->GetSize())
        {
            return this->mpLinearSystem->Solve(currentSolutionOrGuess);
        }
        else
        {
            // When solving on a flagged mesh, the linear system is smaller than the current solution,
            // so we can't use the current solution as an initial guess for the linear solver.
            return this->mpLinearSystem->Solve();
        }
    }
    
    
public:
    AbstractLinearAssembler(unsigned numQuadPoints = 2) :
            AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM, NON_HEART, CONCRETE>(numQuadPoints),
            mMatrixIsConstant(true),
            mLinearSolverRelativeTolerance(1e-6),
            mUseLinearSolverAbsoluteTolerance(false)
    {
    }
    
    /**
     *  Destructor: ensures that the linear solver is thrown away.
     */
    ~AbstractLinearAssembler()
    {
    }
    
    /**
     *  Solve the static pde.
     * 
     *  The mesh, pde and boundary conditions container must be set before Solve() 
     *  is called.
     */
    virtual Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)
    {
        /// \todo move the asserts into PrepareForSolve()
        assert(this->mpMesh!=NULL);

        // have to comment this out because the flagged mesh assembler uses a different
        // bcc and so this is null.
        //assert(this->mpBoundaryConditions!=NULL);
        
        this->PrepareForSolve();
        this->InitialiseForSolve(currentSolutionOrGuess);
        return this->StaticSolve(currentSolutionOrGuess, currentTime);
    }
    
    
    void SetLinearSolverRelativeTolerance(double relativeTolerance)
    {
        assert(this->mpLinearSystem==NULL);
        mLinearSolverRelativeTolerance = relativeTolerance;
    }

    void SetLinearSolverAbsoluteTolerance(double absoluteTolerance)
    {
        assert(this->mpLinearSystem==NULL);
        mUseLinearSolverAbsoluteTolerance = true;
        mLinearSolverAbsoluteTolerance = absoluteTolerance;
    }
    
    /*
    void DebugWithSolution(Vec sol)
    {
        std::cout<<"\n\nWS: This is the matrix:\n";
        mpLinearSystem->DisplayMatrix();
        std::cout<<"\n\nWS: This is the righthand side:\n";
        mpLinearSystem->DisplayRhs();
        std::cout<<"\n\nWS: This is the solution:\n";
        VecView(sol, PETSC_VIEWER_STDOUT_WORLD);
    }
    
    void Debug()
    {
        std::cout<<"\n\nThis is the matrix:\n";
        mpLinearSystem->DisplayMatrix();
        std::cout<<"\n\nThis is the righthand side:\n";
        mpLinearSystem->DisplayRhs();
    }
    */
};

#endif //_ABSTRACTLINEARASSEMBLER_HPP_
