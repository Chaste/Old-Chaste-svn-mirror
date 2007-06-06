#ifndef _ABSTRACTLINEARASSEMBLER_HPP_
#define _ABSTRACTLINEARASSEMBLER_HPP_


#include <vector>
#include <iostream>
#include <petscvec.h>

#include "AbstractAssembler.hpp"
#include "AbstractLinearSolver.hpp"
#include "SimpleLinearSolver.cpp"

/**
 *  AbstractLinearAssembler. See AbstractAssembler for usage.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractLinearAssembler : public AbstractAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{

protected:
    /**
     * The linear solver used to solve the linear system at each time step.
     */
    AbstractLinearSolver *mpLinearSolver;
    bool mWeAllocatedSolverMemory;
    
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
     * Can use a current/initial solution as PETSc template, or base it on the mesh size.
     */
    void InitialiseLinearSystem(Vec currentSolution)
    {
        if (this->mpLinearSystem == NULL)
        {
            if (currentSolution == NULL)
            {
                // Static problem, create linear system using the size
                unsigned size = PROBLEM_DIM * this->mpMesh->GetNumNodes();
                this->mpLinearSystem = new LinearSystem(size);
            }
            else
            {
                // Use the currrent solution (ie the initial solution)
                // as the template in the alternative constructor of
                // LinearSystem. This appears to avoid problems with
                // VecScatter.
                this->mpLinearSystem = new LinearSystem(currentSolution);
            }
        }
    }
    
public:
    AbstractLinearAssembler(unsigned numQuadPoints = 2,
                            double linearSolverRelativeTolerance = 1e-6) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {
        mpLinearSolver = new SimpleLinearSolver(linearSolverRelativeTolerance);
        mWeAllocatedSolverMemory = true;
    }
    
    /**
     *  Destructor: ensures that the LinearSystem is thrown away.
     */
    ~AbstractLinearAssembler()
    {
        if (mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
    }
    
    /**
     * Manually re-set the linear system solver (which by default is a SimpleLinearSolver).
     * 
     * \todo Cover this method!
     */
    virtual void SetLinearSolver(AbstractLinearSolver *pLinearSolver)
    {
        if (mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
        mpLinearSolver = pLinearSolver;
        mWeAllocatedSolverMemory = false;
    }
    
    virtual Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0)=0;
    
    
    
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
