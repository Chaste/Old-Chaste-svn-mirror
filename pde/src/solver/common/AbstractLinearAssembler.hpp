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
template<int ELEMENT_DIM, int SPACE_DIM, int PROBLEM_DIM>
class AbstractLinearAssembler : public virtual AbstractAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{

protected:
    /**
     * The linear solver used to solve the linear system at each time step.
     */
    AbstractLinearSolver *mpLinearSolver;
    bool mWeAllocatedSolverMemory;
 
    
public:
    AbstractLinearAssembler(int numQuadPoints = 2,
                            double linearSolverRelativeTolerance = 1e-6) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(numQuadPoints)
    {
        mpLinearSolver = new SimpleLinearSolver(linearSolverRelativeTolerance);
        mWeAllocatedSolverMemory = true;
        
        this->mpLinearSystem = NULL;
        this->mMatrixIsConstant = false;
        this->mMatrixIsAssembled = false;
        
        this->mProblemIsLinear = true;
    }
    
    
    AbstractLinearAssembler(AbstractBasisFunction<ELEMENT_DIM> *pBasisFunction,
                            AbstractBasisFunction<ELEMENT_DIM-1> *pSurfaceBasisFunction,
                            int numQuadPoints = 2,
                            double linearSolverRelativeTolerance = 1e-6) :
            AbstractAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(pBasisFunction, pSurfaceBasisFunction, numQuadPoints)
    {
        mpLinearSolver = new SimpleLinearSolver(linearSolverRelativeTolerance);
        mWeAllocatedSolverMemory = true;
        
        this->mpLinearSystem = NULL;
        this->mMatrixIsConstant = false;
        this->mMatrixIsAssembled = false;

        this->mProblemIsLinear = true;
    }
    
    /**
     *  Destructor: ensures that the LinearSystem is thrown away.
     */
    ~AbstractLinearAssembler()
    {
        if (this->mpLinearSystem != NULL)
        {
            delete this->mpLinearSystem;
        }
        
        this->mpLinearSystem=NULL;
        
        if(mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
    }
    
    /**
     *  Manually re-set the linear system solver (which by default 
     *  is a SimpleLinearSolver)
     */
    void SetLinearSolver(AbstractLinearSolver *pLinearSolver)
    {
        if(mWeAllocatedSolverMemory)
        {
            delete mpLinearSolver;
        }
        mpLinearSolver = pLinearSolver;
        
        // make sure new solver knows matrix is constant
        if (this->mMatrixIsConstant)
        {
            SetMatrixIsConstant();
        }
    }
    
    virtual Vec Solve()=0;
    
    
    /**
     * Set the boolean mMatrixIsConstant to true to build the matrix only once. 
     */
    void SetMatrixIsConstant()
    {
        this->mMatrixIsConstant = true;
        this->mpLinearSolver->SetMatrixIsConstant();
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
