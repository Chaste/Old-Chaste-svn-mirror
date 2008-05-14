/*

Copyright (C) University of Oxford, 2008

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
    double mLinearSolverTolerance;
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
                // Static problem, create linear system 
                // The following ensures all the unknowns for a particular node
                // are on the same processor
                DistributedVector::SetProblemSize(this->mpMesh->GetNumNodes());
                Vec template_vec = DistributedVector::CreateVec(PROBLEM_DIM);
                
                this->mpLinearSystem = new LinearSystem(template_vec);
                
                VecDestroy(template_vec);
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
                this->mpLinearSystem->SetAbsoluteTolerance(mLinearSolverTolerance);
            }
            else
            {
                this->mpLinearSystem->SetRelativeTolerance(mLinearSolverTolerance);
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

        //this->mpLinearSystem->DisplayRhs();

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
            mLinearSolverTolerance(1e-6),
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
        mUseLinearSolverAbsoluteTolerance = false;
        mLinearSolverTolerance = relativeTolerance;
    }

    void SetLinearSolverAbsoluteTolerance(double absoluteTolerance)
    {
        assert(this->mpLinearSystem==NULL);
        mUseLinearSolverAbsoluteTolerance = true;
        mLinearSolverTolerance = absoluteTolerance;
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
