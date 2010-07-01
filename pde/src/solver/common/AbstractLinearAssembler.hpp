/*

Copyright (C) University of Oxford, 2005-2010

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

    /** Whether the matrix is unchanged each time Solve() is called */
    bool mMatrixIsConstant;

protected:

    /**
     * Hack for dynamic mixin.
     *
     * @param matrixIsConstant defaults to true
     */
    void SetMatrixIsConst(bool matrixIsConstant=true);

    /**
     * Apply Dirichlet boundary conditions to the linear system.
     *
     * @param unusedVector
     * @param applyToMatrix
     */
    void ApplyDirichletConditions(Vec unusedVector, bool applyToMatrix);

    /**
     * Create the linear system object if it hasn't been already.
     * Can use an initial solution as PETSc template, or base it on the mesh size.
     *
     * @param initialSolution an initial guess
     */
    virtual void InitialiseForSolve(Vec initialSolution);

    /**
     * Whether grad_u should be calculated
     */
    bool ProblemIsNonlinear();

    /**
     * Solve a static pde, or a dynamic pde for 1 timestep.
     *
     * The mesh, pde and boundary conditions container must be set first.
     *
     * @param currentSolutionOrGuess  either the current solution (dynamic problem) or
     *     initial guess (static problem) (defaults to NULL)
     * @param currentTime  for a dynamic problem, the current time (defaults to 0.0)
     * @param assembleMatrix  whether to assemble the matrix (defaults to true)
     * @return the solution vector
     */
    virtual Vec StaticSolve(Vec currentSolutionOrGuess=NULL,
                            double currentTime=0.0,
                            bool assembleMatrix=true);

public:

    /**
     * Constructors just call the base class versions.
     *
     * @param numQuadPoints number of quadrature points (defaults to 2)
     */
    AbstractLinearAssembler(unsigned numQuadPoints=2);

    /**
     * Destructor: ensures that the linear solver is thrown away.
     */
    ~AbstractLinearAssembler();

    /**
     * Solve the static pde.
     *
     * The mesh, pde and boundary conditions container must be set before Solve()
     * is called.
     *
     * @param currentSolutionOrGuess either the current solution or initial guess (defaults to NULL)
     * @param currentTime the current time (defaults to 0.0)
     */
    virtual Vec Solve(Vec currentSolutionOrGuess=NULL, double currentTime=0.0);

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


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
void AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::SetMatrixIsConst(bool matrixIsConstant)
{
     mMatrixIsConstant = matrixIsConstant;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
void AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::ApplyDirichletConditions(Vec unusedVector, bool applyToMatrix)
{
    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), applyToMatrix);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
void AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::InitialiseForSolve(Vec initialSolution)
{
    if (this->mpLinearSystem == NULL)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        if (initialSolution == NULL)
        {
            // Static problem, create linear system
            // The following ensures all the unknowns for a particular node
            // are on the same processor
            Vec template_vec = this->mpMesh->GetDistributedVectorFactory()->CreateVec(PROBLEM_DIM);
            
            ///\todo #1216 Choose the row preallocation size more sensibly than just setting it to 54 below.
            this->mpLinearSystem = new LinearSystem(template_vec, 54);

            VecDestroy(template_vec);
        }
        else
        {
            // Use the currrent solution (ie the initial solution)
            // as the template in the alternative constructor of
            // LinearSystem. This is to avoid problems with VecScatter.
            ///\todo #1216 Choose the row preallocation size more sensibly than just setting it to 54 below.
            this->mpLinearSystem = new LinearSystem(initialSolution, 54);
        }

        this->mpLinearSystem->SetMatrixIsConstant(mMatrixIsConstant);
        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
bool AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::ProblemIsNonlinear()
{
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
Vec AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::StaticSolve(Vec currentSolutionOrGuess,
                        double currentTime,
                        bool assembleMatrix)
{
    this->AssembleSystem(true, assembleMatrix, currentSolutionOrGuess, currentTime);
    return this->mpLinearSystem->Solve(currentSolutionOrGuess);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::AbstractLinearAssembler(unsigned numQuadPoints)
    : AbstractStaticAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM, NON_HEART, CONCRETE>(numQuadPoints),
      mMatrixIsConstant(true)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::~AbstractLinearAssembler()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool NON_HEART, class CONCRETE>
Vec AbstractLinearAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, NON_HEART, CONCRETE>::Solve(Vec currentSolutionOrGuess, double currentTime)
{
    // have to comment this out because the flagged mesh assembler uses a different
    // bcc and so this is null.
    //assert(this->mpBoundaryConditions!=NULL);

    this->PrepareForSolve();
    this->InitialiseForSolve(currentSolutionOrGuess);
    return this->StaticSolve(currentSolutionOrGuess, currentTime);
}

#endif //_ABSTRACTLINEARASSEMBLER_HPP_
