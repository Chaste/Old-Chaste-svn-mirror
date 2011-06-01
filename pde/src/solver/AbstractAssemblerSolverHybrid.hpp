
/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTASSEMBLERSOLVERHYBRID_HPP_
#define ABSTRACTASSEMBLERSOLVERHYBRID_HPP_

#include "AbstractFeObjectAssembler.hpp"
#include "AbstractLinearPdeSolver.hpp"

/**
 * A class which inherits from AbstractFeObjectAssembler and
 * implements a method SetupGivenLinearSystem(), which sets up
 * the given linear system using the assembler part of this
 * class, which can be called by SetUpLinearSystem() on a
 * concrete solver.
 *
 * See SimpleLinearEllipticSolver for an example.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractAssemblerSolverHybrid
   : public AbstractFeObjectAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, INTERPOLATION_LEVEL>
{
public:

    /**
     * Constructor.
     *
     * @param pMesh pointer to the mesh
     * @param pBoundaryConditions pointer to the boundary conditions. Can be NULL, to allow concrete assembler-solver
     *        to, say, create standard boundary conditions its constructor, and then set it. If so, the concrete solver
     *        must make sure it calls this->SetApplyNeummanBoundaryConditionsToVector(p_bcc);
     * @param numQuadPoints number of quadrature points in each dimension to use per element (defaults to 2)
     */
    AbstractAssemblerSolverHybrid(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                  BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>* pBoundaryConditions,
                                  unsigned numQuadPoints=2)
        : AbstractFeObjectAssembler<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, true, true, INTERPOLATION_LEVEL>(pMesh,numQuadPoints)
    {
        if (pBoundaryConditions)
        {
            this->SetApplyNeummanBoundaryConditionsToVector(pBoundaryConditions);
        }
    }

    /**
     * Destructor.
     */
    virtual ~AbstractAssemblerSolverHybrid()
    {
    }

    /**
     * Implementation of AbstractLinearPdeSolver::SetupLinearSystem, using the assembler that this class
     * also inherits from. Concrete classes inheriting from both this class and
     * AbstractLinearPdeSolver can then have a one-line implementation of
     * AbstractLinearPdeSolver::SetupLinearSystem which calls this method.
     *
     * @param currentSolution The current solution which can be used in setting up
     *  the linear system if needed (NULL if there isn't a current solution)
     * @param computeMatrix Whether to compute the LHS matrix of the linear system
     *  (mainly for dynamic solves)
     * @param pLinearSystem  The linear system to set up.
     */
    void SetupGivenLinearSystem(Vec currentSolution, bool computeMatrix, LinearSystem* pLinearSystem);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractAssemblerSolverHybrid<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, INTERPOLATION_LEVEL>::SetupGivenLinearSystem(Vec currentSolution, bool computeMatrix, LinearSystem* pLinearSystem)
{
    /*
     * The concrete class should have either passed in boundary conditions
     * into this class's constructor, or passed in NULL and later called
     * this->SetApplyNeummanBoundaryConditionsToVector(p_bcc) with some boundary
     * conditions.
     */
    assert(this->mpBoundaryConditions != NULL);

    assert(pLinearSystem->rGetLhsMatrix() != NULL);
    assert(pLinearSystem->rGetRhsVector() != NULL);

    // Call methods on AbstractFeObjectAssembler
    this->SetMatrixToAssemble(pLinearSystem->rGetLhsMatrix());
    this->SetVectorToAssemble(pLinearSystem->rGetRhsVector(), true);

    if (currentSolution != NULL)
    {
        this->SetCurrentSolution(currentSolution);
    }

    if (computeMatrix)
    {
        this->Assemble();
    }
    else
    {
        this->AssembleVector();
    }

    pLinearSystem->AssembleRhsVector();
    pLinearSystem->AssembleIntermediateLhsMatrix();

    this->mpBoundaryConditions->ApplyDirichletToLinearProblem(*pLinearSystem, true);

    pLinearSystem->AssembleRhsVector();
    pLinearSystem->AssembleFinalLhsMatrix();
}

#endif /*ABSTRACTASSEMBLERSOLVERHYBRID_HPP_*/
