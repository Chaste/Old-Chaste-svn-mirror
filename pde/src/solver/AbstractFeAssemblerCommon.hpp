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

#ifndef ABSTRACTFEASSEMBLERCOMMON_HPP_
#define ABSTRACTFEASSEMBLERCOMMON_HPP_

#include "LinearBasisFunction.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "HeartEventHandler.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "PetscTools.hpp"

/**
 * Enumeration for defining how much interpolation (onto quadrature points) is
 * required by the concrete class.
 *
 * CARDIAC: only interpolates the first component of the unknown (ie the voltage)
 * NORMAL: interpolates the position X and all components of the unknown u
 * NONLINEAR: interpolates X, u and grad(u). Also computes the gradient of the
 *   basis functions when assembling vectors.
 */
typedef enum InterpolationLevel_
{
    CARDIAC = 0,
    NORMAL,
    NONLINEAR
} InterpolationLevel;



/**
 *   A common bass class for AbstractFeObjectAssembler (the main abstract assembler class) and AbstractCableFeObjectAssembler.
 *   See AbstractFeObjectAssembler documentation for info on these assembler classes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractFeAssemblerCommon : boost::noncopyable
{
protected:
    /** The vector to be assembled (only used if CAN_ASSEMBLE_VECTOR == true). */
    Vec mVectorToAssemble;

    /** The matrix to be assembled (only used if CAN_ASSEMBLE_MATRIX == true). */
    Mat mMatrixToAssemble;

    /**
     * Whether to assemble the matrix (an assembler may be able to assemble matrices
     * (CAN_ASSEMBLE_MATRIX==true), but may not want to do so each timestep, hence
     * this second boolean.
     */
    bool mAssembleMatrix;

    /** Whether to assemble the vector. */
    bool mAssembleVector;

    /** Whether to zero the given matrix before assembly, or just add to it. */
    bool mZeroMatrixBeforeAssembly;

    /** Whether to zero the given vector before assembly, or just add to it. */
    bool mZeroVectorBeforeAssembly;

    /** Ownership range of the vector/matrix - lowest component owned. */
    PetscInt mOwnershipRangeLo;

    /** Ownership range of the vector/matrix - highest component owned +1. */
    PetscInt mOwnershipRangeHi;

    /**
     * If the matrix or vector will be dependent on a current solution, say,
     * this is where that information is put.
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;

    /**
     * The main assembly method. Protected, should only be called through Assemble(),
     * AssembleMatrix() or AssembleVector() which set mAssembleMatrix, mAssembleVector
     * accordingly. Pure and therefore is implemented in child classes. Will involve looping
     * over elements (which may be volume, surface or cable elements), and computing
     * integrals and adding them to the vector or matrix
     */
    virtual void DoAssemble()=0;


    /**
     * Useful inline function for getting an entry from the current solution vector.
     *
     * @param nodeIndex node index
     * @param indexOfUnknown index of unknown
     */
    virtual double GetCurrentSolutionOrGuessValue(unsigned nodeIndex, unsigned indexOfUnknown)
    {
        return mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*nodeIndex + indexOfUnknown];
    }

    /**
     * The concrete subclass can overload this and IncrementInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     */
    virtual void ResetInterpolatedQuantities()
    {}

    /**
     * The concrete subclass can overload this and ResetInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     *
     * @param phiI
     * @param pNode pointer to a node
     */
    virtual void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
    {}


public:

    /**
     * Constructor.
     */
    AbstractFeAssemblerCommon();

    /**
     * Set the matrix that needs to be assembled. Requires CAN_ASSEMBLE_MATRIX==true.
     *
     * @param rMatToAssemble Reference to the matrix
     * @param zeroMatrixBeforeAssembly Whether to zero the vector before assembling
     *  (otherwise it is just added to)
     */
    void SetMatrixToAssemble(Mat& rMatToAssemble, bool zeroMatrixBeforeAssembly=true);

    /**
     * Set the vector that needs to be assembled. Requires CAN_ASSEMBLE_VECTOR==true.
     *
     * @param rVecToAssemble Reference to the vector
     * @param zeroVectorBeforeAssembly Whether to zero the vector before assembling
     *  (otherwise it is just added to)
     */
    void SetVectorToAssemble(Vec& rVecToAssemble, bool zeroVectorBeforeAssembly);

    /**
     * Set a current solution vector that will be used in AssembleOnElement and can passed
     * up to ComputeMatrixTerm() or ComputeVectorTerm().
     *
     * @param currentSolution Current solution vector.
     */
    void SetCurrentSolution(Vec currentSolution);

    /**
     * Assemble everything that the class can assemble.
     */
    void Assemble()
    {
        mAssembleMatrix = CAN_ASSEMBLE_MATRIX;
        mAssembleVector = CAN_ASSEMBLE_VECTOR;
        DoAssemble();
    }

    /**
     * Assemble the matrix. Requires CAN_ASSEMBLE_MATRIX==true and ComputeMatrixTerm() to be implemented.
     */
    void AssembleMatrix()
    {
        assert(CAN_ASSEMBLE_MATRIX);
        mAssembleMatrix = true;
        mAssembleVector = false;
        DoAssemble();
    }

    /**
     * Assemble the vector. Requires CAN_ASSEMBLE_VECTOR==true and ComputeVectorTerm() to be implemented.
     */
    void AssembleVector()
    {
        assert(CAN_ASSEMBLE_VECTOR);
        mAssembleMatrix = false;
        mAssembleVector = true;
        DoAssemble();
    }

    /**
     * Destructor.
     */
    virtual ~AbstractFeAssemblerCommon()
    {
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AbstractFeAssemblerCommon()
    : mVectorToAssemble(NULL),
      mMatrixToAssemble(NULL),
      mZeroMatrixBeforeAssembly(true),
      mZeroVectorBeforeAssembly(true)
{
    assert(CAN_ASSEMBLE_VECTOR || CAN_ASSEMBLE_MATRIX);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::SetMatrixToAssemble(Mat& rMatToAssemble, bool zeroMatrixBeforeAssembly)
{
    assert(rMatToAssemble);
    MatGetOwnershipRange(rMatToAssemble, &mOwnershipRangeLo, &mOwnershipRangeHi);

    mMatrixToAssemble = rMatToAssemble;
    mZeroMatrixBeforeAssembly = zeroMatrixBeforeAssembly;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::SetVectorToAssemble(Vec& rVecToAssemble, bool zeroVectorBeforeAssembly)
{
    assert(rVecToAssemble);
    VecGetOwnershipRange(rVecToAssemble, &mOwnershipRangeLo, &mOwnershipRangeHi);

    mVectorToAssemble = rVecToAssemble;
    mZeroVectorBeforeAssembly = zeroVectorBeforeAssembly;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::SetCurrentSolution(Vec currentSolution)
{
    assert(currentSolution != NULL);

    // Replicate the current solution and store so can be used in AssembleOnElement
    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolution);
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    // The AssembleOnElement type methods will determine if a current solution or
    // current guess exists by looking at the size of the replicated vector, so
    // check the size is zero if there isn't a current solution.
    assert(mCurrentSolutionOrGuessReplicated.GetSize() > 0);
}





#endif /* ABSTRACTFEASSEMBLERCOMMON_HPP_ */
