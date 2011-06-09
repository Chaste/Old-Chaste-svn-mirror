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

#include "PetscMatTools.hpp"
#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

void PetscMatTools::SetElement(Mat matrix, PetscInt row, PetscInt col, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    if (row >= lo && row < hi)
    {
        MatSetValue(matrix, row, col, value, INSERT_VALUES);
    }
}

void PetscMatTools::AddToElement(Mat matrix, PetscInt row, PetscInt col, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    if (row >= lo && row < hi)
    {
        MatSetValue(matrix, row, col, value, ADD_VALUES);
    }
}

void PetscMatTools::AssembleFinal(Mat matrix)
{
    MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
}

void PetscMatTools::AssembleIntermediate(Mat matrix)
{
    MatAssemblyBegin(matrix, MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(matrix, MAT_FLUSH_ASSEMBLY);
}

void PetscMatTools::Display(Mat matrix)
{
    MatView(matrix,PETSC_VIEWER_STDOUT_WORLD);
}

void PetscMatTools::SetRow(Mat matrix, PetscInt row, double value)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    if (row >= lo && row < hi)
    {
        PetscInt rows, cols;
        MatGetSize(matrix, &rows, &cols);
        for (PetscInt i=0; i<cols; i++)
        {
            SetElement(matrix, row, i, value);
        }
    }
}

void PetscMatTools::ZeroRowsWithValueOnDiagonal(Mat matrix, std::vector<unsigned>& rRows, double diagonalValue)
{
    AssembleFinal(matrix);

    /*
     * Important! Petsc by default will destroy the sparsity structure for this row and deallocate memory
     * when the row is zeroed, and if there is a next timestep, the memory will have to reallocated when
     * assembly to done again. This can kill performance. The following makes sure the zeroed rows are kept.
     */
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
 #if (PETSC_VERSION_MINOR == 0)
    MatSetOption(matrix, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
 #else
    MatSetOption(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
 #endif
#else
    MatSetOption(matrix, MAT_KEEP_ZEROED_ROWS);
#endif

    PetscInt* rows = new PetscInt[rRows.size()];
    for (unsigned i=0; i<rRows.size(); i++)
    {
        rows[i] = rRows[i];
    }
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    IS is;
    ISCreateGeneral(PETSC_COMM_WORLD, rRows.size(), rows, &is);
    MatZeroRows(matrix, is, &diagonalValue);
    ISDestroy(is);
    /*
     *
[2]PETSC ERROR: MatMissingDiagonal_SeqAIJ() line 1011 in src/mat/impls/aij/seq/aij.c
[2]PETSC ERROR: Petsc has generated inconsistent data!
[2]PETSC ERROR: Matrix is missing diagonal number 15!
[2]PETSC ERROR: MatILUFactorSymbolic_SeqAIJ() line 906 in src/mat/impls/aij/seq/aijfact.c
     *
     *
    // There appears to be a problem with MatZeroRows not setting diagonals correctly
    // While we are supporting PETSc 2.2, we have to do this the slow way

    AssembleFinal(matrix);
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);
    PetscInt size=GetSize(matrix);
    ///assert(rRows.size() == 1);
    for (unsigned index=0; index<rRows.size(); index++)
    {
        PetscInt row = rRows[index];
        if (row >= lo && row < hi)
        {
            std::vector<unsigned> non_zero_cols;
            // This row is local, so zero it.
            for (PetscInt column = 0; column < size; column++)
            {
                if (GetElement(matrix, row, column) != 0.0)
                {
                    non_zero_cols.push_back(column);
                }
            }
            // Separate "gets" from "sets" so that we don't have to keep going into "assembled" mode
            for (unsigned i=0; i<non_zero_cols.size();i++)
            {
                SetElement(matrix, row, non_zero_cols[i], 0.0);
            }
            // Set the diagonal
            SetElement(matrix, row, row, diagonalValue);
        }
        // Everyone communicate after row is finished
        AssembleFinal(matrix);
    }
    */
#else
    MatZeroRows(matrix, rRows.size(), rows, diagonalValue);
#endif
    delete [] rows;
}


void PetscMatTools::ZeroRowsAndColumnsWithValueOnDiagonal(Mat matrix, std::vector<unsigned>& rRowColIndices, double diagonalValue)
{
    AssembleFinal(matrix);

    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);
    std::vector<unsigned>* p_nonzero_rows_per_column = new std::vector<unsigned>[rRowColIndices.size()];

    /*
     * For each column: collect all the row indices corresponding to a non-zero entry.
     * We do all the columns at once, before doing the zeroing, as otherwise a
     * MatAssemblyBegin() & MatAssemblyEnd() would have to be called after every
     * MatSetValues and before the below GetMatrixElement().
     */
    for (unsigned index=0; index<rRowColIndices.size(); index++)
    {
        unsigned column = rRowColIndices[index];

        // Determine which rows in this column are non-zero (and therefore need to be zeroed)
        for (PetscInt row = lo; row < hi; row++)
        {
            if (GetElement(matrix, row, column) != 0.0)
            {
                p_nonzero_rows_per_column[index].push_back(row);
            }
        }
    }

    // Now zero each column in turn
    for (unsigned index=0; index<rRowColIndices.size(); index++)
    {
        // set those rows to be zero by calling MatSetValues
        unsigned size = p_nonzero_rows_per_column[index].size();
        PetscInt* rows = new PetscInt[size];
        PetscInt cols[1];
        double* zeros = new double[size];

        cols[0] = rRowColIndices[index];

        for (unsigned i=0; i<size; i++)
        {
            rows[i] = p_nonzero_rows_per_column[index][i];
            zeros[i] = 0.0;
        }

        MatSetValues(matrix, size, rows, 1, cols, zeros, INSERT_VALUES);
        delete [] rows;
        delete [] zeros;
    }
    delete[] p_nonzero_rows_per_column;

    // Now zero the rows and add the diagonal entries
    ZeroRowsWithValueOnDiagonal(matrix, rRowColIndices, diagonalValue);
}

void PetscMatTools::ZeroColumn(Mat matrix, PetscInt col)
{
    AssembleFinal(matrix);

    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    // Determine which rows in this column are non-zero (and therefore need to be zeroed)
    std::vector<unsigned> nonzero_rows;
    for (PetscInt row = lo; row < hi; row++)
    {
        if (GetElement(matrix, row, col) != 0.0)
        {
            nonzero_rows.push_back(row);
        }
    }

    // Set those rows to be zero by calling MatSetValues
    unsigned size = nonzero_rows.size();
    PetscInt* rows = new PetscInt[size];
    PetscInt cols[1];
    double* zeros = new double[size];

    cols[0] = col;

    for (unsigned i=0; i<size; i++)
    {
        rows[i] = nonzero_rows[i];
        zeros[i] = 0.0;
    }

    MatSetValues(matrix, size, rows, 1, cols, zeros, INSERT_VALUES);
    delete [] rows;
    delete [] zeros;
}

void PetscMatTools::Zero(Mat matrix)
{
    MatZeroEntries(matrix);
}

unsigned PetscMatTools::GetSize(Mat matrix)
{
    PetscInt rows, cols;

    MatGetSize(matrix, &rows, &cols);
    assert(rows == cols);
    return (unsigned) rows;
}

void PetscMatTools::GetOwnershipRange(Mat matrix, PetscInt& lo, PetscInt& hi)
{
    MatGetOwnershipRange(matrix, &lo, &hi);
}

double PetscMatTools::GetElement(Mat matrix, PetscInt row, PetscInt col)
{
    PetscInt lo, hi;
    GetOwnershipRange(matrix, lo, hi);

    assert(lo <= row && row < hi);
    PetscInt row_as_array[1];
    row_as_array[0] = row;
    PetscInt col_as_array[1];
    col_as_array[0] = col;

    double ret_array[1];

    MatGetValues(matrix, 1, row_as_array, 1, col_as_array, ret_array);

    return ret_array[0];
}

void PetscMatTools::SetOption(Mat matrix, MatOption option)
{
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
    MatSetOption(matrix, option, PETSC_TRUE);
#else
    MatSetOption(matrix, option);
#endif
}

Vec PetscMatTools::GetMatrixRowDistributed(Mat matrix, unsigned rowIndex)
{
    /*
     * We need to make sure that lhs_ith_row doesn't ignore off processor entries when assembling,
     * otherwise the VecSetValues call a few lines below will not work as expected.
     */

    PetscInt lo, hi;
    PetscMatTools::GetOwnershipRange(matrix, lo, hi);
    unsigned size = PetscMatTools::GetSize(matrix);

    Vec mat_ith_row = PetscTools::CreateVec(size, hi-lo, false);

    PetscInt num_entries;
    const PetscInt* column_indices;
    const PetscScalar* values;

    bool am_row_owner = (PetscInt)rowIndex >= lo && (PetscInt)rowIndex < hi;

    /*
     * Am I the owner of the row? If so get the non-zero entries and add them lhs_ith_row.
     * In parallel, VecAssembly{Begin,End} will send values to the rest of processors.
     */
    if (am_row_owner)
    {
        MatGetRow(matrix, rowIndex, &num_entries, &column_indices, &values);
        VecSetValues(mat_ith_row, num_entries, column_indices, values, INSERT_VALUES);
    }

    VecAssemblyBegin(mat_ith_row);
    VecAssemblyEnd(mat_ith_row);

    if (am_row_owner)
    {
        MatRestoreRow(matrix, rowIndex, &num_entries, &column_indices, &values);
    }

    return mat_ith_row;
}
