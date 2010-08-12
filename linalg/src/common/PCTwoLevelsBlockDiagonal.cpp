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

#include "PCTwoLevelsBlockDiagonal.hpp"
#include "Exception.hpp"

#include <iostream>

PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonal(KSP& rKspObject, std::vector<PetscInt>& rBathNodes)
{
    PCTwoLevelsBlockDiagonalCreate(rKspObject, rBathNodes);
    PCTwoLevelsBlockDiagonalSetUp();
}

PCTwoLevelsBlockDiagonal::~PCTwoLevelsBlockDiagonal()
{
    MatDestroy(mPCContext.A11_matrix_subblock);
    MatDestroy(mPCContext.A22_B1_matrix_subblock);
    MatDestroy(mPCContext.A22_B2_matrix_subblock);

    PCDestroy(mPCContext.PC_amg_A11);
    PCDestroy(mPCContext.PC_amg_A22_B1);
    PCDestroy(mPCContext.PC_amg_A22_B2);

    VecDestroy(mPCContext.x1_subvector);
    VecDestroy(mPCContext.y1_subvector);

    VecDestroy(mPCContext.x21_subvector);
    VecDestroy(mPCContext.y21_subvector);

    VecDestroy(mPCContext.x22_subvector);
    VecDestroy(mPCContext.y22_subvector);
    
    VecScatterDestroy(mPCContext.A11_scatter_ctx);
    VecScatterDestroy(mPCContext.A22_B1_scatter_ctx);    
    VecScatterDestroy(mPCContext.A22_B2_scatter_ctx);    
}

void PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalCreate(KSP& rKspObject, std::vector<PetscInt>& rBathNodes)
{
    KSPGetPC(rKspObject, &mPetscPCObject);

    Mat system_matrix, dummy;
    MatStructure flag;
    KSPGetOperators(rKspObject, &system_matrix, &dummy, &flag);

    PetscInt num_rows, num_columns;
    MatGetSize(system_matrix, &num_rows, &num_columns);
    assert(num_rows==num_columns);

    PetscInt num_local_rows, num_local_columns;
    MatGetLocalSize(system_matrix, &num_local_rows, &num_local_columns);

    // odd number of rows: impossible in Bidomain.
    // odd number of local rows: impossible if V_m and phi_e for each node are stored in the same processor.
    if ((num_rows%2 != 0) || (num_local_rows%2 != 0))
    {
        TERMINATE("Wrong matrix parallel layout detected in PCLDUFactorisation.");
    }

    // Allocate memory     
    unsigned subvector_num_rows = num_rows/2;    
    unsigned subvector_local_rows = num_local_rows/2;
    
    unsigned subvector_num_rows_tissue = subvector_num_rows - rBathNodes.size();    
    unsigned subvector_local_rows_tissue = subvector_num_rows_tissue; /// \todo: #1082 won't work in parallel

    unsigned subvector_num_rows_bath = rBathNodes.size();    
    unsigned subvector_local_rows_bath = subvector_num_rows_bath; /// \todo: #1082 won't work in parallel
    
    assert(PetscTools::IsSequential());
    
    mPCContext.x1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.x21_subvector = PetscTools::CreateVec(subvector_num_rows_tissue, subvector_local_rows_tissue);
    mPCContext.x22_subvector = PetscTools::CreateVec(subvector_num_rows_bath, subvector_local_rows_bath);
    mPCContext.y1_subvector = PetscTools::CreateVec(subvector_num_rows, subvector_local_rows);
    mPCContext.y21_subvector = PetscTools::CreateVec(subvector_num_rows_tissue, subvector_local_rows_tissue);
    mPCContext.y22_subvector = PetscTools::CreateVec(subvector_num_rows_bath, subvector_local_rows_bath);

    /*
     * Define IS objects that will be used throughout the method.
     */
    IS A11_all_rows;
    ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 2, &A11_all_rows);     
     
    IS A22_all_rows;
    PetscInt A22_size;
    VecGetSize(mPCContext.x1_subvector, &A22_size);
    /// \todo: #1082 assert size(x1) = size(x21) + size(x22)        
    ISCreateStride(PETSC_COMM_WORLD, A22_size, 1, 2, &A22_all_rows);
           
    IS A22_bath_rows;
    PetscInt* phi_e_bath_rows = new PetscInt[rBathNodes.size()];
    for (unsigned index=0; index<rBathNodes.size(); index++)
    {
        phi_e_bath_rows[index] = 2*rBathNodes[index] + 1;
    }
    ISCreateGeneralWithArray(PETSC_COMM_WORLD, rBathNodes.size(), phi_e_bath_rows, &A22_bath_rows);        
    
    IS A22_tissue_rows;
    ISDifference(A22_all_rows, A22_bath_rows, &A22_tissue_rows);            

    // Create scatter contexts
    {
        // Needed by VecScatterCreate in order to find out parallel layout.
        Vec dummy_vec = PetscTools::CreateVec(num_rows, num_local_rows);

        /// \todo: #1082 legacy, no need to use the references
        IS& A11_rows=A11_all_rows;
        IS& A22_B1_rows=A22_tissue_rows;
        IS& A22_B2_rows=A22_bath_rows;        
        
        IS all_vector;    
        ISCreateStride(PETSC_COMM_WORLD, num_rows/2, 0, 1, &all_vector);
        
        IS tissue_vector;    
        ISCreateStride(PETSC_COMM_WORLD, (num_rows/2)-rBathNodes.size(), 0, 1, &tissue_vector);

        IS bath_vector;    
        ISCreateStride(PETSC_COMM_WORLD, rBathNodes.size(), 0, 1, &bath_vector);
        
        VecScatterCreate(dummy_vec, A11_rows, mPCContext.x1_subvector, all_vector, &mPCContext.A11_scatter_ctx);    
        VecScatterCreate(dummy_vec, A22_B1_rows, mPCContext.x21_subvector, tissue_vector, &mPCContext.A22_B1_scatter_ctx);
        VecScatterCreate(dummy_vec, A22_B2_rows, mPCContext.x22_subvector, bath_vector, &mPCContext.A22_B2_scatter_ctx);
        
        ISDestroy(all_vector);
        ISDestroy(tissue_vector);
        ISDestroy(bath_vector);
        
        VecDestroy(dummy_vec);        
    }
            
    // Get matrix sublock A11        
    {
        // Work out local row range for subblock A11 (same as x1 or y1) 
        PetscInt low, high, global_size;
        VecGetOwnershipRange(mPCContext.x1_subvector, &low, &high);
        VecGetSize(mPCContext.x1_subvector, &global_size);        
        assert(global_size == num_rows/2);       
        
        IS A11_local_rows;
        IS& A11_columns=A11_all_rows;
        ISCreateStride(PETSC_COMM_WORLD, high-low, 2*low, 2, &A11_local_rows);
        ISCreateStride(PETSC_COMM_WORLD, global_size, 0, 2, &A11_columns);
    
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1)
        MatGetSubMatrix(system_matrix, A11_local_rows, A11_columns,
			MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A11_local_rows, A11_columns, PETSC_DECIDE, 
			MAT_INITIAL_MATRIX, &mPCContext.A11_matrix_subblock);
#endif
    
        ISDestroy(A11_local_rows);
    }

    // Get matrix sublock A22_B1
    {
//        // Work out local row range for subblock A22 (same as x2 or y2) 
//        PetscInt low, high, global_size;
//        VecGetOwnershipRange(mPCContext.x21_subvector, &low, &high);
//        VecGetSize(mPCContext.x21_subvector, &global_size);        
//        assert(global_size == (num_rows/2) - (PetscInt) rBathNodes.size());       
        
        assert(PetscTools::IsSequential());
        IS& A22_B1_local_rows = A22_tissue_rows; // wrong in parallel, need to give local rows
        IS& A22_B1_columns = A22_tissue_rows;
                
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1)
        MatGetSubMatrix(system_matrix, A22_B1_local_rows, A22_B1_columns,
			MAT_INITIAL_MATRIX, &mPCContext.A22_B1_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A22_B1_local_rows, A22_B1_columns, PETSC_DECIDE, 
			MAT_INITIAL_MATRIX, &mPCContext.A22_B1_matrix_subblock);
#endif
    
    }

    // Get matrix sublock A22_B2
    {
//        // Work out local row range for subblock A22 (same as x2 or y2) 
//        PetscInt low, high, global_size;
//        VecGetOwnershipRange(mPCContext.x21_subvector, &low, &high);
//        VecGetSize(mPCContext.x21_subvector, &global_size);        
//        assert(global_size == (num_rows/2) - (PetscInt) rBathNodes.size());       
        
        assert(PetscTools::IsSequential());
        IS& A22_B2_local_rows = A22_bath_rows; // wrong in parallel, need to give local rows
        IS& A22_B2_columns = A22_bath_rows;
                
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1)
        MatGetSubMatrix(system_matrix, A22_B2_local_rows, A22_B2_columns,
            MAT_INITIAL_MATRIX, &mPCContext.A22_B2_matrix_subblock);
#else
        MatGetSubMatrix(system_matrix, A22_B2_local_rows, A22_B2_columns, PETSC_DECIDE, 
            MAT_INITIAL_MATRIX, &mPCContext.A22_B2_matrix_subblock);
#endif
    
    }
    
    ISDestroy(A11_all_rows);
    ISDestroy(A22_all_rows);
    ISDestroy(A22_bath_rows);
    delete[] phi_e_bath_rows;
    ISDestroy(A22_tissue_rows);

    // Register call-back function and its context
    PCSetType(mPetscPCObject, PCSHELL);
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2)
    PCShellSetApply(mPetscPCObject, PCTwoLevelsBlockDiagonalApply, (void*) &mPCContext);
#else
    // Register PC context so it gets passed to PCTwoLevelsBlockDiagonalApply
    PCShellSetContext(mPetscPCObject, &mPCContext);
    // Register call-back function
    PCShellSetApply(mPetscPCObject, PCTwoLevelsBlockDiagonalApply);
#endif
}

void PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalSetUp()
{
    // These options will get read by PCSetFromOptions
    PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
    PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");
    PetscOptionsSetValue("-pc_hypre_type", "boomeramg");

    // Set up amg preconditioner for block A11
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A11));
    PCSetType(mPCContext.PC_amg_A11, PCBJACOBI);
    PCSetOperators(mPCContext.PC_amg_A11, mPCContext.A11_matrix_subblock, mPCContext.A11_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
    PCSetFromOptions(mPCContext.PC_amg_A11);
    PCSetUp(mPCContext.PC_amg_A11);

    // Set up amg preconditioner for block A22_B1
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22_B1));
    PCSetType(mPCContext.PC_amg_A22_B1, PCBJACOBI);
    PCSetOperators(mPCContext.PC_amg_A22_B1, mPCContext.A22_B1_matrix_subblock, mPCContext.A22_B1_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
    PCSetFromOptions(mPCContext.PC_amg_A22_B1);
    PCSetUp(mPCContext.PC_amg_A22_B1);

    // Set up amg preconditioner for block A22_B2
    PCCreate(PETSC_COMM_WORLD, &(mPCContext.PC_amg_A22_B2));
    PCSetType(mPCContext.PC_amg_A22_B2, PCHYPRE);
    PCSetOperators(mPCContext.PC_amg_A22_B2, mPCContext.A22_B2_matrix_subblock, mPCContext.A22_B2_matrix_subblock, DIFFERENT_NONZERO_PATTERN);//   SAME_PRECONDITIONER);
    PCSetFromOptions(mPCContext.PC_amg_A22_B2);
    PCSetUp(mPCContext.PC_amg_A22_B2);
}

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1) //PETSc 3.1
PetscErrorCode PCTwoLevelsBlockDiagonalApply(PC pc_object, Vec x, Vec y)
{
  void* pc_context;

  PCShellGetContext(pc_object, &pc_context);   
#else
PetscErrorCode PCTwoLevelsBlockDiagonalApply(void* pc_context, Vec x, Vec y)
{
#endif

    // Cast the context pointer to PCTwoLevelsBlockDiagonalContext
    PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalContext* block_diag_context = (PCTwoLevelsBlockDiagonal::PCTwoLevelsBlockDiagonalContext*) pc_context;
    assert(block_diag_context!=NULL);

    /*
     * Scatter x = [x1 x21 x22]'
     */
//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A11_scatter_ctx, x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A11_scatter_ctx, x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A11_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x1_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A11_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B1_scatter_ctx, x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A22_B1_scatter_ctx, x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B1_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x21_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B1_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B2_scatter_ctx, x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(block_diag_context->A22_B2_scatter_ctx, x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B2_scatter_ctx);
    VecScatterEnd(x, block_diag_context->x22_subvector, INSERT_VALUES, SCATTER_FORWARD, block_diag_context->A22_B2_scatter_ctx);
#endif

    /*
     *  y1  = ILU(A11)*x1
     *  y21 = ILU(A22)*x21
     *  y22 = AMG(A22)*x22
     */
    PCApply(block_diag_context->PC_amg_A11, block_diag_context->x1_subvector, block_diag_context->y1_subvector);
    PCApply(block_diag_context->PC_amg_A22_B1, block_diag_context->x21_subvector, block_diag_context->y21_subvector);
    PCApply(block_diag_context->PC_amg_A22_B2, block_diag_context->x22_subvector, block_diag_context->y22_subvector);

    /*
     * Gather y = [y1 y21 y22]'
     */
//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A11_scatter_ctx, block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A11_scatter_ctx);
    VecScatterEnd(block_diag_context->y1_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A11_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B1_scatter_ctx, block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A22_B1_scatter_ctx, block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B1_scatter_ctx);
    VecScatterEnd(block_diag_context->y21_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B1_scatter_ctx);
#endif

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(block_diag_context->A22_B2_scatter_ctx, block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
    VecScatterEnd(block_diag_context->A22_B2_scatter_ctx, block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE);
#else
    VecScatterBegin(block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B2_scatter_ctx);
    VecScatterEnd(block_diag_context->y22_subvector, y, INSERT_VALUES, SCATTER_REVERSE, block_diag_context->A22_B2_scatter_ctx);
#endif

    return 0;
}
