/*

Copyright (C) University of Oxford, 2005-2009

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


#ifndef PCBLOCKDIAGONAL_HPP_
#define PCBLOCKDIAGONAL_HPP_

#include <cassert>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include "PetscTools.hpp"

/**
 *  PETSc will return the control to this function everytime it needs to precondition a vector (i.e. y = inv(M)*x)
 * 
 *  This function needs to be declared global, since I (migb) haven't found a way of defining it inside a class and
 *  be able of passing it by reference.
 * 
 *  @param pc_context preconditioner context struct. Stores preconditioner state (i.e. PC, Mat, and Vec objects used)
 *  @param x unpreconditioned residual.
 *  @param y preconditioned residual. y = inv(M)*x
 */
PetscErrorCode PCBlockDiagonalApply(void *pc_context, Vec x, Vec y);

/**
 *  This class defines a PETSc-compliant purpouse-build preconditioner.
 * 
 *  Let A be a matrix arising in the FEM discretisation of the bidomain equations with the following block structure: 
 * 
 *                 A = (A11  B')
 *                     (B   A22)
 * 
 *  By creating an instance of this class, one will define the following preconditioner: 
 *
 *                 inv(M) = inv( (A11   0)   = (inv(A11)        0) 
 *                               (0   A22) )   (0        inv(A22))
 * 
 *  The inverses are approximate with one cycle of AMG. 
 */
class PCBlockDiagonal
{
public:
    /**
     *  This struct defines the state of the preconditioner (initialised data and objects to be reused)
     */   
    typedef struct{
        Mat A11_matrix_subblock; /**< See \todo */
        Mat A22_matrix_subblock; /**< See \todo */        
        PC  PC_amg_A11; /**< See \todo */
        PC  PC_amg_A22; /**< See \todo */
        
    } PCBlockDiagonalContext;

    PCBlockDiagonalContext mPCContext; /**< As above.  See PCShellSetContext().*/
    PC mPetscPCObject;/**< Actual preconditioner */

public:

    /**
     *  Constructor
     * 
     *  @param ksp_object KSP object where we want to install the block diagonal preconditioner.
     */    
    PCBlockDiagonal(KSP& ksp_object);
    
    ~PCBlockDiagonal();

private:
    /**
     *  Creates all the state data required by the preconditioner
     * 
     *  @param ksp_object KSP object where we want to install the block diagonal preconditioner.
     */ 
    void PCBlockDiagonalCreate(KSP& ksp_object);
    
    /**
     *  Setups preconditioner
     */    
    void PCBlockDiagonalSetUp();
};
#endif /*PCBLOCKDIAGONAL_HPP_*/
