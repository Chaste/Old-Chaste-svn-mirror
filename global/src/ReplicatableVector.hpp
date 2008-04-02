/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef REPLICATABLEVECTOR_HPP_
#define REPLICATABLEVECTOR_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>

class ReplicatableVector
{
private:
    /**
     * The wrapped vector.
     */
    std::vector<double> mData;
    
    VecScatter mToAll;   /**< Variable holding information for replicating a PETSc vector*/
    Vec mReplicated;     /**< Vector to hold concentrated copy of replicated vector*/
    Vec mDistributed;    /**< Vector to hold data before replication*/
    
    void RemovePetscContext();
    
public:
    /**
     * Default constructor.
     * Note that the vector will need to be resized before it can be used.
     */
    ReplicatableVector();
    
    /**
     *  Constructor taking in Petsc vector, which is immediately 
     *  replicated into the internal data
     */
    ReplicatableVector(Vec vec);
    /**
     * Constructor to make a vector of given size.
     */
    ReplicatableVector(unsigned size);
    
    
    /**
     * Default destructor.
     * Remove PETSc context.
     */
    ~ReplicatableVector();
    
    
    /**
     * Return the size of the vector.
     */
    unsigned size(void);
    
    /**
     * Resize the vector.
     * 
     * @param size  The number of elements to allocate memory for.
     */
    void resize(unsigned size);
    
    /**
     * Access the vector.
     */
    double& operator[](unsigned index);
    
    
    /**
     * Replicate this vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes.
     * 
     * @param lo  The start of our ownership range
     * @param hi  One past the end of our ownership range
     */
    void Replicate(unsigned lo, unsigned hi);
    
    /**
     * Replicate the given PETSc vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes, storing it in this object.
     * 
     * Our data vector will automatically be resized to fit the whole PETSc vector.
     * 
     * @param vec  The PETSc vector to replicate.
     */
    void ReplicatePetscVector(Vec vec);
    
};

#endif /*REPLICATABLEVECTOR_HPP_*/
