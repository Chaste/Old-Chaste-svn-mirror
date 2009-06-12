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


#ifndef DISTRIBUTEDVECTOR_HPP_
#define DISTRIBUTEDVECTOR_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <cassert>
#include "DistributedVectorException.hpp"

class DistributedVectorFactory;

/**
 * Gives access to the local portion of a PETSc vector via an iterator.
 *
 * It also provides two nested classes for accessing vectors with particular
 * memory layouts: striped and chunked.
 */
class DistributedVector
{
private:
    friend class TestDistributedVector;
    // Data global to all vectors.
    
    /** The first entry owned by the current processor. */
    static unsigned mLo;
    
    /** One above the last entry owned by the current processor. */
    static unsigned mHi;
    /** The problem size, i.e. the length of the vector of unknowns. */
    static unsigned mGlobalHi;
    /** Whether we've checked that PETSc is initialised. */
    static bool mPetscStatusKnown;

    // Data local to a single vector.
    /** How many processors own parts of this vector. */
    unsigned mNumChunks;
    /** The underlying PETSc vector. */
    Vec mVec;
    /** The local part of the underlying PETSc vector. */
    double *mpVec;
    
    /** 
     * Pointer to the factory that created this DistributedVector. 
     * Gives access to local and global sizes 
     */
    DistributedVectorFactory* mpFactory;

    /**
     * Double check (in debug code) that PETSc has been initialised properly
     */
    static void CheckForPetsc();

public:

    /**
     * Set the problem size specifying distribution over local processor.
     *
     * @param size
     * @param local
     */
    static void SetProblemSizePerProcessor(unsigned size, PetscInt local);

    /**
     * Set the problem size.
     *
     * @param size
     */
    static void SetProblemSize(unsigned size);

    /**
     * Set the problem with an existing PETSc vector -- must have stride=1.
     *
     * @param vec
     */
    static void SetProblemSize(Vec vec);

    /**
     * Get the global problem size.
     */
    static unsigned GetProblemSize();

    /**
     * Test if the given global index is owned by the current process, i.e. is local to it.
     *
     * @param globalIndex
     */
    static bool IsGlobalIndexLocal(unsigned globalIndex);

    /**
     * Create a PETSc vector of the problem size
     */
    static Vec CreateVec();

    /**
     * Create a striped PETSc vector of size: stride * problem size
     *
     * @param stride
     */
    static Vec CreateVec(unsigned stride);

    /**
     * Constructor.
     * This class represents the portion of a distributed PETSc vector on this process.
     *
     * Note that this class does NOT take over responsibility for destroying the Vec.
     *
     * @param vec PETSc vector of which this class shall be a portion
     * @param pFactory pointer to a DistributedVector factory (defaults to NULL)
     */
    DistributedVector(Vec vec, DistributedVectorFactory* pFactory=NULL);

   /**
    * @param globalIndex
    * @return value of distributed vector at globalIndex
    * Do not use if stride>1.
    * For use in tests.
    * Will throw a DistributedVectorException if the specified element is not on this process.
    */
    double& operator[](unsigned globalIndex) throw (DistributedVectorException);

    /**
     * Store elements that have been written to
     * back into the PETSc vector. Call after you have finished writing.
     * It appears that you do not need to call this if you only read from the vector.
     */
    void Restore();

    /**
     * Iterator class allows one to iterate over the elements of the distributed
     * vector on this process.
     */
    class Iterator
    {
    public:
        unsigned Local;  /**< Current index, local to this process. */
        unsigned Global; /**< Current index, global to the whole PETSc vector. */

        /**
         * Compare two indices for inequality.
         *
         * @param other
         */
        bool operator!=(const Iterator& other);

        /** Increment the iterator to the next index. */
        Iterator& operator++();
    };

    /**
     * Provide access to a particular stripe of a striped vector.
     *
     * A striped vector has multiple types of information encoded within a single
     * vector, with a layout like [x_1, y_1, z_1, x_2, y_2, z_2, ... x_n, y_n, z_n].
     * This class provides easy access to, for example, the x values.
     */
    class Stripe
    {
    public:
        unsigned mStride; /**< How many processors own parts of this vector. */
        unsigned mStripe; /**< The number of this stripe within the vector starting from 0. */
        double *mpVec;    /**< The local part of the underlying PETSc vector. */

       /**
        * Constructor.
        *
        * @param parallelVec striped vector
        * @param stripe number of this stripe within the vector starting from 0
        */
        Stripe(DistributedVector parallelVec, unsigned stripe)
        {
            mStride = parallelVec.mNumChunks;
            mStripe = stripe;
            assert(mStripe < mStride);
            mpVec = parallelVec.mpVec;
        }

       /**
        * Access a particular element of the stripe if on this processor.
        * For use in tests. Will throw a DistributedVectorException if
        * the specified element is not on this process.
        *
        * @param globalIndex index within the stripe
        * @return value of striped vector
        */
        double& operator[](unsigned globalIndex) throw (DistributedVectorException)
        {
            if (mLo<=globalIndex && globalIndex <mHi)
            {
                return mpVec[(globalIndex - mLo)*mStride + mStripe];
            }
            throw DistributedVectorException();
        }

        /**
        * @param index
        * @return value of striped distributed vector pointed to by index.
        */
        double& operator[](Iterator index) throw (DistributedVectorException)
        {
            return mpVec[index.Local*mStride+mStripe];
        }

    };

    /**
     * Provide access to a particular chunk of a chunked vector.
     *
     * A chunked vector has multiple types of information encoded within a single
     * vector, with a layout like [x_1, x_2, ..., x_n, y_1, y_2, ... y_n].
     * This class provides easy access to, for example, the x values.
     */
    class Chunk
    {
    public:
        unsigned mOffset; /**< The start of this chunk within the locally-owned part of the vector. */
        double *mpVec;    /**< The local part of the underlying PETSc vector. */

        /**
         * Constructor.
         *
         * @param parallelVec chunked vector
         * @param chunk number of this chunk within the vector starting from 0
         */
        Chunk(DistributedVector parallelVec, unsigned chunk)
        {
            assert(chunk<parallelVec.mNumChunks);
            mOffset = chunk*(parallelVec.mHi - parallelVec.mLo);
            mpVec = parallelVec.mpVec;
        }

       /**
        * Access a particular element of the chunk if on this processor.
        * For use in tests. Will throw a DistributedVectorException if
        * the specified element is not on this process.
        *
        * @param globalIndex index within the chunk
        * @return value of striped vector
        */
        double& operator[](unsigned globalIndex) throw (DistributedVectorException)
        {
            if (mLo<=globalIndex && globalIndex <mHi)
            {
                //localIndex = globalIndex - mLo
                return mpVec[mOffset + globalIndex - mLo];
            }
            throw DistributedVectorException();
         }

        /**
        * @param index
        * @return value of striped distributed vector pointed to by index.
        */
        double& operator[](Iterator index) throw (DistributedVectorException)
        {
            return mpVec[mOffset + index.Local];
        }

    };

    /**
     * @return iterator pointing to the first element of the distributed
     * vector on this process
     */
    Iterator Begin();

    /**
     * @return iterator pointing to one past the last element of the distributed
     * vector on this process
     */
    Iterator End();

    /**
    * @param index
    * @return value of distributed vector pointed to by index.
    * Do not use if stride>1.
    */
    double& operator[](Iterator index) throw (DistributedVectorException);
};



#endif /*DISTRIBUTEDVECTOR_HPP_*/
