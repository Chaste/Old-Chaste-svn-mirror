#ifndef VECTORPORTION_HPP_
#define VECTORPORTION_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>

/***
 * Gives access to the local portion of a PETSc vector via an iterator.
 */

class VectorPortion
{
private:
    unsigned mLo;
    unsigned mHi;
    Vec mVec;
    double *mpVec;
public:
    /***
     * Constructor
     * @param vec PETSc vector of with this class shall be a portion.
     */
    VectorPortion(Vec vec);
    
    /**
     * Store elements that have been written to via the iterator
     * back into the PETSc vector. Call after you have finished writing.
     */
    void Restore();

    /**
     * Iterator class allows one to iterator of then elements in the vector portion.
     */
    class Iterator
    {
    public:
        VectorPortion* mpVectorPortion;
        Iterator(VectorPortion* pVectorPortion): mpVectorPortion(pVectorPortion)
        {
        }

        bool operator==(const Iterator& other)
        {
           return(Global == other.Global);
        }

        bool operator!=(const Iterator& other)
        {
           return(Global != other.Global);
        }
        
        Iterator& operator++()
        {
           if (Local != mpVectorPortion->mHi)
           {
              Local++;
              Global++;
           }
           return(*this);
        }

        /***
         * The local index of the current element.
         */
        unsigned Local;
        
        /**
         * The global index of the current element.
         */
        unsigned Global;
        
        /**
         * Access the element of the vector by dereferencing the iterator.
         */
        double& operator*()
        {
            return mpVectorPortion->mpVec[Local];
        }
    };

    /**
     * @return iterator pointing to the first element in the portion
     */
    Iterator Begin();
    
    /**
     * @return iterator pointing to one past the last element in the portion
     */
    Iterator End();

};



#endif /*VECTORPORTION_HPP_*/
