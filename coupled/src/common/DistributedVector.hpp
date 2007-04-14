#ifndef DISTRIBUTEDVECTOR_HPP_
#define DISTRIBUTEDVECTOR_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>
#include <assert.h>

/***
 * Gives access to the local portion of a PETSc vector via an iterator.
 */

class DistributedVector
{
private:
    static unsigned mLo;
    static unsigned mHi;
    static unsigned mGlobalHi;
    unsigned mStride;
    Vec mVec;
    double *mpVec;
public:
    /**
     * Set the problem size.
     */
    static void SetProblemSize(unsigned size);
    
    /**
     * Set the problem with an existing PERScc vector -- must have stride=1
     */
    static void SetProblemSize(Vec vec);
    
    /*
     * Create a PETSc vector of the problem size
     */
    static Vec CreateVec();
    
    /*
     * Create a striped PETSc vector of size: stride * problem size
     */
    static Vec CreateVec(unsigned stride);

    /***
     * Constructor
     * This class represents the portion of a distributed PETSc vector on this process.
     * @param vec PETSc vector of which this class shall be a portion.
     */
    DistributedVector(Vec vec);
    
    /**
     * Store elements that have been written to
     * back into the PETSc vector. Call after you have finished writing.
     */
    void Restore();

    /**
     * Iterator class allows one to iterator of then elements of the distributed
     * vector on this process
     */    
    class Iterator
    {
    public:
        unsigned Local;
        unsigned Global;
    
        bool operator==(const Iterator& other);

        bool operator!=(const Iterator& other);
        
        Iterator& operator++();
    };
    
    class Stripe
    {
    public:
        unsigned mStride;
        unsigned mStripe;
        double *mpVec;
       /***
        * Constructor
        * @param parallelVec striped vector
        * @param stripe number of this stripe within the vector starting from 0
        */
        
        Stripe(DistributedVector parallelVec, unsigned stripe)
        {
            mStride=parallelVec.mStride;
            mStripe=stripe;
            assert(mStripe<mStride);
            mpVec=parallelVec.mpVec;
        }
        
        /**
        * @param index
        * @return value of striped distributed vector pointed to by index.
        */           
        double& operator[](Iterator index)
        {
            return mpVec[index.Local*mStride+mStripe];
        }
        
    };

    /**
     * @return iterator pointing to the first element of the distributed
     * vector on this process 
     */
    static Iterator Begin();
    
    /**
     * @return iterator pointing to one past the last element of the distributed
     * vector on this process
     */
    static Iterator End();
 
    /**
    * @param index
    * @return value of distributed vector pointed to by index.
    * Do not use if stride>1.
    */   
    double& operator[](Iterator index);



};



#endif /*DISTRIBUTEDVECTOR_HPP_*/
