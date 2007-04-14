#ifndef PARALLELVECTOR_HPP_
#define PARALLELVECTOR_HPP_

#include <vector>
#include <petscvec.h>
#include <iostream>

/***
 * Gives access to the local portion of a PETSc vector via an iterator.
 */

class ParallelVector
{
private:
    static unsigned mLo;
    static unsigned mHi;
    Vec mVec;
    double *mpVec;
public:
    /***
     * Constructor
     * @param vec PETSc vector of with this class shall be a portion.
     */
    ParallelVector(Vec vec);
    
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
            Local++;
            Global++;
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
    };

    /**
     * @return iterator pointing to the first element in the portion
     */
    static Iterator Begin();
    
    /**
     * @return iterator pointing to one past the last element in the portion
     */
    static Iterator End();
    
    double& operator[](Iterator index);

};



#endif /*PARALLELVECTOR_HPP_*/
