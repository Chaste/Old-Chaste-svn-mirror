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
