#ifndef REPLICATABLEVECTOR_HPP_
#define REPLICATABLEVECTOR_HPP_

#include <vector>
#include <petscvec.h>

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
    
    void RemovePetscContext()
    {
        if (mToAll != NULL)
        {
            VecScatterDestroy(mToAll);
            mToAll=NULL;
        }
        
        if (mReplicated != NULL)
        {
            VecDestroy(mReplicated);
            mReplicated=NULL;
        }
        
        if (mDistributed != NULL)
        {
            VecDestroy(mDistributed);
            mDistributed=NULL;
        }
    }
    
public:
    /**
     * Default constructor.
     * Note that the vector will need to be resized before it can be used.
     */
    ReplicatableVector()
    {
        mToAll=NULL;
        mReplicated=NULL;
        mDistributed=NULL;
    }
    /**
     * Default destructor.
     * Remove PETSc context.
     */
    ~ReplicatableVector()
    {
        RemovePetscContext();
    }
    
    /**
     * Constructor to make a vector of given size.
     */
    ReplicatableVector(unsigned size)
    {
        mToAll=NULL;
        mReplicated=NULL;
        mDistributed=NULL;
        resize(size);
    }
    
    /**
     * Return the size of the vector.
     */
    unsigned size(void)
    {
        return mData.size();
    }
    
    /**
     * Resize the vector.
     * 
     * @param size  The number of elements to allocate memory for.
     */
    void resize(unsigned size)
    {
        //PETSc stuff will be out of date
        RemovePetscContext();
        mData.resize(size);
    }
    
    /**
     * Access the vector.
     */
    double& operator[](unsigned index)
    {
        return mData[index];
    }
    
    
    /**
     * Replicate this vector over all processes.
     * 
     * Each process knows its local part of the vector.  This method shares that knowledge
     * across all the processes.
     * 
     * @param lo  The start of our ownership range
     * @param hi  One past the end of our ownership range
     */
    void Replicate(unsigned lo, unsigned hi)
    {
        //Copy information into a PetSC vector
        if (mDistributed==NULL)
        {
            VecCreateMPI(PETSC_COMM_WORLD, hi-lo, this->size(), &mDistributed);
        }
        
        double *p_distributed;
        VecGetArray(mDistributed, &p_distributed);
        for (unsigned global_index=lo; global_index<hi; global_index++)
        {
            p_distributed[ (global_index-lo) ]= mData[global_index];
        }
        VecAssemblyBegin(mDistributed);
        VecAssemblyEnd(mDistributed);
        
        //Now do the real replication
        ReplicatePetscVector(mDistributed);
    }
    
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
    void ReplicatePetscVector(Vec vec)
    {
        //If the size has changed then we'll need to make a new context
        int size;
        VecGetSize(vec, &size);
        if ((int) this->size() != size)
        {
            resize(size);
        }
        if (mReplicated == NULL)
        {
            //This creates mReplicated (the scatter context) and mReplicated (to store values)
            VecScatterCreateToAll(vec, &mToAll, &mReplicated);
        }
        
        //Replicate the data
        VecScatterBegin(vec, mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
        VecScatterEnd(vec,   mReplicated, INSERT_VALUES, SCATTER_FORWARD, mToAll);
        
        //Information is now in mReplicated PETSc vector
        //Copy into mData
        double *p_replicated;
        VecGetArray(mReplicated, &p_replicated);
        for (int i=0; i<size; i++)
        {
            mData[i]=p_replicated[i];
        }
    }
    
};

#endif /*REPLICATABLEVECTOR_HPP_*/
