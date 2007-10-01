#ifndef CELLWISEDATA_HPP_
#define CELLWISEDATA_HPP_

#include "CancerParameters.hpp"
#include "Crypt.cpp"
#include <vector>

/**
 *  A singleton object for storing data certain cell cycle models
 *  will need to know about. Eg, storing nutrient concentrations
 *  computed via some pde in a place a nutrient based cell 
 *  cycle model can get to it
 */
template<unsigned DIM>
class CellwiseData
{
private:
    /* the single instance of the singleton object */
    static CellwiseData* mpInstance;

    /* a reference to a crypt so a cell's node can be found */
    Crypt<DIM>* mpCrypt;

    /*< allocated memory for mData object */
    bool mAllocatedMemory;

    /*< number of variables per node to be stored */
    unsigned mNumberOfVariables;

    /*< store of the data */
    std::vector<double> mData;
    
protected:
    /**
     *  Protected constuctor. Not to be called, use Instance() instead
     */
    CellwiseData();
    
public:
    /**
     *  Get an instance of the object
     */
    static CellwiseData* Instance();
    
    virtual ~CellwiseData();
    
    /** 
     *  Destroy the current instance. Should be called at the end of a 
     *  simulation.
     */
    static void Destroy();
    
    /**
     *  Get the value for particular cell and given variable number (defaults
     *  to zero)
     */
    double GetValue(MeinekeCryptCell* pCell, unsigned variableNumber=0);
    
    /**
     *  Set the value for particular node and given variable number (defaults
     *  to zero)
     */
    void SetValue(double value, Node<DIM>* pNode, unsigned variableNumber=0);
    
    /**
     *  Set the crypt. Must be called before GetValue(). This calls 
     *  crypt.Initialise()
     */
    void SetCrypt(Crypt<DIM>& rCrypt);
    
    /**
     *  Set the number of variables to be stored per cell. The constructor
     *  assumes 1 variable so only really needs to be called if num_vars > 1
     */
    void SetNumNodesAndVars(unsigned numNodes, unsigned numVars);
    
    /**
     *  Is the instance in existence and fully set up
     */
    bool IsSetUp();
    
    /**
     * Reallocate size of mData. Needed because of growth/death.
     */
    void ReallocateMemory(unsigned numNodes);
};

#endif /*CELLWISEDATA_HPP_*/
