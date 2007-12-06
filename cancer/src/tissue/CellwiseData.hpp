#ifndef CELLWISEDATA_HPP_
#define CELLWISEDATA_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>
#include "Tissue.cpp"


/**
 *  A singleton object for storing data that certain cell cycle models
 *  need to know about, e.g. nutrient concentrations computed via some PDE 
 *  for use in nutrient based cell cycle models.
 */
template<unsigned DIM>
class CellwiseData
{
    friend class TestCellwiseData;
    
private:
    /* the single instance of the singleton object */
    static CellwiseData* mpInstance;

    /* a reference to a Tissue so a cell's node can be found */
    Tissue<DIM>* mpTissue;

    /*< allocated memory for mData object */
    bool mAllocatedMemory;

    /*< number of variables per node to be stored */
    unsigned mNumberOfVariables;

    /*< store of the data */
    std::vector<double> mData;  
      
    std::vector<double> mConstantDataForTesting;
    bool mUseConstantDataForTesting;    
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
std::cout << " Entering serialize method \n" << std::flush;
        archive & mpTissue;
        archive & mAllocatedMemory;
        archive & mNumberOfVariables;
        archive & mData;
        archive & mConstantDataForTesting;
        archive & mUseConstantDataForTesting;
    }
    
    
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
    double GetValue(TissueCell* pCell, unsigned variableNumber=0);
    
    /**
     *  Set the value for particular node and given variable number (defaults
     *  to zero)
     */
    void SetValue(double value, Node<DIM>* pNode, unsigned variableNumber=0);
    
    /**
     *  Set the Tissue. Must be called before GetValue(). This calls 
     *  Tissue.Initialise()
     */
    void SetTissue(Tissue<DIM>& rTissue);
    
    /**
     *  Set the number of variables to be stored per cell. The constructor
     *  assumes 1 variable so only really needs to be called if num_vars > 1
     */
    void SetNumNodesAndVars(unsigned numNodes, unsigned numVars);
    
    /**
     *  Force the data to return given values for all cells (only for testing)
     */
    void SetConstantDataForTesting(std::vector<double> values)
    {        
        mConstantDataForTesting = values;
        mUseConstantDataForTesting = true;
    }
        
    /**
     *  Is the instance in existence and fully set up
     */
    bool IsSetUp();
    
    /**
     *  Reallocate size of mData. Needed because of growth/death. Reallocates
     *  according to the number of nodes in the mesh in the Tissue member variable
     */
    void ReallocateMemory();
    
};

#endif /*CELLWISEDATA_HPP_*/
