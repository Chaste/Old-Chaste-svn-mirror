#ifndef _ABSTRACTMESHREADER_CPP_
#define _ABSTRACTMESHREADER_CPP_
#include "AbstractMeshReader.hpp"
#include "Exception.hpp"

// AbstractMeshReader.cpp


/**
 * Reads an input file fileName, removes comments (indicated by a #) and blank
 * lines and returns a vector of strings. Each string corresponds to one line
 * of the input file.
 *
 *
 */


template<int ELEMENT_DIM, int SPACE_DIM>
std::vector<std::string> AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetRawDataFromFile(std::string fileName)
{
    // Open raw data file
    
    std::vector<std::string> RawDataFromFile;
    std::ifstream dataFile(fileName.c_str());
    
    
    
    /**
     * Checks that input file has been opened correctly. If not throws an 
     * exception that should be caught by the user.
     * 
     */
    
    if (!dataFile.is_open())
    {
        EXCEPTION("Could not open data file "+fileName+" .");
    }
    
    
    
    // Read each line in turn
    
    std::string RawLineFromFile;
    getline(dataFile, RawLineFromFile);
    
    while (dataFile)
    {
    
        // Remove comments
        
        unsigned hashLocation=RawLineFromFile.find('#',0);
        if (hashLocation >= 0)
        {
            RawLineFromFile=RawLineFromFile.substr(0,hashLocation);
        }
        
        
        // Remove blank lines
        
        unsigned notBlankLocation=RawLineFromFile.find_first_not_of(" \t",0);
        if (notBlankLocation >= 0)
        {
            RawDataFromFile.push_back(RawLineFromFile);
        }
        
        
        // Move onto next line
        
        getline(dataFile, RawLineFromFile);
    }
    
    
    dataFile.close(); // Closes the data file
    
    return(RawDataFromFile);
}



/**
 * Returns the maximum node index. Used in testing to check that output nodes
 * are always indexed from zero even if they are input indexed from one.
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
int AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMaxNodeIndex()
{
    //Initialize an interator for the vector of nodes
    std::vector<std::vector<unsigned> >::iterator the_iterator;
    
    int max_node_index = -1; // Must be negative -- can't go unsigned yet
    
    for (the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
    {
        std::vector<unsigned> indices = *the_iterator; // the_iterator points at each line in turn
        
        for (unsigned i = 0; i < ELEMENT_DIM+1; i++)
        {
            if ((int) indices[i] >  max_node_index)
            {
                max_node_index = indices[i];
            }
        }
    }
    
    return max_node_index;
}



/**
 * Returns the minimum node index. Used in testing to check that output nodes
 * are always indexed from zero even if they are input indexed from one.
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
unsigned AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMinNodeIndex()
{
    //Initialize an interator for the vector of nodes
    std::vector<std::vector<unsigned> >::iterator the_iterator;
    
    unsigned min_node_index = 1000000; // A large integer
    
    for (the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
    {
        std::vector<unsigned> indices = *the_iterator; // the_iterator points at each line in turn
        
        for (unsigned i = 0; i < ELEMENT_DIM+1; i++)
        {
            if (indices[i] < min_node_index)
            {
                min_node_index = indices[i];
            }
        }
    }
    
    return min_node_index;
}



/**
 * Returns a vector of the coordinates of each node in turn, starting with
 * node 0 the first time it is called followed by nodes 1, 2, ... , mNumNodes-1.
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
std::vector<double> AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    /**
     * Checks that there are still some nodes left to read. If not throws an
     * exception that must be caught by the user.
     * 
     */
    
    if (mpNodeIterator == mNodeData.end())
    {
        EXCEPTION("All nodes already got");
    }
    
    std::vector<double> next_node = *mpNodeIterator;
    
    mpNodeIterator++;
    
    return next_node;
}



/**
 * Returns a vector of the nodes of each element in turn, starting with
 * element 0 the first time it is called followed by elements 1, 2, ... ,
 * mNumElements-1.
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
std::vector<unsigned> AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    /**
     * Checks that there are still some elements left to read. If not throws an
     * exception that must be caught by the user.
     * 
     */
    
    if (mpElementIterator == mElementData.end())
    {
        EXCEPTION("All elements already got");
    }
    
    std::vector<unsigned> next_element = *mpElementIterator;
    
    mpElementIterator++;
    
    return next_element;
}

/**
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
void AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    mpElementIterator = mElementData.begin();
    mpFaceIterator = mFaceData.begin();
    mpNodeIterator = mNodeData.begin();
}


/**
 * Returns a vector of the nodes of each face in turn, starting with face 0 the
 * first time it is called followed by faces 1, 2, ... , mNumFaces-1.
 *
 * Is a synonum of GetNextEdge(). The two functions can be used interchangeably,
 * i.e. they use the same iterator.
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
std::vector<unsigned> AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFace()
{
    /**
     * Checks that there are still some faces left to read. If not throws an
     * exception that must be caught by the user.
     * 
     */
    
    if (mpFaceIterator == mFaceData.end())
    {
        EXCEPTION("All faces (or edges) already got");
    }
    
    std::vector<unsigned> next_face = *mpFaceIterator;
    
    mpFaceIterator++;
    
    return next_face;
}



/**
 * Returns a vector of the nodes of each edge in turn, starting with edge 0 the
 * first time it is called followed by edges 1, 2, ... , mNumFaces-1.
 *
 * Is a synonym of GetNextFace(). The two functions can be used interchangeably,
 * i.e. they use the same iterator.
 *
 */

template<int ELEMENT_DIM, int SPACE_DIM>
std::vector<unsigned> AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextEdge()
{
    // Call GetNextFace()
    return GetNextFace();
}
#endif //_ABSTRACTMESHREADER_CPP_

