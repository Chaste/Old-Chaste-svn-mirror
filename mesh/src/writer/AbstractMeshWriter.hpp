#ifndef _ABSTRACTMESHWRITER_HPP_
#define _ABSTRACTMESHWRITER_HPP_

#include "ConformingTetrahedralMesh.cpp"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractMeshReader.cpp"


template<int ELEMENT_DIM, int SPACE_DIM>
class AbstractMeshWriter
{
protected:
    OutputFileHandler *mpOutputFileHandler; /**< Output file handler */
    std::string mBaseName; /**< Base name for the input files */
    
    std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
    std::vector< std::vector<unsigned> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
    std::vector< std::vector<unsigned> > mBoundaryFaceData; /**< Is an array of the nodes on each boundary face ((i,j)th entry is the jth node of face i) */
    
    std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
    std::vector< std::vector<unsigned> >::iterator mpElementIterator; /**< Is an iterator for the element data */
    std::vector< std::vector<unsigned> >::iterator mpBoundaryFaceIterator; /**< Is an iterator for the boundary face data */
    
    bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
    bool mWriteMetaFile;
public:
    /** Constructor */
    AbstractMeshWriter(const std::string &rDirectory,
                       const std::string &rBaseName)
            : mBaseName(rBaseName)
    {
        mpOutputFileHandler = new OutputFileHandler(rDirectory);
    }
    /** Destructor */
    virtual ~AbstractMeshWriter()
    {
        delete mpOutputFileHandler;
    }
    std::string GetOutputDirectory(void);
    
    void SetNextNode(std::vector<double> nextNode);
    void SetNextElement(std::vector<unsigned> nextElement);
    void SetNextBoundaryFace(std::vector<unsigned> nextFace);
    void SetNextBoundaryEdge(std::vector<unsigned> nextEdge);
    virtual void WriteFiles()=0;
    unsigned GetNumNodes()
    {
        return mNodeData.size();
    }
    unsigned GetNumElements()
    {
        return mElementData.size();
    }
    unsigned GetNumBoundaryFaces()
    {
        return mBoundaryFaceData.size();
    }
    unsigned GetNumBoundaryEdges()
    {
        return mBoundaryFaceData.size();
    }
    void WriteFilesUsingMesh(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);
    void WriteFilesUsingMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);
};

#endif //_ABSTRACTMESHWRITER_HPP_
