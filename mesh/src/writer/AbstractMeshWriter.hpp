#ifndef _ABSTRACTMESHWRITER_HPP_
#define _ABSTRACTMESHWRITER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Exception.hpp"
#include "OutputFileHandler.hpp"

class AbstractMeshWriter 
{
	protected:
        OutputFileHandler *mpOutputFileHandler; /**< Output file handler */
        std::string mBaseName; /**< Base name for the input files */

		unsigned int mDimension; /**< Is the dimension the mesh*/
	
		std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
		std::vector< std::vector<int> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
		std::vector< std::vector<int> > mBoundaryFaceData; /**< Is an array of the nodes on each boundary face ((i,j)th entry is the jth node of face i) */		
		
		std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
		std::vector< std::vector<int> >::iterator mpElementIterator; /**< Is an iterator for the element data */
		std::vector< std::vector<int> >::iterator mpBoundaryFaceIterator; /**< Is an iterator for the boundary face data */		
	
		bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
		bool mWriteMetaFile; 
	public:	
        /** Constructor */
		AbstractMeshWriter(const std::string &rDirectory, 
                           const std::string &rBaseName, 
                           const unsigned int &rDimension)
            : mBaseName(rBaseName),
              mDimension(rDimension)
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
    	void SetNextElement(std::vector<int> nextElement);
    	void SetNextBoundaryFace(std::vector<int> nextFace);
    	void SetNextBoundaryEdge(std::vector<int> nextEdge);
    	virtual void WriteFiles(){};
    	int GetNumNodes(){return mNodeData.size();}
    	int GetNumElements(){return mElementData.size();}
    	int GetNumBoundaryFaces(){return mBoundaryFaceData.size();}
    	int GetNumBoundaryEdges(){return mBoundaryFaceData.size();}
};

#endif //_ABSTRACTMESHWRITER_HPP_
