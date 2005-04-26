// AbstractMeshReadWrite.hpp

/**
 * Abstract mesh reader/writer class.
 * "AbstractMeshReadWrite" contains <strong>all</strong>
 * data necessary to either read or write a mesh.
 */
 
 #ifndef _ABSTRACTMESHREADWRITE_HPP_
#define _ABSTRACTMESHREADWRITE_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Exception.hpp"

class AbstractMeshReadWrite 
{
	protected:
		int mNumElements; /**< Is the number of elements in the mesh*/
		int mNumNodes; /**< Is the number of nodes in the mesh*/
		int mNumFaces; /**< Is the number of faces (edges in 2-d) in the mesh*/
		int mNumBoundaryFaces; /**< Is the number of boundary faces (edges in 2-d) in the mesh*/
		int mDimension; /**< Is the dimension the mesh*/
		
		int mNumNodeAttributes; /**< Is the number of attributes stored at each node */
		int mMaxNodeBdyMarker; /**< Is the maximum node boundary marker */
		int mNumElementNodes; /** Is the number of nodes per element*/
		int mNumElementAttributes; /**< Is the number of attributes stored for each element */
		int mMaxFaceBdyMarker; /**< Is the maximum face (or edge) boundary marker */
		
		std::string mPathBaseName; /**< Path to the directory where the input files are stored */
		std::vector<std::string> mNodeRawData;  /**< Contents of node input file with comments removed */
		std::vector<std::string> mElementRawData;  /**< Contents of element input file with comments removed */
		std::vector<std::string> mFaceRawData;  /**< Contents of face (or edge) input file with comments removed */
	
		std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
		std::vector< std::vector<int> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
		std::vector< std::vector<int> > mFaceData; /**< Is an array of the nodes in each face ((i,j)th entry is the jth node of face i) */
		std::vector< std::vector<int> > mBoundaryFaceData; /**< Is an array of the nodes on each boundary face ((i,j)th entry is the jth node of face i) */		
		
		std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
		std::vector< std::vector<int> >::iterator mpElementIterator; /**< Is an iterator for the element data */
		std::vector< std::vector<int> >::iterator mpFaceIterator; /**< Is an iterator for the face data */
		std::vector< std::vector<int> >::iterator mpBoundaryFaceIterator; /**< Is an iterator for the boundary face data */		
	
		bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */
	public:	
		AbstractMeshReadWrite() /**< Constructor */
		{
			mNumElements = 0;
			mNumNodes = 0;
			mNumFaces = 0;
			mNumBoundaryFaces = 0;			
			mDimension = 0;
			
			mNumNodeAttributes = 0;
			mMaxNodeBdyMarker = 0;
			mNumElementNodes = 0;
			mNumElementAttributes = 0;
			mMaxFaceBdyMarker = 0;
			
			// We have initialized all numeric variables to zero
			
			mPathBaseName = "";
			mIndexFromZero = false; // Initially assume that nodes are not numbered from zero
		}
};
#endif //_ABSTRACTMESHREADWRITE_HPP_
