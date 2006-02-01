#ifndef _ABSTRACTMESHREADER_HPP_
#define _ABSTRACTMESHREADER_HPP_

/**
 * Abstract mesh reader class. Reads output generated by a mesh generator
 * and converts it to a standard format for use in constructing a finite
 * element mesh structure.
 * 
 * A derived class TrianglesMeshReader exists for reading meshes generated
 * by Triangles (in 2-d) and TetGen (in 3-d).
 * 
 * A derived class MemfemMeshReader reads 3D data from the Tulane University code
 * 
 * A derived class FemlabMeshReader reads 2D data from Femlab or Matlab PDEToolbox
 * 
 */


#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>


class AbstractMeshReader 
{
	protected:
		int mDimension; /**< Is the dimension the mesh*/
		
		int mNumNodeAttributes; /**< Is the number of attributes stored at each node */
		int mMaxNodeBdyMarker; /**< Is the maximum node boundary marker */
		int mNumElementNodes; /** Is the number of nodes per element*/
		int mNumElementAttributes; /**< Is the number of attributes stored for each element */
		int mMaxFaceBdyMarker; /**< Is the maximum face (or edge) boundary marker */
		
		std::vector<std::string> mNodeRawData;  /**< Contents of node input file with comments removed */
		std::vector<std::string> mElementRawData;  /**< Contents of element input file with comments removed */
		std::vector<std::string> mFaceRawData;  /**< Contents of face (or edge) input file with comments removed */
	
		std::vector< std::vector<double> > mNodeData; /**< Is an array of node coordinates ((i,j)th entry is the jth coordinate of node i)*/
		std::vector< std::vector<int> > mElementData; /**< Is an array of the nodes in each element ((i,j)th entry is the jth node of element i) */
		std::vector< std::vector<int> > mFaceData; /**< Is an array of the nodes in each face ((i,j)th entry is the jth node of face i) */
		
		std::vector< std::vector<double> >::iterator mpNodeIterator; /**< Is an iterator for the node data */
		std::vector< std::vector<int> >::iterator mpElementIterator; /**< Is an iterator for the element data */
		std::vector< std::vector<int> >::iterator mpFaceIterator; /**< Is an iterator for the face data */
	
		bool mIndexFromZero; /**< True if input data is numbered from zero, false otherwise */

		std::vector<std::string> GetRawDataFromFile(std::string fileName); /**< Reads an input file fileName, removes comments (indicated by a #) and blank lines */


	public:	
		AbstractMeshReader() /**< Constructor */
		{
			mDimension = 0;
			
			mNumNodeAttributes = 0;
			mMaxNodeBdyMarker = 0;
			mNumElementNodes = 0;
			mNumElementAttributes = 0;
			mMaxFaceBdyMarker = 0;
			
			// We have initialized all numeric variables to zero
			
			mIndexFromZero = false; // Initially assume that nodes are not numbered from zero
		}
        virtual ~AbstractMeshReader(){}

		
		int GetNumElements() const {return mElementData.size();} /**< Returns the number of elements in the mesh */
		int GetNumNodes() const {return mNodeData.size();} /**< Returns the number of nodes in the mesh */
		int GetNumFaces() const {return mFaceData.size();} /**< Returns the number of faces in the mesh (synonym of GetNumEdges()) */
		int GetNumEdges() const {return mFaceData.size();}	/**< Returns the number of edges in the mesh (synonym of GetNumFaces()) */
		int GetDimension() const {return mDimension;} /**< Returns the dimension of the system */
		
		int GetMaxNodeIndex(); /**< Returns the maximum node index */
		int GetMinNodeIndex(); /**< Returns the minimum node index */
		
		std::vector<double> GetNextNode(); /**< Returns a vector of the coordinates of each node in turn */
		void Reset(); /**< Resets pointers to beginning*/
		std::vector<int> GetNextElement(); /**< Returns a vector of the nodes of each element in turn */
		std::vector<int> GetNextEdge(); /**< Returns a vector of the nodes of each edge in turn (synonym of GetNextFace()) */
		std::vector<int> GetNextFace(); /**< Returns a vector of the nodes of each face in turn (synonym of GetNextEdge()) */
};

#endif //_ABSTRACTMESHREADER_HPP_
