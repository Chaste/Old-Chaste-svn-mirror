#ifndef _CONFORMINGTETRAHEDRALMESH_HPP_
#define _CONFORMINGTETRAHEDRALMESH_HPP_

#include <vector>
#include "AbstractMeshReader.hpp"
#include "Element.hpp"
#include "Node.hpp"

/**
 * This is a test version of the class using an AbstractMeshReader to construct
 * the mesh.
 * 
 * The AddElement and AddNode methods are now redundant, but still used for
 * tests, so I've left them in (for now).
 * 
 * \todo
 * Work still needs to be done with boundary nodes & elements?
 */

template<int ELEMENT_DIM, int SPACE_DIM>
class ConformingTetrahedralMesh 
{
public:
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> >::const_iterator MeshIterator;
    typedef typename std::vector<const Element<ELEMENT_DIM-1, SPACE_DIM>*>::const_iterator BoundaryElementIterator;
    typedef typename std::vector<const Node<SPACE_DIM>*>::const_iterator BoundaryNodeIterator;
private:

	// Note that since these are vectors of objects, not pointers, push_back
	// will copy the objects.
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> > mElements;    
    std::vector<Node<SPACE_DIM> > mNodes;
    
    std::vector<const Element<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;
    std::vector<const Node<SPACE_DIM> *> mBoundaryNodes;
    
    MeshIterator mpConstIter;
    BoundaryElementIterator mpBoundaryElementIter;
    BoundaryNodeIterator mpBoundaryNodeIter;

	// Methods used for testing before we had MeshReader.
    void AddElement(Element<ELEMENT_DIM, SPACE_DIM> newElement);
    void AddNode(Node<SPACE_DIM> newNode);
    /**
     * \todo
     *  Tests should all be changed to use MeshReaders.
     */
    void AddSurfaceElement(const Element<ELEMENT_DIM-1, SPACE_DIM> *pNewElement);
    /**
     * \todo
     * Create our list of boundary nodes from the node list.
     * This should be removed, and tests updated to use MeshReader.
     */
    void GenerateBoundaryNodeList(void)
    {
    	for (int i=0; i<mNodes.size(); i++)
    	{
			if (mNodes[i].IsBoundaryNode())
		    {
		    	mBoundaryNodes.push_back(&mNodes[i]);
		    }
    	}
    }
    friend class TestConformingTetrahedralMesh;
    
public:
    
	ConformingTetrahedralMesh();
    ConformingTetrahedralMesh(long numElements);
    
    void ConstructFromMeshReader(AbstractMeshReader &rMeshReader);
            
    const Node<SPACE_DIM> *GetNodeAt(long index) const;

    long GetNumNodes();
    long GetNumElements();
    long GetNumBoundaryNodes();
    long GetNumBoundaryElements();
    
    MeshIterator GetFirstElement()
    {
         mpConstIter = mElements.begin();
         return mpConstIter;
    }
    MeshIterator GetLastElement()
    {
         mpConstIter = mElements.end();
         return mpConstIter;
    }
    
    BoundaryElementIterator GetFirstBoundaryElement()
    {	
    	mpBoundaryElementIter = mBoundaryElements.begin();
    	return mpBoundaryElementIter;
    }
    BoundaryElementIterator GetLastBoundaryElement()
    {
    	mpBoundaryElementIter = mBoundaryElements.end();
    	return mpBoundaryElementIter;
    }

    BoundaryNodeIterator GetFirstBoundaryNode()
    {
    	mpBoundaryNodeIter = mBoundaryNodes.begin();
    	return mpBoundaryNodeIter;
    }
    BoundaryNodeIterator GetLastBoundaryNode()
    {
    	mpBoundaryNodeIter = mBoundaryNodes.end();
    	return mpBoundaryNodeIter;
    }
    
};

#endif //_CONFORMINGTETRAHEDRALMESH_HPP_
