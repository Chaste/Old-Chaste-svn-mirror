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
 * Work still needs to be done with boundary nodes & elements.
 */

template<int ELEMENT_DIM, int SPACE_DIM>
class ConformingTetrahedralMesh 
{
public:
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> >::const_iterator MeshIterator;
    typedef typename std::vector<const Element<ELEMENT_DIM-1, SPACE_DIM>*>::const_iterator BoundaryElementIterator;
private:

	// Note that since these are vectors of objects, not pointers, push_back
	// will copy the objects.
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> > mElements;    
    std::vector<Node<SPACE_DIM> > mNodes;
    
    std::vector<const Element<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;
    std::vector<Node<SPACE_DIM> *> mBoundaryNodes;
    
    MeshIterator mpConstIter;
    BoundaryElementIterator mpBoundaryElementIter;
    
public:
    
	ConformingTetrahedralMesh();
    ConformingTetrahedralMesh(long numElements);
    
    void ConstructFromMeshReader(AbstractMeshReader &rMeshReader);
    
    void AddElement(Element<ELEMENT_DIM, SPACE_DIM> newElement);
    void AddNode(Node<SPACE_DIM> newNode);
    
    // Until we have a more permanent solution...
    void AddSurfaceElement(const Element<ELEMENT_DIM-1, SPACE_DIM> *pNewElement);
    
    //void AddElement(Element<ELEMENT_DIM, SPACE_DIM> newElement,std::vector<int> boundaryElementIndices);
    //void AddNode(Node<SPACE_DIM> newNode ,bool isBoundaryNode);
        
    const Node<SPACE_DIM> *GetNodeAt(long index) const;

    long GetNumNodes();
    long GetNumElements();
    
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
    
};

#endif //_CONFORMINGTETRAHEDRALMESH_HPP_
