#ifndef _CONFORMINGTETRAHEDRALMESH_HPP_
#define _CONFORMINGTETRAHEDRALMESH_HPP_

#include <vector>
#include "AbstractMeshReader.hpp"
#include "Element.cpp"
#include "Node.hpp"
#include "AbstractMaterial.hpp"

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
	
	int mNumCornerNodes;
	
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
//    /**
//     * \todo
//     *  Tests should all be changed to use MeshReaders.
//     */
//    void AddSurfaceElement(const Element<ELEMENT_DIM-1, SPACE_DIM> *pNewElement);
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
    
    ~ConformingTetrahedralMesh();
    
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM> &rMeshReader, int orderOfBasisFunctions=1);
    
    void RescaleMeshFromBoundaryNode(Point<1> updatedPoint, int boundaryNodeIndex);
            
    Node<SPACE_DIM> *GetNodeAt(long index);

    long GetNumNodes();
    long GetNumElements();
    long GetNumBoundaryNodes();
    long GetNumBoundaryElements();
    long GetNumCornerNodes();
    long GetNumAllNodes();
    
    /**
     * Return a pointer to the first element in the mesh.
     */
    MeshIterator GetElementIteratorBegin()
    {
         mpConstIter = mElements.begin();
         return mpConstIter;
    }
    /**
     * Return a pointer to *one past* the last element in the mesh
     * (for consistency with STL iterators).
     */
    MeshIterator GetElementIteratorEnd()
    {
         mpConstIter = mElements.end();
         return mpConstIter;
    }
    
    /**
     * Return a pointer to the first boundary element in the mesh.
     */
    
    BoundaryElementIterator GetBoundaryElementIteratorBegin()
    {	
    	mpBoundaryElementIter = mBoundaryElements.begin();
    	return mpBoundaryElementIter;
    }
    /**
     * Return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd()
    {
    	mpBoundaryElementIter = mBoundaryElements.end();
    	return mpBoundaryElementIter;
    }

    /**
     * Return a pointer to the first boundary node in the mesh.
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorBegin()
    {
    	mpBoundaryNodeIter = mBoundaryNodes.begin();
    	return mpBoundaryNodeIter;
    }
    /**
     * Return a pointer to *one past* the last boundary node in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryNodeIterator GetBoundaryNodeIteratorEnd()
    {
    	mpBoundaryNodeIter = mBoundaryNodes.end();
    	return mpBoundaryNodeIter;
    }
    
    void SetMaterialToElement(int index, AbstractMaterial<SPACE_DIM>* pMaterial)
    {
    	mElements[index].SetMaterial(pMaterial);
    }
    
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(int index)
    {
        return &(mElements[index]);
    }
    
    void SetNode(unsigned index, Point<SPACE_DIM> point, bool verify=true);

  
    
};

#endif //_CONFORMINGTETRAHEDRALMESH_HPP_
