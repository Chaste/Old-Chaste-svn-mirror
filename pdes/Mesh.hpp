#ifndef _MESH_HPP_
#define _MESH_HPP_

#include <vector>
#include "Element.hpp"
#include "Node.hpp"
#include "AbstractMeshReader.hpp"

/**
 * This is a test version of the class using an AbstractMeshReader to construct
 * the mesh.
 * 
 * In particular, the AddElement and AddNode methods will
 * be redundant, and should be removed (or possibly made private).
 * 
 * Work still needs to be done with boundary nodes & elements.
 */

template<int ELEMENT_DIM, int SPACE_DIM>
class Mesh 
{
public:
	/** Type of an iterator over mesh elements. */
    typedef typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> >::const_iterator MeshIterator;

private:

    std::vector<Element<ELEMENT_DIM, SPACE_DIM> > mElements;    
    std::vector<Node<SPACE_DIM> > mNodes;
    
    std::vector<Element<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;
    std::vector<Node<SPACE_DIM> *> mBoundaryNodes;
    
    MeshIterator mpConstIter;
    
public:
    
	Mesh();
	
	void ConstructFromMeshReader(AbstractMeshReader &rMeshReader);
    
    //void AddElement(Element<ELEMENT_DIM, SPACE_DIM>& rNewElement);
    //void AddNode(Node<SPACE_DIM>& rNewNode);
    
    //void AddElement(Element<ELEMENT_DIM, SPACE_DIM> newElement,std::vector<int> boundaryElementIndices);
    //void AddNode(Node<SPACE_DIM> newNode ,bool isBoundaryNode);
        
    const Node<SPACE_DIM>& GetNodeAt(long index) const;

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
};

#endif //_MESH_HPP_
