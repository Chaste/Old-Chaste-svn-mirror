#ifndef _REFINEDTETRAHEDRALMESH_CPP_
#define _REFINEDTETRAHEDRALMESH_CPP_

#include "ConformingTetrahedralMesh.cpp"
#include "NodeMap.hpp"

template<int ELEMENT_DIM, int SPACE_DIM>
class RefinedTetrahedralMesh : 
                  public ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{
private:
    ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *mpFineMesh;
    NodeMap *mpNodeMap;
    std::vector <std::set <Element <ELEMENT_DIM,SPACE_DIM>* > > mCoarseFineElementsMap;

public:
    RefinedTetrahedralMesh() : ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>()
    {
        mpNodeMap = NULL; 
    }
    
    ~RefinedTetrahedralMesh()
    {
    	delete mpNodeMap;
    }
    void SetFineMesh(ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pFineMesh)
    {
        mpFineMesh=pFineMesh;
        
        assert(this->GetNumNodes() > 0);
        mpNodeMap = new NodeMap(this->GetNumNodes());
        
        std::vector<Node<SPACE_DIM> * > coarse_nodes;
        for (unsigned i=0;i<this->GetNumNodes();i++)
        {
            coarse_nodes.push_back(this->GetNode(i));
        }
        std::sort(coarse_nodes.begin(), coarse_nodes.end(), CompareNodesLex());
        
        std::vector<Node<SPACE_DIM> * > fine_nodes;
        for (unsigned i=0;i<mpFineMesh->GetNumNodes();i++)
        {
            fine_nodes.push_back(mpFineMesh->GetNode(i));
        }
        std::sort(fine_nodes.begin(), fine_nodes.end(), CompareNodesLex());
        
       
        unsigned fine_mesh_index=0;
        //.. update node map
        for (unsigned coarse_mesh_index=0;coarse_mesh_index<this->GetNumNodes();coarse_mesh_index++)
        {
            while (!EqualNodes(fine_nodes[fine_mesh_index], coarse_nodes[coarse_mesh_index]) )
            {
                fine_mesh_index++;
                if (fine_mesh_index==fine_nodes.size())
                {
                    EXCEPTION("Coarse mesh node doesn't have a partner in the fine mesh.");
                }
            }
            //Same node, set map
            mpNodeMap->SetNewIndex(coarse_mesh_index, fine_mesh_index);
        } 
            
        mCoarseFineElementsMap.resize(this->GetNumElements());
        typename ConformingTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator i_fine_element;
        for (i_fine_element=mpFineMesh->GetElementIteratorBegin();
             i_fine_element< mpFineMesh->GetElementIteratorEnd();
             i_fine_element++)
        {
            bool coarse_elements_found=false;
            for (unsigned node_local_index = 0; node_local_index <= ELEMENT_DIM; node_local_index++)
            {
                Point<SPACE_DIM> vertex = (*i_fine_element)->GetNode(node_local_index)->GetPoint();
                try
                {
                    unsigned coarse_element_index = this->GetContainingElementIndex(vertex, true);
                    coarse_elements_found=true;
                    mCoarseFineElementsMap[coarse_element_index].insert(*i_fine_element);
                }
                catch (Exception &e)
                {
                    // vertex must coincide with coarse node
                } 
            }
            if (!coarse_elements_found)
            {
                // fine element must coincide with coarse element
                Point<SPACE_DIM> centroid = Point<SPACE_DIM>((*i_fine_element)->CalculateCentroid());
                try
                {
                    unsigned coarse_element_index = this->GetContainingElementIndex(centroid, true);
                    mCoarseFineElementsMap[coarse_element_index].insert(*i_fine_element);
                }
                catch (Exception &e)
                {
                    EXCEPTION("Fine mesh contains an element which does not overlap any coarse mesh element");
                }
            }    
        }
        
        
    }
    
    std::set< Element<ELEMENT_DIM, SPACE_DIM>* > GetFineElementsForCoarseElementIndex(unsigned coarse_element_index)
    {
        return mCoarseFineElementsMap[coarse_element_index];
    }
    
    NodeMap& rGetCoarseFineNodeMap()
    {
        return *mpNodeMap;
    }
    
    bool EqualNodes (const Node<SPACE_DIM>* pNode1, const Node<SPACE_DIM>* pNode2)
    {
       return (norm_2(pNode1->rGetLocation() - pNode2->rGetLocation()) < DBL_EPSILON*10);
       
    
        
    }
    class CompareNodesLex : public std::binary_function<Node<SPACE_DIM>*, Node<SPACE_DIM> *, bool>
    {
    public:
        bool operator () (const Node<SPACE_DIM>* pNode1, const Node<SPACE_DIM>* pNode2)
        {
            //Test if node1 is strictly less than node2.
            //Lexigraphical ordering tests the highest dimension first.
            unsigned dimension=SPACE_DIM;
            do 
            {
                dimension--;
                
                if (pNode1->rGetLocation()[dimension] < pNode2->rGetLocation()[dimension])
                {
                    return true;
                }
                if (pNode1->rGetLocation()[dimension] > pNode2->rGetLocation()[dimension])
                {
                    return false;
                }
            
            
            }
            while (dimension>0);
            
            //Otherwise the nodes are colocated  
            return false;
         }
    };
};




#endif // _REFINEDTETRAHEDRALMESH_CPP_
