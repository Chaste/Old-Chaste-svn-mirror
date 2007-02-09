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

public:
    
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
        for (int i=0;i<this->GetNumNodes();i++)
        {
            coarse_nodes.push_back(this->GetNode(i));
        }
        std::sort(coarse_nodes.begin(), coarse_nodes.end(), CompareNodesLex());
        
        std::vector<Node<SPACE_DIM> * > fine_nodes;
        for (int i=0;i<mpFineMesh->GetNumNodes();i++)
        {
            fine_nodes.push_back(mpFineMesh->GetNode(i));
        }
        std::sort(fine_nodes.begin(), fine_nodes.end(), CompareNodesLex());
        
       
        int fine_mesh_index=0;
        //.. update node map
        for (int coarse_mesh_index=0;coarse_mesh_index<this->GetNumNodes();coarse_mesh_index++)
        {
            while (!EqualNodes(fine_nodes[fine_mesh_index], coarse_nodes[coarse_mesh_index]) )
            {
                fine_mesh_index++;
                if (fine_mesh_index==(int)fine_nodes.size())
                {
                    EXCEPTION("Coarse mesh node doesn't have a partner in the fine mesh.");
                }
            }
            //Same node, set map
            mpNodeMap->SetNewIndex(coarse_mesh_index, fine_mesh_index);
        } 
            
        
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
            if (pNode1->rGetLocation()[2] < pNode2->rGetLocation()[2])
            {
                return true;
            }
            if (pNode1->rGetLocation()[2] > pNode2->rGetLocation()[2])
            {
                return false;
            }
            if (pNode1->rGetLocation()[1] < pNode2->rGetLocation()[1])
            {
                return true;
            }
            if (pNode1->rGetLocation()[1] > pNode2->rGetLocation()[1])
            {
                return false;
            }
            if (pNode1->rGetLocation()[0] < pNode2->rGetLocation()[0])
            {
                return true;
            }
            return false;
            
        }
    };
};




#endif // _REFINEDTETRAHEDRALMESH_CPP_
