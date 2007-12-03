#ifndef DISCRETESYSTEMFORCECALCULATOR_HPP_
#define DISCRETESYSTEMFORCECALCULATOR_HPP_

#include "Meineke2001SpringSystem.hpp"


class DiscreteSystemForceCalculator
{
public:
    Meineke2001SpringSystem<2>& mrMeinekeSpringSystem; 
    
    std::vector<double> CalculateFtAndFn(unsigned index, double theta)
    {
        ConformingTetrahedralMesh<2,2>& r_mesh = mrMeinekeSpringSystem.GetTissue().rGetMesh();
        Node<2>* p_node = r_mesh.GetNode(index);
        
        std::set<unsigned> neighbouring_node_indices;
        
        for (Node<2>::ContainingElementIterator it = p_node->ContainingElementsBegin();
             it != p_node->ContainingElementsEnd();
             ++it)
        {
            Element<2,2>* p_element = r_mesh.GetElement(*it);
            for(unsigned i=0; i<p_element->GetNumNodes(); i++)
            {
                unsigned node_index = p_element->GetNodeGlobalIndex(i);
                if(node_index!=index)
                {
                    neighbouring_node_indices.insert(node_index);
                }
            }
        }
        
        double total_tangential_force = 0.0;
        double total_normal_force = 0.0;
        
        for(std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
            iter != neighbouring_node_indices.end();
            iter++)
        {
            double alpha = r_mesh.GetAngleBetweenNodes(index, *iter);

            if(sin(alpha-theta)>0)
            {
                c_vector<double,2> force_between_nodes = mrMeinekeSpringSystem.CalculateForceBetweenNodes(index,*iter);

                c_vector<double,2> unit_vec_betweeen_nodes;
                unit_vec_betweeen_nodes[0] = cos(alpha);
                unit_vec_betweeen_nodes[1] = sin(alpha);

                double plusminus_norm_force = inner_prod(force_between_nodes,unit_vec_betweeen_nodes);         
                total_tangential_force += plusminus_norm_force * cos(alpha-theta);
                total_normal_force += plusminus_norm_force * sin(alpha-theta);
            }
        }

        std::vector<double> ret(2);
        ret[0] = total_tangential_force;
        ret[1] = total_normal_force;
        
        return ret;
    }


public:
    DiscreteSystemForceCalculator(Meineke2001SpringSystem<2>& rMeinekeSpringSystem)
        : mrMeinekeSpringSystem(rMeinekeSpringSystem)
    {
    }
};
#endif /*DISCRETESYSTEMFORCECALCULATOR_HPP_*/
