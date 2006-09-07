#ifndef ELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_
#define ELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_

#include "BoundaryConditionsContainer.hpp"


/** 
 *  ElasticityBoundaryConditionsContainer
 * 
 *  A boundary conditions container with some added functionality for use in 
 *  linear and finite elasticity problems.
 */
template <int DIM>
class ElasticityBoundaryConditionsContainer : public BoundaryConditionsContainer<DIM,DIM,DIM>
{

public :
    /** 
     *  Fix a node (ie specify zero displacement on that node). The node
     *  must be a boundary node.
     */
    void FixNode(const Node<DIM>* pBoundaryNode)
    {
        SetDisplacement(pBoundaryNode, zero_vector<double>(DIM));
    }

    /** 
     *  Prescribe a given displacement for node (ie specify zero displacement on that node). 
     *  The node must be a boundary node.
     */
    void SetDisplacement(const Node<DIM>* pNode, c_vector<double, DIM> displacement)
    {
        for(unsigned space_index=0; space_index<DIM; space_index++)
        {
            ConstBoundaryCondition<DIM>* p_boundary_condition 
              = new ConstBoundaryCondition<DIM>( displacement(space_index) );

            AddDirichletBoundaryCondition(pNode, p_boundary_condition, space_index);
        }
    }
    
    
    ElasticityBoundaryConditionsContainer(int numNodes) 
      : BoundaryConditionsContainer<DIM,DIM,DIM>(numNodes)
    {
    }
};
#endif /*ELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_*/
