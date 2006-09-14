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
     *  Prescribe a given displacement for node. The node must be a boundary node.
     */
    void SetDisplacement(const Node<DIM>* pNode, c_vector<double,DIM> displacement)
    {
        for(unsigned space_index=0; space_index<DIM; space_index++)
        {
            ConstBoundaryCondition<DIM>* p_boundary_condition 
              = new ConstBoundaryCondition<DIM>( displacement(space_index) );

            this->AddDirichletBoundaryCondition(pNode, p_boundary_condition, space_index);
        }
    }

    
    /** 
     *  Prescribe zero traction (force per surface area) on a boundary element. The 
     *  precise boundary condition on the equations is then: "sigma_{ij} n_j = 0" where 
     *  sigma_{ij} is the stress, and n the unit normal.
     */
    void SetZeroTraction(const BoundaryElement<DIM-1,DIM>* pElement)
    {
        SetTraction(pElement, zero_vector<double>(DIM));
    }


    /** 
     *  Prescribe a given traction (force per surface area) on a boundary element. The 
     *  precise boundary condition on the equations is then: "sigma_{ij} n_j = t_i" where 
     *  sigma_{ij} is the stress, n the unit normal and t the given traction.
     */
    void SetTraction(const BoundaryElement<DIM-1,DIM>* pElement, c_vector<double,DIM> traction)
    {
        for(unsigned space_index=0; space_index<DIM; space_index++)
        {
            ConstBoundaryCondition<DIM>* p_boundary_condition
              = new ConstBoundaryCondition<DIM>( traction(space_index) );
              
            this->AddNeumannBoundaryCondition(pElement, p_boundary_condition, space_index);
        }
    }
        

    /** 
     *  Constructor just calls base class constructor
     */
    ElasticityBoundaryConditionsContainer(int numNodes) 
      : BoundaryConditionsContainer<DIM,DIM,DIM>(numNodes)
    {
    }
};
#endif /*ELASTICITYBOUNDARYCONDITIONSCONTAINER_HPP_*/
