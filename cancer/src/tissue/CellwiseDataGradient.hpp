#ifndef CELLWISEDATAGRADIENT_HPP_
#define CELLWISEDATAGRADIENT_HPP_

#include "CellwiseData.hpp"
#include "LinearBasisFunction.cpp"

/**
 *  A class for calculating the gradients of the CellwiseData.
 */
template<unsigned DIM>
class CellwiseDataGradient
{
private:
    /**
     *  The final gradients at the nodes
     */
    std::vector<c_vector<double, DIM> > mGradients;
    
public:
    
    /**
     *  Compute the gradients at the nodes
     *  This is done by averaging the gradients at all the containing (non-ghost)
     *  elements for that node. Note that the gradients are piecewise constant- 
     *  constant in each element
     */ 
    void SetupGradients()
    {
        Tissue<DIM>& r_tissue = CellwiseData<DIM>::Instance()->rGetTissue();
        ConformingTetrahedralMesh<DIM,DIM>& r_mesh = r_tissue.rGetMesh();

        // initialise gradients size        
        unsigned num_nodes = r_tissue.rGetMesh().GetNumNodes();
        mGradients.resize(num_nodes, zero_vector<double>(DIM));
        
        // the constant gradients at each element
        std::vector<c_vector<double, DIM> > gradients_on_elements;
        unsigned num_elements = r_mesh.GetNumElements();
        gradients_on_elements.resize(num_elements, zero_vector<double>(DIM));
        
        // the number of elements containing a given node (excl ghost elements)
        std::vector<unsigned> num_real_elems_for_node(num_nodes, 0);
        
        for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
        {
            Element<DIM,DIM>& r_elem = *(r_mesh.GetElement(elem_index));

            // calculate the basis functions at any point (eg zero) in the element            
            const c_matrix<double, DIM, DIM> *p_inverse_jacobian = r_elem.GetInverseJacobian();
            const ChastePoint<DIM> zero_point; 
            c_matrix<double, DIM, DIM+1> grad_phi 
              = LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, *p_inverse_jacobian);
        
            bool is_ghost_element = false;
        
            for (unsigned node_index=0; node_index<DIM+1; node_index++)
            {
                unsigned node_global_index = r_elem.GetNodeGlobalIndex(node_index);

                // check whether ghost element
                if( r_tissue.rGetGhostNodes()[node_global_index]==true )
                {
                    is_ghost_element = true;
                    break;
                }

                // if no ghost element, get nutrient conc
                TissueCell& r_cell = r_tissue.rGetCellAtNodeIndex(node_global_index);
                double nutrient_concentration = CellwiseData<DIM>::Instance()->GetValue(&r_cell,0);
                
                // interpolate gradient
                for (unsigned i=0; i<DIM; i++)
                {
                    gradients_on_elements[elem_index](i) += nutrient_concentration* grad_phi(i, node_index);  
                }                
            }
            
            // add gradient at element to gradient at node
            if(!is_ghost_element)
            {
                for (unsigned node_index=0; node_index<DIM+1; node_index++)
                {
                    unsigned node_global_index = r_elem.GetNodeGlobalIndex(node_index);
                    mGradients[node_global_index] += gradients_on_elements[elem_index];
                    num_real_elems_for_node[node_global_index]++;
                }
            }
        }
        
        // divide to obtain average gradient
        for (typename Tissue<DIM>::Iterator cell_iter = r_tissue.Begin();
             cell_iter != r_tissue.End();
             ++cell_iter)
        {
            unsigned node_global_index = cell_iter->GetNodeIndex();

            //// if this fails the node is real node which is not in any real element 
            //assert(num_real_elems_for_node[node_global_index]>0); 
            mGradients[node_global_index] /= num_real_elems_for_node[node_global_index];
        }
    }
    
    
    /**
     *  Get the gradient at each node. Not set up for ghost nodes
     */    
    c_vector<double, DIM>& rGetGradient(unsigned nodeIndex)
    {
        assert( !(CellwiseData<DIM>::Instance()->rGetTissue().rGetGhostNodes()[nodeIndex]) );
        return mGradients[nodeIndex];
    }
};


#endif /*CELLWISEDATAGRADIENT_HPP_*/
