#ifndef DISTANCEMAPCALCULATOR_HPP_
#define DISTANCEMAPCALCULATOR_HPP_

#include <vector>
#include <queue>

/**
 *  This class defines  
 * 
 * 
 */

template<unsigned SPACE_DIM>
class DistanceMapCalculator
{
private:
    ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>& mrMesh;
    unsigned mNumNodes;
    
    /**
     *  Computes the euclidean distance of two given points
     * 
     *  @param pointA First point 
     *  @param pointB Second point
     */    
    inline double EuclideanDistanceTwoPoints(const c_vector<double, SPACE_DIM>& pointA, const c_vector<double, SPACE_DIM>& pointB) const
    {
        double dist=0.0;
        
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            dist+=(pointA(dim)-pointB(dim)) * (pointA(dim)-pointB(dim));
        }
        
        return sqrt(dist);
    }

    /**
     *  Given a cartesian distance, computes the associated euclidean distance
     * 
     *  @param cartDistance Cartesian distance of a given point
     */
    inline double CartToEucliDistance(c_vector<double, SPACE_DIM>& cartDistance) const
    {
        double dist=0.0;
        
        for (unsigned dim=0; dim<SPACE_DIM; dim++)
        {
            dist+=cartDistance(dim) * cartDistance(dim);
        }
        
        return sqrt(dist);
    }


public:
    DistanceMapCalculator(ConformingTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh):
        mrMesh(rMesh)
    {
        mNumNodes = mrMesh.GetNumNodes();
    }
    
    /**
     *  Generates a distance map of all the nodes of the mesh to the given surface
     * 
     *  @param rOriginSurface set of node indexes defining the surface
     *  @param rNodeDistances distance map computed
     */
    void ComputeDistanceMap(std::vector<unsigned>& rOriginSurface, std::vector<double>& rNodeDistances)
    {
        if (rNodeDistances.size() != mNumNodes)
        {
            rNodeDistances.reserve(mNumNodes);
        }   
        
        /*
         * Matrix of distances along each dimension (initialised to +inf)
         */
        std::vector< c_vector<double, SPACE_DIM> >  cart_distances(mNumNodes);
        for (unsigned index=0; index<mNumNodes; index++)
        {
            for (unsigned dim=0; dim<SPACE_DIM; dim++)
            {
                cart_distances[index][dim] = DBL_MAX;
            }                     
        } 
        
        /*
         * Queue of nodes to be processed (initialised with the nodes defining the surface)
         */
        std::queue<unsigned> active_nodes;
        for (unsigned index=0; index<rOriginSurface.size(); index++)
        {
            active_nodes.push(rOriginSurface[index]);
            
            for (unsigned dim=0; dim<SPACE_DIM; dim++)
            {
                cart_distances[rOriginSurface[index]][dim] = 0u;
            }     
        }
        
        while (!active_nodes.empty())
        {
            // Get a pointer to the next node in the queue
            unsigned current_node_index = active_nodes.front();
            Node<SPACE_DIM>* p_current_node = mrMesh.GetNode(current_node_index);
            
            // Loop over the elements containing the given node
            for(typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
                element_iterator != p_current_node->ContainingElementsEnd();
                ++element_iterator)
            {
                // Get a pointer to the container element
                Element<SPACE_DIM,SPACE_DIM>* p_containing_element = mrMesh.GetElement(*element_iterator);

               // Loop over the nodes of the element
               for(unsigned node_local_index=0; 
                   node_local_index<p_containing_element->GetNumNodes();
                   node_local_index++)
               {
                    Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
                    unsigned neighbour_node_index = p_neighbour_node->GetIndex();
                
                    // Avoid revisiting the active node
                    if(neighbour_node_index != current_node_index)
                    {
                        // Test if we have found a shorter path from the origin surface to the current neighbour through current node
                        if ((CartToEucliDistance(cart_distances[current_node_index]) 
                              + EuclideanDistanceTwoPoints(p_current_node->rGetLocation(), p_neighbour_node->rGetLocation()))
                            < CartToEucliDistance(cart_distances[neighbour_node_index]))
                        {
                            cart_distances[neighbour_node_index] = cart_distances[current_node_index] 
                                                                   + (p_current_node->rGetLocation() - p_neighbour_node->rGetLocation());
                            active_nodes.push(neighbour_node_index); 
                        }                                                                                             
                    }                 
               }
            }            
            
            active_nodes.pop();
        }

        for (unsigned index=0; index<mNumNodes; index++)
        {
            rNodeDistances[index] = CartToEucliDistance(cart_distances[index]);             
        }
        
    }
    
};

#endif /*DISTANCEMAPCALCULATOR_HPP_*/
