/*

Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "VertexElement.hpp"
#include "RandomNumberGenerator.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index, std::vector<Node<SPACE_DIM>*> nodes)
    : AbstractElement<ELEMENT_DIM, SPACE_DIM>(index, nodes)
{
    RegisterWithNodes();
    mVertexElementArea = DOUBLE_UNSET;
    mVertexElementPerimeter = DOUBLE_UNSET;
    mElementModified = true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::~VertexElement()
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::RegisterWithNodes()
{
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        this->mNodes[i]->AddElement(this->mIndex);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;
    
    // Update nodes in the element so they know they are not contained by it
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        this->mNodes[i]->RemoveElement(this->mIndex);
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::ResetIndex(unsigned index)
{
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
       this->mNodes[i]->RemoveElement(this->mIndex);
    }
    this->mIndex = index;
    RegisterWithNodes();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Remove it from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Update the node at this location
    this->mNodes[rIndex] = pNode;

    // Add element to this node
    this->mNodes[rIndex]->AddElement(this->mIndex);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>   
void VertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteNode(const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Remove element from the node at this location
    this->mNodes[rIndex]->RemoveElement(this->mIndex);

    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin( ) + rIndex);

    // Flag that element has changed and we need to recalculate area and perimeter
    mElementModified = true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode)
{
    assert(rIndex < this->mNodes.size());

    // Add pNode to rIndex+1 element of mNodes pushing the others up
    this->mNodes.insert( this->mNodes.begin( ) + rIndex+1,  pNode);

    // Add element to this node
    this->mNodes[rIndex+1]->AddElement(this->mIndex);

    // Flag that element has changed and we need to recalculate area and perimeter
    mElementModified = true;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>    
void VertexElement<ELEMENT_DIM, SPACE_DIM>::CalculateVertexElementAreaAndPerimeter()
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> current_node;
    c_vector<double, SPACE_DIM> anticlockwise_node; 

    double temp_vertex_element_area = 0;
    double temp_vertex_element_perimeter = 0;
    unsigned number_of_nodes = this->GetNumNodes();

    for (unsigned i=0; i<number_of_nodes; i++)
    {
        // Find locations of current node and anticlockwise node
        current_node = this->GetNodeLocation(i);
        anticlockwise_node = this->GetNodeLocation((i+1)%number_of_nodes);

        /// \todo will need to change length calculation to something like GetVectorFromAtoB (see #825)

        temp_vertex_element_area += 0.5*(current_node[0]*anticlockwise_node[1] 
                                          - anticlockwise_node[0]*current_node[1]);

        temp_vertex_element_perimeter += norm_2(current_node-anticlockwise_node);
    }

    mVertexElementArea = temp_vertex_element_area;
    mVertexElementPerimeter = temp_vertex_element_perimeter;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>    
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetArea()
{
//    if ((mVertexElementArea == DOUBLE_UNSET)||(mElementModified==true))
//    {
//        this->CalculateVertexElementAreaAndPerimeter();
//        mElementModified = false;
//    }
/// \todo the above code was commented as this method must recalculate the element area
///       whenever it is called by the force law (see #861)
    this->CalculateVertexElementAreaAndPerimeter();
    return mVertexElementArea;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::GetAreaGradientAtNode(unsigned localIndex)
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> area_gradient;

    unsigned next_index = (localIndex+1)%(this->GetNumNodes());
	unsigned previous_index = (this->GetNumNodes()+localIndex-1)%(this->GetNumNodes()); // As localIndex-1 can be -ve which breaks %      
    
    c_vector<double, SPACE_DIM> previous_node_location = this->GetNode(previous_index)->rGetLocation();
    c_vector<double, SPACE_DIM> next_node_location = this->GetNode(next_index)->rGetLocation();

    area_gradient[0] = 0.5*(next_node_location[1] - previous_node_location[1]);
    area_gradient[1] = 0.5*(previous_node_location[0] - next_node_location[0]);

    return area_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::GetPreviousEdgeGradientAtNode(unsigned localIndex)
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> previous_edge_gradient;

    unsigned previous_index = (this->GetNumNodes()+localIndex-1)%(this->GetNumNodes()); // As localIndex-1 can be -ve which breaks %

    c_vector<double, SPACE_DIM> current_node_location = this->GetNode(localIndex)->rGetLocation();
    c_vector<double, SPACE_DIM> previous_node_location = this->GetNode(previous_index)->rGetLocation();

    double previous_edge_length = norm_2(current_node_location - previous_node_location);

    if (previous_edge_length < 1e-12) /// \todo magic number - replace with DBL_EPSILON?
    {   
        std::cout << "\n Should Not Reach "<< std::flush;
        previous_edge_gradient[0] = 0.0;
        previous_edge_gradient[1] = 0.0;
    }
    else
    {   
        previous_edge_gradient[0] = (current_node_location[0] - previous_node_location[0])/previous_edge_length;
        previous_edge_gradient[1] = (current_node_location[1] - previous_node_location[1])/previous_edge_length;
    }       
    return previous_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNextEdgeGradientAtNode(unsigned localIndex)
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> next_edge_gradient;

    unsigned next_index = (localIndex+1)%(this->GetNumNodes());

    c_vector<double, SPACE_DIM> current_node_location = this->GetNode(localIndex)->rGetLocation();
    c_vector<double, SPACE_DIM> next_node_location = this->GetNode(next_index)->rGetLocation();

    double next_edge_length = norm_2(next_node_location - current_node_location);

    if (next_edge_length < 1e-12) /// \todo magic number - replace with DBL_EPSILON?
    {   
        std::cout << "\n Should Not Reach "<< std::flush;
        next_edge_gradient[0] = 0.0;
        next_edge_gradient[1] = 0.0;
    }
    else
    {   
        next_edge_gradient[0] = (current_node_location[0] - next_node_location[0])/next_edge_length;
        next_edge_gradient[1] = (current_node_location[1] - next_node_location[1])/next_edge_length;
    }
    return next_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::GetPerimeterGradientAtNode(unsigned localIndex)
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM==2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> previous_edge_gradient = GetPreviousEdgeGradientAtNode(localIndex);
    c_vector<double, SPACE_DIM> next_edge_gradient = GetNextEdgeGradientAtNode(localIndex);

    return previous_edge_gradient + next_edge_gradient;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetPerimeter()
{
//    if ((mVertexElementArea == DOUBLE_UNSET)||(mElementModified==true))
//    {
//        this->CalculateVertexElementAreaAndPerimeter();
//        mElementModified = false;
//    }
/// \todo the above code was commented as this method must recalculate the element perimeter
///       whenever it is called by the force law (see #861)
    this->CalculateVertexElementAreaAndPerimeter();
    return mVertexElementPerimeter;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
c_vector<double, 3> VertexElement<ELEMENT_DIM, SPACE_DIM>::CalculateMoments()
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
#undef COVERAGE_IGNORE

    c_vector<double, 3> moments = zero_vector<double>(3);

    unsigned node_1;
    unsigned node_2;
    unsigned num_nodes = this->GetNumNodes();

    for (unsigned i=0; i<num_nodes; i++)
    {
        node_1 = i;
        node_2 = (i+1)%num_nodes;

        c_vector<double, 2> pos_1 = this->mNodes[node_1]->rGetLocation();
        c_vector<double, 2> pos_2 = this->mNodes[node_2]->rGetLocation();

        // Ixx
        moments(0) += (pos_2(0)-pos_1(0))*(  pos_1(1)*pos_1(1)*pos_1(1)
                                           + pos_1(1)*pos_1(1)*pos_2(1)
                                           + pos_1(1)*pos_2(1)*pos_2(1)
                                           + pos_2(1)*pos_2(1)*pos_2(1));

        // Iyy
        moments(1) += (pos_2(1)-pos_1(1))*(  pos_1(0)*pos_1(0)*pos_1(0)
                                           + pos_1(0)*pos_1(0)*pos_2(0)
                                           + pos_1(0)*pos_2(0)*pos_2(0)
                                           + pos_2(0)*pos_2(0)*pos_2(0));

        // Ixy
        moments(2) +=   pos_1(0)*pos_1(0)*pos_2(1)*(pos_1(1)*2 + pos_2(1))
                      - pos_2(0)*pos_2(0)*pos_1(1)*(pos_1(1) + pos_2(1)*2)
                      + 2*pos_1(0)*pos_2(0)*(pos_2(1)*pos_2(1) - pos_1(1)*pos_1(1));
    }

    moments(0) /= -12;
    moments(1) /= 12;
    moments(2) /= 24;

    return moments;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::CalculateCentroid()
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);
    c_vector<double, SPACE_DIM> current_node;
    c_vector<double, SPACE_DIM> anticlockwise_node; 

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    unsigned num_nodes = this->GetNumNodes();

    for (unsigned i=0; i<num_nodes; i++)
    {
        // Find locations of current node and anticlockwise node
        current_node = this->GetNodeLocation(i);
        anticlockwise_node = this->GetNodeLocation((i+1)%num_nodes);

        /// \todo will need to change length calculation to something like GetVectorFromAtoB (see #825)

        temp_centroid_x += (current_node[0]+anticlockwise_node[0])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);
        temp_centroid_y += (current_node[1]+anticlockwise_node[1])*(current_node[0]*anticlockwise_node[1]-current_node[1]*anticlockwise_node[0]);               
    }

    double vertex_area = GetArea();
    double centroid_coefficient = 1.0/6.0/vertex_area;

    centroid(0) = centroid_coefficient*temp_centroid_x;
    centroid(1) = centroid_coefficient*temp_centroid_y;

    return centroid;        
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VertexElement<ELEMENT_DIM, SPACE_DIM>::CalculateShortAxis()
{
#define COVERAGE_IGNORE
    assert(SPACE_DIM == 2);
#undef COVERAGE_IGNORE

    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);
    c_vector<double, 3> moments = CalculateMoments();

    double largest_eigenvalue, discriminant;            

    discriminant = sqrt((moments(0) - moments(1))*(moments(0) - moments(1)) + 4.0*moments(2)*moments(2));

    // This is always the largest eigenvalue as both eigenvalues are real as it is a 
    // symmetric matrix
    largest_eigenvalue = ((moments(0) + moments(1)) + discriminant)*0.5;       

    if (fabs(discriminant) < 1e-10)
    {
        // Return a random unit vector
        short_axis(0) = RandomNumberGenerator::Instance()->ranf();
        short_axis(1) = sqrt(1.0-short_axis(0)*short_axis(0));
    }

    else
    {                      
        if (moments(2) == 0.0)
        {
             short_axis(0) = 0.0;
             short_axis(1) = 1.0;                       
        }
        else
        {
            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - largest_eigenvalue)/moments(2);

            double length_short_axis = norm_2(short_axis);

            short_axis /= length_short_axis;
        }     
    }
    return short_axis;        
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNodeLocalIndex(unsigned globalIndex)
{
    unsigned local_index= UINT_MAX;
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        if (this->GetNodeGlobalIndex(i) == globalIndex)
        {
            local_index = i;
        }
    }
    return local_index;  
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VertexElement<1,1>;
template class VertexElement<1,2>;
template class VertexElement<2,2>;
template class VertexElement<2,3>;
template class VertexElement<3,3>;
