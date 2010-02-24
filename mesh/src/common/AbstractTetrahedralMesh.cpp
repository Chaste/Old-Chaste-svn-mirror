/*

Copyright (C) University of Oxford, 2005-2010

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

#include "AbstractTetrahedralMesh.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
{
    assert(hi >= lo);
    for (unsigned element_index=0; element_index<mElements.size(); element_index++)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[element_index];
        p_element->SetOwnership(false);
        for (unsigned local_node_index=0; local_node_index< p_element->GetNumNodes(); local_node_index++)
        {
            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
            if (lo<=global_node_index && global_node_index<hi)
            {
                p_element->SetOwnership(true);
                break;
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::AbstractTetrahedralMesh()
    : mMeshIsLinear(true)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::~AbstractTetrahedralMesh()
{
    // Iterate over elements and free the memory
    for (unsigned i=0; i<mElements.size(); i++)
    {
        delete mElements[i];
    }
    // Iterate over boundary elements and free the memory
    for (unsigned i=0; i<mBoundaryElements.size(); i++)
    {
        delete mBoundaryElements[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllBoundaryElements() const
{
    return mBoundaryElements.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryElements() const
{
    return mBoundaryElements.size();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    unsigned local_index = SolveElementMapping(index);
    return mElements[local_index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElement(unsigned index) const
{
    unsigned local_index = SolveBoundaryElementMapping(index);
    return mBoundaryElements[local_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorBegin() const
{
    return mBoundaryElements.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetBoundaryElementIteratorEnd() const
{
    return mBoundaryElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(
        unsigned elementIndex,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
        double& rJacobianDeterminant,
        c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    mElements[SolveElementMapping(elementIndex)]->CalculateInverseJacobian(rJacobian, rJacobianDeterminant, rInverseJacobian);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(
        unsigned elementIndex,
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double& rJacobianDeterminant) const
{
    mBoundaryElements[SolveBoundaryElementMapping(elementIndex)]->CalculateWeightedDirection(rWeightedDirection, rJacobianDeterminant );
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructLinearMesh(unsigned width)
{
    assert(ELEMENT_DIM == 1);

    for (unsigned node_index=0; node_index<=width; node_index++)
    {
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index, node_index==0 || node_index==width, node_index);
        this->mNodes.push_back(p_node); // create node
        if (node_index==0) // create left boundary node and boundary element
        {
            this->mBoundaryNodes.push_back(p_node);
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(0, p_node) );
        }
        if (node_index==width) // create right boundary node and boundary element
        {
            this->mBoundaryNodes.push_back(p_node);
            this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(1, p_node) );
        }
        if (node_index>0) // create element
        {
            std::vector<Node<SPACE_DIM>*> nodes;
            nodes.push_back(this->mNodes[node_index-1]);
            nodes.push_back(this->mNodes[node_index]);
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(node_index-1, nodes) );
        }
    }

    this->RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
{
    assert(SPACE_DIM == 2);
    assert(ELEMENT_DIM == 2);

    //Construct the nodes
    unsigned node_index=0;
    for (unsigned j=0; j<height+1; j++)
    {
        for (unsigned i=0; i<width+1; i++)
        {
            bool is_boundary=false;
            if (i==0 || j==0 || i==width || j==height)
            {
                is_boundary=true;
            }
            //Check in place for parallel
            assert(node_index==(width+1)*(j) + i);
            Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index++, is_boundary, i, j);
            this->mNodes.push_back(p_node);
            if (is_boundary)
            {
                this->mBoundaryNodes.push_back(p_node);
            }
        }
    }

    //Construct the boundary elements
    unsigned belem_index=0;
    //Top
    for (unsigned i=0; i<width; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[height*(width+1)+i]);
        nodes.push_back(this->mNodes[height*(width+1)+i+1]);
        assert(belem_index==i);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Right
    for (unsigned j=1; j<=height; j++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[(width+1)*(j+1)-1]);
        nodes.push_back(this->mNodes[(width+1)*j-1]);
        assert(belem_index==width+j-1);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Bottom
    for (unsigned i=0; i<width; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[i+1]);
        nodes.push_back(this->mNodes[i]);
        assert(belem_index==width+height+i);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Left
    for (unsigned j=0; j<height; j++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[(width+1)*(j+1)]);
        nodes.push_back(this->mNodes[(width+1)*(j)]);
        assert(belem_index==2*width+height+j);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }

    //Construct the elements
    unsigned elem_index = 0;
    for (unsigned j=0; j<height; j++)
    {
        for (unsigned i=0; i<width; i++)
        {
            unsigned parity=(i+(height-j))%2;//Note that parity is measured from the top-left (not bottom left) for historical reasons
            unsigned nw=(j+1)*(width+1)+i; //ne=nw+1
            unsigned sw=(j)*(width+1)+i;   //se=sw+1
            std::vector<Node<SPACE_DIM>*> upper_nodes;
            upper_nodes.push_back(this->mNodes[nw]);
            upper_nodes.push_back(this->mNodes[nw+1]);
            if (stagger==false  || parity == 1)
            {
                upper_nodes.push_back(this->mNodes[sw+1]);
            }
            else
            {
                upper_nodes.push_back(this->mNodes[sw]);
            }
            assert(elem_index==2*(j*width+i));
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,upper_nodes));
            std::vector<Node<SPACE_DIM>*> lower_nodes;
            lower_nodes.push_back(this->mNodes[sw+1]);
            lower_nodes.push_back(this->mNodes[sw]);
            if (stagger==false  ||parity == 1)
            {
                lower_nodes.push_back(this->mNodes[nw]);
            }
            else
            {
                lower_nodes.push_back(this->mNodes[nw+1]);
            }
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,lower_nodes));
        }
    }

    this->RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
        unsigned height,
        unsigned depth)
{
    assert(SPACE_DIM == 3);
    assert(ELEMENT_DIM == 3);
    //Construct the nodes

    unsigned node_index = 0;
    for (unsigned k=0; k<depth+1; k++)
    {
        for (unsigned j=0; j<height+1; j++)
        {
            for (unsigned i=0; i<width+1; i++)
            {
                bool is_boundary = false;
                if (i==0 || j==0 || k==0 || i==width || j==height || k==depth)
                {
                    is_boundary = true;
                }
                assert(node_index == (k*(height+1)+j)*(width+1)+i);
                Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(node_index++, is_boundary, i, j, k);

                this->mNodes.push_back(p_node);
                if (is_boundary)
                {
                    this->mBoundaryNodes.push_back(p_node);
                }
            }
        }
    }

    // Construct the elements

    unsigned elem_index = 0;
    unsigned belem_index = 0;
    unsigned element_nodes[6][4] = {{0, 1, 5, 7}, {0, 1, 3, 7},
                                        {0, 2, 3, 7}, {0, 2, 6, 7},
                                        {0, 4, 6, 7}, {0, 4, 5, 7}};
/* Alternative tessellation - (gerardus)
 * Note that our method (above) has a bias that all tetrahedra share a
 * common edge (the diagonal 0 - 7).  In the following method the cube is
 * split along the "face diagonal" 1-2-5-6 into two prisms.  This also has a bias.
 *
    unsigned element_nodes[6][4] = {{ 0, 6, 5, 4},
                                    { 0, 2, 6, 1},
                                    { 0, 1, 6, 5},
                                    { 1, 2, 3, 7},
                                    { 1, 2, 6, 7},
                                    { 1, 6, 7, 5 }};
*/
    std::vector<Node<SPACE_DIM>*> tetrahedra_nodes;

    for (unsigned k=0; k<depth; k++)
    {
        if (k!=0)
        {
            // height*width squares on upper face, k layers of 2*height+2*width square aroun
            assert(belem_index ==   2*(height*width+k*2*(height+width)) );
        }
        for (unsigned j=0; j<height; j++)
        {
            for (unsigned i=0; i<width; i++)
            {
                // Compute the nodes' index
                unsigned global_node_indices[8];
                unsigned local_node_index = 0;

                for (unsigned z = 0; z < 2; z++)
                {
                    for (unsigned y = 0; y < 2; y++)
                    {
                        for (unsigned x = 0; x < 2; x++)
                        {
                            global_node_indices[local_node_index] = i+x+(width+1)*(j+y+(height+1)*(k+z));

                            local_node_index++;
                        }
                    }
                }

                for (unsigned m = 0; m < 6; m++)
                {
                    // Tetrahedra #m

                    tetrahedra_nodes.clear();

                    for (unsigned n = 0; n < 4; n++)
                    {
                        tetrahedra_nodes.push_back(this->mNodes[global_node_indices[element_nodes[m][n]]]);
                    }

                    assert(elem_index == 6 * ((k*height+j)*width+i)+m );
                    this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++, tetrahedra_nodes));
                }

                //Are we at a boundary?
                std::vector<Node<SPACE_DIM>*> triangle_nodes;

                if (i == 0) //low face at x==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (i == width-1) //high face at x=width
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == 0) //low face at y==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (j == height-1) //high face at y=height
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == 0) //low face at z==0
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[2]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[0]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[1]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[3]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
                if (k == depth-1) //high face at z=depth
                {
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[5]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                    triangle_nodes.clear();
                    triangle_nodes.push_back(this->mNodes[global_node_indices[4]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[6]]);
                    triangle_nodes.push_back(this->mNodes[global_node_indices[7]]);
                    this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,triangle_nodes));
                }
            }//i
        }//j
    }//k

    this->RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex ) 
{
        unsigned tie_break_index = this->GetBoundaryElement(faceIndex)->GetNodeGlobalIndex(0);

        unsigned hi = this->GetDistributedVectorFactory()->GetHigh();
        unsigned lo = this->GetDistributedVectorFactory()->GetLow();
        //if it is in my range
        if (tie_break_index>=lo && tie_break_index<hi)
        {
            return true;
        }
        else
        {
            return false;
        }

}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractTetrahedralMesh<1,1>;
template class AbstractTetrahedralMesh<1,2>;
template class AbstractTetrahedralMesh<1,3>;
template class AbstractTetrahedralMesh<2,2>;
template class AbstractTetrahedralMesh<2,3>;
template class AbstractTetrahedralMesh<3,3>;
