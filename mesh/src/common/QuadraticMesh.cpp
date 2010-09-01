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

#include "QuadraticMesh.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"

//Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen.h"
#undef REAL
#undef VOID


template<unsigned DIM>
QuadraticMesh<DIM>::QuadraticMesh(double spaceStep, double width, double height, double depth)
{
    this->ConstructRegularSlabMesh(spaceStep, width, height, depth);
}

template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructLinearMesh(unsigned numElemX)
{
    this->mNodes.resize(2*numElemX+1);
    mNumVertices = numElemX+1;
    
    // create the left-most node
    Node<DIM>* p_edge_node = new Node<DIM>(0, true, 0.0);
    this->mNodes[0] = p_edge_node; // create nodes
    this->mBoundaryNodes.push_back(p_edge_node);
    this->mBoundaryElements.push_back(new BoundaryElement<DIM-1,DIM>(0, p_edge_node) );
    
    for (unsigned element_index=0; element_index<numElemX; element_index++)
    {
        unsigned right_node_index = element_index+1;
        unsigned mid_node_index = mNumVertices + element_index;
        
        double x_value_right_x = right_node_index;
        double x_value_mid_node = x_value_right_x-0.5;
        
        bool is_boundary = (element_index+1==numElemX);
        Node<DIM>* p_right_node = new Node<DIM>(right_node_index, is_boundary, x_value_right_x); 
        Node<DIM>* p_mid_node   = new Node<DIM>(mid_node_index, false, x_value_mid_node);

        this->mNodes[right_node_index] = p_right_node;
        this->mNodes[mid_node_index] = p_mid_node;
        
        if (element_index+1==numElemX) // right boundary
        {
            this->mBoundaryNodes.push_back(p_right_node);
            this->mBoundaryElements.push_back(new BoundaryElement<DIM-1,DIM>(1, p_right_node) );
        }

        std::vector<Node<DIM>*> nodes;
        nodes.push_back(this->mNodes[right_node_index-1]);
        nodes.push_back(this->mNodes[right_node_index]);
        nodes.push_back(this->mNodes[mid_node_index]);
        this->mElements.push_back(new Element<DIM,DIM>(element_index, nodes) );
    }
    
    this->RefreshMesh();
}




template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructRectangularMesh(unsigned numElemX, unsigned numElemY, bool unused)
{
    assert(DIM==2);

    assert(numElemX>0);
    assert(numElemY>0);
    assert(unused);
    
    this->mMeshIsLinear=false;
    unsigned num_nodes=(numElemX+1)*(numElemY+1);
    struct triangulateio mesher_input;
    mesher_input.pointlist = (double *) malloc( num_nodes * DIM * sizeof(double));
    mesher_input.numberofpoints = num_nodes;
    mesher_input.numberofpointattributes = 0;
    mesher_input.pointmarkerlist = NULL;
    mesher_input.numberofsegments = 0;
    mesher_input.numberofholes = 0;
    mesher_input.numberofregions = 0;

    unsigned new_index = 0;
    for (unsigned j=0; j<=numElemY; j++)
    {
        double y = j;
        for (unsigned i=0; i<=numElemX; i++)
        {
            double x = i;

            mesher_input.pointlist[DIM*new_index] = x;
            mesher_input.pointlist[DIM*new_index + 1] = y;
            new_index++;
        }
    }

    // Make structure for output
    struct triangulateio mesher_output;
    mesher_output.pointlist = NULL;
    mesher_output.pointattributelist = (double *) NULL;
    mesher_output.pointmarkerlist = (int *) NULL;
    mesher_output.trianglelist = (int *) NULL;
    mesher_output.triangleattributelist = (double *) NULL;
    mesher_output.edgelist = (int *) NULL;
    mesher_output.edgemarkerlist = (int *) NULL;

    // Library call
    triangulate((char*)"Qzeo2", &mesher_input, &mesher_output, NULL);

    assert(mesher_output.numberofcorners == (DIM+1)*(DIM+2)/2);//Nodes per element (including internals, one per edge)

    // Construct the nodes
    for (unsigned node_index=0; node_index<(unsigned)mesher_output.numberofpoints; node_index++)
    {
        if (mesher_output.pointmarkerlist[node_index] == 1)
        {
            // Boundary node
            Node<DIM>* p_node = new Node<DIM>(node_index, true,
              mesher_output.pointlist[node_index * DIM],
              mesher_output.pointlist[node_index * DIM+1]);
            this->mNodes.push_back(p_node);
            this->mBoundaryNodes.push_back(p_node);
            ///\todo #1545 Make this code more like the 3D code (where we don't have the edge markers)
        }
        else
        {
            this->mNodes.push_back(new Node<DIM>(node_index, false,
              mesher_output.pointlist[node_index * DIM],
              mesher_output.pointlist[node_index * DIM+1]));
        }
    }

    mIsInternalNode.resize(this->GetNumNodes(), true);

    // Construct the elements
    this->mElements.reserve(mesher_output.numberoftriangles);
    for (unsigned element_index=0; element_index<(unsigned)mesher_output.numberoftriangles; element_index++)
    {
        std::vector<Node<DIM>*> nodes;
        //First (DIM+1) are the vertices
        for (unsigned j=0; j<DIM+1; j++)
        {
            unsigned global_node_index = mesher_output.trianglelist[element_index*mesher_output.numberofcorners + j];
            assert(global_node_index < this->mNodes.size());
            nodes.push_back(this->mNodes[global_node_index]);
            mIsInternalNode[global_node_index]=false;
        }
        //Construct with just the vertices
        this->mElements.push_back(new Element<DIM, DIM>(element_index, nodes));
        //Add the internals
        for (unsigned j=DIM+1; j<(unsigned)mesher_output.numberofcorners; j++)
        {
            unsigned global_node_index = mesher_output.trianglelist[element_index*mesher_output.numberofcorners + j];
            assert(global_node_index < this->mNodes.size());
            this->mElements[element_index]->AddNode( this->mNodes[global_node_index] );
            this->mNodes[global_node_index]->AddElement(element_index);
        }
    }
    //Check that the nodes go "Vertex, Vertex, ...., Internal, Internal..."
    bool vertices_mode = true;
    mNumVertices = 0u;
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        if (mIsInternalNode[i]==false)
        {
            mNumVertices++;
            assert(vertices_mode);//If this trips, then the nodes were not in the expected order -- investigate the library call to triangle
        }
        if ( (vertices_mode == true)  && (mIsInternalNode[i]==true) )
        {
            vertices_mode = false;
        }
    }
    unsigned next_boundary_element_index = 0;
    for (unsigned boundary_element_index=0; boundary_element_index<(unsigned)mesher_output.numberofedges; boundary_element_index++)
    {
        if (mesher_output.edgemarkerlist[boundary_element_index] == 1)
        {
            std::vector<Node<DIM>*> nodes;
            for (unsigned j=0; j<DIM; j++)
            {
                unsigned global_node_index=mesher_output.edgelist[boundary_element_index*DIM + j];
                assert(global_node_index < this->mNodes.size());
                nodes.push_back(this->mNodes[global_node_index]);
            }
            this->mBoundaryElements.push_back(new BoundaryElement<DIM-1, DIM>(next_boundary_element_index++, nodes));
        }
    }

    this->RefreshJacobianCachedData();

    AddNodesToBoundaryElements(NULL);

    free(mesher_input.pointlist);

    free(mesher_output.pointlist);
    free(mesher_output.pointattributelist);
    free(mesher_output.pointmarkerlist);
    free(mesher_output.trianglelist);
    free(mesher_output.triangleattributelist);
    free(mesher_output.edgelist);
    free(mesher_output.edgemarkerlist);
}



template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructCuboid(unsigned numElemX, unsigned numElemY, unsigned numElemZ)
{
    assert(DIM==3);


    assert(numElemX>0);
    assert(numElemY>0);
    assert(numElemZ>0);
    this->mMeshIsLinear=false;
    unsigned num_nodes=(numElemX+1)*(numElemY+1)*(numElemZ+1);
    
    struct tetgen::tetgenio mesher_input;
    mesher_input.pointlist =  new double[num_nodes * DIM];
    mesher_input.numberofpoints = num_nodes;
    unsigned new_index = 0;
    for (unsigned k=0; k<=numElemZ; k++)
    {
        double z = k;
        for (unsigned j=0; j<=numElemY; j++)
        {
            double y = j;
            for (unsigned i=0; i<=numElemX; i++)
            {
                double x = i;
                mesher_input.pointlist[DIM*new_index] = x;
                mesher_input.pointlist[DIM*new_index + 1] = y;
                mesher_input.pointlist[DIM*new_index + 2] = z;
                new_index++;
            }
        }
    }

    
    // Library call
    struct tetgen::tetgenio mesher_output;
    tetgen::tetrahedralize((char*)"Qo2", &mesher_input, &mesher_output, NULL);
    
    assert(mesher_output.numberofcorners == (DIM+1)*(DIM+2)/2);//Nodes per element (including internals, one per edge)

    // Construct the nodes
    for (unsigned node_index=0; node_index<(unsigned)mesher_output.numberofpoints; node_index++)
    {
        this->mNodes.push_back(new Node<DIM>(node_index, false,
          mesher_output.pointlist[node_index * DIM],
          mesher_output.pointlist[node_index * DIM+1],
          mesher_output.pointlist[node_index * DIM+2]));
    }
    

    mIsInternalNode.resize(this->GetNumNodes(), true);

   // Construct the elements
    this->mElements.reserve(mesher_output.numberoftetrahedra);
    for (unsigned element_index=0; element_index<(unsigned)mesher_output.numberoftetrahedra; element_index++)
    {
        std::vector<Node<DIM>*> nodes;
        //First (DIM+1) are the vertices
        for (unsigned j=0; j<DIM+1; j++)
        {
            unsigned global_node_index = mesher_output.tetrahedronlist[element_index*mesher_output.numberofcorners + j];
            assert(global_node_index < this->mNodes.size());
            nodes.push_back(this->mNodes[global_node_index]);
            mIsInternalNode[global_node_index]=false;
        }
        //Construct with just the vertices
        this->mElements.push_back(new Element<DIM, DIM>(element_index, nodes));
        //Add the internals
        for (unsigned j=DIM+1; j<(unsigned)mesher_output.numberofcorners; j++)
        {
            unsigned global_node_index = mesher_output.tetrahedronlist[element_index*mesher_output.numberofcorners + j];
            assert(global_node_index < this->mNodes.size());
            this->mElements[element_index]->AddNode( this->mNodes[global_node_index] );
            this->mNodes[global_node_index]->AddElement(element_index);
        }
    }
    
    //Check that the nodes go "Vertex, Vertex, ...., Internal, Internal..."
    bool vertices_mode = true;
    mNumVertices = 0u;
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        if (mIsInternalNode[i]==false)
        {
            mNumVertices++;
            assert(vertices_mode);//If this trips, then the nodes were not in the expected order -- investigate the library call to triangle
        }
        if ( (vertices_mode == true)  && (mIsInternalNode[i]==true) )
        {
            vertices_mode = false;
        }
    }
    
    unsigned next_boundary_element_index = 0;
    for (unsigned boundary_element_index=0; boundary_element_index<(unsigned)mesher_output.numberoftrifaces; boundary_element_index++)
    {
        std::vector<Node<DIM>*> nodes;
        for (unsigned j=0; j<DIM; j++)
        {
            unsigned global_node_index=mesher_output.trifacelist[boundary_element_index*DIM + j];
            assert(global_node_index < this->mNodes.size());
            nodes.push_back(this->mNodes[global_node_index]);
            if (!nodes[j]->IsBoundaryNode())
            {
                nodes[j]->SetAsBoundaryNode();
                this->mBoundaryNodes.push_back(nodes[j]);
            }
        }
        this->mBoundaryElements.push_back(new BoundaryElement<DIM-1, DIM>(next_boundary_element_index++, nodes));
    }
    
    this->RefreshJacobianCachedData();
    
    AddNodesToBoundaryElements(NULL);

}


template<unsigned DIM>
unsigned QuadraticMesh<DIM>::GetNumVertices()
{
    return mNumVertices;
}


template<unsigned DIM>
void QuadraticMesh<DIM>::RunMesherAndReadMesh(std::string outputDir,
                                              std::string fileStem)
{
    assert(DIM == 3);

    if(PetscTools::AmMaster())
    {
        std::string args = "-Qo2";

        std::string command =  "tetgen " + args + " " + outputDir
                           + "/" + fileStem + ".node" + " > /dev/null";

        int return_value = system(command.c_str());

        if (return_value != 0)
        {
            #define COVERAGE_IGNORE
            EXCEPTION("Remeshing (by calling tetgen) failed.  Do you have it in your path?\nThe quadratic mesh relies on functionality from tetgen (http://tetgen.berlios.de/).");
            #undef COVERAGE_IGNORE
        }

        // move the output files to the chaste directory
        command =   "mv " + outputDir + "/"
                  + fileStem + ".1.* .";

        // NOTE: we don't check whether the return value here is zero, because if CHASTE_TESTOUTPUT
        // is "." (ie if it hasn't been exported), then the mv will fail (source and destination files
        // are the same), but this isn't a problem.
        return_value = system(command.c_str());
    }
    PetscTools::Barrier();

    // load
    TrianglesMeshReader<DIM,DIM> mesh_reader(fileStem + ".1", 2, 1, false); // false as tetgen/triangle has been used and therefore boundary elems will be linear;
    ConstructFromMeshReader(mesh_reader);
    ///\todo: Could use the '-nn' flag when calling tetgen and then face file would have containing element info and second false
    // could be a true instead. Currently though there is the intermediate step of having to delete manually the attribute values
    // column after using this flag.

    if(PetscTools::AmMaster())
    {
        // delete the temporary files
        std::string command = "rm -f " + outputDir + "/" + fileStem + ".node";
        EXPECT0(system, command);
        EXPECT0(system, "rm -f " + fileStem + ".1.node");
        EXPECT0(system, "rm -f " + fileStem + ".1.ele");
        EXPECT0(system, "rm -f " + fileStem + ".1.face");
    }
}


template<unsigned DIM>
void QuadraticMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rAbsMeshReader)
{
    TrianglesMeshReader<DIM, DIM>* p_mesh_reader=dynamic_cast<TrianglesMeshReader<DIM, DIM>*>(&rAbsMeshReader);
    assert(p_mesh_reader != NULL);


    if (p_mesh_reader->GetOrderOfElements() == 1)
    {
        EXCEPTION("Supplied mesh reader is reading a linear mesh into quadratic mesh");
    }

    TetrahedralMesh<DIM,DIM>::ConstructFromMeshReader(*p_mesh_reader);
    assert(this->GetNumBoundaryElements()>0);

    // set up the information on whether a node is an internal node or not (if not,
    // it'll be a vertex)
    mIsInternalNode.resize(this->GetNumNodes(), true);
    for (unsigned elem_index=0; elem_index<this->GetNumElements(); elem_index++)
    {
        for (unsigned i=0; i<DIM+1 /*num vertices*/; i++)
        {
            unsigned node_index = this->GetElement(elem_index)->GetNodeGlobalIndex(i);
            mIsInternalNode[ node_index ] = false;
        }
    }
    // count the number of vertices, and also check all vertices come before the
    // rest of the nodes (as this is assumed in other parts of the code)
    mNumVertices = 0;
    bool vertices_mode = true;
    for (unsigned i=0; i<this->GetNumNodes(); i++)
    {
        if (mIsInternalNode[i]==false)
        {
            mNumVertices++;
        }
        if ((vertices_mode == false)  && (mIsInternalNode[i]==false ) )
        {
            //Covered in the 1D case by special mesh file data
            EXCEPTION("The quadratic mesh doesn't appear to have all vertices before the rest of the nodes");
        }
        if ( (vertices_mode == true)  && (mIsInternalNode[i]==true) )
        {
            vertices_mode = false;
        }
    }

    p_mesh_reader->Reset();


    // add the extra nodes (1 extra node in 1D, 3 in 2D, 6 in 3D) to the element
    // data.
    for (unsigned i=0; i<this->GetNumElements(); i++)
    {
        std::vector<unsigned> nodes = p_mesh_reader->GetNextElementData().NodeIndices;
        assert(nodes.size()==(DIM+1)*(DIM+2)/2);
        assert(this->GetElement(i)->GetNumNodes()==DIM+1); // element is initially linear
        // add extra nodes to make it a quad element
        for (unsigned j=DIM+1; j<(DIM+1)*(DIM+2)/2; j++)
        {
            this->GetElement(i)->AddNode( this->GetNode(nodes[j]) );
            this->GetNode(nodes[j])->AddElement(this->GetElement(i)->GetIndex());
        }
    }
    if (DIM > 1)
    {
        // if  OrderOfBoundaryElements is 2 it can read in the extra nodes for each boundary element, other have to compute them.
        if (p_mesh_reader->GetOrderOfBoundaryElements() == 2u)
        {
            p_mesh_reader->Reset();
            for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
                  = this->GetBoundaryElementIteratorBegin();
                 iter != this->GetBoundaryElementIteratorEnd();
                 ++iter)
            {
                std::vector<unsigned> nodes = p_mesh_reader->GetNextFaceData().NodeIndices;

                assert((*iter)->GetNumNodes()==DIM); // so far just the vertices
                assert(nodes.size()==DIM*(DIM+1)/2); // the reader should have got 6 nodes (3d) for each face

                for (unsigned j=DIM; j<DIM*(DIM+1)/2; j++)
                {
                    (*iter)->AddNode( this->GetNode(nodes[j]) );
                }
            }
        }
        else
        {
            AddNodesToBoundaryElements(p_mesh_reader);
        }
    }

    // Check each boundary element has a quadratic number of nodes
#ifndef NDEBUG
    unsigned expected_num_nodes = DIM*(DIM+1)/2;
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
          = this->GetBoundaryElementIteratorBegin();
          iter != this->GetBoundaryElementIteratorEnd();
          ++iter)
    {
        assert((*iter)->GetNumNodes()==expected_num_nodes);
    }
#endif
}

template<unsigned DIM>
void QuadraticMesh<DIM>::AddNodesToBoundaryElements(TrianglesMeshReader<DIM,DIM>* pMeshReader)
 {
    // Loop over all boundary elements, find the equivalent face from all
    // the elements, and add the extra nodes to the boundary element
    bool boundary_element_file_has_containing_element_info=false;

    if (pMeshReader)
    {
        boundary_element_file_has_containing_element_info=pMeshReader->GetReadContainingElementOfBoundaryElement();
    }

    if (DIM>1)
    {
        if(boundary_element_file_has_containing_element_info)
        {
            pMeshReader->Reset();
        }

        //unsigned counter = 0;
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
               = this->GetBoundaryElementIteratorBegin();
             iter != this->GetBoundaryElementIteratorEnd();
             ++iter)
        {
            //std::cout << "\rAddNodesToBoundaryElements: " << counter++ << " of " << this->GetNumBoundaryElements() << std::flush;

            // collect the nodes of this boundary element in a set
            std::set<unsigned> boundary_element_node_indices;
            for (unsigned i=0; i<DIM; i++)
            {
                boundary_element_node_indices.insert( (*iter)->GetNodeGlobalIndex(i) );
            }

            bool found_this_boundary_element = false;

            // Loop over elements, then loop over each face of that element, and see if it matches
            // this boundary element.
            // Note, if we know what elem it should be in (boundary_element_file_has_containing_element_info==true)
            // we will reset elem_index immediately (below)
            for (unsigned elem_index=0; elem_index<this->GetNumElements(); elem_index++)
            {
                // we know what elem it should be in
                if(boundary_element_file_has_containing_element_info)
                {
                    elem_index = pMeshReader->GetNextFaceData().ContainingElement;
                }

                Element<DIM,DIM>* p_element = this->GetElement(elem_index);

                // for each element, loop over faces (the opposites to a node)
                for (unsigned face=0; face<DIM+1; face++)
                {
                    // collect the node indices for this face
                    std::set<unsigned> node_indices;
                    for (unsigned local_node_index=0; local_node_index<DIM+1; local_node_index++)
                    {
                        if (local_node_index!=face)
                        {
                            node_indices.insert( p_element->GetNodeGlobalIndex(local_node_index) );
                        }
                    }

                    assert(node_indices.size()==DIM);

                    // see if this face matches the boundary element,
                    // and call AddExtraBoundaryNodes() if so
                    if (node_indices==boundary_element_node_indices)
                    {
                        AddExtraBoundaryNodes(*iter, p_element, face);

                        found_this_boundary_element = true;
                        break;
                    }
                }

                // if the containing element info was given, we should certainly have found the
                // face first time.
                if(boundary_element_file_has_containing_element_info && !found_this_boundary_element)
                {
                    #define COVERAGE_IGNORE
                    std::cout << (*iter)->GetIndex() << " " <<  pMeshReader->GetNextFaceData().ContainingElement << "\n";
                    std::stringstream ss;
                    ss << "Boundary element " << (*iter)->GetIndex()
                       << "wasn't found in the containing element given for it "
                       << elem_index;

                    EXCEPTION(ss.str());
                    #undef COVERAGE_IGNORE
                }

                if (found_this_boundary_element)
                {
                    break;
                }
            }

            if (!found_this_boundary_element)
            {
                #define COVERAGE_IGNORE
                EXCEPTION("Unable to find a face of an element which matches one of the boundary elements");
                #undef COVERAGE_IGNORE
            }
        }
    }
}


template<unsigned DIM>
void QuadraticMesh<DIM>::AddNodeToBoundaryElement(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                                  Element<DIM,DIM>* pElement,
                                                  unsigned internalNode)
{
    assert(DIM>1);
    assert(internalNode >= DIM+1);
    assert(internalNode < (DIM+1)*(DIM+2)/2);
    Node<DIM>* p_internal_node = pElement->GetNode(internalNode);

    // add node to the boundary node list
    if (!p_internal_node->IsBoundaryNode())
    {
        p_internal_node->SetAsBoundaryNode();
        this->mBoundaryNodes.push_back(p_internal_node);
    }

    pBoundaryElement->AddNode( p_internal_node );
}


template<unsigned DIM>
void QuadraticMesh<DIM>::AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                               Element<DIM,DIM>* pElement,
                                               unsigned nodeIndexOppositeToFace)
{
    assert(DIM!=1);
    if (DIM==2)
    {
        assert(nodeIndexOppositeToFace<3);
        // the single internal node of the elements face will be numbered 'face+3'
        AddNodeToBoundaryElement(pBoundaryElement, pElement, nodeIndexOppositeToFace+3);
    }
    else
    {
        assert(DIM==3);

        unsigned b_elem_n0 = pBoundaryElement->GetNodeGlobalIndex(0);
        unsigned b_elem_n1 = pBoundaryElement->GetNodeGlobalIndex(1);

        unsigned offset;
        bool reverse;

        if (nodeIndexOppositeToFace==0)
        {
            // face opposite to node 0 = {1,2,3}, with corresponding internals {9,8,5}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 1, 2, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 9, 8, 5, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==1)
        {
            // face opposite to node 1 = {2,0,3}, with corresponding internals {7,9,6}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 2, 0, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 7, 9, 6, offset, reverse);
        }
        else if (nodeIndexOppositeToFace==2)
        {
            // face opposite to node 2 = {0,1,3}, with corresponding internals {8,7,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 3, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 8, 7, 4, offset, reverse);
        }
        else
        {
            assert(nodeIndexOppositeToFace==3);
            // face opposite to node 3 = {0,1,2}, with corresponding internals {5,6,4}
            HelperMethod1(b_elem_n0, b_elem_n1, pElement, 0, 1, 2, offset, reverse);
            HelperMethod2(pBoundaryElement, pElement, 5, 6, 4, offset, reverse);
        }
    }
}

template<unsigned DIM>
void QuadraticMesh<DIM>::WriteBoundaryElementFile(std::string directory, std::string fileName)
{
    OutputFileHandler handler(directory, false);
    out_stream p_file = handler.OpenOutputFile(fileName);

    unsigned expected_num_nodes;
    assert(DIM > 1);
    if (DIM == 2)
    {
        expected_num_nodes = 3;
    }
    else if (DIM == 3)
    {
        expected_num_nodes = 6;
    }

    unsigned num_elements = 0;

    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
          = this->GetBoundaryElementIteratorBegin();
          iter != this->GetBoundaryElementIteratorEnd();
          ++iter)
    {
        assert((*iter)->GetNumNodes()==expected_num_nodes);
        num_elements++;
    }

    *p_file << num_elements << " 0\n";

    unsigned counter = 0;
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
          = this->GetBoundaryElementIteratorBegin();
          iter != this->GetBoundaryElementIteratorEnd();
          ++iter)
    {
        *p_file << counter++ << " ";
        for (unsigned i=0; i<(*iter)->GetNumNodes(); i++)
        {
            *p_file << (*iter)->GetNodeGlobalIndex(i) << " ";
        }
        *p_file << "\n";
    }

    p_file->close();
}

///////////////////////////////////////////////////////////////////////////////
// two unpleasant helper methods for AddExtraBoundaryNodes()
///////////////////////////////////////////////////////////////////////////////

#define COVERAGE_IGNORE /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void QuadraticMesh<DIM>::HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                                       Element<DIM,DIM>* pElement,
                                       unsigned node0, unsigned node1, unsigned node2,
                                       unsigned& rOffset,
                                       bool& rReverse)
{
    if (pElement->GetNodeGlobalIndex(node0)==boundaryElemNode0)
    {
        rOffset = 0;
        if (pElement->GetNodeGlobalIndex(node1)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node2)==boundaryElemNode1);
            rReverse = true;
        }
    }
    else if (pElement->GetNodeGlobalIndex(node1)==boundaryElemNode0)
    {
        rOffset = 1;
        if (pElement->GetNodeGlobalIndex(node2)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node0)==boundaryElemNode1);
            rReverse = true;
        }
    }
    else
    {
        assert(pElement->GetNodeGlobalIndex(node2)==boundaryElemNode0);
        rOffset = 2;
        if (pElement->GetNodeGlobalIndex(node0)==boundaryElemNode1)
        {
            rReverse = false;
        }
        else
        {
            assert(pElement->GetNodeGlobalIndex(node1)==boundaryElemNode1);
            rReverse = true;
        }
    }
}
#undef COVERAGE_IGNORE /// \todo These helper methods aren't properly covered


#define COVERAGE_IGNORE /// \todo These helper methods aren't properly covered
template<unsigned DIM>
void QuadraticMesh<DIM>::HelperMethod2(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                       Element<DIM,DIM>* pElement,
                                       unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                                       unsigned offset,
                                       bool reverse)
{
    if (offset==1)
    {
        unsigned temp = internalNode0;
        internalNode0 = internalNode1;
        internalNode1 = internalNode2;
        internalNode2 = temp;
    }
    else if (offset == 2)
    {
        unsigned temp = internalNode0;
        internalNode0 = internalNode2;
        internalNode2 = internalNode1;
        internalNode1 = temp;
    }

    if (reverse)
    {
        unsigned temp = internalNode1;
        internalNode1 = internalNode2;
        internalNode2 = temp;
    }

    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode0);
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode1);
    AddNodeToBoundaryElement(pBoundaryElement, pElement, internalNode2);
}
#undef COVERAGE_IGNORE /// \todo These helper methods aren't properly covered


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class QuadraticMesh<1>;
template class QuadraticMesh<2>;
template class QuadraticMesh<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh)
