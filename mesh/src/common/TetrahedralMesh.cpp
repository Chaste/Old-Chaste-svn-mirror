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

#include "TetrahedralMesh.hpp"

#include <iostream>
#include <cassert>
#include <sstream>
#include <map>

#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "Exception.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"
#include "RandomNumberGenerator.hpp"

/////////////////////////////////////////////////////////////////////////////////////
//   IMPLEMENTATION
/////////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralMesh()
{
    Clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructFromMeshReader(
    AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    this->mMeshFileBaseName = rMeshReader.GetMeshFileBaseName();

    // Record number of corner nodes
    unsigned num_nodes = rMeshReader.GetNumNodes();

    // Reserve memory for nodes, so we don't have problems with pointers stored in
    // elements becoming invalid.
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    //typename std::map<std::pair<unsigned,unsigned>,unsigned>::const_iterator iterator;
    //std::map<std::pair<unsigned,unsigned>,unsigned> internal_nodes_map;

    // Add corner nodes
    std::vector<double> coords;
    for (unsigned i=0; i < num_nodes; i++)
    {
        coords = rMeshReader.GetNextNode();
        this->mNodes.push_back(new Node<SPACE_DIM>(i, coords, false));
    }

    //unsigned new_node_index = mNumCornerNodes;

    rMeshReader.Reset();
    // Add elements
    //new_node_index = mNumCornerNodes;
    this->mElements.reserve(rMeshReader.GetNumElements());

    for (unsigned element_index=0; element_index < (unsigned) rMeshReader.GetNumElements(); element_index++)
    {
        ElementData element_data = rMeshReader.GetNextElementData();
        std::vector<Node<SPACE_DIM>*> nodes;

        // NOTE: currently just reading element vertices from mesh reader - even if it
        // does contain information about internal nodes (ie for quadratics) this is
        // ignored here and used elsewhere: ie don't do this:
        //   unsigned nodes_size = node_indices.size();

        for (unsigned j=0; j<ELEMENT_DIM+1; j++) // num vertices=ELEMENT_DIM+1, may not be equal to nodes_size.
        {
            assert(element_data.NodeIndices[j] <  this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        Element<ELEMENT_DIM,SPACE_DIM>* p_element = new Element<ELEMENT_DIM,SPACE_DIM>(element_index, nodes);
        this->mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetRegion(attribute_value);
        }
    }

    // Add boundary elements and nodes
    for (unsigned face_index=0; face_index<(unsigned)rMeshReader.GetNumFaces(); face_index++)
    {
        ElementData face_data = rMeshReader.GetNextFaceData();
        std::vector<unsigned> node_indices = face_data.NodeIndices;

        // NOTE: as above just read boundary element *vertices* from mesh reader - even if
        // it is a quadratic mesh with internal elements, the extra nodes are ignored here 
        // and used elsewhere: ie, we don't do this:
        //   unsigned nodes_size = node_indices.size();

        std::vector<Node<SPACE_DIM>*> nodes;
        for (unsigned node_index=0; node_index<ELEMENT_DIM; node_index++) // node_index from 0 to DIM-1, not 0 to node.size()-1
        {
            assert(node_indices[node_index] < this->mNodes.size());
            // Add Node pointer to list for creating an element
            nodes.push_back(this->mNodes[node_indices[node_index]]);
        }

        // This is a boundary face
        // Ensure all its nodes are marked as boundary nodes
        
        assert(nodes.size()==ELEMENT_DIM); // just taken vertices of boundary node from 
        for (unsigned j=0; j<nodes.size(); j++)
        {
            if (!nodes[j]->IsBoundaryNode())
            {
                nodes[j]->SetAsBoundaryNode();
                this->mBoundaryNodes.push_back(nodes[j]);
            }
            //Register the index that this bounday element will have
            //with the node
            nodes[j]->AddBoundaryElement(face_index);
        }

        // The added elements will be deleted in our destructor
        BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* p_boundary_element = new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(face_index, nodes);
        this->mBoundaryElements.push_back(p_boundary_element);

        if (rMeshReader.GetNumFaceAttributes() > 0)
        {
            assert(rMeshReader.GetNumFaceAttributes() == 1);
            unsigned attribute_value = face_data.AttributeValue;
            p_boundary_element->SetRegion(attribute_value);
        }
    }

    RefreshJacobianCachedData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ReadNodesPerProcessorFile(const std::string& rNodesPerProcessorFile)
{
    std::vector<unsigned> nodes_per_processor_vec;

    std::ifstream file_stream(rNodesPerProcessorFile.c_str());
    if (file_stream.is_open())
    {
        while (file_stream)
        {
            unsigned nodes_per_processor;
            file_stream >> nodes_per_processor;

            if (file_stream)
            {
                nodes_per_processor_vec.push_back(nodes_per_processor);
            }
        }
    }
    else
    {
        EXCEPTION("Unable to read nodes per processor file " + rNodesPerProcessorFile);
    }

    unsigned sum = 0;
    for (unsigned i=0; i<nodes_per_processor_vec.size(); i++)
    {
        sum += nodes_per_processor_vec[i];
    }

    if (sum != this->GetNumNodes())
    {
        std::stringstream string_stream;
        string_stream << "Sum of nodes per processor, " << sum
                     << ", not equal to number of nodes in mesh, " << this->GetNumNodes();
        EXCEPTION(string_stream.str());
    }
    
    unsigned num_owned=nodes_per_processor_vec[PetscTools::GetMyRank()];
    
    if (nodes_per_processor_vec.size() != PetscTools::GetNumProcs())
    {
        EXCEPTION("Number of processes doesn't match the size of the nodes-per-processor file");
    }
    delete this->mpDistributedVectorFactory;
    this->mpDistributedVectorFactory=new DistributedVectorFactory(this->GetNumNodes(), num_owned);
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CheckIsConforming()
{
    /* Each face of each element is a set of node indices
     * We form a set of these in order to get their parity:
     *   all faces which appear once are inserted into the set
     *   all faces which appear twice are inserted and then removed from the set
     *   we're assuming that faces never appear more than twice
     */
    std::set< std::set<unsigned> > odd_parity_faces;
    
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        for (unsigned face_index=0; face_index<=ELEMENT_DIM; face_index++)
        {
            std::set<unsigned> face_info;
            for (unsigned node_index=0; node_index<=ELEMENT_DIM; node_index++)
            {
                //Leave one index out each time
                if (node_index != face_index)
                {
                    face_info.insert(iter->GetNodeGlobalIndex(node_index));
                }
            }
            //Face is now formed - attempt to find it
            std::set< std::set<unsigned> >::iterator find_face=odd_parity_faces.find(face_info);
            if( find_face != odd_parity_faces.end())
            {
                //Face was in set, so it now has even parity.
                //Remove it via the iterator
                odd_parity_faces.erase(find_face);
            }
            else
            {
                //Face is not in set so it now has odd parity.  Insert it 
                odd_parity_faces.insert(face_info);
            }
            
        }
    }
    /* At this point the odd parity faces should be the same as the
     * boundary elements.  We could check this explicitly or we
     * could just count them.
     */
    return( odd_parity_faces.size() == this->GetNumBoundaryElements());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetVolume()
{
    double mesh_volume = 0.0;

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        mesh_volume += iter->GetVolume(mElementJacobianDeterminants[iter->GetIndex()]);
    }

    return mesh_volume;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceArea()
{
    //ELEMENT_DIM-1 is the dimension of the boundary element
    assert(ELEMENT_DIM >= 1);
    const unsigned bound_element_dim = ELEMENT_DIM-1;
    assert(bound_element_dim < 3);
    if ( bound_element_dim == 0)
    {
        return 0.0;
    }

    double mesh_surface = 0.0;
    typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator it = this->GetBoundaryElementIteratorBegin();

    while (it != this->GetBoundaryElementIteratorEnd())
    {
        mesh_surface += mBoundaryElementJacobianDeterminants[(*it)->GetIndex()];
        it++;
    }

    if ( bound_element_dim == 2)
    {
        mesh_surface /= 2.0;
    }

    return mesh_surface;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(
    const double xMovement,
    const double yMovement,
    const double zMovement)
{
    c_vector<double, SPACE_DIM> displacement;

    switch (SPACE_DIM)
    {
        case 3:
            displacement[2] = zMovement;
        case 2:
            displacement[1] = yMovement;
        case 1:
            displacement[0] = xMovement;
    }

    Translate(displacement);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Translate(c_vector<double, SPACE_DIM> transVec)
{
    unsigned num_nodes = this->GetNumAllNodes();

    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = this->mNodes[i]->rGetModifiableLocation();
        r_location += transVec;
    }

    RefreshMesh();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(
    c_matrix<double , SPACE_DIM, SPACE_DIM> rotationMatrix)
{
    unsigned num_nodes = this->GetNumAllNodes();

    for (unsigned i=0; i<num_nodes; i++)
    {
        c_vector<double, SPACE_DIM>& r_location = this->mNodes[i]->rGetModifiableLocation();
        r_location = prod(rotationMatrix, r_location);
    }

    RefreshMesh();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(c_vector<double,3> axis, double angle)
{
    assert(SPACE_DIM == 3);
    double norm = norm_2(axis);
    c_vector<double,3> unit_axis=axis/norm;

    c_matrix<double, SPACE_DIM,SPACE_DIM> rotation_matrix;

    double c = cos(angle);
    double s = sin(angle);

    rotation_matrix(0,0) = unit_axis(0)*unit_axis(0)+c*(1-unit_axis(0)*unit_axis(0));
    rotation_matrix(0,1) = unit_axis(0)*unit_axis(1)*(1-c) - unit_axis(2)*s;
    rotation_matrix(1,0) = unit_axis(0)*unit_axis(1)*(1-c) + unit_axis(2)*s;
    rotation_matrix(1,1) = unit_axis(1)*unit_axis(1)+c*(1-unit_axis(1)*unit_axis(1));
    rotation_matrix(0,2) = unit_axis(0)*unit_axis(2)*(1-c)+unit_axis(1)*s;
    rotation_matrix(1,2) = unit_axis(1)*unit_axis(2)*(1-c)-unit_axis(0)*s;
    rotation_matrix(2,0) = unit_axis(0)*unit_axis(2)*(1-c)-unit_axis(1)*s;
    rotation_matrix(2,1) = unit_axis(1)*unit_axis(2)*(1-c)+unit_axis(0)*s;
    rotation_matrix(2,2) = unit_axis(2)*unit_axis(2)+c*(1-unit_axis(2)*unit_axis(2));

    Rotate(rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateX(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> x_rotation_matrix=identity_matrix<double>(SPACE_DIM);

    x_rotation_matrix(1,1) = cos(theta);
    x_rotation_matrix(1,2) = sin(theta);
    x_rotation_matrix(2,1) = -sin(theta);
    x_rotation_matrix(2,2) = cos(theta);
    Rotate(x_rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateY(const double theta)
{
    if (SPACE_DIM != 3)
    {
        EXCEPTION("This rotation is only valid in 3D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> y_rotation_matrix=identity_matrix<double>(SPACE_DIM);

    y_rotation_matrix(0,0) = cos(theta);
    y_rotation_matrix(0,2) = -sin(theta);
    y_rotation_matrix(2,0) = sin(theta);
    y_rotation_matrix(2,2) = cos(theta);


    Rotate(y_rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RotateZ(const double theta)
{
    if (SPACE_DIM < 2)
    {
        EXCEPTION("This rotation is not valid in less than 2D");
    }
    c_matrix<double , SPACE_DIM, SPACE_DIM> z_rotation_matrix=identity_matrix<double>(SPACE_DIM);


    z_rotation_matrix(0,0) = cos(theta);
    z_rotation_matrix(0,1) = sin(theta);
    z_rotation_matrix(1,0) = -sin(theta);
    z_rotation_matrix(1,1) = cos(theta);

    Rotate(z_rotation_matrix);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Rotate(double theta)
{
    RotateZ(theta);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes()
{
    RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();

    // Working from the back, each node is swapped with a random node that precedes it in the array
    for (unsigned index=this->mNodes.size()-1; index>0; index--)
    {
        unsigned  other=p_rng->randMod(index+1); //includes the possibility of rolling "index"
        // Swap index and other
        Node<SPACE_DIM> *temp=this->mNodes[index];
        this->mNodes[index]=this->mNodes[other];
        this->mNodes[other]=temp;
    }

    // Update indices
    for (unsigned index=0; index<this->mNodes.size(); index++)
    {
        this->mNodes[index]->SetIndex(index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodes(std::vector<unsigned>& perm)
{
    // Let's not do this if there are any deleted nodes
    assert( this->GetNumAllNodes() == this->GetNumNodes());

    assert(perm.size() == this->mNodes.size());

    // Copy the node pointers
    std::vector< Node<SPACE_DIM>* > copy_m_nodes;
    copy_m_nodes.assign(this->mNodes.begin(), this->mNodes.end());

    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        assert(perm[i] < this->mNodes.size());
        this->mNodes[ perm[i] ] = copy_m_nodes[i];
    }

    // Update indices
    for (unsigned index=0; index<this->mNodes.size(); index++)
    {
        this->mNodes[index]->SetIndex(index);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::PermuteNodesWithMetisBinaries(unsigned numProcs)
{
    assert(ELEMENT_DIM==2 || ELEMENT_DIM==3);
    assert( this->GetNumAllElements() == this->GetNumElements());
    assert( this->GetNumAllNodes() == this->GetNumNodes());

    // Open a file for the elements
    OutputFileHandler handler("");

    // Filenames
    std::string basename = "metis.mesh";
    std::stringstream output_file;
    output_file << basename << ".npart." << numProcs;
    std::string nodes_per_proc_file = basename + ".nodesperproc";

    // Only the master process should do IO and call METIS
    if (handler.IsMaster())
    {
        out_stream metis_file = handler.OpenOutputFile(basename);

        (*metis_file)<<this->GetNumElements()<<"\t";
        if (ELEMENT_DIM==2)
        {
            (*metis_file)<<1<<"\n"; //1 is Metis speak for triangles
        }
        else
        {
            (*metis_file)<<2<<"\n"; //2 is Metis speak for tetrahedra
        }

        for (unsigned i=0; i<this->GetNumElements(); i++)
        {
            for (unsigned j=0; j<ELEMENT_DIM+1; j++)
            {
                //Note the +1 since Metis wants meshes indexed from 1
                (*metis_file)<<this->mElements[i]->GetNode(j)->GetIndex() + 1<<"\t";
            }
            (*metis_file)<<"\n";
        }
        metis_file->close();

        /*
         *  Call METIS binary to perform the partitioning.
         *  It will output a file called metis.mesh.npart.numProcs
         */
        std::stringstream permute_command;
        permute_command <<  "./bin/partdmesh "
                        <<  handler.GetOutputDirectoryFullPath("")
                        <<  basename << " "
                        <<  numProcs
                        <<  " > /dev/null";

        // METIS doesn't return 0 after a successful execution
        IGNORE_RET(system, permute_command.str());

        /*
         *  Create a file with the number of nodes per partition
         */
        // Make sure it doesn't exist, since values will be appended with >>
        std::stringstream clear_command;
        clear_command << "rm -f "
                      << handler.GetOutputDirectoryFullPath("")
                      << nodes_per_proc_file
                      << " > /dev/null";
        EXPECT0(system, clear_command.str());

        // Loop over the partition number (i.e. processor number) and count how many nodes
        for (unsigned proc_index=0; proc_index<numProcs; proc_index++)
        {
            std::stringstream count_command;
            count_command << "grep "
                          << proc_index << " "
                          << handler.GetOutputDirectoryFullPath("")
                          << output_file.str()
                          << " | wc -l >> "
                          << handler.GetOutputDirectoryFullPath("")
                          << nodes_per_proc_file;

            EXPECT0(system, count_command.str());
        }

    }

    // Wait for the permutation to be available
    PetscTools::Barrier();

    /*
     *  Read partition file back into a vector.
     */
    std::vector<unsigned> partition(this->GetNumNodes());
    std::vector<unsigned> offset(numProcs,0u);

    std::ifstream partition_stream;
    std::string full_path = handler.GetOutputDirectoryFullPath("")
                            + output_file.str();

    partition_stream.open(full_path.c_str());
    assert(partition_stream.is_open());

    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        unsigned part_read;

        partition_stream >> part_read;

        partition[node_index] = part_read;
        for (unsigned proc=part_read+1; proc<numProcs; proc++)
        {
            offset[proc]++;
        }
    }
    partition_stream.close();

    /*
     *  Create the permutation vector based on Metis output
     */
    std::vector<unsigned> permutation(this->GetNumNodes(), UINT_MAX);
    std::vector<unsigned> count(numProcs,0u);

    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        unsigned part = partition[node_index];
        // Permutation defined like: new index for node node_index is "offset[part] + count[part]"
        permutation [ node_index ] = offset[part] + count[part];

        count[part]++;
    }

    PermuteNodes(permutation);

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructLinearMesh(unsigned width)
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

    RefreshJacobianCachedData();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructRectangularMesh(unsigned width, unsigned height, bool stagger)
{
    assert(SPACE_DIM == 2);
    assert(ELEMENT_DIM == 2);

    //Construct the nodes
    unsigned node_index=0;
    for (int j=(int)height; j>=0; j--) //j must be signed for this loop to terminate
    {
        for (unsigned i=0; i<width+1; i++)
        {
            bool is_boundary=false;
            if (i==0 || j==0 || i==width || j==(int)height)
            {
                is_boundary=true;
            }
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
        nodes.push_back(this->mNodes[i]);
        nodes.push_back(this->mNodes[i+1]);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Right
    for (unsigned i=1; i<height+1; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[(width+1)*i-1]);
        nodes.push_back(this->mNodes[(width+1)*(i+1)-1]);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Bottom
    for (unsigned i=0; i<width; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[height*(width+1)+i+1]);
        nodes.push_back(this->mNodes[height*(width+1)+i]);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }
    //Left
    for (unsigned i=0; i<height; i++)
    {
        std::vector<Node<SPACE_DIM>*> nodes;
        nodes.push_back(this->mNodes[(width+1)*(i+1)]);
        nodes.push_back(this->mNodes[(width+1)*(i)]);
        this->mBoundaryElements.push_back(new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(belem_index++,nodes));
    }

    //Construct the elements
    unsigned elem_index = 0;
    for (unsigned j=0; j<height; j++)
    {
        for (unsigned i=0; i<width; i++)
        {
            unsigned parity=(i+j)%2;
            std::vector<Node<SPACE_DIM>*> upper_nodes;
            upper_nodes.push_back(this->mNodes[j*(width+1)+i]);
            upper_nodes.push_back(this->mNodes[j*(width+1)+i+1]);
            if (stagger==false  || parity == 0)
            {
                upper_nodes.push_back(this->mNodes[(j+1)*(width+1)+i+1]);
            }
            else
            {
                upper_nodes.push_back(this->mNodes[(j+1)*(width+1)+i]);
            }
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,upper_nodes));
            std::vector<Node<SPACE_DIM>*> lower_nodes;
            lower_nodes.push_back(this->mNodes[(j+1)*(width+1)+i+1]);
            lower_nodes.push_back(this->mNodes[(j+1)*(width+1)+i]);
            if (stagger==false  ||parity == 0)
            {
                lower_nodes.push_back(this->mNodes[j*(width+1)+i]);
            }
            else
            {
                lower_nodes.push_back(this->mNodes[j*(width+1)+i+1]);
            }
            this->mElements.push_back(new Element<ELEMENT_DIM,SPACE_DIM>(elem_index++,lower_nodes));
        }
    }

    RefreshJacobianCachedData();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndex(ChastePoint<SPACE_DIM> testPoint, bool strict, std::set<unsigned> testElements)
{
    for (std::set<unsigned>::iterator iter=testElements.begin(); iter!=testElements.end(); iter++)
    {
        assert(*iter<this->GetNumElements());
        ///\todo What if the element is deleted?
        if (this->mElements[*iter]->IncludesPoint(testPoint, strict))
        {
            return *iter;
        }
    }

    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        ///\todo What if the element is deleted?
        if (this->mElements[i]->IncludesPoint(testPoint, strict))
        {
            return i;
        }
    }

    // If it's in none of the elements, then throw
    std::stringstream ss; 
    ss << "Point [";
    for(unsigned j=0; (int)j<(int)SPACE_DIM-1; j++)
    {
        ss << testPoint[j] << ",";
    }
    ss << testPoint[SPACE_DIM-1] << "] is not in mesh - all elements tested";
    EXCEPTION(ss.str());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndexWithInitialGuess(ChastePoint<SPACE_DIM> testPoint, unsigned startingElementGuess, bool strict)
{
    assert(startingElementGuess<this->GetNumElements());
    
    // let m=startingElementGuess, N=num_elem-1
    // We search from in this order: m, m+1, m+2, .. , N, 0, 1, .., m-1.

    unsigned i = startingElementGuess;
    bool reached_end = false;

    while(!reached_end)
    {
        ///\todo What if the element is deleted?
        if (this->mElements[i]->IncludesPoint(testPoint, strict))
        {
            return i;
        }

        // increment
        i++;
        if(i==this->GetNumElements())
        {
            i=0;
        }

        // back to the beginning yet?
        if(i==startingElementGuess)
        {
            reached_end = true;
        }
    }

    // If it's in none of the elements, then throw
    std::stringstream ss; 
    ss << "Point [";    
    for(unsigned j=0; (int)j<(int)SPACE_DIM-1; j++)
    {
        ss << testPoint[j] << ",";
    }
    ss << testPoint[SPACE_DIM-1] << "] is not in mesh - all elements tested";
    EXCEPTION(ss.str());
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetNearestElementIndex(ChastePoint<SPACE_DIM> testPoint)
{
    ///\todo This ought to return a set of all elements that contain the point (if the point is a node in the mesh then it's contained in multiple elements)
    ///\todo Polling every element is unnecessary.  We ought to start from a likely place and hill climb

    double max_min_weight = -INFINITY;
    unsigned closest_index = 0;
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        ///\todo What if the element is deleted?
        c_vector<double, ELEMENT_DIM+1> weight=this->mElements[i]->CalculateInterpolationWeights(testPoint);
        double neg_weight_sum=0.0;
        for (unsigned j=0; j<=ELEMENT_DIM; j++)
        {
            if (weight[j] < 0.0)
            {
                neg_weight_sum += weight[j];
            }
        }
        if (neg_weight_sum > max_min_weight)
        {
            max_min_weight = neg_weight_sum;
            closest_index=i;
        }

    }
    return closest_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(ChastePoint<SPACE_DIM> testPoint)
{
    std::vector<unsigned> element_indices;
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        ///\todo What if the element is deleted?
        if (this->mElements[i]->IncludesPoint(testPoint))
        {
            element_indices.push_back(i);
        }
    }
    return element_indices;
}

//template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SetElementOwnerships(unsigned lo, unsigned hi)
//{
//    assert(hi>=lo);
//    for (unsigned element_index=0; element_index<this->mElements.size(); element_index++)
//    {
//        Element<ELEMENT_DIM, SPACE_DIM>* p_element=this->mElements[element_index];
//        p_element->SetOwnership(false);
//        for (unsigned local_node_index=0; local_node_index< p_element->GetNumNodes(); local_node_index++)
//        {
//            unsigned global_node_index = p_element->GetNodeGlobalIndex(local_node_index);
//            if (lo<=global_node_index && global_node_index<hi)
//            {
//                p_element->SetOwnership(true);
//                break;
//            }
//        }
//
//    }
//}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ConstructCuboid(unsigned width,
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

    RefreshJacobianCachedData();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // three loops, just like the destructor. note we don't delete boundary nodes.
    for (unsigned i=0; i<this->mBoundaryElements.size(); i++)
    {
        delete this->mBoundaryElements[i];
    }
    for (unsigned i=0; i<this->mElements.size(); i++)
    {
        delete this->mElements[i];
    }
    for (unsigned i=0; i<this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }

    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryElements.clear();
    this->mBoundaryNodes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::CalculateBoundaryOfFlaggedRegion()
{
    // A set of nodes which lie on the face, size 3 in 2D, size 4 in 3D
    typedef std::set<unsigned> FaceNodes;

    // Face maps to true the first time it is encountered, and false subsequent
    // times. Thus, faces mapping to true at the end are boundary faces
    std::map<FaceNodes,bool> face_on_boundary;

    // Loop over all elements
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        if (iter->IsFlagged())
        {
            // To get faces, initially start with all nodes
            std::set<unsigned> all_nodes;
            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                all_nodes.insert(iter->GetNodeGlobalIndex(i));
            }

            // Remove one node in turn to obtain each face
            for (unsigned i=0; i<ELEMENT_DIM+1; i++)
            {
                FaceNodes face_nodes = all_nodes;
                face_nodes.erase(iter->GetNodeGlobalIndex(i));

                // Search the map of faces to see if it contains this face
                std::map<FaceNodes,bool>::iterator it = face_on_boundary.find(face_nodes);

                if (it == face_on_boundary.end())
                {
                    // Face not found, add and assume on boundary
                    face_on_boundary[face_nodes]=true;
                }
                else
                {
                    // Face found in map, so not on boundary
                    it->second = false;
                }
            }
        }
    }

    // Boundary nodes to be returned
    std::set<unsigned> boundary_of_flagged_region;

    // Get all faces in the map
    std::map<FaceNodes,bool>::iterator it=face_on_boundary.begin();
    while (it!=face_on_boundary.end())
    {
        // If the face maps to true it is on the boundary
        if (it->second==true)
        {
            // Get all nodes in the face and put in set to be returned
            boundary_of_flagged_region.insert(it->first.begin(),it->first.end());
        }
        it++;
    }

    return boundary_of_flagged_region;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetAngleBetweenNodes(unsigned indexA, unsigned indexB)
{
    assert(SPACE_DIM == 2);
    assert(SPACE_DIM == ELEMENT_DIM);

    double x_diff = this->mNodes[indexB]->rGetLocation()[0] - this->mNodes[indexA]->rGetLocation()[0];
    double y_diff = this->mNodes[indexB]->rGetLocation()[1] - this->mNodes[indexA]->rGetLocation()[1];

    if (x_diff==0)
    {
        if (y_diff>0)
        {
            return M_PI/2.0;
        }
        else if (y_diff<0)
        {
            return -M_PI/2.0;
        }
        else
        {
            EXCEPTION("Tried to compute polar angle of (0,0)");
        }
    }

    double angle = atan2(y_diff,x_diff);
    return angle;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::UnflagAllElements()
{
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        iter->Unflag();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::FlagElementsNotContainingNodes(std::set<unsigned> nodesList)
{
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        bool found_node = false;

        for (unsigned i=0; i<iter->GetNumNodes(); i++)
        {
            unsigned node_index = iter->GetNodeGlobalIndex(i);

            std::set<unsigned>::iterator set_iter = nodesList.find(node_index);
            if (set_iter != nodesList.end())
            {
                found_node = true;
            }
        }

        if (!found_node)
        {
            iter->Flag();
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
//                          edge iterator class                             //
//////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeA()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeALocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Node<SPACE_DIM>* TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::GetNodeB()
{
    assert((*this) != mrMesh.EdgesEnd());
    Element<ELEMENT_DIM,SPACE_DIM>* p_element = mrMesh.GetElement(mElemIndex);
    return p_element->GetNode(mNodeBLocalIndex);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator!=(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& rOther)
{
    return (mElemIndex != rOther.mElemIndex ||
            mNodeALocalIndex != rOther.mNodeALocalIndex ||
            mNodeBLocalIndex != rOther.mNodeBLocalIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator& TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::operator++()
{
    std::set<unsigned> current_node_pair;
    std::set<std::set<unsigned> >::iterator set_iter;

    do
    {
        // Advance to the next edge in the mesh.
        // Node indices are incremented modulo #nodes_per_elem
        mNodeBLocalIndex = (mNodeBLocalIndex + 1) % (ELEMENT_DIM+1);
        if (mNodeBLocalIndex == mNodeALocalIndex)
        {
            mNodeALocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM+1);
            mNodeBLocalIndex = (mNodeALocalIndex + 1) % (ELEMENT_DIM+1);
        }

        if (mNodeALocalIndex == 0 && mNodeBLocalIndex == 1) // advance to next element...
        {
            mElemIndex++;
            // ...skipping deleted ones
            while (mElemIndex!=mrMesh.GetNumAllElements() && mrMesh.GetElement(mElemIndex)->IsDeleted())
            {
                mElemIndex++;
            }
        }

        if (mElemIndex != mrMesh.GetNumAllElements())
        {
            unsigned node_a_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeALocalIndex);
            unsigned node_b_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeBLocalIndex);

            // Check we haven't seen it before
            current_node_pair.clear();
            current_node_pair.insert(node_a_global_index);
            current_node_pair.insert(node_b_global_index);
            set_iter = mEdgesVisited.find(current_node_pair);
        }
    }
    while (*this != mrMesh.EdgesEnd() && set_iter != mEdgesVisited.end());
    mEdgesVisited.insert(current_node_pair);

    return (*this);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator::EdgeIterator(TetrahedralMesh& rMesh, unsigned elemIndex)
    : mrMesh(rMesh),
      mElemIndex(elemIndex),
      mNodeALocalIndex(0),
      mNodeBLocalIndex(1)
{
    if (elemIndex==mrMesh.GetNumAllElements())
    {
        return;
    }

    mEdgesVisited.clear();

    // add the current node pair to the store
    std::set<unsigned> current_node_pair;
    unsigned node_a_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeALocalIndex);
    unsigned node_b_global_index = mrMesh.GetElement(mElemIndex)->GetNodeGlobalIndex(mNodeBLocalIndex);
    current_node_pair.insert(node_a_global_index);
    current_node_pair.insert(node_b_global_index);

    mEdgesVisited.insert(current_node_pair);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesBegin()
{
    unsigned first_element_index=0;
    while (first_element_index!=this->GetNumAllElements() && this->GetElement(first_element_index)->IsDeleted())
    {
        first_element_index++;
    }
    return EdgeIterator(*this, first_element_index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgeIterator TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::EdgesEnd()
{
    return EdgeIterator(*this, this->GetNumAllElements());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshMesh()
{
    RefreshJacobianCachedData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    assert(index < this->mBoundaryElements.size() );
    return index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianCachedData()
{
    // Make sure we have enough space
    this->mElementJacobians.resize(this->GetNumAllElements());
    this->mElementInverseJacobians.resize(this->GetNumAllElements());
    
    if (ELEMENT_DIM < SPACE_DIM)
    {
        this->mElementWeightedDirections.resize(this->GetNumAllElements());
    }

    this->mBoundaryElementWeightedDirections.resize(this->GetNumAllBoundaryElements());

    this->mElementJacobianDeterminants.resize(this->GetNumAllElements());
    this->mBoundaryElementJacobianDeterminants.resize(this->GetNumAllBoundaryElements());

    // Update caches
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
         iter != this->GetElementIteratorEnd();
         ++iter)
    {
        unsigned index = iter->GetIndex();
        iter->CalculateInverseJacobian(this->mElementJacobians[index], this->mElementJacobianDeterminants[index], this->mElementInverseJacobians[index]);
    }
        
    if (ELEMENT_DIM < SPACE_DIM)
    {
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = this->GetElementIteratorBegin();
             iter != this->GetElementIteratorEnd();
             ++iter)
        {
             unsigned index = iter->GetIndex();
             iter->CalculateWeightedDirection(this->mElementWeightedDirections[index], this->mElementJacobianDeterminants[index]);
        }
    }

    for ( typename TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator itb = this->GetBoundaryElementIteratorBegin();
          itb != this->GetBoundaryElementIteratorEnd();
          itb++)
    {
        unsigned index = (*itb)->GetIndex();
        (*itb)->CalculateWeightedDirection(this->mBoundaryElementWeightedDirections[index], this->mBoundaryElementJacobianDeterminants[index]);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double& rJacobianDeterminant) const
{
    assert(ELEMENT_DIM <= SPACE_DIM);
    assert(elementIndex < this->mElementJacobians.size());
    rJacobian = this->mElementJacobians[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    assert(ELEMENT_DIM <= SPACE_DIM);
    assert(elementIndex < this->mElementInverseJacobians.size());
    rInverseJacobian = this->mElementInverseJacobians[elementIndex];
    rJacobian = this->mElementJacobians[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    assert(ELEMENT_DIM < SPACE_DIM);
    assert(elementIndex < this->mElementWeightedDirections.size());
    rWeightedDirection = this->mElementWeightedDirections[elementIndex];
    rJacobianDeterminant = this->mElementJacobianDeterminants[elementIndex];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    assert(elementIndex < this->mBoundaryElementWeightedDirections.size());
    rWeightedDirection = this->mBoundaryElementWeightedDirections[elementIndex];
    rJacobianDeterminant = this->mBoundaryElementJacobianDeterminants[elementIndex];
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class TetrahedralMesh<1,1>;
template class TetrahedralMesh<1,2>;
template class TetrahedralMesh<1,3>;
template class TetrahedralMesh<2,2>;
template class TetrahedralMesh<2,3>;
template class TetrahedralMesh<3,3>;
