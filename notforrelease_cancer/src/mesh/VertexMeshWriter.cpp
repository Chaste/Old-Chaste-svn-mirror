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

#include "VertexMeshWriter.hpp"
#include "NodeMap.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::VertexMeshWriter(const std::string& rDirectory,
                                                           const std::string& rBaseName,
                                                           const bool clearOutputDir)
    : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir)
{
#ifdef CHASTE_VTK
     // Dubious, since we shouldn't yet know what any details of the mesh are.
     mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif //CHASTE_VTK
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~VertexMeshWriter()
{
 #ifdef CHASTE_VTK
     // Dubious, since we shouldn't yet know what any details of the mesh are.
     mpVtkUnstructedMesh->Delete(); // Reference counted
#endif //CHASTE_VTK
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp)
{
#ifdef CHASTE_VTK
    //Make the Vtk mesh
    assert(SPACE_DIM==3 || SPACE_DIM == 2);
    vtkPoints *p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned node_num=0; node_num<rMesh.GetNumNodes(); node_num++)
    {
        c_vector<double, SPACE_DIM> position = rMesh.GetNode(node_num)->rGetLocation();
        if (SPACE_DIM==2)
        {
            p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
        }
        else
        {
            p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
        }
    }

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); //Reference counted
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator iter = rMesh.GetElementIteratorBegin();
             iter != rMesh.GetElementIteratorEnd();
             ++iter)
    {
        vtkCell *p_cell;
        if (SPACE_DIM == 2)
        {
            p_cell = vtkPolygon::New();
        }
        else
        {
            p_cell = vtkConvexPointSet::New();
        }
        vtkIdList *p_cell_id_list = p_cell->GetPointIds();
        p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());
        for (unsigned j = 0; j < iter->GetNumNodes(); ++j)
        {
            p_cell_id_list->SetId(j, iter->GetNodeGlobalIndex(j));
        }
        mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
        p_cell->Delete(); //Reference counted
    }

    //Vtk mesh is now made
    assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
    vtkXMLUnstructuredGridWriter *p_writer = vtkXMLUnstructuredGridWriter::New();
    p_writer->SetInput(mpVtkUnstructedMesh);
    //Uninitialised stuff arises (see #1079), but you can remove
    //valgrind problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    p_writer->SetCompressor(NULL);
    // **** REMOVE WITH CAUTION *****

    std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName;
    if (stamp != "")
    {
        vtk_file_name += "_" + stamp;
    }
    vtk_file_name += ".vtu";

    p_writer->SetFileName(vtk_file_name.c_str());
    //p_writer->PrintSelf(std::cout, vtkIndent());
    p_writer->Write();
    p_writer->Delete(); //Reference counted
#endif //CHASTE_VTK
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray *p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData *p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); //Reference counted
#endif //CHASTE_VTK
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray *p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkPointData *p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); //Reference counted
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(const VertexMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    NodeMap node_map(rMesh.GetNumAllNodes());
    unsigned new_index = 0;
    for (unsigned i=0; i<rMesh.GetNumAllNodes(); i++)
    {
        Node<SPACE_DIM> *p_node = rMesh.GetNode(i);

        if (p_node->IsDeleted() == false)
        {
            std::vector<double> coords(SPACE_DIM);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                coords[j] = p_node->GetPoint()[j];
            }
            this->SetNextNode(coords);
            node_map.SetNewIndex(i, new_index++);
        }
        else
        {
            node_map.SetDeleted(i);
        }
    }
    assert(new_index==rMesh.GetNumNodes());

    for (unsigned i=0; i<rMesh.GetNumAllElements(); i++) // can't use an element interator as the mesh is const
    {
        VertexElement<ELEMENT_DIM, SPACE_DIM> *p_element = rMesh.GetElement(i);
        if (!p_element->IsDeleted())
        {
            std::vector<unsigned> indices(p_element->GetNumNodes()+1);
            indices[0] = p_element->GetNumNodes(); // first entry is the num nodes contained in this element

            for (unsigned j=1; j<indices.size(); j++)
            {
                unsigned old_index = p_element->GetNodeGlobalIndex(j-1);
                indices[j] = node_map.GetNewIndex(old_index);
            }
            this->SetNextElement(indices);
        }
    }

    WriteFiles();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    std::string comment = "#Generated by Chaste vertex mesh file writer";

    // Write node file
    std::string node_file_name = this->mBaseName + ".node";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    // Write the node header
    unsigned num_attr = 0;
    unsigned max_bdy_marker = 0;
    unsigned num_nodes = this->GetNumNodes();

    *p_node_file << num_nodes << "\t";
    *p_node_file << SPACE_DIM << "\t";
    *p_node_file << num_attr << "\t";
    *p_node_file << max_bdy_marker << "\n";
    *p_node_file << std::setprecision(6);

    // Write each node's data
    unsigned default_marker = 0;
    for (unsigned item_num=0; item_num<num_nodes; item_num++)
    {
        std::vector<double> current_item = this->mNodeData[item_num];
        *p_node_file << item_num;
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *p_node_file << "\t" << current_item[i];
        }
        *p_node_file << "\t" << default_marker << "\n";

    }
    *p_node_file << comment << "\n";
    p_node_file->close();

    // Write element file
    std::string element_file_name = this->mBaseName + ".cell";
    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

    // Write the element header
    unsigned num_elements = this->GetNumElements();

    *p_element_file << num_elements << "\t";
    *p_element_file << num_attr << "\n";

    // Write each element's data
    /// \todo need to think about how best to do this in 3D (see #866)
    for (unsigned item_num=0; item_num<num_elements; item_num++)
    {
        std::vector<unsigned> current_item = this->mElementData[item_num];
        *p_element_file << item_num;
        for (unsigned i=0; i<current_item.size(); i++)
        {
            *p_element_file << "\t" << current_item[i];
        }
        *p_element_file << "\n";
    }
    *p_element_file << comment << "\n";
    p_element_file->close();
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexMeshWriter<1,1>;
template class VertexMeshWriter<1,2>;
template class VertexMeshWriter<1,3>;
template class VertexMeshWriter<2,2>;
template class VertexMeshWriter<2,3>;
template class VertexMeshWriter<3,3>;
