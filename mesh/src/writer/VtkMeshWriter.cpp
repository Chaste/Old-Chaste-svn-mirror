/*

Copyright (C) University of Oxford, 2005-2011

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

#include "VtkMeshWriter.hpp"
#include "DistributedTetrahedralMesh.hpp"

#ifdef CHASTE_VTK
///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshWriter<ELEMENT_DIM, SPACE_DIM>::VtkMeshWriter(const std::string& rDirectory,
                     const std::string& rBaseName,
                     const bool& rCleanDirectory)
    : AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, rCleanDirectory),
      mWriteParallelFiles(false)
{
    this->mIndexFromZero = true;

    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::~VtkMeshWriter()
{
    mpVtkUnstructedMesh->Delete(); // Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::MakeVtkMesh()
{
    assert(SPACE_DIM==3 || SPACE_DIM == 2);
    assert(SPACE_DIM==ELEMENT_DIM);
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->GetNextNode(); //this->mNodeData[item_num];
        if (SPACE_DIM==2)
        {
            current_item.push_back(0.0);//For z-coordinate
        }
        assert(current_item.size() == 3);
        p_pts->InsertPoint(item_num, current_item[0], current_item[1], current_item[2]);
    }

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); //Reference counted
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        std::vector<unsigned> current_element = this->GetNextElement().NodeIndices; // this->mElementData[item_num];
        assert(current_element.size() == ELEMENT_DIM + 1);
        vtkCell* p_cell=NULL;
        if (SPACE_DIM == 3)
        {
            p_cell = vtkTetra::New();
        }
        if (SPACE_DIM == 2)
        {
            p_cell = vtkTriangle::New();
        }
        vtkIdList* p_cell_id_list = p_cell->GetPointIds();
        for (unsigned j = 0; j < ELEMENT_DIM+1; ++j)
        {
            p_cell_id_list->SetId(j, current_element[j]);
        }
        mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
        p_cell->Delete(); //Reference counted
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddProvenance(std::string fileName)
{
    std::string comment = "<!-- " + ChasteBuildInfo::GetProvenanceString() + "-->";
    
    out_stream p_vtu_file = this->mpOutputFileHandler->OpenOutputFile(fileName, std::ios::out | std::ios::app);

    *p_vtu_file << "\n" << comment << "\n";    
    p_vtu_file->close();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteFiles()
{
    // Using separate scope here to make sure file is properly closed before re-opening it to add provenance info.
    {
        MakeVtkMesh();
        assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
        vtkXMLUnstructuredGridWriter* p_writer = vtkXMLUnstructuredGridWriter::New();
        p_writer->SetInput(mpVtkUnstructedMesh);
        //Uninitialised stuff arises (see #1079), but you can remove
        //valgrind problems by removing compression:
        // **** REMOVE WITH CAUTION *****
        p_writer->SetCompressor(NULL);
        // **** REMOVE WITH CAUTION *****
        std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName+".vtu";
        p_writer->SetFileName(vtk_file_name.c_str());
        //p_writer->PrintSelf(std::cout, vtkIndent());
        p_writer->Write();
        p_writer->Delete(); //Reference counted
    }    

    AddProvenance(this->mBaseName + ".vtu");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddCellData(std::string dataName, std::vector<c_vector<double, SPACE_DIM> > dataPayload)
{
    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(3);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            p_vectors->InsertNextValue(dataPayload[i][j]);
        }
        //When SPACE_DIM<3, then pad 
        for (unsigned j=SPACE_DIM; j<3; j++)
        {
            p_vectors->InsertNextValue(0.0);
        }
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); //Reference counted

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::AddPointData(std::string dataName, std::vector<c_vector<double, SPACE_DIM> > dataPayload)
{
    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(3);
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            p_vectors->InsertNextValue(dataPayload[i][j]);
        }
        //When SPACE_DIM<3, then pad 
        for (unsigned j=SPACE_DIM; j<3; j++)
        {
            p_vectors->InsertNextValue(0.0);
        }
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_vectors);
    p_vectors->Delete(); //Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM,SPACE_DIM>::SetParallelFiles()
{
    mWriteParallelFiles = true;
}

///\todo #1322 Mesh should be const
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(
      AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
      bool keepOriginalElementIndexing)
{
    //Have we got a parallel mesh?
    this->mpDistributedMesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* >(&rMesh);
    
    if ( PetscTools::IsSequential() || this->mpDistributedMesh == NULL )
    {
        AbstractTetrahedralMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteFilesUsingMesh( rMesh,keepOriginalElementIndexing );
    }
    else
    {
        //Make the local mesh into a VtkMesh
        assert(SPACE_DIM==3 || SPACE_DIM == 2);
        assert(SPACE_DIM==ELEMENT_DIM);
        vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
        p_pts->GetData()->SetName("Vertex positions");
        
        //Node index that we are writing to VTK (index into mNodes and mHaloNodes as if they were concatenated)
        unsigned index = 0;
        //Map a global node index into a local index (into mNodes and mHaloNodes as if they were concatenated)
        std::map<unsigned, unsigned> global_to_local_index_map;
        
        // Owned nodes
        for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
             node_iter != rMesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double, SPACE_DIM> current_item = node_iter->rGetLocation();
            p_pts->InsertNextPoint(current_item[0], current_item[1], (SPACE_DIM==3)?current_item[2]:0.0);
            global_to_local_index_map[node_iter->GetIndex()] = index;
            index++;
        }
        
        // Halo nodes
        for(typename DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::HaloNodeIterator halo_iter=this->mpDistributedMesh->GetHaloNodeIteratorBegin(); 
                halo_iter != this->mpDistributedMesh->GetHaloNodeIteratorEnd();
                ++halo_iter)
        {
            c_vector<double, SPACE_DIM> current_item = (*halo_iter)->rGetLocation();
            p_pts->InsertNextPoint(current_item[0], current_item[1], (SPACE_DIM==3)?current_item[2]:0.0);
            global_to_local_index_map[(*halo_iter)->GetIndex()] = index;
            index++;
        }       
        
        mpVtkUnstructedMesh->SetPoints(p_pts);
        p_pts->Delete(); //Reference counted
        
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
             elem_iter != rMesh.GetElementIteratorEnd();
             ++elem_iter)
        {
        
            vtkCell* p_cell=NULL;
            if (SPACE_DIM == 3)
            {
                p_cell = vtkTetra::New();
            }
            if (SPACE_DIM == 2)
            {
                p_cell = vtkTriangle::New();
            }
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            for (unsigned j = 0; j < ELEMENT_DIM+1; ++j)
            {
                unsigned global_node_index = elem_iter->GetNodeGlobalIndex(j);
                p_cell_id_list->SetId(j, global_to_local_index_map[global_node_index]);
            }
            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); //Reference counted
        }
        //This block is to guard the mesh writers (vtkXMLPUnstructuredGridWriter) so that they
        //go out of scope, flush buffers and close files 
        {
            assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
            vtkXMLPUnstructuredGridWriter* p_writer = vtkXMLPUnstructuredGridWriter::New();

            p_writer->SetDataModeToBinary();
 
            p_writer->SetNumberOfPieces(PetscTools::GetNumProcs());
            //p_writer->SetGhostLevel(-1);
            p_writer->SetStartPiece(PetscTools::GetMyRank());
            p_writer->SetEndPiece(PetscTools::GetMyRank());


            p_writer->SetInput(mpVtkUnstructedMesh);
            //Uninitialised stuff arises (see #1079), but you can remove
            //valgrind problems by removing compression:
            // **** REMOVE WITH CAUTION *****
            p_writer->SetCompressor(NULL);
            // **** REMOVE WITH CAUTION *****
            std::string pvtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName+ ".pvtu";
            p_writer->SetFileName(pvtk_file_name.c_str());
            //p_writer->PrintSelf(std::cout, vtkIndent());
            p_writer->Write();
            p_writer->Delete(); //Reference counted
        }
        
        //Add provenance to the individual files
        std::stringstream filepath;
        filepath << this->mBaseName << "_" << PetscTools::GetMyRank() << ".vtu";
        AddProvenance(filepath.str());  
        /// Add to the main file \todo #1494 Do we need a barrier?
        if (PetscTools::AmMaster())
        {
            AddProvenance(this->mBaseName+ ".pvtu");
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class VtkMeshWriter<1,1>;
template class VtkMeshWriter<1,2>;
template class VtkMeshWriter<1,3>;
template class VtkMeshWriter<2,2>; //Actually used
template class VtkMeshWriter<2,3>;
template class VtkMeshWriter<3,3>; //Actually used

#endif //CHASTE_VTK
