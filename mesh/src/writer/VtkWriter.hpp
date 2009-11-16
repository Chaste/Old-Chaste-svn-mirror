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

#ifndef VTKWRITER_HPP_
#define VTKWRITER_HPP_

#ifdef CHASTE_VTK
//Requires  "sudo aptitude install libvtk5-dev" or similar
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkDataCompressor.h>
#include "AbstractTetrahedralMeshWriter.hpp"

/**
 *  VtkWriter
 *
 *  Writes a mesh in VTK .vtu format (that's an XML-based, data compressed unstructured mesh)
 *
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VtkWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>
{

//Requires  "sudo aptitude install libvtk5-dev" or similar

private:
    vtkUnstructuredGrid* mpVtkUnstructedMesh;

    void MakeVtkMesh();
public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    VtkWriter(const std::string& rDirectory, const std::string& rBaseName, const bool& rCleanDirectory=true);

    /**
     * Write mesh data to files.
     */
    void WriteFiles();

    void AddCellData(std::string name, std::vector<double> data);
    void AddPointData(std::string name, std::vector<double> data);

    /**
     * Destructor.
     */
    virtual ~VtkWriter();
};


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkWriter<ELEMENT_DIM, SPACE_DIM>::VtkWriter(const std::string& rDirectory,
                     const std::string& rBaseName,
                     const bool& rCleanDirectory)
    : AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, rCleanDirectory)
{
    this->mIndexFromZero = true;

    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VtkWriter<ELEMENT_DIM,SPACE_DIM>::~VtkWriter()
{
    mpVtkUnstructedMesh->Delete(); // Reference counted
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkWriter<ELEMENT_DIM,SPACE_DIM>::MakeVtkMesh()
{
    assert(SPACE_DIM==3 || SPACE_DIM == 2);
    assert(SPACE_DIM==ELEMENT_DIM);
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->mNodeData[item_num];
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
        std::vector<unsigned> current_element = this->mElementData[item_num];
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
void VtkWriter<ELEMENT_DIM,SPACE_DIM>::WriteFiles()
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VtkWriter<ELEMENT_DIM,SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
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
void VtkWriter<ELEMENT_DIM,SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
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
#endif //CHASTE_VTK

#endif /*VTKWRITER_HPP_*/
