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
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "AbstractMeshWriter.hpp"



/**
 *  VtkWriter
 * 
 *  Writes a mesh in VTK  vtu format (that's an XMl-based, data compressed unstructured mesh)
 * 
 */
class VtkWriter : public AbstractMeshWriter<3,3>
{

#ifdef CHASTE_VTK 
//Requires  "sudo aptitude install libvtk5-dev" or similar 
    
private:
    vtkUnstructuredGrid* mpVtkUnstructedMesh;

    void MakeVtkMesh();
#endif //CHASTE_VTK
public:
    VtkWriter(const std::string &rDirectory,
                const std::string &rBaseName,
                const bool &rCleanDirectory=true);
    void WriteFiles();
    virtual ~VtkWriter()
    {}
};


VtkWriter::VtkWriter(const std::string &rDirectory,
                         const std::string &rBaseName,
                         const bool &rCleanDirectory)
        : AbstractMeshWriter<3,3>(rDirectory, rBaseName, rCleanDirectory)
{
    this->mIndexFromZero=true;
    MakeVtkMesh();
}


void VtkWriter::MakeVtkMesh()
{
    vtkPoints *p_pts = vtkPoints::New(VTK_DOUBLE);
    //p_pts->SetDataTypeToDouble(); 
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->mNodeData[item_num];
        p_pts->InsertPoint(item_num, current_item.data());
    }
  
    mpVtkUnstructedMesh=vtkUnstructuredGrid::New();
    //mpVtkUnstructedMesh->Allocate(rMesh.GetNumNodes(), rMesh.GetNumNodes());
    mpVtkUnstructedMesh->SetPoints(p_pts); 
        
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        std::vector<unsigned> current_element = this->mElementData[item_num];
        vtkTetra * p_tetra = vtkTetra::New();
        vtkIdList * p_tetra_id_list = p_tetra->GetPointIds();
        for (int j = 0; j < 4; ++j)
        {
            p_tetra_id_list->SetId(j, current_element[j]);
        }
        mpVtkUnstructedMesh->InsertNextCell(p_tetra->GetCellType(), p_tetra_id_list);
    }
}    


void VtkWriter::WriteFiles()
{
    vtkXMLUnstructuredGridWriter *p_writer = vtkXMLUnstructuredGridWriter::New();
    p_writer->SetInput(mpVtkUnstructedMesh);
    std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName+".vtu";
    p_writer->SetFileName(vtk_file_name.c_str());
    p_writer->Write();
}
#endif //CHASTE_VTK 

#endif /*VTKWRITER_HPP_*/
