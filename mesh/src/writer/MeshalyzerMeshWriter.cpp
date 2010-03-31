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

#include "MeshalyzerMeshWriter.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshalyzerMeshWriter<ELEMENT_DIM, SPACE_DIM>::MeshalyzerMeshWriter(const std::string &rDirectory,
        const std::string &rBaseName,
        const bool &rCleanDirectory,
        const bool &rSetCoolGraphics)
        : AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, rCleanDirectory)
{
   /* if (ELEMENT_DIM != SPACE_DIM)
    {
        EXCEPTION("ELEMENT_DIM must be equal to SPACE_DIM");
    }*/

    if (rSetCoolGraphics)
    {
        this->mIndexFromZero = false;
        this->mWriteMetaFile = true;
    }
    else
    {
        this->mIndexFromZero = true;
        this->mWriteMetaFile = false;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MeshalyzerMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();
    
    //Write node file
    std::string node_file_name = this->mBaseName + ".pts";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    //Write the node header
    unsigned num_nodes = this->GetNumNodes();
    *p_node_file << num_nodes << "\n";

    // Write each node's data
    for (unsigned item_num=0; item_num<num_nodes; item_num++)
    {
        std::vector<double> current_item = this->GetNextNode(); //this->mNodeData[item_num];
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *p_node_file << current_item[i] << "\t";
        }
        if (SPACE_DIM==2)
        {
            *p_node_file << 0 << "\t";
        }
        if (SPACE_DIM==1)
        {
            *p_node_file << 0 << "\t" << 0 << "\t";
        }
        *p_node_file << "\n";

    }
    *p_node_file << comment;
    p_node_file->close();

    // Write element file
    std::string element_file_name;

    if (ELEMENT_DIM == 3)
    {
        element_file_name = this->mBaseName + ".tetras";
    }
    else if (ELEMENT_DIM == 2)
    {
        element_file_name = this->mBaseName + ".tri";
    }
    else //ELEMENT_DIM == 1
    {
        element_file_name = this->mBaseName + ".cnnx";
    }

    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

    // Write the element header
    unsigned num_elements = this->GetNumElements();

    *p_element_file << num_elements << "\n";

    // Write each element's data
    unsigned nodes_per_element = ELEMENT_DIM+1;
    for (unsigned item_num=0; item_num<num_elements; item_num++)
    {
        ElementData element_data = this->GetNextElement();

        std::vector<unsigned> current_item = element_data.NodeIndices;
        for (unsigned i=0; i<nodes_per_element; i++)
        {
            if (this->mIndexFromZero)
            {
                *p_element_file << current_item[i] << "\t";
            }
            else
            {
                *p_element_file << current_item[i]+1 << "\t";
            }
        }

        *p_element_file << element_data.AttributeValue << "\n";
    }
    *p_element_file << comment;
    p_element_file->close();

    if (ELEMENT_DIM==3)
    {
        // Write boundary face file
        std::string face_file_name = this->mBaseName + ".tri";
        out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name);

        // Write the boundary face header
        unsigned num_faces = this->GetNumBoundaryFaces();

        *p_face_file << num_faces << "\n";

        // Write each face's data
        double material_property = 0.0;
        for (unsigned item_num=0; item_num<num_faces; item_num++)
        {
            std::vector<unsigned> current_item = this->mBoundaryFaceData[item_num];
            for (unsigned i=0; i<ELEMENT_DIM; i++)
            {
                if (this->mIndexFromZero)
                {
                    *p_face_file << current_item[i] << "\t";
                }
                else
                {
                    *p_face_file << current_item[i]+1 <<"\t";
                }
            }
            *p_face_file << material_property << "\n";
        }
        *p_face_file << comment;
        p_face_file->close();

        if (this->mWriteMetaFile)
        {
            std::string meta_file_name = this->mBaseName + ".cg_in";
            out_stream p_meta_file = this->mpOutputFileHandler->OpenOutputFile(meta_file_name);

            *p_meta_file << "1\n" << "0\n";
            *p_meta_file << face_file_name <<"\n";
            *p_meta_file << comment;
            p_meta_file->close();
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MeshalyzerMeshWriter<ELEMENT_DIM, SPACE_DIM>::~MeshalyzerMeshWriter()
{
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class MeshalyzerMeshWriter<1,1>;
template class MeshalyzerMeshWriter<1,2>;
template class MeshalyzerMeshWriter<1,3>;
template class MeshalyzerMeshWriter<2,2>;
template class MeshalyzerMeshWriter<2,3>;
template class MeshalyzerMeshWriter<3,3>;
