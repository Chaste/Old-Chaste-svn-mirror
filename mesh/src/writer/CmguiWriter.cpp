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
#include "Exception.hpp"
#include "CmguiWriter.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

CmguiWriter::CmguiWriter(const std::string &rDirectory,
                         const std::string &rBaseName,
                         const bool &rCleanDirectory)
        : AbstractTetrahedralMeshWriter<3,3>(rDirectory, rBaseName, rCleanDirectory)
{
    this->mIndexFromZero = false;
}

void CmguiWriter::WriteFiles()
{
    //////////////////////////
    // Write the exnode file
    //////////////////////////
    std::string node_file_name = this->mBaseName + ".exnode";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);
    
    // Write the node header
   * p_node_file << "Group name: " << this->mBaseName << "\n";
   * p_node_file << CmguiNodeFileHeader;

    // Write each node's data
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->mNodeData[item_num];

       * p_node_file << "Node:\t" << item_num+1 << "\t";
        for (unsigned i=0; i<3; i++)
        {
           * p_node_file << current_item[i] << "\t";
        }

       * p_node_file << "\n";
    }
    p_node_file->close();

    //////////////////////////
    // Write the exlem file
    //////////////////////////
    std::string elem_file_name = this->mBaseName + ".exelem";
    out_stream p_elem_file = this->mpOutputFileHandler->OpenOutputFile(elem_file_name);

    // Write the elem header
   * p_elem_file << "Group name: " << this->mBaseName << "\n";
   * p_elem_file << CmguiElementFileHeader;

    //now we need to figure out how many additional fields we have 
    unsigned number_of_fields = mAdditionalFieldNames.size();
    std::stringstream string_of_number_of_fields;
    //we write the number of additional fields + 1 because the coordinates field gets written anyway
    string_of_number_of_fields << number_of_fields+1;
    //and write accordingly the total number of fields
   * p_elem_file << " #Fields="<<string_of_number_of_fields.str()<<"\n"; 
    
    //first field (the coordinates field is fixed and always there
   * p_elem_file << CmguiCoordinatesFileHeader;
    
    //now write the specification for each additional field
    for (unsigned i = 0; i <  number_of_fields; i++)
    {
        //unsigned to string
        std::stringstream i_string;
        i_string << i+2;
       * p_elem_file<<i_string.str()<<")  "<<mAdditionalFieldNames[i]<<" ,";
       * p_elem_file << CmguiAdditonalFieldHeader;
    }
    
    // Write each elements's data
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        std::vector<unsigned> current_element = this->mElementData[item_num];

       * p_elem_file << "Element:\t" << item_num+1 << " 0 0 Nodes:\t";
        for (unsigned i=0; i<4; i++)
        {
           * p_elem_file << current_element[i]+1 << "\t";
        }

       * p_elem_file << "\n";
    }
    p_elem_file->close();
}

void CmguiWriter::SetAdditionalFieldNames(std::vector<std::string>& rFieldNames)
{
    mAdditionalFieldNames = rFieldNames;
}

