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
#include "Exception.hpp"
#include "CmguiMeshWriter.hpp"
#include "Version.hpp"
#include <boost/shared_ptr.hpp>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::CmguiMeshWriter(const std::string &rDirectory,
                                                        const std::string &rBaseName,
                                                        bool cleanDirectory)
        : AbstractTetrahedralMeshWriter<ELEMENT_DIM,SPACE_DIM>(rDirectory, rBaseName, cleanDirectory)
{
    this->mIndexFromZero = false;
    mGroupName = this->mBaseName;
    
    switch(ELEMENT_DIM)
    {
        case 1:
        {
            mElementFileHeader = CmguiElementFileHeader1D;
            mCoordinatesFileHeader = CmguiCoordinatesFileHeader1D;
            mAdditionalFieldHeader = CmguiAdditionalFieldHeader1D;
            break;
        }
        case 2:
        {
            mElementFileHeader = CmguiElementFileHeader2D;
            mCoordinatesFileHeader = CmguiCoordinatesFileHeader2D;
            mAdditionalFieldHeader = CmguiAdditionalFieldHeader2D;
            break;
        }
        case 3:
        {
            mElementFileHeader = CmguiElementFileHeader3D;
            mCoordinatesFileHeader = CmguiCoordinatesFileHeader3D;
            mAdditionalFieldHeader = CmguiAdditionalFieldHeader3D;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }    
    

    mNumNodesPerElement = ELEMENT_DIM+1;
    mReordering.resize(mNumNodesPerElement);
    for(unsigned i=0; i<mNumNodesPerElement; i++)
    {
        mReordering[i] = i;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteFiles()
{
    //////////////////////////
    // Write the exnode file
    //////////////////////////
    std::string node_file_name = this->mBaseName + ".exnode";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    WriteNodeFileHeader(p_node_file);

    // Write each node's data
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->GetNextNode();

        *p_node_file << "Node:\t" << item_num+1 << "\t";
        for (unsigned i=0; i<ELEMENT_DIM; i++)
        {
            *p_node_file << current_item[i] << "\t";
        }

        *p_node_file << "\n";
    }
    p_node_file->close();

    //////////////////////////
    // Write the exlem file
    //////////////////////////
    
    // If nobody defined region names we default to the same name as the file.
    if (mRegionNames.size() == 0)
    {
        mRegionNames.push_back(this->mBaseName);
    }

    // Array with file descriptors for each of regions    
    std::vector<boost::shared_ptr<std::ofstream> > p_elem_file;
    
    p_elem_file.resize(mRegionNames.size());
    
    for (unsigned region_index=0; region_index<mRegionNames.size(); region_index++)     
    {
        std::string elem_file_name = mRegionNames[region_index] + ".exelem";
        out_stream p_tmp_file = this->mpOutputFileHandler->OpenOutputFile(elem_file_name);
        p_elem_file[region_index] = p_tmp_file;
        
        // Write the elem header
        
        //write provenance info
        std::string comment = "! " + ChasteBuildInfo::GetProvenanceString();
        *p_elem_file[region_index] << comment;
    
        *p_elem_file[region_index] << "Group name: " << mGroupName << "\n";
        *p_elem_file[region_index] << mElementFileHeader;
    
    
        //now we need to figure out how many additional fields we have
        unsigned number_of_fields = mAdditionalFieldNames.size();
        std::stringstream string_of_number_of_fields;
        //we write the number of additional fields + 1 because the coordinates field gets written anyway
        string_of_number_of_fields << number_of_fields+1;
        //and write accordingly the total number of fields
        *p_elem_file[region_index] << " #Fields="<<string_of_number_of_fields.str()<<"\n";
    
        //first field (the coordinates field is fixed and alwys there)
        *p_elem_file[region_index] << mCoordinatesFileHeader;
    
    
        //now write the specification for each additional field
        for (unsigned i = 0; i <  number_of_fields; i++)
        {
            //unsigned to string
            std::stringstream i_string;
            i_string << i+2;
            *p_elem_file[region_index]<<i_string.str()<<")  "<<mAdditionalFieldNames[i]<<" ,";
            *p_elem_file[region_index] << mAdditionalFieldHeader;
        }
    }

    // Write each elements's data        
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        ElementData elem =this->GetNextElement();
        std::vector<unsigned> current_element = elem.NodeIndices;
        
        /// \todo: EXCEPTION maybe...
        assert(elem.AttributeValue < mRegionNames.size());
           
        *p_elem_file[elem.AttributeValue] << "Element:\t" << item_num+1 << " 0 0 Nodes:\t";
        for (unsigned i=0; i<mNumNodesPerElement; i++)
        {
            *p_elem_file[elem.AttributeValue] << current_element[mReordering[i]]+1 << "\t";
        }

        *p_elem_file[elem.AttributeValue] << "\n";
            
    }
    
    for (unsigned region_index=0; region_index<mRegionNames.size(); region_index++)     
    {
        p_elem_file[region_index]->close();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::SetAdditionalFieldNames(std::vector<std::string>& rFieldNames)
{
    mAdditionalFieldNames = rFieldNames;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::SetRegionNames(std::vector<std::string>& rRegionNames)
{    
    mRegionNames = rRegionNames;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::WriteNodeFileHeader(out_stream& rpNodeFile)
{
    //write provenance info
    std::string comment = "! " + ChasteBuildInfo::GetProvenanceString();
    *rpNodeFile << comment;
    
    // Write the node header
    *rpNodeFile << "Group name: " << this->mGroupName << "\n";
    switch (SPACE_DIM)
    {
        case 1:
        {
            *rpNodeFile << CmguiNodeFileHeader1D;
            break;
        }
        case 2:
        {
            *rpNodeFile << CmguiNodeFileHeader2D;
            break;
        }
        case 3:
        {
            *rpNodeFile << CmguiNodeFileHeader3D;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM>::CompareCmguiFiles(std::string& rPath1, std::string& rPath2)
{
    std::string compare_command = "diff --ignore-matching-lines=\"! \" ";
    compare_command += rPath1;
    compare_command += " ";
    compare_command += rPath2;

    //Compare the new test file with one from the repository
    if (system(compare_command.c_str()) == 0)
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

template class CmguiMeshWriter<1,1>;
template class CmguiMeshWriter<1,2>;
template class CmguiMeshWriter<1,3>;
template class CmguiMeshWriter<2,2>;
template class CmguiMeshWriter<2,3>;
template class CmguiMeshWriter<3,3>;
