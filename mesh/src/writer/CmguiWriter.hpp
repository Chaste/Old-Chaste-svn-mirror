/*

Copyright (C) University of Oxford, 2008

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

#ifndef CMGUIWRITER_HPP_
#define CMGUIWRITER_HPP_

#include "AbstractMeshWriter.hpp"


static const char CmguiNodeFileHeader[] = " #Fields=1\n\
 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  Value index= 1, #Derivatives= 0\n\
   y.  Value index= 2, #Derivatives= 0\n\
   z.  Value index= 3, #Derivatives= 0\n";
   
static const char CmguiElementFileHeader[] = "Shape.  Dimension=3, simplex(2;3)*simplex*simplex\n\
 #Scale factor sets= 0\n\
 #Nodes= 4\n\
 #Fields=1\n\
 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n";


/**
 *  CmguiWriter
 * 
 *  Writes a mesh in Cmgui (the visualisation frontend of CMISS) format. Creates an exnode
 *  file and a exelem file. Note that the lines and faces are not written in the exelem
 *  file, so to load the data in Cmgui, you must use 'generate_faces_and_lines', i.e.
 * 
 *  gfx read node <base_file>
 *  gfx read elem <base_file> generate_faces_and_lines
 *  gfx cr win
 */
class CmguiWriter : public AbstractMeshWriter<3,3>
{
public:
    CmguiWriter(const std::string &rDirectory,
                const std::string &rBaseName,
                const bool &rCleanDirectory=true);
    void WriteFiles();
    virtual ~CmguiWriter()
    {}
};


CmguiWriter::CmguiWriter(const std::string &rDirectory,
                         const std::string &rBaseName,
                         const bool &rCleanDirectory)
        : AbstractMeshWriter<3,3>(rDirectory, rBaseName, rCleanDirectory)
{
    this->mIndexFromZero=false;
}

void CmguiWriter::WriteFiles()
{
    //////////////////////////
    // Write the exnode file
    //////////////////////////
    std::string node_file_name = this->mBaseName+".exnode";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    //Write the node header
    *p_node_file << "Group name: " << this->mBaseName << "\n";
    *p_node_file << CmguiNodeFileHeader;

    //Write each node's data
    for (unsigned item_num=0; item_num<this->GetNumNodes(); item_num++)
    {
        std::vector<double> current_item = this->mNodeData[item_num];

        *p_node_file << "Node:\t" << item_num+1 << "\t"; 
        for (unsigned i=0;i<3;i++)
        {
            *p_node_file << current_item[i] << "\t";
        }

        *p_node_file << "\n";
    }
    p_node_file->close();

    //////////////////////////
    // Write the exlem file
    //////////////////////////
    std::string elem_file_name = this->mBaseName+".exelem";
    out_stream p_elem_file = this->mpOutputFileHandler->OpenOutputFile(elem_file_name);

    //Write the elem header
    *p_elem_file << "Group name: " << this->mBaseName << "\n";
    *p_elem_file << CmguiElementFileHeader;

    //Write each elements's data
    for (unsigned item_num=0; item_num<this->GetNumElements(); item_num++)
    {
        std::vector<unsigned> current_element = this->mElementData[item_num];

        *p_elem_file << "Element:\t" << item_num+1 << " 0 0 Nodes:\t"; 
        for (unsigned i=0; i<4; i++)
        {
            *p_elem_file << current_element[i]+1 << "\t";
        }

        *p_elem_file << "\n";
    }
    p_elem_file->close();
    
}

#endif /*CMGUIWRITER_HPP_*/
