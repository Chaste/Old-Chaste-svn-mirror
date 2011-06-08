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

#include "FibreWriter.hpp"
#include "Version.hpp"

template<unsigned DIM>
FibreWriter<DIM>::FibreWriter(const std::string& rDirectory,
                              const std::string& rBaseName,
                              const bool clearOutputDir)
    : mBaseName(rBaseName),
      mFileIsBinary(false)               
{
    mpOutputFileHandler= new OutputFileHandler(rDirectory, clearOutputDir); 
}

template<unsigned DIM>
FibreWriter<DIM>::~FibreWriter()
{
    delete mpOutputFileHandler;
}

template<unsigned DIM>
void FibreWriter<DIM>::WriteAllAxi(const std::vector< c_vector<double, DIM> >& fibres)
{
    // Write axi file
    std::string axi_file_name = this->mBaseName + ".axi";
    out_stream p_axi_file = this->mpOutputFileHandler->OpenOutputFile(axi_file_name);
  
    //Header line says how many fibres are coming...
    *p_axi_file << fibres.size();
    //... and whether the fibres are binary
    if (this->mFileIsBinary)
    {
        *p_axi_file << "\tBIN\n";
    }
    else
    {
        *p_axi_file << "\n";
    }
    //Now give the fibre directions
    for (unsigned i=0; i<fibres.size();i++ )
    {
        if (this->mFileIsBinary)
        {
            p_axi_file->write((char*)&fibres[i][0], DIM*sizeof(double));
        }
        else
        {
            for(unsigned j=0; j<DIM; j++)
            {
                *p_axi_file << fibres[i][j] << "\t";
            }
            *p_axi_file <<"\n";
        }
    }
    *p_axi_file << "#\n# " <<  ChasteBuildInfo::GetProvenanceString();
    p_axi_file->close();
}


template<unsigned DIM>
void FibreWriter<DIM>::SetWriteFileAsBinary()
{
    mFileIsBinary = true;   
}

template class FibreWriter<1>;
template class FibreWriter<2>;
template class FibreWriter<3>;





