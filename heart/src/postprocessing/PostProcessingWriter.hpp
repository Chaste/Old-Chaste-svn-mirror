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

#ifndef POSTPROCESSINGWRITER_HPP_
#define POSTPROCESSINGWRITER_HPP_

#include "PetscTools.hpp"
#include "Hdf5DataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include <iostream>

class PostProcessingWriter
{
  private:
    PropagationPropertiesCalculator* mpCalculator;   
    unsigned mNumberOfNodes;
    std::string mOutputDirectory;
    
  public:
    PostProcessingWriter(Hdf5DataReader* pDataReader, std::string outputDirectory)
    {
         mpCalculator = new PropagationPropertiesCalculator(pDataReader);
         mNumberOfNodes = pDataReader->GetNumberOfRows();
         mOutputDirectory = outputDirectory;
    }
    
    ~PostProcessingWriter()
    {
        delete mpCalculator;
    }
    
    void WriteApdMapFile(double threshold, double repolarisationPercentage)
    {
        if(PetscTools::AmMaster())
        {
            out_stream p_file=out_stream(NULL);
            OutputFileHandler output_file_handler(mOutputDirectory, false);
            p_file = output_file_handler.OpenOutputFile("ApdMap.dat");
            for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
            { 
                std::vector<double> apds;
                try
                {
                    apds = mpCalculator->CalculateAllActionPotentialDurations(repolarisationPercentage, node_index, threshold);
                }
                catch(Exception& e)
                {
                    
                    apds.push_back(0);
                }
                for (unsigned i = 0; i < apds.size(); i++)
                {
                    *p_file << apds[i] << "\t";
                }
                *p_file << std::endl;
            }
            p_file->close();
        }
    }
};

#endif /*POSTPROCESSINGWRITER_HPP_*/
