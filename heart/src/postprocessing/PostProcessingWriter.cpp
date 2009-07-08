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

#include "UblasCustomFunctions.hpp"
#include "HeartConfig.hpp"
#include "PostProcessingWriter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include <iostream>

PostProcessingWriter::PostProcessingWriter(Hdf5DataReader* pDataReader)
{
    mpCalculator = new PropagationPropertiesCalculator(pDataReader);
    mNumberOfNodes = pDataReader->GetNumberOfRows();
}
    
PostProcessingWriter::~PostProcessingWriter()
{
    delete mpCalculator;
}


void PostProcessingWriter::WriteApdMapFile(double threshold, double repolarisationPercentage)
{
    OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
    if(PetscTools::AmMaster())
    {
        out_stream p_file=out_stream(NULL);
        std::stringstream stream;
        stream << repolarisationPercentage;
        p_file = output_file_handler.OpenOutputFile("Apd" + stream.str() + "Map.dat");
        for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
        { 
            std::vector<double> apds;
            try
            {
                apds = mpCalculator->CalculateAllActionPotentialDurations(repolarisationPercentage, node_index, threshold);
                assert(apds.size() != 0);
            }
            catch(Exception& e)
            {                    
                apds.push_back(0);
                assert(apds.size() == 1);
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


void PostProcessingWriter::WriteUpstrokeTimeMap(double threshold)
{
    if(PetscTools::AmMaster())
    {
        out_stream p_file=out_stream(NULL);
        OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
        p_file = output_file_handler.OpenOutputFile("UpstrokeTimeMap.dat");
        for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
        { 
            std::vector<double> upstroke_times;
            upstroke_times = mpCalculator->CalculateUpstrokeTimes(node_index, threshold);
            assert(upstroke_times.size()!=0); ///\todo Fix (see above)
            for (unsigned i = 0; i < upstroke_times.size(); i++)
            {
                *p_file << upstroke_times[i] << "\t";
            }
            *p_file << std::endl;
        }
        p_file->close();
    }
}

void PostProcessingWriter::WriteMaxUpstrokeVelocityMap(double threshold)
{
    if(PetscTools::AmMaster())
    {
        out_stream p_file=out_stream(NULL);
        OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
        p_file = output_file_handler.OpenOutputFile("MaxUpstrokeVelocityMap.dat");
        for (unsigned node_index = 0; node_index < mNumberOfNodes; node_index++)
        { 
            std::vector<double> upstroke_velocities;
            upstroke_velocities = mpCalculator->CalculateAllMaximumUpstrokeVelocities(node_index, threshold);
            assert(upstroke_velocities.size()!=0); ///\todo Fix (see above)
            for (unsigned i = 0; i < upstroke_velocities.size(); i++)
            {
                *p_file << upstroke_velocities[i] << "\t";
            }
            *p_file << std::endl;
         }
         p_file->close();
    }
}

void PostProcessingWriter::WriteConductionVelocityMap(unsigned originNode, std::vector<double> distancesFromOriginNode)
{
    if(PetscTools::AmMaster())
    {
        out_stream p_file=out_stream(NULL);
        OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
        
        std::stringstream filename;
        filename << "ConductionVelocityFromNode" << originNode << ".dat";               
        p_file = output_file_handler.OpenOutputFile(filename.str());
        for (unsigned dest_node = 0; dest_node < mNumberOfNodes; dest_node++)
        { 
            std::vector<double> conduction_velocities;
            conduction_velocities = mpCalculator->CalculateAllConductionVelocities(originNode, dest_node, distancesFromOriginNode[dest_node]);
            assert(conduction_velocities.size()!=0);
            for (unsigned i = 0; i < conduction_velocities.size(); i++)
            {
                *p_file << conduction_velocities[i] << "\t";
            }
            *p_file << std::endl;
         }
         p_file->close();
    }        
}

