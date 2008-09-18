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


#ifndef HDF5TOMESHALYZERCONVERTER_HPP_
#define HDF5TOMESHALYZERCONVERTER_HPP_

#include "Hdf5DataReader.hpp"
#include "PetscTools.hpp"
#include "HeartConfig.hpp"

/**
 *  This class converts from Hdf5 format to meshalyzer format, ie, for
 *  voltage, one file, which looks like
 *
 *  V_node_0_time_0
 *  ..
 *  V_node_N_time_0
 *  V_node_0_time_1
 *  ..
 *  V_node_N_time_1
 *  V_node_0_time_2
 *  ..
 *  V_node_N_time_M
 *
 *  The files that are written are <base_name>_V.dat or <base_name>_Phi_e.dat,
 *  where <base_name> is the base name of the original .h5 file. The new files
 *  are written in the same directory as the .h5 file. All paths are relative
 *  to the CHASTE_TEST_OUTPUT directory.
 */
class Hdf5ToMeshalyzerConverter
{
private:
    Hdf5DataReader* mpReader;
    std::string mOutputDirectory;
    std::string mFileBaseName;

    /** A helper method which takes in a string, which must be 'V' or 'Phi_e'
     *  and reads the data corresponding to that string, writing it out in
     *  meshalyzer format.
     */
    void Write(std::string type)
    {
        assert(type=="V" || type=="Phi_e");
        
        out_stream p_file=out_stream(NULL);
        if (PetscTools::AmMaster())
        {
            //Note that we don't want the child processes to create
            //a fresh directory if it doesn't already exist
            OutputFileHandler output_file_handler(mOutputDirectory, false);
            p_file = output_file_handler.OpenOutputFile(mFileBaseName + "_" + type + ".dat");
        }

        unsigned num_nodes = mpReader->GetNumberOfRows();
        unsigned num_timesteps = mpReader->GetUnlimitedDimensionValues().size();

        DistributedVector::SetProblemSize(num_nodes);
        Vec data = DistributedVector::CreateVec();
        for (unsigned time_step=0; time_step<num_timesteps; time_step++)
        {
            mpReader->GetVariableOverNodes(data, type, time_step);
            ReplicatableVector repl_data(data);

            assert(repl_data.size()==num_nodes);

            if(PetscTools::AmMaster())
            {
                for(unsigned i=0; i<num_nodes; i++)
                {
                    *p_file << repl_data[i] << "\n";
                }
            }
        }
        VecDestroy(data);
        if(PetscTools::AmMaster())
        {
            p_file->close();
        }
    }


public:
    /** Constructor, which does the conversion.
     *  @param inputDirectory The input directory, relative to CHASTE_TEST_OUTPUT,
     *  where the .h5 file has been written
     *  @param outputDirectory  The output directory, relative to CHASTE_TEST_OUTPUT, where the output will be place
     *  @param fileBaseName The base name of the data file.
     */
    Hdf5ToMeshalyzerConverter(std::string inputDirectory,
                              std::string outputDirectory,
                              std::string fileBaseName)
    {
        // store dir and filenames, and create a reader
        mOutputDirectory = outputDirectory;
        mFileBaseName = fileBaseName;
        mpReader = new Hdf5DataReader(inputDirectory, mFileBaseName);

        // check the data file read has one or two variables (ie V; or V and PhiE)
        std::vector<std::string> variable_names = mpReader->GetVariableNames();
        if((variable_names.size()==0) || (variable_names.size()>2))
        {
            delete mpReader;
            EXCEPTION("Data has zero or more than two variables - doesn't appear to be mono or bidomain");
        }

        // if one variable, a monodomain problem
        if(variable_names.size()==1)
        {
            if(variable_names[0]!="V")
            {
                delete mpReader;
                EXCEPTION("One variable, but it is not called 'V'");
            }

            Write("V");
        }

        // if two variable, a bidomain problem
        if(variable_names.size()==2)
        {
            if(variable_names[0]!="V" || variable_names[1]!="Phi_e")
            {
                delete mpReader;
                EXCEPTION("Two variables, but they are not called 'V' and 'Phi_e'");
            }

            Write("V");
            Write("Phi_e");
        }
        
        MPI_Barrier(PETSC_COMM_WORLD);
       if (PetscTools::AmMaster())
        {
            //Note that we don't want the child processes to create
            //a fresh directory if it doesn't already exist
            OutputFileHandler output_file_handler(mOutputDirectory, false);
            out_stream p_file = output_file_handler.OpenOutputFile(mFileBaseName + "_times.info");
            unsigned num_timesteps = mpReader->GetUnlimitedDimensionValues().size();
            *p_file << "Number of timesteps "<<num_timesteps<<"\n";
            *p_file << "timestep "<<HeartConfig::Instance()->GetPrintingTimeStep()<<"\n";
            double first_timestep=mpReader->GetUnlimitedDimensionValues().front();
            *p_file << "First timestep "<<first_timestep<<"\n";
            double last_timestep=mpReader->GetUnlimitedDimensionValues().back();
            *p_file << "Last timestep "<<last_timestep<<"\n";
            
            p_file->close();
            
        }



    }

    ~Hdf5ToMeshalyzerConverter()
    {
        delete mpReader;
    }
};

#endif /*HDF5TOMESHALYZERCONVERTER_HPP_*/
