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


#include <vector>

#include "UblasCustomFunctions.hpp"
#include "HeartConfig.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"

void Hdf5ToMeshalyzerConverter::Write(std::string type)
{
    assert(type=="V" || type=="Phi_e");

    out_stream p_file=out_stream(NULL);
    OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
    if (PetscTools::AmMaster())
    {
        p_file = output_file_handler.OpenOutputFile(mFileBaseName + "_" + type + ".dat");
    }

    unsigned num_nodes = mpReader->GetNumberOfRows();
    unsigned num_timesteps = mpReader->GetUnlimitedDimensionValues().size();

    DistributedVectorFactory factory(num_nodes);

    Vec data = factory.CreateVec();
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        mpReader->GetVariableOverNodes(data, type, time_step);
        ReplicatableVector repl_data(data);

        assert(repl_data.GetSize()==num_nodes);

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


Hdf5ToMeshalyzerConverter::Hdf5ToMeshalyzerConverter(std::string inputDirectory,
                          std::string fileBaseName)
{
    // store dir and filenames, and create a reader
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

    OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);
    if (PetscTools::AmMaster())
    {
        //Note that we don't want the child processes to write info files

        out_stream p_file = output_file_handler.OpenOutputFile(mFileBaseName + "_times.info");
        unsigned num_timesteps = mpReader->GetUnlimitedDimensionValues().size();
       * p_file << "Number of timesteps "<<num_timesteps<<"\n";
       * p_file << "timestep "<<HeartConfig::Instance()->GetPrintingTimeStep()<<"\n";
        double first_timestep=mpReader->GetUnlimitedDimensionValues().front();
       * p_file << "First timestep "<<first_timestep<<"\n";
        double last_timestep=mpReader->GetUnlimitedDimensionValues().back();
       * p_file << "Last timestep "<<last_timestep<<"\n";

        p_file->close();

    }

    PetscTools::Barrier();
}

Hdf5ToMeshalyzerConverter::~Hdf5ToMeshalyzerConverter()
{
    delete mpReader;
}
