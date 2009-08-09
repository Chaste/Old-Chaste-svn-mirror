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
#include "Hdf5ToCmguiConverter.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"

void Hdf5ToCmguiConverter::Write(std::string type)
{
    assert(type=="Mono" || type=="Bi");
    out_stream p_file=out_stream(NULL);
    OutputFileHandler output_file_handler(HeartConfig::Instance()->GetOutputDirectory() + "/output", false);

    unsigned num_nodes = mpReader->GetNumberOfRows();
    unsigned num_timesteps = mpReader->GetUnlimitedDimensionValues().size();

    DistributedVectorFactory factory(num_nodes);

    Vec data = factory.CreateVec();//for V
    Vec data_phie = factory.CreateVec();//for phi_e
    
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        //create the file for this time step
        std::stringstream time_step_string;
        //unsigned to string
        time_step_string << time_step;
        if (PetscTools::AmMaster())
        {
            p_file = output_file_handler.OpenOutputFile(mFileBaseName + "_" + time_step_string.str() + ".exnode");
        }
        
        //read the data for this time step
        mpReader->GetVariableOverNodes(data, "V", time_step);
        ReplicatableVector repl_data(data);
        assert(repl_data.size()==num_nodes);
        
        //get the data for phie only if needed
        ReplicatableVector repl_data_phie;
        if (type=="Bi")
        {
            repl_data_phie.resize(num_nodes);
            mpReader->GetVariableOverNodes(data_phie, "Phi_e", time_step);
            repl_data_phie.ReplicatePetscVector(data_phie);          
        }
       
        if(PetscTools::AmMaster())
        {
            //The header first
            *p_file << "Group name: " << this->mFileBaseName << "\n";
            
            //we need two fields for bidomain and one only for monodomain
            if(type=="Mono")
            {
                *p_file << "#Fields=1" << "\n" << " 1) " << "V , field, rectangular cartesian, #Components=1" << "\n" << "x.  Value index=1, #Derivatives=0, #Versions=1"<<"\n";
            }
            else
            {
                *p_file << "#Fields=2" << "\n" << " 1) " << "V , field, rectangular cartesian, #Components=1" << "\n" << "x.  Value index=1, #Derivatives=0, #Versions=1"<<"\n";
                //the details of the second field
                *p_file << "\n" << " 2) " << "Phie , field, rectangular cartesian, #Components=1" << "\n" << "x.  Value index=1, #Derivatives=0, #Versions=1"<<"\n";
            }
            
            //write the data 
            for(unsigned i=0; i<num_nodes; i++)
            {
                //cmgui counts nodes from 1
                *p_file << "Node: "<< i+1 << "\n" << repl_data[i] << "\n";
                //if it is a bidomain simulation, write the data for phie
                if (type=="Bi")
                {
                    *p_file <<  repl_data_phie[i] << "\n";
                }
            }
        }
    }
    VecDestroy(data);
    VecDestroy(data_phie);

    if(PetscTools::AmMaster())
    {
        p_file->close();
    }
}

Hdf5ToCmguiConverter::Hdf5ToCmguiConverter(std::string inputDirectory,
                          std::string fileBaseName)
{
    // store dir and filenames, and create the reader
    mFileBaseName = fileBaseName;
    mpReader = new Hdf5DataReader(inputDirectory, mFileBaseName);
    
    // check the data file read has one or two variables (ie V; or V and PhiE)
    std::vector<std::string> variable_names = mpReader->GetVariableNames();
    if((variable_names.size()==0) || (variable_names.size()>2))
    {
        delete mpReader;
        EXCEPTION("Data has zero or more than two variables - doesn't appear to be mono or bidomain");
    }

    // if one variable, it is a monodomain problem
    if(variable_names.size()==1)
    {
        if(variable_names[0]!="V")
        {
            delete mpReader;
            EXCEPTION("One variable, but it is not called 'V'");
        }

        Write("Mono");
    }

    // if two variables, it is a bidomain problem
    if(variable_names.size()==2)
    {
        if(variable_names[0]!="V" || variable_names[1]!="Phi_e")
        {
            delete mpReader;
            EXCEPTION("Two variables, but they are not called 'V' and 'Phi_e'");
        }

        Write("Bi");
    }

    MPI_Barrier(PETSC_COMM_WORLD);


}

Hdf5ToCmguiConverter::~Hdf5ToCmguiConverter()
{
    delete mpReader;
}
