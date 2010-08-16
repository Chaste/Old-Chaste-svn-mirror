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

#include "Hdf5ToMeshalyzerConverter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "GenericMeshReader.hpp"
#include "UblasCustomFunctions.hpp"
#include "HeartConfig.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "Version.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM>::Write(std::string type)
{
    //assert(type=="V" || type=="Phi_e");

    out_stream p_file=out_stream(NULL);
    if (PetscTools::AmMaster())
    {
        p_file = this->mpOutputFileHandler->OpenOutputFile(this->mFileBaseName + "_" + type + ".dat");
    }

    unsigned num_nodes = this->mpReader->GetNumberOfRows();
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    DistributedVectorFactory factory(num_nodes);

    Vec data = factory.CreateVec();
    ReplicatableVector repl_data(num_nodes);
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        this->mpReader->GetVariableOverNodes(data, type, time_step);
        repl_data.ReplicatePetscVector(data);
        
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
    	std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();
    	*p_file << comment;
        p_file->close();
    }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM>::Hdf5ToMeshalyzerConverter(std::string inputDirectory,
                          std::string fileBaseName,
                          AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> *pMesh) :
                     AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(inputDirectory, fileBaseName, pMesh, "output")
{

    Write(this->mpReader->GetVariableNames()[0]);//V
    if (this->mNumVariables == 2)//if it is a bidomain problem
    {
        Write(this->mpReader->GetVariableNames()[1]);//Phi_e
    }
    if (this->mNumVariables == 3)//for extended bidomain problems
    {
    	Write(this->mpReader->GetVariableNames()[1]);//V for the second cell
        Write(this->mpReader->GetVariableNames()[2]);//then Phi_e
    }

    //Write mesh in a suitable form for meshalyzer
    std::string output_directory =  HeartConfig::Instance()->GetOutputDirectory() + "/output";
    MeshalyzerMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(output_directory, HeartConfig::Instance()->GetOutputFilenamePrefix()+"_mesh", false);
    //Normal case is that the in-memory mesh is converted
    if (HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering() == false )
    {
        // The second argument tells the writer to not follow original element ordering for performance reasons.
        mesh_writer.WriteFilesUsingMesh(*(this->mpMesh), false);
    }
    else
    {
        //In this case we expect the mesh to have been read in from file
        ///\todo What if the mesh has been scaled, translated or rotated?
        //Note that the next line will throw if the mesh has not been read from file
        std::string original_file=this->mpMesh->GetMeshFileBaseName();
        GenericMeshReader<ELEMENT_DIM, SPACE_DIM> original_mesh_reader(original_file);
        mesh_writer.WriteFilesUsingMeshReader(original_mesh_reader);
    }
    PetscTools::Barrier("Hdf5ToMeshalyzerConverter");
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Hdf5ToMeshalyzerConverter<1,1>;
template class Hdf5ToMeshalyzerConverter<1,2>;
template class Hdf5ToMeshalyzerConverter<2,2>;
template class Hdf5ToMeshalyzerConverter<1,3>;
template class Hdf5ToMeshalyzerConverter<2,3>;
template class Hdf5ToMeshalyzerConverter<3,3>;
