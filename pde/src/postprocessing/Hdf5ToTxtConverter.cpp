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

#include "UblasCustomFunctions.hpp"
#include "Hdf5ToTxtConverter.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "GenericMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "Warnings.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToTxtConverter<ELEMENT_DIM, SPACE_DIM>::Hdf5ToTxtConverter(std::string inputDirectory,
                                                               std::string fileBaseName,
                                                               AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(inputDirectory, fileBaseName, pMesh, "txt_output")
{
    // Make sure that we are never trying to write from an incomplete data HDF5 file
    assert(this->mpReader->GetNumberOfRows() == pMesh->GetNumNodes());

    std::string output_directory = inputDirectory + "/" + this->mRelativeSubdirectory;
    OutputFileHandler handler(output_directory);

    DistributedVectorFactory* p_factory = pMesh->GetDistributedVectorFactory();
    Vec data = p_factory->CreateVec();

    unsigned num_nodes = pMesh->GetNumNodes();
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    // Loop over time steps
    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        // Create a .txt file for this time step
        std::stringstream file_name;
        file_name << fileBaseName << "_" << time_step << ".txt";
        out_stream p_file = handler.OpenOutputFile(file_name.str());

        ///\todo this seems like a very inefficient way to do this! (#1841)

        // Loop over nodes
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            // Get node location
            c_vector<double, SPACE_DIM> node_location = pMesh->GetNode(node_index)->rGetLocation();
            *p_file << node_location[0];
            for (unsigned dim=1; dim<SPACE_DIM; dim++)
            {
                *p_file << " " << node_location[dim];
            }

            // Loop over variables
            for (unsigned var_index=0; var_index<this->mNumVariables; var_index++)
            {
                // Get variable at this time step from HDF5 archive
                std::string variable_name = this->mpReader->GetVariableNames()[var_index];
                this->mpReader->GetVariableOverNodes(data, variable_name, time_step);
                ReplicatableVector repl_data(data);

                *p_file << " " << repl_data[node_index];
            }
            *p_file << "\n";
        }
    }

    // Tidy up
    VecDestroy(data);
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Hdf5ToTxtConverter<1,1>;
template class Hdf5ToTxtConverter<1,2>;
template class Hdf5ToTxtConverter<2,2>;
template class Hdf5ToTxtConverter<1,3>;
template class Hdf5ToTxtConverter<2,3>;
template class Hdf5ToTxtConverter<3,3>;
