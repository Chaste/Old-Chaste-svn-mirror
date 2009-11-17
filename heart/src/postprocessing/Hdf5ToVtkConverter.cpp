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
#include "Hdf5ToVtkConverter.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "VtkWriter.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToVtkConverter<ELEMENT_DIM, SPACE_DIM>::Hdf5ToVtkConverter(std::string inputDirectory,
                          std::string fileBaseName,
                          AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh) :
                    AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(inputDirectory, fileBaseName, pMesh, "vtk_output")
{   
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar

    VtkWriter<ELEMENT_DIM,SPACE_DIM> vtk_writer(HeartConfig::Instance()->GetOutputDirectory() + "/vtk_output", fileBaseName, false);
    
    unsigned num_nodes = this->mpReader->GetNumberOfRows();
    DistributedVectorFactory factory(num_nodes);
    ///\todo Can we get this as a std::vector?
    Vec data = factory.CreateVec();//for V
    
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();
    
    // Loop over time steps
    for (unsigned time_step=0; time_step<num_timesteps; time_step++) //num_timesteps; time_step++)
    {
        this->mpReader->GetVariableOverNodes(data, "V", time_step); // Gets V at this time step from HDF5 archive
        ReplicatableVector repl_data(data);
        std::vector<double> v_for_vtk;
        v_for_vtk.resize(num_nodes);
        for (unsigned index=0; index<num_nodes; index++)
        {
            v_for_vtk[index]  = repl_data[index];
        }
        
        std::ostringstream V_point_data_name;
        V_point_data_name << "V_" << std::setw(6) << std::setfill('0') << time_step; 
        
        if (PetscTools::AmMaster())
        {
            // Add V into the node "point" data
            vtk_writer.AddPointData(V_point_data_name.str(), v_for_vtk);  
        }
        if(this->mNumVariables==2)
        {
            this->mpReader->GetVariableOverNodes(data, "Phi_e", time_step); // Gets Phi at this time step from HDF5 archive
            ReplicatableVector repl_phi(data);
            std::vector<double> phi_for_vtk;
            phi_for_vtk.resize(num_nodes);
            for (unsigned index=0; index<num_nodes; index++)
            {
                phi_for_vtk[index]  = repl_phi[index];
            }
            
            std::ostringstream Phi_point_data_name;
            Phi_point_data_name << "Phi_e_" << std::setw(6) << std::setfill('0') << time_step; 
        
            if (PetscTools::AmMaster())
            {
                // Add Phi into the node "point" data
                vtk_writer.AddPointData(Phi_point_data_name.str(), phi_for_vtk);
            }              
        }  
    }
    VecDestroy(data);
    vtk_writer.WriteFilesUsingMesh( *(this->mpMesh) );
#endif //CHASTE_VTK

}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Hdf5ToVtkConverter<1,1>;
template class Hdf5ToVtkConverter<1,2>;
template class Hdf5ToVtkConverter<2,2>;
template class Hdf5ToVtkConverter<1,3>;
template class Hdf5ToVtkConverter<2,3>;
template class Hdf5ToVtkConverter<3,3>;
