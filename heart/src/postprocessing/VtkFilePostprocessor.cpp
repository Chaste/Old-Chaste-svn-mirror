/*

Copyright (C) Fujitsu Laboratories of Europe, 2009

*/

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


#ifdef CHASTE_VTK

#include "VtkFilePostprocessor.hpp"

#include <set>
#include <utility>

#include "VtkMeshReader.hpp"
#include "Exception.hpp"

VtkFilePostprocessor::VtkFilePostprocessor()
{

}

void VtkFilePostprocessor::CalculateEdgeLengthStatistics( std::string filePath,
                                                          double& rMinLength,
                                                          double& rMaxLength,
                                                          double& rMeanLength )
{
    VtkMeshReader<3,3> mesh_reader( filePath );
    vtkUnstructuredGrid* p_vtk_unstructured_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();

    rMinLength = DBL_MAX;
    rMaxLength = -DBL_MAX;
    rMeanLength = -DBL_MAX;

    double sum_of_edge_lengths = 0.0;
    int num_edges = 0;

    // Figure out the edges.
    std::set<std::pair<int, int> > edges;

    // Form node-node list.
    for(unsigned element = 0; element < (unsigned) p_vtk_unstructured_grid->GetNumberOfCells(); element++)
    {
        int k, l;
        vtkTetra *tetra = (vtkTetra *)p_vtk_unstructured_grid->GetCell(element);
        for(unsigned i=0; i<4; i++)
        {
          for(unsigned j=i+1;j<4;j++)
          {
            k = tetra->GetPointId(i);
            l = tetra->GetPointId(j);
            edges.insert(std::pair<int, int>(std::min(k, l), std::max(k, l)));
          }
        }
    }

    // Calculate the lengths of the edges and update min, max, sum
    for (std::set<std::pair<int, int> >::iterator edge = edges.begin(); edge != edges.end(); edge++)
    {
        num_edges++;

        double x1, y1, z1, x2, y2, z2, dx, dy, dz, this_edge_length;

        x1 = p_vtk_unstructured_grid->GetPoint(edge->first)[0];
        y1 = p_vtk_unstructured_grid->GetPoint(edge->first)[1];
        z1 = p_vtk_unstructured_grid->GetPoint(edge->first)[2];

        x2 = p_vtk_unstructured_grid->GetPoint(edge->second)[0];
        y2 = p_vtk_unstructured_grid->GetPoint(edge->second)[1];
        z2 = p_vtk_unstructured_grid->GetPoint(edge->second)[2];

        dx = x1 - x2;
        dy = y1 - y2;
        dz = z1 - z2;

        this_edge_length = sqrt(dx*dx + dy*dy + dz*dz);

        if (this_edge_length > rMaxLength)
        {
            rMaxLength = this_edge_length;
        }

        if (this_edge_length < rMinLength)
        {
            rMinLength = this_edge_length;
        }

        sum_of_edge_lengths += this_edge_length;
    }

    rMeanLength = sum_of_edge_lengths/num_edges;

}

unsigned VtkFilePostprocessor::LocateCoordinatesInElement( std::string filePath,
                                                           const std::vector<double> coords,
                                                           std::vector<double>& rWeights )
{
    VtkMeshReader<3,3> mesh_reader( filePath );
    vtkUnstructuredGrid* p_vtk_unstructured_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();

    unsigned element = LocateCoordinatesInElement( p_vtk_unstructured_grid, coords, rWeights );

    return element;
}

unsigned VtkFilePostprocessor::LocateCoordinatesInElement( vtkUnstructuredGrid* p_vtk_unstructured_grid,
                                                           const std::vector<double> coords,
                                                           std::vector<double>& rWeights )
{
    unsigned candidate_element = 0;

    int sub_id;
    double global_coords[3], local_coords[3], local_weights[4];
    double dist2;

    global_coords[0] = coords[0];
    global_coords[1] = coords[1];
    global_coords[2] = coords[2];

    while ( candidate_element < (unsigned) p_vtk_unstructured_grid->GetNumberOfCells() )
    {
        vtkTetra *tetra = (vtkTetra *)p_vtk_unstructured_grid->GetCell(candidate_element);

        int inside = tetra->EvaluatePosition( global_coords, NULL, sub_id, local_coords, dist2, local_weights );

        if ( inside == 1 && local_coords[0] > -1e-10 && local_coords[1] > -1e-10 && local_coords[2] > -1e-10 )
        {
            rWeights[0] = local_weights[0];
            rWeights[1] = local_weights[1];
            rWeights[2] = local_weights[2];
            rWeights[3] = local_weights[3];

            return candidate_element;
        }

        candidate_element++;
    }

    // Replace this with exception
    return INT_MAX;
}

std::vector<double> VtkFilePostprocessor::GenerateActionPotentialAtPointOfStaticMesh(
                                std::string filePrefix,
                                unsigned hiFileIndex,
                                std::vector<double> coords )
{
    vtkUnstructuredGrid* p_vtk_unstructured_grid;
    std::vector<double> rWeights(4), action_potential_at_coords(hiFileIndex+1);
    unsigned element = 0u;
    double vm_at_coords = 0.0;
    std::string data_name = "Vm";

    for (unsigned file = 0; file <= hiFileIndex; file++)
    {
        std::ostringstream next_file_name;
        next_file_name.str("");
        next_file_name << filePrefix << std::setw(4) << std::setfill('0') << file << ".vtu";
        VtkMeshReader<3,3> mesh_reader( next_file_name.str() );
        p_vtk_unstructured_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();

        if (file == 0)
        {
            element = LocateCoordinatesInElement( p_vtk_unstructured_grid, coords, rWeights );
        }

        vtkPointData *p_point_data = p_vtk_unstructured_grid->GetPointData();
        if ( !p_point_data->HasArray(data_name.c_str()) )
        {
            EXCEPTION("No point data '" + data_name + "'");
        }
        vtkDataArray *p_scalars = p_point_data->GetArray( data_name.c_str() );

        vm_at_coords = 0.0;
        for (unsigned node = 0; node < 4; node++)
        {
            vm_at_coords += rWeights[node] * p_scalars->GetTuple( p_vtk_unstructured_grid->GetCell(element)->GetPointId(node) )[0];
        }
        action_potential_at_coords[file] = vm_at_coords;

        mesh_reader.Initialize();
    }

    return action_potential_at_coords;
}

std::vector<double> VtkFilePostprocessor::GenerateActionPotentialAtPointOfAdaptingMesh(
                                std::string filePrefix,
                                unsigned hiFileIndex,
                                std::vector<double> coords )
{
    vtkUnstructuredGrid* p_vtk_unstructured_grid;
    std::vector<double> rWeights(4), action_potential_at_coords(hiFileIndex+1);
    unsigned element;
    double vm_at_coords = 0.0;
    std::string data_name = "Vm";

    for (unsigned file = 0; file <= hiFileIndex; file++)
    {
        std::ostringstream next_file_name;
        next_file_name.str("");
        next_file_name << filePrefix << std::setw(4) << std::setfill('0') << file << ".vtu";
        VtkMeshReader<3,3> mesh_reader( next_file_name.str() );
        p_vtk_unstructured_grid = mesh_reader.OutputMeshAsVtkUnstructuredGrid();
        element = LocateCoordinatesInElement( p_vtk_unstructured_grid, coords, rWeights );

        vtkPointData *p_point_data = p_vtk_unstructured_grid->GetPointData();
        if ( !p_point_data->HasArray(data_name.c_str()) )
        {
            EXCEPTION("No point data '" + data_name + "'");
        }
        vtkDataArray *p_scalars = p_point_data->GetArray( data_name.c_str() );

        vm_at_coords = 0.0;
        for (unsigned node = 0; node < 4; node++)
        {
            vm_at_coords += rWeights[node] * p_scalars->GetTuple( p_vtk_unstructured_grid->GetCell(element)->GetPointId(node) )[0];
        }
        action_potential_at_coords[file] = vm_at_coords;

        mesh_reader.Initialize();
    }

    return action_potential_at_coords;
}

#endif
