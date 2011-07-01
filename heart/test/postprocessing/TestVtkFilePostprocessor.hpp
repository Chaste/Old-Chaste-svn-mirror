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


#ifndef TESTVTKFILEPOSTPROCESSOR_HPP_
#define TESTVTKFILEPOSTPROCESSOR_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <set>

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkTriangle.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellDerivatives.h>
#include <vtkDataSetToUnstructuredGridFilter.h>
#include <vtkCellDataToPointData.h>
#include <vtkExtractVectorComponents.h>
#include <vtkPointData.h>
#endif //CHASTE_VTK

#include "VtkFilePostprocessor.hpp"
#include "OutputFileHandler.hpp"

class TestVtkFilePostprocessor : public CxxTest::TestSuite
{
//Requires  "sudo aptitude install libvtk5-dev" or similar
public:

    void TestCalculateMinMaxAndAverageEdgeLength(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkFilePostprocessor postprocessor;
        double min_edge_length, max_edge_length, mean_edge_length;
        std::string file_path;

        file_path = "mesh/test/data/cube_2mm_12_elements.vtu";
        postprocessor.CalculateEdgeLengthStatistics( file_path, min_edge_length, max_edge_length, mean_edge_length );
        TS_ASSERT_DELTA( min_edge_length,  9.990700e-02, 1e-6 );
        TS_ASSERT_DELTA( max_edge_length,  2.828427e-01, 1e-6 );
        TS_ASSERT_DELTA( mean_edge_length, 1.860265e-01, 1e-6 );
#endif // CHASTE_VTK
    }

    void TestLocateCoordinatesInElement(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkFilePostprocessor postprocessor;

        std::string file_path = "mesh/test/data/cube_2mm_12_elements.vtu";

        std::vector<double> coords(3), weights(4);
        coords[0] = 0.1;
        coords[1] = 0.1;
        coords[2] = 0.1;

        unsigned element = postprocessor.LocateCoordinatesInElement( file_path, coords, weights );

        TS_ASSERT_EQUALS( element, 11u );
        TS_ASSERT_DELTA( weights[0], 0.499470000, 1e-6 );
        TS_ASSERT_DELTA( weights[1], 0.499735000, 1e-6 );
        TS_ASSERT_DELTA( weights[2], 0.000264618, 1e-9 );
        TS_ASSERT_DELTA( weights[3], 0.000529729, 1e-9 );
#endif // CHASTE_VTK
    }

    void TestGenerateActionPotentialAtPointOfStaticMesh(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkFilePostprocessor postprocessor;

        std::string file_path_head = "notforrelease/test/data/coarse_slab_neumann_no_adapt";
        int hi_file_index = 2;
        std::vector<double> coords(3), weights(4);
        coords[0] = 0.05;
        coords[1] = 0.05;
        coords[2] = 0.05;

        std::vector<double> action_potential =
                        postprocessor.GenerateActionPotentialAtPointOfStaticMesh( file_path_head,
                                                                                  hi_file_index,
                                                                                  coords );

        TS_ASSERT_DELTA( action_potential[0], -83.8530, 1e-4 );
        TS_ASSERT_DELTA( action_potential[1], -83.8567, 1e-4 );
        TS_ASSERT_DELTA( action_potential[2], -83.7268, 1e-4 );

        OutputFileHandler output_file_handler("TestVtkFilePostprocessor", false);
        std::string file_name = "no_adapt_action_potential.txt";
        out_stream p_node_file = output_file_handler.OpenOutputFile(file_name);

        for (unsigned index = 0; index < action_potential.size(); index++)
        {
            *p_node_file << action_potential[index] << "\n";
        }
#endif // CHASTE_VTK
    }

    void TestGenerateActionPotentialAtPointOfAdaptingMesh(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkFilePostprocessor postprocessor;

        std::string file_path_head = "notforrelease/test/data/coarse_slab_neumann";
        int hi_file_index = 2;
        std::vector<double> coords(3), weights(4);
        coords[0] = 0.05;
        coords[1] = 0.05;
        coords[2] = 0.05;

        std::vector<double> action_potential =
                        postprocessor.GenerateActionPotentialAtPointOfAdaptingMesh( file_path_head,
                                                                                    hi_file_index,
                                                                                    coords );

        TS_ASSERT_DELTA( action_potential[0], -83.8530, 1e-4 );
        TS_ASSERT_DELTA( action_potential[1], -83.8567, 1e-4  );
        TS_ASSERT_DELTA( action_potential[2], -83.4991, 1e-4  );

        OutputFileHandler output_file_handler("TestVtkFilePostprocessor", false);
        std::string file_name = "adaptive_action_potential.txt";
        out_stream p_node_file = output_file_handler.OpenOutputFile(file_name);

        for (unsigned index = 0; index < action_potential.size(); index++)
        {
            *p_node_file << action_potential[index] << "\n";
        }
#endif // CHASTE_VTK
    }

    void TestExceptionIfNoVoltageStoredInMesh(void) throw(Exception)
    {
#ifdef CHASTE_VTK
        VtkFilePostprocessor postprocessor;

        std::string file_path_head = "notforrelease/test/data/twin_flow";
        int hi_file_index = 0;
        std::vector<double> coords(3), weights(4);
        coords[0] = 0.05;
        coords[1] = 0.05;
        coords[2] = 0.05;

        TS_ASSERT_THROWS_THIS( std::vector<double> action_potential =
                                postprocessor.GenerateActionPotentialAtPointOfStaticMesh(
                                file_path_head, hi_file_index,  coords ),
                                "No point data 'Vm'");

        TS_ASSERT_THROWS_THIS( std::vector<double> action_potential =
                                postprocessor.GenerateActionPotentialAtPointOfAdaptingMesh(
                                file_path_head, hi_file_index,  coords ),
                                "No point data 'Vm'");
#endif // CHASTE_VTK
    }

};

#endif /*TESTVTKFILEPOSTPROCESSOR_HPP_*/
