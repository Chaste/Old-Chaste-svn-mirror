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
#ifndef TESTSTREETERFIBREGENERATORNIGHTLY_HPP_
#define TESTSTREETERFIBREGENERATORNIGHTLY_HPP_

#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "MemfemMeshReader.hpp"
#include "../../global/test/NumericFileComparison.hpp"


class TestStreeterFibreGeneratorNightly : public CxxTest::TestSuite
{
public:
    void TestSimpleOrthotropic() throw (Exception)
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "point50.ortho", true);

        OutputFileHandler handler("streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/point50_heart_mesh/point50.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }

    };

#endif /*TESTSTREETERFIBREGENERATORNIGHTLY_HPP_*/
