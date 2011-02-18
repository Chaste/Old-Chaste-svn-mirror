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
#ifndef TESTSTREETERFIBREGENERATOR_HPP_
#define TESTSTREETERFIBREGENERATOR_HPP_

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "NumericFileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestStreeterFibreGenerator : public CxxTest::TestSuite
{
public:

    void TestSimpleOrthotropic() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);

        fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "box_heart.ortho", true);

        OutputFileHandler handler("shorter_streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "box_heart.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/box_shaped_heart/box_heart.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }

    void TestSimpleOrthotropicNotDistributed() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);

        fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "box_heart_not_dist.ortho", true);

        OutputFileHandler handler("shorter_streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "box_heart_not_dist.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/box_shaped_heart/box_heart.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }
    
    void todoTestDownSampled() throw (Exception)
    {
        
        TrianglesMeshReader<3,3> mesh_reader("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
        std::string epi_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.epi";
        std::string rv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.rv";
        std::string lv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.lv";
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);

        fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "downsampled.ortho");

    }
    
    void TestExceptions()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);

        // No surfaces defined
        TS_ASSERT_THROWS_THIS(fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "file.fibres"),
                "Files defining the heart surfaces not set");

        // Wrong surface filename
        TS_ASSERT_THROWS_THIS(fibre_generator.SetSurfaceFiles("wrong_name", "wrong_name", "wrong_name", false),
                "Wrong surface definition file name wrong_name");
    }
    void TestOrientationException()
    {
        TrianglesMeshReader<3,3> mesh_reader("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
        std::string epi_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.epi";
        std::string rv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.rv";
        std::string lv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.lv";
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);

        TS_ASSERT_THROWS_THIS(fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "downsampled.ortho"),
            "Ventricular surfaces overlap too much in the y-axis");
        
    }
};

#endif /*TESTSTREETERFIBREGENERATOR_HPP_*/
