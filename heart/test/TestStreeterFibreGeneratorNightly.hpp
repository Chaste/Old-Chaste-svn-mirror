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
#ifndef TESTSTREETERFIBREGENERATORNIGHTLY_HPP_
#define TESTSTREETERFIBREGENERATORNIGHTLY_HPP_

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "MemfemMeshReader.hpp"
#include "NumericFileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestStreeterFibreGeneratorNightly : public CxxTest::TestSuite
{
public:
    void TestSimpleOrthotropic() throw (Exception)
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter_parallel", "point50.ortho", true);

        OutputFileHandler handler("streeter_parallel", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50.ortho";
        std::string wall_file = handler.GetOutputDirectoryFullPath() + "wall_thickness.data";

        NumericFileComparison comp_ortho(fibre_file,"heart/test/data/point50_heart_mesh/point50.ortho");
        TS_ASSERT(comp_ortho.CompareFiles(1e-11));
        NumericFileComparison comp_wall(wall_file,"heart/test/data/point50_heart_mesh/wall_thickness.data");
        TS_ASSERT(comp_wall.CompareFiles(1e-11));
    }

    void TestSimpleOrthotropicNotDistributed() throw (Exception)
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "point50_not_dist.ortho", true);

        OutputFileHandler handler("streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50_not_dist.ortho";
        std::string wall_file = handler.GetOutputDirectoryFullPath() + "wall_thickness.data";

        NumericFileComparison comp_ortho(fibre_file,"heart/test/data/point50_heart_mesh/point50.ortho");
        TS_ASSERT(comp_ortho.CompareFiles(1e-11));
        NumericFileComparison comp_wall(wall_file,"heart/test/data/point50_heart_mesh/wall_thickness.data");
        TS_ASSERT(comp_wall.CompareFiles(1e-11));
    }
        void TestDownSampledRabbit() throw (Exception)
    {
        
        TrianglesMeshReader<3,3> mesh_reader("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
        std::string epi_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.epi";
        std::string rv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.rv";
        std::string lv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.lv";
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "downsampled.ortho");

        OutputFileHandler handler("streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "downsampled.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/fibre_tests/downsampled.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
        
        //Output to VTK.  Output isn't tested, so a non-VTK compilation ought to work.
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "downsampled_fibres", false);
        FileFinder file("heart/test/data/fibre_tests/downsampled.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);
        TS_ASSERT_EQUALS(fibres.size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(second.size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(third.size(), mesh.GetNumElements());
        writer.AddCellData("OrthoFibres", fibres);
        writer.AddCellData("OrthoSecond", second);
        writer.AddCellData("OrthoThird", third);
        writer.WriteFilesUsingMesh(mesh);
    
    }
    

};

#endif /*TESTSTREETERFIBREGENERATORNIGHTLY_HPP_*/
