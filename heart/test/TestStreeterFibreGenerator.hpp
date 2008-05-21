/*

Copyright (C) University of Oxford, 2008

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

#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "../../global/test/NumericFileComparison.hpp"


class TestStreeterFibreGenerator : public CxxTest::TestSuite
{
public:
    void TestSimpleOrthotropic() throw (Exception)
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";        
                
        ConformingTetrahedralMesh<3,3> mesh;               
        mesh.ConstructFromMeshReader(mesh_reader);
        
        StreeterFibreGenerator<3> fibre_generator(mesh);        
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file);
                        
        fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "ortho.fibres", true);

        OutputFileHandler handler("streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "ortho.fibres";
        
        //TS_ASSERT_EQUALS(system(("ndiff  -abserr 1e-11 " + fibre_file + " heart/test/data/streeter_point50_heart_mesh.ortho").c_str()), 0);        
       
        NumericFileComparison comp(fibre_file,"heart/test/data/streeter_point50_heart_mesh.ortho");
        TS_ASSERT(comp.CompareFiles());
    }    
    
    void TestExceptions()
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
                
        ConformingTetrahedralMesh<3,3> mesh;               
        mesh.ConstructFromMeshReader(mesh_reader);
        
        StreeterFibreGenerator<3> fibre_generator(mesh);        

        // No surfaces defined
        TS_ASSERT_THROWS_ANYTHING(fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "file.fibres"));
        
        // Wrong surface filename
        fibre_generator.SetSurfaceFiles("wrong_name", "wrong_name", "wrong_name");        
        TS_ASSERT_THROWS_ANYTHING(fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "file.fibres"));
        
        // Wrong surface format
        std::string wrong_face_file = "heart/test/data/point50_heart_mesh/wrong_format.tri";      
        fibre_generator.SetSurfaceFiles(wrong_face_file, wrong_face_file, wrong_face_file);
                
        TS_ASSERT_THROWS_ANYTHING(fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "file.fibres"));
        
        
    }
};

#endif /*TESTSTREETERFIBREGENERATOR_HPP_*/
