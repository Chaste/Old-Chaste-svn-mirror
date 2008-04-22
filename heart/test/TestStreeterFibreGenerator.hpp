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
                
        TS_ASSERT_THROWS_ANYTHING(fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "ortho.fibres"));
        
        
    }
};

#endif /*TESTSTREETERFIBREGENERATOR_HPP_*/
