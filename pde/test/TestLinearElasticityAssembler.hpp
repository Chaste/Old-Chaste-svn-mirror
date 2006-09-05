#ifndef TESTLINEARELASTICITYASSEMBLER_HPP_
#define TESTLINEARELASTICITYASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include <petsc.h>
#include "TrianglesMeshReader.cpp"
#include "LinearElasticityAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestLinearElasticityAssembler : public CxxTest::TestSuite
{
public:
    void Test3dExample(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        double E = 10;
        double nu = 0.3;
        c_vector<double,3> g;
        g(0) = -9.81;
        g(1) = 0; 
        g(2) = 0;
        
        BoundaryConditionsContainer<3,3,3> bcc(mesh.GetNumNodes());
        
        LinearElasticityAssembler<3> assembler(&mesh,&bcc);
        assembler.SetParameters(E,nu,g);

        Vec result = assembler.Solve(); // doesnt do anything as ComputeRhsTerm etc not yet filled it
        ReplicatableVector result_repl(result);
        
        for(int i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(result_repl[i], 0.0, 1e-6);
        }
    }
};


#endif /*TESTLINEARELASTICITYASSEMBLER_HPP_*/
