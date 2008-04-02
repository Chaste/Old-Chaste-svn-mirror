/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TESTVOLUMECALCULATOR_HPP_
#define TESTVOLUMECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "ConformingTetrahedralMesh.cpp"
#include "TrianglesMeshReader.cpp"
#include <cmath>

#include <vector>

class TestVolumeCalculator : public CxxTest::TestSuite
{
public:

    void CheckVolume(char* meshFile, double expectedVolume, double errorTolerance)
    {
        TrianglesMeshReader<3,3> meshReader(meshFile);
        ConformingTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(meshReader);
        
        TS_ASSERT_DELTA(mesh.CalculateMeshVolume()/expectedVolume, 1.0, errorTolerance);
    }
    
    void TestCube(void)
    {
        CheckVolume("mesh/test/data/cube_136_elements", 1.0, 1e-6);
    }
    
    void TestCube2(void)
    {
        CheckVolume("mesh/test/data/cube_1626_elements", 1.0, 1e-6);
    }
    
    void TestCylinderWithHole(void)
    {
        CheckVolume("mesh/test/data/cylinder_with_hole_840_elements", (15*5*5-4.0/3.0)*M_PI, 5e-2);
    }
    
    void TestCalculate1DLengthIn1DSpace()
    {
        // Calculate length of non-uniform mesh
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_4_non_uniform_elements");
        
        ConformingTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        double mesh_length = mesh.CalculateMeshVolume();
        
        TS_ASSERT_DELTA(mesh_length, 1.0, 1.0e-5);
    }
    
    
    void Test2DMeshArea()
    {
        // Tests area calculation for a square 2D mesh
        TrianglesMeshReader<2,2> square_mesh_reader("mesh/test/data/2D_0_to_1mm_200_elements");
        
        ConformingTetrahedralMesh<2,2> mesh_square;
        
        mesh_square.ConstructFromMeshReader(square_mesh_reader);
        
        double mesh_square_area = mesh_square.CalculateMeshVolume();
        
        TS_ASSERT_DELTA(mesh_square_area, 0.01, 1e-6);
        
        // Tests area calculation for a annular 2D mesh
        TrianglesMeshReader<2,2> annuluar_mesh_reader("mesh/test/data/annulus_256_elements");
        
        ConformingTetrahedralMesh<2,2> mesh_annulus;
        
        mesh_annulus.ConstructFromMeshReader(annuluar_mesh_reader);
        
        double mesh_annulus_area = mesh_annulus.CalculateMeshVolume();
        
        TS_ASSERT_DELTA(mesh_annulus_area, 8.0*M_PI, 2e-1);
    }
    
    
    void TestCubeSurfaceArea(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get a iterator over the boundary elements of the 3D mesh
        ConformingTetrahedralMesh<3, 3>::BoundaryElementIterator boundary_iter =
            mesh.GetBoundaryElementIteratorBegin();
            
        double area=0.0;
        while (boundary_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            // There are 16 regular triangles on each face of the cube
            // since determinant is double the area we are expecting (1/8).
            TS_ASSERT_DELTA((*boundary_iter)->GetJacobianDeterminant(), 1/8.0, 1e-6);
            area+=(*boundary_iter)->GetJacobianDeterminant()/2.0;
            boundary_iter++;
            
        }
        // It is a unit cube, 6 faces.
        TS_ASSERT_DELTA(area, 6.0, 1e-6);
    }
    
    
    void TestCubeSurfaceAreaMoreElements(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        
        ConformingTetrahedralMesh<3,3> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Get a iterator over the boundary elements of the 3D mesh
        ConformingTetrahedralMesh<3, 3>::BoundaryElementIterator boundary_iter =
            mesh.GetBoundaryElementIteratorBegin();
            
        double area=0.0;
        while (boundary_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            area+=(*boundary_iter)->GetJacobianDeterminant()/2.0;
            boundary_iter++;
            
        }
        // It is a unit cube, 6 faces.
        TS_ASSERT_DELTA(area, 6.0, 1e-6);
    }
    
    void TestDiskPerimeter(void)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh;
        
        mesh.ConstructFromMeshReader(mesh_reader);
        // Get a iterator over the boundary elements of the 2D mesh
        ConformingTetrahedralMesh<2,2>::BoundaryElementIterator boundary_iter =
            mesh.GetBoundaryElementIteratorBegin();
            
        double perimeter=0.0;
        while (boundary_iter != mesh.GetBoundaryElementIteratorEnd())
        {
            perimeter+=(*boundary_iter)->GetJacobianDeterminant();
            boundary_iter++;
            
        }
        // It is a unit disk radius 1, perimeter 2Pi (approximated by triangles)
        TS_ASSERT_DELTA(perimeter, 2.0*M_PI, 2e-3);
    }
    
    //Repeat the tests above
    void TestSurfaceWithMethod(void)
    {
        TrianglesMeshReader<3,3> mesh_reader_1("mesh/test/data/cube_136_elements");
        ConformingTetrahedralMesh<3,3> mesh_1;
        mesh_1.ConstructFromMeshReader(mesh_reader_1);
        TS_ASSERT_DELTA(mesh_1.CalculateMeshSurface(),6.0,1e-6);
        TrianglesMeshReader<3,3> mesh_reader_2("mesh/test/data/cube_1626_elements");
        ConformingTetrahedralMesh<3,3> mesh_2;
        mesh_2.ConstructFromMeshReader(mesh_reader_2);
        TS_ASSERT_DELTA(mesh_2.CalculateMeshSurface(),6.0,1e-6);
        
        TrianglesMeshReader<2,2> mesh_reader_3("mesh/test/data/disk_984_elements");
        ConformingTetrahedralMesh<2,2> mesh_3;
        mesh_3.ConstructFromMeshReader(mesh_reader_3);
        TS_ASSERT_DELTA(mesh_3.CalculateMeshSurface(), 2*M_PI, 2e-3);
        
        TrianglesMeshReader<1,1> mesh_reader_4("mesh/test/data/1D_0_to_1_200_elements");
        ConformingTetrahedralMesh<1,1> mesh_4;
        mesh_4.ConstructFromMeshReader(mesh_reader_4);
        TS_ASSERT_DELTA(mesh_4.CalculateMeshSurface(),0.0,1e-6);
    }
    
};

#endif /*TESTVOLUMECALCULATOR_HPP_*/
