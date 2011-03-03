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


#ifndef TESTFIBREREADER_HPP_
#define TESTFIBREREADER_HPP_

#include <cxxtest/TestSuite.h>

#include "FibreReader.hpp"
#include "HeartFileFinder.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"

// simple helper function
template<unsigned DIM>
double UblasMatrixInfinityNorm(c_matrix<double,DIM,DIM> mat)
{
    double ret = fabs(mat(0,0));
    for(unsigned i=0; i<DIM; i++)
    {
        for(unsigned j=0; j<DIM; j++)
        {
            if( fabs(mat(i,j)) > ret )
            {
                ret = fabs(mat(i,j));
            }
        }
    }
    return ret;
}



class TestFibreReader : public CxxTest::TestSuite
{
public:
    void TestOrthoReaderSetup()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/random_fibres.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader(file_finder, ORTHO);
        
        TS_ASSERT_EQUALS(fibre_reader.GetNumLinesOfData(), 4u);

        c_matrix<double, 2, 2> fibre_matrix;

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        c_matrix<double, 2, 2> correct_matrix = identity_matrix<double>(2,2);
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        correct_matrix(1,1) = -1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        // this one isn't orthogonal - the false prevents this being checked
        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix, false);
        correct_matrix(0,1) = 1.0;
        correct_matrix(1,0) = 1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);
        
        // Non-symmetrical test case, standard rotation matrix
        // [cos(theta)  sin(theta)]
        // [-sin(theta) cos(theta)]
        correct_matrix(0,0) =  1.0/sqrt(2);  // fibre0
        correct_matrix(1,0) = -1.0/sqrt(2);  // fibre1
        correct_matrix(0,1) =  1.0/sqrt(2);  // sheet0
        correct_matrix(1,1) =  1.0/sqrt(2);  // sheet1
        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        // next matrix is not orthogonal, here we make sure this is checked
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "not orthogonal")

        // called too many times
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "End of file")
    }

    void TestAxiReaderSetup()
    {
        FileFinder file_finder("heart/test/data/fibre_tests/random_fibres.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader(file_finder, AXISYM);

        TS_ASSERT_EQUALS(fibre_reader.GetNumLinesOfData(), 3u);

        c_vector<double, 2> fibre_vector;

        fibre_reader.GetNextFibreVector(fibre_vector);
        TS_ASSERT_DELTA(fibre_vector(0), 0, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), 1, 1e-9);

        fibre_reader.GetNextFibreVector(fibre_vector);
        TS_ASSERT_DELTA(fibre_vector(0), 0.6, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), -0.8, 1e-9);

        // next vector is not normalised, the false below prevents this being checked 
        fibre_reader.GetNextFibreVector(fibre_vector, false);
        TS_ASSERT_DELTA(fibre_vector(0), 2, 1e-9);
        TS_ASSERT_DELTA(fibre_vector(1), 6, 1e-9);

        // next vector is not normalised, here we make sure this is checked
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetNextFibreVector(fibre_vector), "not normalised")

        // called too many times
        TS_ASSERT_THROWS_CONTAINS(fibre_reader.GetNextFibreVector(fibre_vector), "End of file")
    }
    
    void TestFibretoVtk()
    {
#ifdef CHASTE_VTK
        //See TestConductivityTensors
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "simple_fibres", false);

        {
            FileFinder file("heart/test/data/fibre_tests/SimpleAxisymmetric.axi", RelativeTo::ChasteSourceRoot);
            FibreReader<3> fibre_reader(file, AXISYM);
            std::vector< c_vector<double, 3> > fibres;
            fibre_reader.GetAllAxi(fibres);
            TS_ASSERT_EQUALS(fibres.size(), mesh.GetNumElements());
            writer.AddCellData("AxiFibres", fibres);
        }
        {
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
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
        }
        writer.WriteFilesUsingMesh(mesh);
        //Check that it has been written
        OutputFileHandler handler("TestVtkMeshWriter", false); 
        std::ifstream vtk_file;
        std::string command = handler.GetOutputDirectoryFullPath()+"/simple_fibres.vtu";
        vtk_file.open(command.c_str());
        TS_ASSERT(vtk_file.is_open());
        vtk_file.close();
#endif
    }

    void TestFibreReaderExceptions()
    {
        c_matrix<double, 2, 2> fibre_matrix;

        // file doesn't exist
        FileFinder finder0("heart/test/data/fibre_tests/dgfsdgjdf.ortho", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS( FibreReader<2> fibre_reader(finder0,ORTHO), "Failed to open fibre file");

        // line for first element is incomplete
        FileFinder finder1("heart/test/data/fibre_tests/bad_ortho1.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader1(finder1, ORTHO);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader1.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "A line is incomplete in");

        // line for third element is missing
        FileFinder finder2("heart/test/data/fibre_tests/bad_ortho2.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader2(finder2, ORTHO);
        fibre_reader2.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        fibre_reader2.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader2.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "End of file");

        // line for second element has too many entries
        FileFinder finder3("heart/test/data/fibre_tests/bad_ortho3.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader3(finder3, ORTHO);
        fibre_reader3.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader3.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "Too many entries in a line in");

        // first line doesn't give the number of lines of data
        FileFinder finder4("heart/test/data/fibre_tests/bad_ortho4.ortho", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS( FibreReader<2> fibre_reader(finder4,ORTHO), "First (non comment) line of the fibre orientation file should contain the number of lines");

        // Wrong method call, can't read an 'orthotropic vector'
        c_vector<double, 2> fibre_vector;
        FileFinder finder5("heart/test/data/fibre_tests/random_fibres.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader5(finder5, ORTHO);
        TS_ASSERT_THROWS_THIS(fibre_reader5.GetNextFibreVector(fibre_vector), "Use GetNextFibreSheetAndNormalMatrix when reading orthotropic fibres");
        std::vector<c_vector<double,2> > v1;
        TS_ASSERT_THROWS_THIS(fibre_reader5.GetAllAxi(v1), "Use GetAllOrtho when reading orthotropic fibres");
        // wrong method call, can't read an 'axisymmetric matrix'
        FileFinder finder6("heart/test/data/fibre_tests/random_fibres.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader6(finder6, AXISYM);
        TS_ASSERT_THROWS_THIS(fibre_reader6.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "Use GetNextFibreVector when reading axisymmetric fibres");
        std::vector<c_vector<double,2> > v2;
        std::vector<c_vector<double,2> > v3;
        TS_ASSERT_THROWS_THIS(fibre_reader6.GetAllOrtho(v1, v2, v3), "Use GetAllAxi when reading axisymmetric fibres");

        // Incomplete axi data
        FileFinder finder7("heart/test/data/fibre_tests/bad_axi.axi", RelativeTo::ChasteSourceRoot);
        FibreReader<2> fibre_reader7(finder7, AXISYM);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader7.GetNextFibreVector(fibre_vector), "A line is incomplete in");
    }
};


#endif /*TESTFIBREREADER_HPP_*/
