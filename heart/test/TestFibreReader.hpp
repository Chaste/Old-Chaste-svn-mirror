/*

Copyright (C) University of Oxford, 2005-2010

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
    void TestFibreReaderSetup()
    {
        cp::path_type path("heart/test/data/random_fibres.ortho");
        path.relative_to(cp::relative_to_type::chaste_source_root);
        HeartFileFinder file_finder(path);
        
        FibreReader<2> fibre_reader(file_finder);
        
        TS_ASSERT_EQUALS(fibre_reader.mNumElements, 3u);

        c_matrix<double, 2, 2> fibre_matrix;

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        c_matrix<double, 2, 2> correct_matrix = identity_matrix<double>(2,2);
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        correct_matrix(1,1) = -1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        correct_matrix(0,1) = 1.0;
        correct_matrix(1,0) = 1.0;
        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);
    }

    void TestFibreReaderExceptions()
    {
        c_matrix<double, 2, 2> fibre_matrix;

        // file doesn't exist
        cp::path_type  bad_path("heart/test/data/dgfsdgjdf.ortho");
        bad_path.relative_to(cp::relative_to_type::chaste_source_root);
        HeartFileFinder bad_file_finder(bad_path);
        TS_ASSERT_THROWS_CONTAINS( FibreReader<2> fibre_reader(bad_file_finder), "Failed to open");

        // line for first element is incomplele
        cp::path_type  path1("heart/test/data/bad_ortho1.ortho");
        path1.relative_to(cp::relative_to_type::chaste_source_root);
        HeartFileFinder finder1(path1);
        FibreReader<2> fibre_reader1(finder1);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader1.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "A line is incomplete in");

        // line for third element is missing
        cp::path_type  path2("heart/test/data/bad_ortho2.ortho");
        path2.relative_to(cp::relative_to_type::chaste_source_root);
        HeartFileFinder finder2(path2);
        FibreReader<2> fibre_reader2(finder2);
        fibre_reader2.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        fibre_reader2.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader2.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "Fibre orientation file contains less");

        // line for second element has too many entries
        cp::path_type  path3("heart/test/data/bad_ortho3.ortho");
        path3.relative_to(cp::relative_to_type::chaste_source_root);
        HeartFileFinder finder3(path3);
        FibreReader<2> fibre_reader3(finder3);
        fibre_reader3.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        TS_ASSERT_THROWS_CONTAINS(fibre_reader3.GetNextFibreSheetAndNormalMatrix(fibre_matrix), "Too many entries in a line in");

        // first line doesn't give the number of elements
        cp::path_type  path4("heart/test/data/bad_ortho4.ortho");
        path4.relative_to(cp::relative_to_type::chaste_source_root);
        HeartFileFinder finder4(path4);
        TS_ASSERT_THROWS_CONTAINS( FibreReader<2> fibre_reader(finder4), "First (non comment) line of the fibre orientation file should contain the number of elements");
    }
};


#endif /*TESTFIBREREADER_HPP_*/
