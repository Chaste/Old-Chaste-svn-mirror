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
        
        c_matrix<double, 2, 2> fibre_matrix;

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        c_matrix<double, 2, 2> correct_matrix = identity_matrix<double>(2,2);
//        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        correct_matrix(1,1) = -1.0;
//        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);

        fibre_reader.GetNextFibreSheetAndNormalMatrix(fibre_matrix);
        correct_matrix(0,1) = 1.0;
        correct_matrix(1,0) = 1.0;
//        TS_ASSERT_DELTA(UblasMatrixInfinityNorm<2>(fibre_matrix-correct_matrix), 0, 1e-9);
    }
};

#endif /*TESTFIBREREADER_HPP_*/
