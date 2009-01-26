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


#ifndef TESTDEALIILINEARSYSTEM_HPP_
#define TESTDEALIILINEARSYSTEM_HPP_


#include <cxxtest/TestSuite.h>
#include "DealiiLinearSystem.hpp"
#include "TetrahedralMesh.hpp"

class TestDealiiLinearSystem : public CxxTest::TestSuite
{
public :
    void TestMe() throw(Exception)
    {
        unsigned size = 101;
        double h = 1.0/(size-1);

        DealiiLinearSystem system(size);

        system.mLhsMatrix.add(0, 0, 1.0);
        system.mLhsMatrix.add(size-1, size-1, 1.0);

        for(unsigned i=1; i<size-1; i++)
        {
            system.mLhsMatrix.add(i, i-1,-1.0/h);
            system.mLhsMatrix.add(i, i,   2.0/h);
            system.mLhsMatrix.add(i, i+1,-1.0/h);

            system.mRhsVector(i) = h;
        }
        
        system.Solve();
        
        for(unsigned i=0; i<size; i++)
        {
            double x = h*i;
            double u = 0.5*x*(1-x);
            TS_ASSERT_DELTA(system.rGetLhsVector()(i), u, 1e-12); // finite element solution is exact at the nodes
        }
    }
    
    void TestMoreOfMe() throw(Exception)
    {
        unsigned size = 11;
        double h = 1.0/10;

        DealiiLinearSystem system(size);

        system.mRhsVector(4) = exp(1);
        system.ZeroRhsVector();

        system.mLhsMatrix.add(2, 1, 3.1415);
        system.ZeroLhsMatrix();

        unsigned num_elem = size-1;
        for(unsigned i=0; i<num_elem; i++)
        {
            unsigned nodes[2];
            nodes[0] = i;
            nodes[1] = i+1;

            c_vector<double,2> b;
            b(0) = h/2;
            b(1) = h/2;
        
            c_matrix<double,2,2> a;
            a(0,0) = 1/h;
            a(0,1) = a(1,0) = -1/h;
            a(1,1) = 1/h;
            
            ///\todo This code (added in r4597) has never compiled
         //   system.AddLhsMultipleValues(nodes,a);
         //   system.AddRhsMultipleValues(nodes,b);
        }

        system.ZeroMatrixRow(0);
        system.SetMatrixElement(0,0,1.0);
        system.SetRhsVectorElement(0,0.0);

        system.ZeroMatrixRow(size-1);
        system.SetMatrixElement(size-1,size-1,1.0);
        system.SetRhsVectorElement(size-1,0.0);

        system.Solve();
        
        for(unsigned i=0; i<size; i++)
        {
            double x = h*i;
            double u = 0.5*x*(1-x);
            TS_ASSERT_DELTA(system.rGetLhsVector()(i), u, 1e-12); // finite element solution is exact at the nodes
        }
    }
    

    void TestWithMesh() throw(Exception)
    {
        QuadraticMesh<2> quad_mesh("mesh/test/data/square_128_elements_quadratic");

        // difficult to test - have printed and looked at sparsity pattern, looks correct
        DealiiLinearSystem system(quad_mesh);
    }
};


#endif /*TESTDEALIILINEARSYSTEM_HPP_*/
